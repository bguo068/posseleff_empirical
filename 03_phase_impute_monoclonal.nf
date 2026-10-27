nextflow.enable.types = true

// Parameters used when this script is run on its own (see the entry workflow
// at the bottom). When included by main.nf, main.nf's params block is used.
params {
    beagle_options_monoclonal: String = ''
    beagle_jar: Path = "${projectDir}/lib/beagle.08Feb22.fa4.jar"
    gmap: Path = "${projectDir}/assets/Pf3D7_all.txt"
    deploid_nopanel_options: String = "-rate 8 -burn 0.67 -k 4"
}

record PopVcf {
    grp: String
    vcf: Path
    fws: Path
}

record MonoclonalParams {
    deploid_nopanel_options: String
    beagle_jar: Path
    gmap: Path
    beagle_options: String
}

// All files produced for one population (group). `imputed_vcf` contains the
// phased and imputed monoclonal samples.
record MonoclonalPop {
    grp: String
    plaf: Path
    sample_vcfs: Set<Path>
    deploid_vcfs: List<Path>
    deploid_csis: List<Path>
    collate_vcfs: List<Path>
    collate_csis: List<Path>
    merged_vcf: Path
    imputed_vcf: Path
}

// This process serves two purposes: (1) calculating the population-level allele
// frequencies and (2) splitting the multi-sample VCF file into multiple
// single-sample VCF files.
//
// Note: If there are too many samples in the input VCF file, bcftools +split
// may not be able to split all of them. This is because it opens a separate
// VCF.gz file for each sample, which can cause an error. To avoid this, it is
// better to use the following method:
// - Create a text file 'sample_list.txt' with the names of all the samples in
// the VCF file, and then
// - run the following command:  split -l 399 sample_list.txt group_split_
process BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF {
    tag "${grp}"

    input:
    record(grp: String, fws: Path, vcf: Path)

    output:
    record(grp: grp, plaf: file("${grp}_plaf.tsv"), sample_vcfs: files("single_sample_vcfs/*.vcf.gz"))

    script:
    """
    calc_plaf_from_vcf.py --vcf ${vcf} --out ${grp}_plaf.tsv

    awk '\$2>=0.95 {print \$1}' ${fws} > monoclonal_samples.txt
    split -l 400 monoclonal_samples.txt split_
    bcftools index ${vcf}
    mkdir single_sample_vcfs
    for SPLIT in split_*
    do
        bcftools +split ${vcf} -S \${SPLIT} -Oz -o single_sample_vcfs
    done
    """

    stub:
    """
    touch ${grp}_plaf.tsv
    mkdir single_sample_vcfs
    touch single_sample_vcfs/${grp}_{0..100}.vcf.gz
    """
}

// This process is used to phase a single monoclonal sample without the use of a
// panel, and it involves the following steps:
// (1) Cleaning up sites with missing genotypes.
// (2) Filtering out high leverage SNPs (see Zhu et al. 2018 paper).
// (3) Subsetting the PLAF file so that it contains the same SNP sites as the
// filtered single-sample VCF file.
// (4) Running dEploid with the -noPanel parameter.
// (5) Extracting haploid genomes from the dEploid result files.
// (6) Cleaning up and converting the haploid genomes into pseudo-homozygous
// diploid samples.

process DEPLOID_PHASE_MONOCLONAL {
    tag "${grp}:${sample}"

    input:
    record(grp: String, grp_size: Integer, sample: String, plaf: Path, vcf: Path)
    deploid_options: String

    output:
    record(
        grp: grp,
        grp_size: grp_size,
        sample: sample,
        vcf: file("${grp}:${sample}_out.vcf.gz"),
        csi: file("${grp}:${sample}_out.vcf.gz.csi")
    )

    script:
    def prefix = "${grp}:${sample}"
    """
    set -uexo pipefail
    # remove missing sites if nopanel; do want gt infor at missing site to be used for making panel
    bcftools view -i "F_MISSING==0" ${vcf} -Oz -o tmp.vcf.gz

    # find regions to exclude
    deploid_filter_high_leverage_snps.py --vcf tmp.vcf.gz --sites pos_to_exclude.txt

    # make sure plaf have the same number of sites as in vcf
    deploid_subset_plaf.py --in_plaf ${plaf} --in_vcf tmp.vcf.gz --out_plaf plaf.subset.tmp.txt

    dEploid -vcf tmp.vcf.gz -plaf plaf.subset.tmp.txt -noPanel -vcfOut -o ${prefix} \
        ${deploid_options} -exclude pos_to_exclude.txt
    rm -f tmp.vcf.gz

    # generate a .tsv file, *valid.txt file, which.txt, new_sample_name.txt file
    deploid_get_hap.py --sample ${sample} --hap ${prefix}.hap --log ${prefix}.log

    VALID="`cat *~valid.txt`"

    sed '/^##DEploid/d' ${prefix}.vcf > tmp.vcf

    # select haplotype and rename sample name
    bcftools view tmp.vcf -S which.txt | bcftools reheader -s new_sample_name.txt | bcftools view -Oz -o tmp2.vcf.gz

    # convert haploid to homozygous diploid
    gzip -dc tmp2.vcf.gz \
        | awk -v OFS='\t' -v FS='\t' '/^#/{print \$0} /^[^#]/{\$10=\$10 "|" \$10; print \$0}' \
        | bgzip -c > tmp3.vcf.gz

    # convert to bcftools format
    bcftools view tmp3.vcf.gz -Oz -o ${prefix}_out.vcf.gz
    bcftools index ${prefix}_out.vcf.gz

    # clean up
    rm tmp.vcf tmp2.vcf.gz tmp3.vcf.gz
    """

    stub:
    def prefix = "${grp}:${sample}"
    """
    touch  ${prefix}_out.vcf.gz ${prefix}_out.vcf.gz.csi
    """
}

// The following two processes are used to merge many single-sample phased VCF
// files. The reason for splitting the merge into two processes is that merging
// them all together in a single process can cause a failure in bcftools due to
// opening too many files at once. These two processes allow the merging to be
// done at two levels, so that each level only opens a relatively small number
// of single-sample VCF files.
//
// NOTE: stdin is the list of VCF file names, in the same order as vcfs. This
// avoids explicitly pasting all file names on the command line.

// Level 1 merge

process BCFTOOLS_MERGE_SUBSUBPOP {
    tag "${grp}_${subgrp}"

    input:
    record(grp: String, subgrp: Integer, n: Integer, vcfs: List<Path>, csis: List<Path>, file_list: String)

    stage:
    stdin file_list

    output:
    record(
        grp: grp,
        subgrp: subgrp,
        n: n,
        vcf: file("${grp}_${subgrp}_merged.vcf.gz"),
        csi: file("${grp}_${subgrp}_merged.vcf.gz.csi")
    )

    script:
    if( vcfs.size() == 1 )  // in case there is only one sample per group
        """
        cat - > file_list.txt
        mv ${vcfs[0]} ${grp}_${subgrp}_merged.vcf.gz
        bcftools index  ${grp}_${subgrp}_merged.vcf.gz
        """
    else
        """
        cat - > file_list.txt
        bcftools merge -l file_list.txt -Oz -o ${grp}_${subgrp}_merged.vcf.gz
        bcftools index ${grp}_${subgrp}_merged.vcf.gz
        """

    stub:
    """
    touch  ${grp}_${subgrp}_merged.vcf{.gz,.gz.csi}
    """
}

// Level 2 merge
// See comments in the above process
process BCFTOOLS_MERGE_FILETER_SUBPOP {
    tag "${grp}"

    input:
    record(grp: String, vcfs: List<Path>, csis: List<Path>, file_list: String)

    stage:
    stdin file_list

    output:
    record(grp: grp, vcf: file("${grp}_merged.vcf.gz"))

    script:
    if( vcfs.size() == 1 )  // in case there is only one subgroup per group
        """
        cat - > file_list.txt; bcftools view -i "F_MISSING<0.3" ${vcfs[0]} -Oz -o ${grp}_merged.vcf.gz
        """
    else
        """
        cat - > file_list.txt; bcftools merge -l file_list.txt | bcftools view -i "F_MISSING<0.3" -Oz -o ${grp}_merged.vcf.gz
        """

    stub:
    """
    touch  ${grp}_merged.vcf.gz
    """
}


// This process calls the Beagle software to impute the missing genotypes of
// samples (without using a panel) within the VCF files that contain all phased
// single-sample VCF files.
process BEAGLE_IMPUTE {
    tag "${grp}"

    input:
    record(grp: String, vcf: Path)
    beagle_jar: Path
    gmap: Path
    beagle_options: String

    output:
    record(grp: grp, vcf: file("${grp}_imputed.vcf.gz"))

    script:
    """
    java -Xmx${task.memory.toGiga()}g -jar ${beagle_jar} gt=${vcf} map=${gmap} out=${grp}_imputed nthreads=${task.cpus} \\
        ${beagle_options}
    """

    stub:
    """
    touch  ${grp}_imputed.vcf.gz
    """
}


workflow WF_PHASE_IMPUTE_MONOCLONAL {
    take:
    pops: Channel<PopVcf>
    opts: MonoclonalParams

    main:
    ch_split = BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF(pops)

    // one record per single-sample VCF
    ch_vcf_plaf = ch_split.flatMap { r ->
        r.sample_vcfs.collect { vcf ->
            record(grp: r.grp, grp_size: r.sample_vcfs.size(), sample: vcf.simpleName, plaf: r.plaf, vcf: vcf)
        }
    }

    ch_phased_vcf = DEPLOID_PHASE_MONOCLONAL(ch_vcf_plaf, opts.deploid_nopanel_options)

    // group by population (emitted as soon as all samples of the group are
    // phased), then sort by sample name so the merge order is reproducible
    ch_grouped_phased_vcf = ch_phased_vcf
        .map { r -> tuple(r.grp, r.grp_size, r) }
        .groupBy()
        .map { grp, rs ->
            def sorted = rs.toSorted { a, b -> a.sample <=> b.sample }
            record(grp: grp, vcfs: sorted.collect { r -> r.vcf }.toList(), csis: sorted.collect { r -> r.csi }.toList())
        }

    // collate into chunks of at most 200 VCFs
    ch_collated_vcf = ch_grouped_phased_vcf.flatMap { r ->
        def vcf_chunks = r.vcfs.collate(200)
        def csi_chunks = r.csis.collate(200)
        vcf_chunks.withIndex().collect { vcfs, subgrp ->
            record(
                grp: r.grp,
                subgrp: subgrp,
                n: vcf_chunks.size(),
                vcfs: vcfs,
                csis: csi_chunks[subgrp],
                file_list: vcfs.collect { p -> p.name }.join('\n') + '\n'
            )
        }
    }

    // Merge level 1
    ch_merge_level1 = BCFTOOLS_MERGE_SUBSUBPOP(ch_collated_vcf)

    ch_grouped_level1 = ch_merge_level1
        .map { r -> tuple(r.grp, r.n, r) }
        .groupBy()
        .map { grp, rs ->
            def sorted = rs.toSorted { a, b -> a.subgrp <=> b.subgrp }
            def vcfs = sorted.collect { r -> r.vcf }.toList()
            record(
                grp: grp,
                vcfs: vcfs,
                csis: sorted.collect { r -> r.csi }.toList(),
                file_list: vcfs.collect { p -> p.name }.join('\n') + '\n'
            )
        }

    // Merge level 2
    ch_merge_level2 = BCFTOOLS_MERGE_FILETER_SUBPOP(ch_grouped_level1)

    ch_imputed = BEAGLE_IMPUTE(ch_merge_level2, opts.beagle_jar, opts.gmap, opts.beagle_options)

    // Gather every file produced for a group into a single record
    ch_pops = ch_split
        .join(ch_grouped_phased_vcf.map { r -> record(grp: r.grp, deploid_vcfs: r.vcfs, deploid_csis: r.csis) }, by: 'grp')
        .join(ch_grouped_level1.map { r -> record(grp: r.grp, collate_vcfs: r.vcfs, collate_csis: r.csis) }, by: 'grp')
        .join(ch_merge_level2.map { r -> record(grp: r.grp, merged_vcf: r.vcf) }, by: 'grp')
        .join(ch_imputed.map { r -> record(grp: r.grp, imputed_vcf: r.vcf) }, by: 'grp')

    emit:
    pops: Channel<MonoclonalPop> = ch_pops
}

// --------------------------- PUBLISHING ---------------------------------------------
// Where each output file goes, relative to the output directory. These are
// shared by the output blocks of this script and main.nf.

def monoclonal_pop_publish_paths(r: MonoclonalPop) -> Map {
    def d = '03_phase_impute_monoclonal'
    return [
        (r.plaf):         "${d}/01_single_sample_vcf_and_group_plaf/",
        (r.sample_vcfs):  "${d}/01_single_sample_vcf_and_group_plaf/",
        (r.deploid_vcfs): "${d}/02_deploid_monoclonal/",
        (r.deploid_csis): "${d}/02_deploid_monoclonal/",
        (r.collate_vcfs): "${d}/03_merge_collate/",
        (r.collate_csis): "${d}/03_merge_collate/",
        (r.merged_vcf):   "${d}/03_merge_filter/",
        (r.imputed_vcf):  "${d}/04_beagle_imputed/",
    ]
}

// --------------------- FOR TESTING -------------------------------
// will be ignored when included as subworkflow script
workflow {
    main:
    ch_fws = channel.fromPath("results/02_pop_fws/04_moimix_fws/*.fws.txt")
        .map { p -> record(grp: p.simpleName, fws: p) }
    ch_subpop_vcf = channel.fromPath("results/02_pop_fws/02_subpop_vcf/*.vcf.gz")
        .map { p -> record(grp: p.simpleName, vcf: p) }
    ch_pops = ch_subpop_vcf.join(ch_fws, by: 'grp')

    def opts = record(
        deploid_nopanel_options: params.deploid_nopanel_options,
        beagle_jar: params.beagle_jar,
        gmap: params.gmap,
        beagle_options: params.beagle_options_monoclonal,
    )
    mono = WF_PHASE_IMPUTE_MONOCLONAL(ch_pops, opts)

    publish:
    phase_impute_monoclonal = mono
}

output {
    phase_impute_monoclonal: Channel<MonoclonalPop> {
        path { r -> monoclonal_pop_publish_paths(r) }
    }
}

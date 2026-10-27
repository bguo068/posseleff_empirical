nextflow.enable.types = true

// Parameters used when this script is run on its own (see the entry workflow
// at the bottom). When included by main.nf, main.nf's params block is used.
params {
    beagle_options_polyclonal: String = 'overlap=3.0' // default is 2.0, cause error for some vcf files
    deploid_panel_size: Integer = 16
    beagle_jar: Path = "${projectDir}/lib/beagle.08Feb22.fa4.jar"
    gmap: Path = "${projectDir}/assets/Pf3D7_all.txt"
    deploid_panel_options: String = "-nSample 250 -rate 8 -burn 0.67 -k 4 -ibd"
}

record PopVcf {
    grp: String
    vcf: Path
    fws: Path
}

// Phased and imputed monoclonal samples of a group (from WF_PHASE_IMPUTE_MONOCLONAL)
record MonoclonalImputed {
    grp: String
    imputed_vcf: Path
}

record PolyclonalParams {
    deploid_panel_size: Integer
    deploid_panel_options: String
    beagle_jar: Path
    gmap: Path
    beagle_options: String
}

// Panel and single-sample VCFs of a group
record PolyclonalSplit {
    grp: String
    panel: Path
    plaf: Path
    sample_vcfs: Set<Path>
}

// All dEploid-IBD results of a group, including samples without a dominant clone
record PolyclonalDeploid {
    grp: String
    deploid_vcfs: List<Path>
    deploid_csis: List<Path>
}

// Merge and imputation results of a group with at least one dominant clone
record PolyclonalMerged {
    grp: String
    collate_vcfs: List<Path>
    collate_csis: List<Path>
    merged_vcf: Path
    imputed_vcf: Path
}

// A group without any polyclonal sample that has a dominant clone. Its final
// VCF is the monoclonal VCF.
record MonoclonalOnly {
    grp: String
    vcf: Path
}

record GroupVcf {
    grp: String
    vcf: Path
}

// This process is used to identify the most diverse samples from the
// phased & imputed monoclonal samples. These samples are then extracted and
// converted into a table format that is compatible with the dEploid-IBD
// program.
//
// Note that the `vcf` file contains both monoclonal and polyclonal samples,
// and it also includes Format/AD fields which are used for finding
// representatives. The intermediate file, `mono_vcf`, contains only
// monoclonal samples, with the sample names having a ~xxx suffix, and with GT
// imputed/phased and the AD field removed.
process ALLEL_GENERATE_PANEL {
    tag "${grp}"

    input:
    record(grp: String, fws: Path, vcf: Path, mono_vcf: Path)
    panel_size: Integer

    output:
    record(grp: grp, panel: file("${grp}_panel.tsv"))

    script:
    """
    awk '\$2>=0.95 {print \$1}' ${fws} > clonal_samples.txt
    deploid_find_representives.py --vcf ${vcf} --clonal_samples clonal_samples.txt --n_samples ${panel_size} \\
        --out ${grp}_panel_samples.txt
    bcftools index ${vcf}
    bcftools view -S ${grp}_panel_samples.txt ${vcf} -Oz -o panel.vcf.gz
    deploid_panel_vcf2table.py --vcf ${mono_vcf} --samples ${grp}_panel_samples.txt --panel ${grp}_panel.tsv
    """

    stub:
    """
    touch ${grp}_panel.tsv
    """
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
//
// NOTE: some group might not have polyclonal samples. Then `sample_vcfs` is
// empty and `no_polyclonal` is set.
process BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF {
    tag "${grp}"

    input:
    record(grp: String, fws: Path, vcf: Path)

    output:
    record(
        grp: grp,
        plaf: file("${grp}_plaf.tsv"),
        sample_vcfs: files("single_sample_vcfs/*.vcf.gz", optional: true),
        no_polyclonal: file('GROUP_HAS_NO_POLYCLONAL_SAMPLE', optional: true)
    )

    script:
    """
    calc_plaf_from_vcf.py --vcf ${vcf} --out ${grp}_plaf.tsv
    awk '\$2<0.95 {print \$1}' ${fws} > polyclonal_samples.txt
    mkdir single_sample_vcfs
    bcftools index ${vcf}

    # NOTE: some group might not have polyclonal samples
    if [ `cat  polyclonal_samples.txt | wc -l` == 0 ]; then
        touch GROUP_HAS_NO_POLYCLONAL_SAMPLE
    else
        split -l 400 polyclonal_samples.txt split_
        for SPLIT in split_*
        do
            bcftools +split ${vcf} -S \${SPLIT} -Oz -o single_sample_vcfs
        done
    fi
    """

    stub:
    """
    touch ${grp}_plaf.tsv
    mkdir single_sample_vcfs
    if [ ${grp} == "SAM" ]; then
       touch GROUP_HAS_NO_POLYCLONAL_SAMPLE
    else
        touch single_sample_vcfs/${grp}_{0..100}.vcf.gz
    fi
    """
}

// This process is used to phase a single polyclonal sample by using a panel,
//  and it involves the following steps:
// (1) Cleaning up sites not present in the panel
// (2) Filtering out high leverage SNPs (see Zhu et al. 2018 paper).
// (3) Subsetting the PLAF file so that it contains the same SNP sites as the
// filtered single-sample VCF file.
// (4) Running dEploid-IBD use the subset of diverse monoclonal samples (from
// the same population) as the panel.
// (5) Extracting haploid genomes from the dEploid result files.
// Note: the *valid.txt file is used to indicate whether this polyclonal sample
// contains a dominant haploid clone
// (6) Cleaning up and converting the haploid genomes into pseudo-homozygous
// diploid samples.

process DEPLOID_PHASE_POLYCLONAL {
    tag "${grp}:${sample}"

    input:
    record(grp: String, grp_size: Integer, sample: String, plaf: Path, vcf: Path, panel: Path)
    deploid_options: String

    output:
    record(
        grp: grp,
        grp_size: grp_size,
        sample: sample,
        vcf: file("${grp}:${sample}_out.vcf.gz"),
        csi: file("${grp}:${sample}_out.vcf.gz.csi"),
        valid: env('VALID')
    )

    script:
    def prefix = "${grp}:${sample}"
    """
    set -uexo pipefail
    # remove sites not in panel; otherwise dEploid will complain
    cat ${panel} | sed 1d | cut -f 1-2 > tmp_sites_to_keep.txt
    bcftools index ${vcf}
    bcftools view -R tmp_sites_to_keep.txt  ${vcf} -Oz -o tmp.vcf.gz
    # delete the index of the input vcf; it is not needed any more
    rm -f ${vcf}.csi

    # find regions to exclude
    deploid_filter_high_leverage_snps.py --vcf tmp.vcf.gz --sites pos_to_exclude.txt

    # make sure plaf have the same number of sites as in vcf
    deploid_subset_plaf.py --in_plaf ${plaf} --in_vcf tmp.vcf.gz --out_plaf plaf.subset.tmp.txt

    dEploid -vcf tmp.vcf.gz -plaf plaf.subset.tmp.txt -panel ${panel} -vcfOut -o ${prefix} \
        ${deploid_options} -exclude pos_to_exclude.txt
    rm -f tmp.vcf.gz

    # generate a .tsv file, *valid.txt file, which.txt, new_sample_name.txt file
    deploid_get_hap.py --sample ${sample} --hap ${prefix}.hap --log ${prefix}.log

    # valid criteria in get_hap.py: if second/first < 1/3 and first > 0.7:
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
    VALID=true
    """
}

// The following two processes are used to merge many single-sample phased VCF
// files. The reason for splitting the merge into two processes is that merging
// them all together in a single process can cause a failure in bcftools due to
// opening too many files at once. These two processes allow the merging to be
// done at two levels, so that each level only opens a relatively small number
// of single-sample VCF files.
//
// Note the level 2 merge also include the monoclonal samples of the sample group.
// After this process, both the monoclonal samples and the polyclonal samples with
// dominant haplotype are merged.
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

// Level 2 merge + monoclonal samples
// See comments in the above process
// NOTE: the last of `vcfs` is the monoclonal VCF, which has no index yet.
process BCFTOOLS_MERGE_WITH_MONOCLONALS_FILETER_SUBPOP {
    tag "${grp}"

    input:
    record(grp: String, vcfs: List<Path>, csis: List<Path>, file_list: String)

    stage:
    stdin file_list

    output:
    record(grp: grp, vcf: file("${grp}_merged.vcf.gz"))

    script:
    if( vcfs.size() == 1 )  // in case there is only one vcf per group
        """
        cat - > file_list.txt; bcftools view -i "F_MISSING<0.3" ${vcfs[0]} -Oz -o ${grp}_merged.vcf.gz
        """
    else
        """
        cat - > file_list.txt
        bcftools index -f ${vcfs.last()}
        bcftools merge -l file_list.txt | bcftools view -i "F_MISSING<0.3" -Oz -o ${grp}_merged.vcf.gz
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

workflow WF_PHASE_IMPUTE_POLYCLONAL {
    take:
    pops: Channel<PopVcf>
    mono: Channel<MonoclonalImputed>
    opts: PolyclonalParams

    main:
    ch_mono_vcf = mono.map { r -> record(grp: r.grp, mono_vcf: r.imputed_vcf) }

    ch_panel = ALLEL_GENERATE_PANEL(pops.join(ch_mono_vcf, by: 'grp'), opts.deploid_panel_size)

    ch_split = BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF(pops)

    // one record per single-sample VCF, joined with the panel of its group
    ch_vcf_plaf_panel = ch_split
        .filter { r -> !r.sample_vcfs.isEmpty() }
        .flatMap { r ->
            r.sample_vcfs.collect { vcf ->
                record(grp: r.grp, grp_size: r.sample_vcfs.size(), sample: vcf.simpleName, plaf: r.plaf, vcf: vcf)
            }
        }
        .join(ch_panel, by: 'grp')

    ch_phased_vcf = DEPLOID_PHASE_POLYCLONAL(ch_vcf_plaf_panel, opts.deploid_panel_options)

    // grouping, sorting and filtering (keep only samples with a dominant clone)
    ch_grouped_phased_vcf = ch_phased_vcf
        .map { r -> tuple(r.grp, r.grp_size, r) }
        .groupBy()
        .map { grp, rs ->
            def sorted = rs.toSorted { a, b -> a.sample <=> b.sample }
            def valid = sorted.findAll { r -> r.valid == 'true' }
            record(
                grp: grp,
                deploid_vcfs: sorted.collect { r -> r.vcf }.toList(),
                deploid_csis: sorted.collect { r -> r.csi }.toList(),
                vcfs: valid.collect { r -> r.vcf }.toList(),
                csis: valid.collect { r -> r.csi }.toList()
            )
        }

    ch_predominant_clones = ch_grouped_phased_vcf.filter { r -> r.vcfs.size() > 0 }

    // Groups whose final VCF is the monoclonal VCF:
    // 1. groups without polyclonal samples, which did not make it to the
    //    deconvolution of polyclonal samples
    // 2. groups with polyclonal samples, none of which has a dominant clone
    ch_mono_only = ch_split.filter { r -> r.no_polyclonal != null }.map { r -> r.grp }
        .mix(ch_grouped_phased_vcf.filter { r -> r.vcfs.size() == 0 }.map { r -> r.grp })
        .map { grp -> record(grp: grp) }
        .join(ch_mono_vcf, by: 'grp')
        .map { r -> record(grp: r.grp, vcf: r.mono_vcf) }

    // collate into chunks of at most 200 VCFs
    ch_collated_vcf = ch_predominant_clones.flatMap { r ->
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
            record(grp: grp, vcfs: sorted.collect { r -> r.vcf }.toList(), csis: sorted.collect { r -> r.csi }.toList())
        }

    // add the monoclonal VCF of the group as the last VCF to merge
    ch_level1_with_monoclonals = ch_grouped_level1
        .join(ch_mono_vcf, by: 'grp')
        .map { r ->
            def vcfs = r.vcfs + [r.mono_vcf]
            record(grp: r.grp, vcfs: vcfs, csis: r.csis, file_list: vcfs.collect { p -> p.name }.join('\n') + '\n')
        }

    // Merge level 2
    ch_merge_level2 = BCFTOOLS_MERGE_WITH_MONOCLONALS_FILETER_SUBPOP(ch_level1_with_monoclonals)

    ch_imputed = BEAGLE_IMPUTE(ch_merge_level2, opts.beagle_jar, opts.gmap, opts.beagle_options)

    ch_phased_imputed_polyclonal_monoclonal_vcf = ch_imputed.mix(ch_mono_only)

    // Gather the files produced for each group into records for publishing
    ch_split_panel = ch_split
        .map { r -> record(grp: r.grp, plaf: r.plaf, sample_vcfs: r.sample_vcfs) }
        .join(ch_panel, by: 'grp')

    ch_deploid = ch_grouped_phased_vcf
        .map { r -> record(grp: r.grp, deploid_vcfs: r.deploid_vcfs, deploid_csis: r.deploid_csis) }

    ch_merged = ch_grouped_level1.map { r -> record(grp: r.grp, collate_vcfs: r.vcfs, collate_csis: r.csis) }
        .join(ch_merge_level2.map { r -> record(grp: r.grp, merged_vcf: r.vcf) }, by: 'grp')
        .join(ch_imputed.map { r -> record(grp: r.grp, imputed_vcf: r.vcf) }, by: 'grp')

    emit:
    phased: Channel<GroupVcf> = ch_phased_imputed_polyclonal_monoclonal_vcf
    split: Channel<PolyclonalSplit> = ch_split_panel
    deploid: Channel<PolyclonalDeploid> = ch_deploid
    merged: Channel<PolyclonalMerged> = ch_merged
    mono_only: Channel<MonoclonalOnly> = ch_mono_only
}

// --------------------------- PUBLISHING ---------------------------------------------
// Where each output file goes, relative to the output directory. These are
// shared by the output blocks of this script and main.nf.

def polyclonal_split_publish_paths(r: PolyclonalSplit) -> Map {
    def d = '04_phase_impute_polyclonal'
    return [
        (r.panel):       "${d}/01_deploid_panel_from_clonal_samples/",
        (r.plaf):        "${d}/02_single_sample_vcf_and_group_plaf/",
        (r.sample_vcfs): "${d}/02_single_sample_vcf_and_group_plaf/",
    ]
}

def polyclonal_deploid_publish_paths(r: PolyclonalDeploid) -> Map {
    def d = '04_phase_impute_polyclonal'
    return [
        (r.deploid_vcfs): "${d}/03_deploid_polyclonal/",
        (r.deploid_csis): "${d}/03_deploid_polyclonal/",
    ]
}

def polyclonal_merged_publish_paths(r: PolyclonalMerged) -> Map {
    def d = '04_phase_impute_polyclonal'
    return [
        (r.collate_vcfs): "${d}/04_merge_collate/",
        (r.collate_csis): "${d}/04_merge_collate/",
        (r.merged_vcf):   "${d}/04_merge_filter/",
        (r.imputed_vcf):  "${d}/05_phased_imputed_polyclonal_monoclonal/",
    ]
}

// The monoclonal VCF is published next to the imputed VCFs of the other groups
def monoclonal_only_publish_paths(r: MonoclonalOnly) -> Map {
    def d = '04_phase_impute_polyclonal'
    return [
        (r.vcf): "${d}/05_phased_imputed_polyclonal_monoclonal/${r.vcf.simpleName}_exported.vcf.gz",
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

    ch_mono = channel.fromPath("results/03_phase_impute_monoclonal/04_beagle_imputed/*.vcf.gz")
        .map { p -> record(grp: p.simpleName.replaceAll(/_imputed/, ''), imputed_vcf: p) }

    def opts = record(
        deploid_panel_size: params.deploid_panel_size,
        deploid_panel_options: params.deploid_panel_options,
        beagle_jar: params.beagle_jar,
        gmap: params.gmap,
        beagle_options: params.beagle_options_polyclonal,
    )
    poly = WF_PHASE_IMPUTE_POLYCLONAL(ch_pops, ch_mono, opts)

    publish:
    poly_split = poly.split
    poly_deploid = poly.deploid
    poly_merged = poly.merged
    poly_mono_only = poly.mono_only
}

output {
    poly_split: Channel<PolyclonalSplit> {
        path { r -> polyclonal_split_publish_paths(r) }
    }
    poly_deploid: Channel<PolyclonalDeploid> {
        path { r -> polyclonal_deploid_publish_paths(r) }
    }
    poly_merged: Channel<PolyclonalMerged> {
        path { r -> polyclonal_merged_publish_paths(r) }
    }
    poly_mono_only: Channel<MonoclonalOnly> {
        path { r -> monoclonal_only_publish_paths(r) }
        // the monoclonal VCF is an input of this script when run on its own
        includeInputs true
    }
}

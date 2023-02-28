nextflow.enable.dsl=2

params.beagle_options = ''
params.beagle_jar = "${projectDir}/lib/beagle.08Feb22.fa4.jar"
params.gmap = "/autofs/projects-t3/toconnor_grp/bing.guo/ref/recombination_map/Pf3D7_all.txt"
params.deploid_nopanel_options = "-rate 8 -burn 0.67 -k 4"
def publish_dir = "${params.outdir ?: launchDir}/03_phase_impute_monoclonal"

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
    tag "$grp"
    publishDir "${publish_dir}/01_single_sample_vcf_and_group_plaf", mode: 'symlink'
    input: tuple val(grp), path(fws_file), path(grp_vcf) 
    output: tuple val(grp), path("*_plaf.tsv"), path("single_sample_vcfs/*.vcf.gz")
    script: 
        """
        calc_plaf_from_vcf.py --vcf ${grp_vcf} --out ${grp}_plaf.tsv

        awk 'NR>1 && \$2>=0.95 {print \$1}' ${fws_file} > monoclonal_samples.txt
        split -l 400 monoclonal_samples.txt split_
        bcftools index ${grp_vcf}
        mkdir single_sample_vcfs
        for SPLIT in split_*
        do
            bcftools +split ${grp_vcf} -S \${SPLIT} -Oz -o single_sample_vcfs
        done
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
    publishDir "${publish_dir}/02_deploid_monoclonal", mode: 'symlink'
    input: tuple val(grp), val(grp_size), val(sample), path(plaf_fn), path(single_sample_vcf)
    output: tuple val(grp), val(grp_size), val(sample), path("*_out.vcf.gz"), path("*_out.vcf.gz.csi")
    script:
        def prefix = "${grp}:${sample}"
        """
        set -uexo pipefail
        # remove missing sites if nopanel; do want gt infor at missing site to be used for making panel 
        bcftools view -i "F_MISSING==0" ${single_sample_vcf} -Oz -o tmp.vcf.gz

        # find regions to exclude
        deploid_filter_high_leverage_snps.py --vcf tmp.vcf.gz --sites pos_to_exclude.txt

        # make sure plaf have the same number of sites as in vcf
        deploid_subset_plaf.py --in_plaf ${plaf_fn} --in_vcf tmp.vcf.gz --out_plaf plaf.subset.tmp.txt

        dEploid -vcf tmp.vcf.gz -plaf plaf.subset.tmp.txt -noPanel -vcfOut -o ${prefix} \
            ${params.deploid_nopanel_options} -exclude pos_to_exclude.txt
        rm -f tmp.vcf.gz 

        # generate a .tsv file, *valid.txt file, which.txt, new_sample_name.txt file
        deploid_get_hap.py --sample ${sample} --hap ${prefix}.hap --log ${prefix}.log

        VALID="`cat *~valid.txt`"

        sed '/^##DEploid/d' ${prefix}.vcf > tmp.vcf

        # select haplotype and rename sample name
        bcftools view tmp.vcf -S which.txt | bcftools reheader -s new_sample_name.txt | bcftools view -Oz -o tmp2.vcf.gz

        # convert haploid to homozygous diploid
        zcat tmp2.vcf.gz \
            | awk -v OFS='\t' -v FS='\t' '/^#/{print \$0} /^[^#]/{\$10=\$10 "|" \$10; print \$0}' \
            | bgzip -c > tmp3.vcf.gz

        # convert to bcftools format
        bcftools view tmp3.vcf.gz -Oz -o ${prefix}_out.vcf.gz
        bcftools index ${prefix}_out.vcf.gz 

        # clean up
        rm tmp.vcf tmp2.vcf.gz tmp3.vcf.gz
    """
}

// The following two processes are used to merge many single-sample phased VCF
// files. The reason for splitting the merge into two processes is that merging
// them all together in a single process can cause a failure in bcftools due to
// opening too many files at once. These two processes allow the merging to be
// done at two levels, so that each level only opens a relatively small number
// of single-sample VCF files.

// Level 1 merge

process BCFTOOLS_MERGE_SUBSUBPOP {
    tag "${grp}_${subgrp}"
    publishDir "${publish_dir}/03_merge_collate", mode: 'symlink'
    input: tuple val(grp), val(subgrp), val(n), path(vcfs), path(vcf_idxs), stdin 
    output: tuple val(grp), val(subgrp), val(n), path("${grp}_${subgrp}_merged.vcf.gz"), path("${grp}_${subgrp}_merged.vcf.gz.csi")
    script: 
        def is_single_path = vcfs instanceof nextflow.processor.TaskPath  // in case there is only some sample per group
        if(is_single_path)
            """ 
            cat - > file_list.txt
            mv ${vcfs} ${grp}_${subgrp}_merged.vcf.gz
            bcftools index  ${grp}_${subgrp}_merged.vcf.gz
            """
        else
            """
            cat - > file_list.txt
            bcftools merge -l file_list.txt -Oz -o ${grp}_${subgrp}_merged.vcf.gz
            bcftools index ${grp}_${subgrp}_merged.vcf.gz
            """
    // NOTE: stdin for make file contain list of sample in the same order as vcfs
    // avoid explicitly paste all sample names on command line
}


// Level 2 merge
// See comments in the above process
process BCFTOOLS_MERGE_FILETER_SUBPOP {
    tag "${grp}"
    publishDir "${publish_dir}/03_merge_filter", mode: 'symlink'
    input: tuple val(grp), path(vcfs), path(vcf_idxs), stdin 
    output: tuple val(grp), path("${grp}_merged.vcf.gz")
    script: 
        def is_single_path = vcfs instanceof nextflow.processor.TaskPath  // in case there is only some sample per group
        if(is_single_path)
            """ cat - > file_list.txt; bcftools view -i "F_MISSING<0.3" ${vcfs} -Oz -o ${grp}_merged.vcf.gz """
        else
            """ cat - > file_list.txt; bcftools merge -l file_list.txt | bcftools view -i "F_MISSING<0.3" -Oz -o ${grp}_merged.vcf.gz """
}


// This process calls the Beagle software to impute the missing genotypes of
// samples (without using a panel) within the VCF files that contain all phased
// single-sample VCF files.
process BEAGLE_IMPUTE {
    tag "${grp}"
    publishDir "${publish_dir}/04_beagle_imputed", mode: 'symlink'
    input: tuple val(grp), path(merged_vcf)
    output: tuple val(grp), path("*_imputed.vcf.gz")
    script: 
    """ java -Xmx${task.memory.giga}g -jar ${params.beagle_jar} gt=${merged_vcf} map=${params.gmap} out=${grp}_imputed nthreads=${task.cpus} \\
        ${params.beagle_options} 
    """
}


workflow WF_PHASE_IMPUTE_MONOCLONAL {
    take:
        ch_fws          //Format: tuple val(grp), path(grp_fws_txt)
        ch_subpop_vcf   //Format: tuple val(grp), path("${grp}.vcf.gz") 
    main:
        ch_fws_vcf = ch_fws.combine(ch_subpop_vcf, by: 0)
        BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF(ch_fws_vcf)

        ch_svcfs_plaf = BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF.out
        ch_vcf_plaf = ch_svcfs_plaf
            .flatMap{grp, plaf, vcfs -> 
                def grp_size = 0
                if(vcfs instanceof java.util.ArrayList){          // conditionally flatten out single vcf files
                    grp_size = vcfs.size()
                    return vcfs.collect{vcf-> [grp, grp_size, vcf.getSimpleName(), plaf, vcf]}
                }
                else {
                    grp_size = 1
                    return [[grp, grp_size, vcfs.getSimpleName(), plaf, vcfs]]
                }
            }
        DEPLOID_PHASE_MONOCLONAL(ch_vcf_plaf)
        ch_phased_vcf = DEPLOID_PHASE_MONOCLONAL.out

        // NOTE did not use groupTuple size argument as the groupResult seems to be randome not good for resume
        ch_grouped_phased_vcf = ch_phased_vcf
            .map {grp, grp_size, sample, vcf, vcf_idx -> [groupKey(grp, grp_size), [sample, vcf, vcf_idx]]}
            .groupTuple(by:0)
            .map {grp, ll-> 
                ll.sort{a, b->a[0]<=>b[0]};                     // sort
                def sorted_vcfs = ll.collect{it -> it[1]}       // extract
                def sorted_vcf_idxs =ll.collect{it -> it[2]}    // extract
                return [grp, sorted_vcfs, sorted_vcf_idxs]
                }

        ch_collated_vcf = ch_grouped_phased_vcf.flatMap{grp, vcfs, idx->
            def a = vcfs.collate(200)
            def b = idx.collate(200)
            def n = a.size()
            def subgrps = (0..(n-1))
            return [subgrps, a, b].transpose().collect{subgrp, vcf_lst, idx_lst->
                def stdin_vcf_str = vcf_lst.join('\n') + '\n'
                return [grp, subgrp, n, vcf_lst, idx_lst, stdin_vcf_str]
            }
        }
            // .view{[it[0], it[1].collect{x->x.getSimpleName()}, it[2].collect{x->x.getSimpleName()}]}

        // Merge level 1
        BCFTOOLS_MERGE_SUBSUBPOP(ch_collated_vcf)

        ch_merge_level1 = BCFTOOLS_MERGE_SUBSUBPOP.out.map {grp, subgrp, n, m_vcf, m_idx ->
                [ groupKey(grp, n), [subgrp, m_vcf, m_idx]]}
            .groupTuple()
            .map {grp, ll-> 
                ll.sort{a, b->a[0]<=>b[0]};                     // sort
                def sorted_vcfs = ll.collect{it -> it[1]}       // extract
                def sorted_vcf_idxs =ll.collect{it -> it[2]}    // extract
                def stdin_vcf_list_str = sorted_vcfs.collect{p->p.getName()}.join('\n') + '\n'
                return [grp, sorted_vcfs, sorted_vcf_idxs, stdin_vcf_list_str]
                }
        // Merge level 2
        BCFTOOLS_MERGE_FILETER_SUBPOP(ch_merge_level1)

        ch_merge_level2 = BCFTOOLS_MERGE_FILETER_SUBPOP.out
        BEAGLE_IMPUTE(ch_merge_level2)

        ch_phased_imputed_vcf = BEAGLE_IMPUTE.out

    emit:
        ch_phased_imputed_vcf // Format: tuple val(grp), path("*_imputed.vcf.gz")

}

// --------------------- FOR TESTING -------------------------------
// will be ignored when included as subworkflow script
workflow {
    ch_fws = channel.fromPath("results/02_pop_fws/04_moimix_fws/*.fws.txt")
        .map {fws_p-> def grp=fws_p.getSimpleName(); return [grp, fws_p]}
    ch_subpop_vcf = channel.fromPath("results/02_pop_fws/02_subpop_vcf/*.vcf.gz")
        .map {subpop_vcf_p -> def grp=subpop_vcf_p.getSimpleName(); return [grp, subpop_vcf_p]}

    WF_PHASE_IMPUTE_MONOCLONAL(ch_fws, ch_subpop_vcf)
    
}

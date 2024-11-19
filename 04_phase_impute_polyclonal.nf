nextflow.enable.dsl=2

params.hapibd_options = ''
params.beagle_options = 'overlap=3.0' // default is 2.0, cause error for some vcf files
params.deploid_panel_size = 16 
params.beagle_jar = "${projectDir}/lib/beagle.08Feb22.fa4.jar"
params.gmap = "/autofs/projects-t3/toconnor_grp/bing.guo/ref/recombination_map/Pf3D7_all.txt"
params.deploid_panel_options = "-nSample 250 -rate 8 -burn 0.67 -k 4 -ibd"
def publish_dir = "${params.outdir ?: launchDir}/04_phase_impute_polyclonal"

// This process is used to identify the most diverse samples from the
// phased & imputed monoclonal samples. These samples are then extracted and
// converted into a table format that is compatible with the dEploid-IBD
// program.
//
// Note that the `grp_vcf` file contains both monoclonal and polyclonal samples,
// and it also includes Format/AD fields which are used for finding
// representatives. The intermediate file, `grp_monomputed_vcf`, contains only
// monoclonal samples, with the sample names having a ~xxx suffix, and with GT
// imputed/phased and the AD field removed.
process ALLEL_GENERATE_PANEL {
    tag "$grp"
    publishDir "${publish_dir}/01_deploid_panel_from_clonal_samples", mode: 'symlink'
    input: tuple val(grp),  path(fws_fn), path(grp_vcf), path(grp_monomputed_vcf)
    output: tuple val(grp), path("*_panel.tsv")
    script: 
        """
        awk 'NR>1 && \$2>=0.95 {print \$1}' ${fws_fn} > clonal_samples.txt
        deploid_find_representives.py --vcf ${grp_vcf} --clonal_samples clonal_samples.txt --n_samples ${params.deploid_panel_size} \\
            --out1 ${grp}_panel_samples_1.txt --out2 ${grp}_panel_samples_2.txt
        bcftools index ${grp_vcf}
        bcftools view -S ${grp}_panel_samples_2.txt ${grp_vcf} -Oz -o panel2.vcf.gz
        deploid_panel_vcf2table.py --vcf ${grp_monomputed_vcf} --samples ${grp}_panel_samples_2.txt --panel ${grp}_panel.tsv
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
process BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF {
    tag "$grp"
    publishDir "${publish_dir}/02_single_sample_vcf_and_group_plaf", mode: 'symlink'
    input: tuple val(grp), path(fws_file), path(grp_vcf) 
    output: 
        tuple val(grp), path("*_plaf.tsv"), path("single_sample_vcfs/*.vcf.gz"), emit: single_sample_vcf optional true
        tuple val(grp), path(grp_vcf), path('GROUP_HAS_NO_POLYCLONAL_SAMPLE'), emit: mono_only_grp_vcf optional true
    script: 
        """
        calc_plaf_from_vcf.py --vcf ${grp_vcf} --out ${grp}_plaf.tsv
        awk 'NR>1 && \$2<0.95 {print \$1}' ${fws_file} > polyclonal_samples.txt
        mkdir single_sample_vcfs
        bcftools index ${grp_vcf} 
        
        # NOTE: some group might not have polyclonal samples
        if [ `cat  polyclonal_samples.txt | wc -l` == 0 ]; then 
            touch GROUP_HAS_NO_POLYCLONAL_SAMPLE
        else 
            split -l 400 polyclonal_samples.txt split_
            for SPLIT in split_*
            do
                bcftools +split ${grp_vcf} -S \${SPLIT} -Oz -o single_sample_vcfs
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
    publishDir "${publish_dir}/03_deploid_polyclonal", mode: 'symlink'
    input: tuple val(grp), val(grp_size), val(sample), path(plaf_fn), path(single_sample_vcf), path(panel_fn)
    output: tuple val(grp), val(grp_size), val(sample), path("*_out.vcf.gz"), path("*_out.vcf.gz.csi"), env(VALID)
    script:
        def prefix = "${grp}:${sample}"
        """
        set -uexo pipefail
        # remove sites not in panel; otherwise dEploid will complain
        cat $panel_fn | sed 1d | cut -f 1-2 > tmp_sites_to_keep.txt
        bcftools index ${single_sample_vcf}
        bcftools view -R tmp_sites_to_keep.txt  ${single_sample_vcf} -Oz -o tmp.vcf.gz
        # delete the index to avoid output matching two csi files
        rm -f ${single_sample_vcf}.csi

        # find regions to exclude
        deploid_filter_high_leverage_snps.py --vcf tmp.vcf.gz --sites pos_to_exclude.txt

        # make sure plaf have the same number of sites as in vcf
        deploid_subset_plaf.py --in_plaf ${plaf_fn} --in_vcf tmp.vcf.gz --out_plaf plaf.subset.tmp.txt

        dEploid -vcf tmp.vcf.gz -plaf plaf.subset.tmp.txt -panel ${panel_fn} -vcfOut -o ${prefix} \
            ${params.deploid_panel_options} -exclude pos_to_exclude.txt
        rm -f tmp.vcf.gz 

        # generate a .tsv file, *valid.txt file, which.txt, new_sample_name.txt file
        deploid_get_hap.py --sample ${sample} --hap ${prefix}.hap --log ${prefix}.log

        # valid criteria in get_hap.py: if second/first < 1/3 and first > 0.7:
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
// Note the level 1 merge also include the monoclonal samples of the sample group.
// After this process, both the monoclonal samples and the polyclonal samples with
// dominant haplotype are merged. 

// Level 1 merge

process BCFTOOLS_MERGE_SUBSUBPOP {
    tag "${grp}_${subgrp}"
    publishDir "${publish_dir}/04_merge_collate", mode: 'symlink'
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
    stub:
    """
    touch  ${grp}_${subgrp}_merged.vcf{.gz,.gz.csi}
    """
}

// Level 2 merge + monoclonal samples
// See comments in the above process
process BCFTOOLS_MERGE_WITH_MONOCLONALS_FILETER_SUBPOP {
    tag "${grp}"
    publishDir "${publish_dir}/04_merge_filter", mode: 'symlink'
    input: tuple val(grp), path(vcfs), path(vcf_idxs), stdin 
    output: tuple val(grp), path("${grp}_merged.vcf.gz")
    script: 
        def is_single_path = vcfs instanceof nextflow.processor.TaskPath  // in case there is only one sample per group
        if(is_single_path)
            """ cat - > file_list.txt; bcftools view -i "F_MISSING<0.3" ${vcfs} -Oz -o ${grp}_merged.vcf.gz """
        else
            """ cat - > file_list.txt
            bcftools index -f ${vcfs.last()}
             bcftools merge -l file_list.txt | bcftools view -i "F_MISSING<0.3" -Oz -o ${grp}_merged.vcf.gz """
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
    publishDir "${publish_dir}/05_phased_imputed_polyclonal_monoclonal", mode: 'symlink'
    input: tuple val(grp), path(merged_vcf)
    output: tuple val(grp), path("*_imputed.vcf.gz")
    script: 
    """ java -Xmx${task.memory.giga}g -jar ${params.beagle_jar} gt=${merged_vcf} map=${params.gmap} out=${grp}_imputed nthreads=${task.cpus} \\
        ${params.beagle_options} 
    """
    stub:
    """
    touch  ${grp}_imputed.vcf.gz
    """
}

// An empty process is used to organize result files.
// NOTE: mono_only_vcf was not already linked to result directory. 
// This process will do that.
process PUBLISH_MONO_ONLY_VCF {
    tag "${grp}"
    publishDir "${publish_dir}/05_phased_imputed_polyclonal_monoclonal", mode: 'symlink'
    input: tuple val(grp), path(mono_only_vcf)
    output: tuple val(grp), path("*_exported.vcf.gz")
    script: 
    """
    cp ${mono_only_vcf} ${mono_only_vcf.getSimpleName() + "_exported.vcf.gz"}
    """
}

workflow WF_PHASE_IMPUTE_POLYCLONAL {
    take:
        ch_fws          //Format: tuple val(grp), path(grp_fws_txt)
        ch_subpop_vcf   //Format: tuple val(grp), path("${grp}.vcf.gz") 
        ch_subpop_clonal_phased_imputed_vcf // Format: tuple val(grp), path("*_imputed.vcf.gz")
    main:
        ch_fws_vcf = ch_fws.combine(ch_subpop_vcf, by: 0)
        ch_fws_vcf_monoimputedvcf = ch_fws_vcf.combine(ch_subpop_clonal_phased_imputed_vcf, by: 0)
        

        ALLEL_GENERATE_PANEL(ch_fws_vcf_monoimputedvcf)
        ch_subpop_panel = ALLEL_GENERATE_PANEL.out

        BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF(ch_fws_vcf)


        ch_svcfs_plaf = BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF.out.single_sample_vcf
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
        // join vcf_plaf_with_panel
        ch_vcf_plaf_panel = ch_vcf_plaf.combine(ch_subpop_panel, by: 0)
            //.view{grp, sz, sample, plaf, vcf, panel -> [grp, sz, sample, plaf.getSimpleName(), vcf.getSimpleName(), panel.getSimpleName()]}

        DEPLOID_PHASE_POLYCLONAL(ch_vcf_plaf_panel)
        ch_phased_vcf = DEPLOID_PHASE_POLYCLONAL.out

        // This is used to circument the 'DEPLOID_PHASE_POLYCLONAL' process
        //// ch_phased_vcf = channel.fromPath('/autofs/chib/toconnor_grp/bing/20220314_analyze_joint_call_pf/runs/run_real/bk_deploid_vcf_map')
        ////     .splitCsv(sep: '\t')
        ////     .map{grp, grp_size, sample, vcf, vcf_idx, valid->
        ////         [grp, grp_size.toInteger(), sample, file(vcf), file(vcf_idx), valid]
        ////     }

            
        // grouping, sorting and filtering (keep only samples with a dominant clone)
        ch_grouped_phased_vcf = ch_phased_vcf
            .map {grp, grp_size, sample, vcf, vcf_idx, valid -> [groupKey(grp, grp_size), [sample, vcf, vcf_idx, valid]]}
            .groupTuple(by:0)
            .map {grp, ll-> 
                ll.sort{a, b->a[0]<=>b[0]};                     // sort
                ll_valid = ll.findAll {sample, vcf, vcf_idx, valid -> valid=='true'}
                def sorted_vcfs = ll_valid.collect{it -> it[1]}       // extract
                def sorted_vcf_idxs =ll_valid.collect{it -> it[2]}    // extract
                return [grp, sorted_vcfs, sorted_vcf_idxs]
                }
            .branch {grp, vcfs, idx ->
                predominant_clones: vcfs.size() > 0
                grp_without_valid_clones: vcfs.size() == 0
            }
             // .view{[it[0], it[1].collect{x->x.getSimpleName()}, it[2].collect{x->x.getSimpleName()}]}

        // make a channel to save grp vcf of two cases
        // 1. when the group has no polyclonal samples
        // grps that did not make it to deconvolution for polyclonal
        ch_grps1 = BCFTOOLS_SPLIT_PREP_SINGLE_SAMPLE_VCF_AND_CALC_PLAF.out.mono_only_grp_vcf
            .map{grp, vcf, placeholder -> grp} // get grp name
        // 2. when the grp has polyclonal samples but do not have predominant clones
        // grps that make it to the deconvolution for polyclonal but
        // there is no polyclonal samples having dominant clones
        ch_grps2 = ch_grouped_phased_vcf.grp_without_valid_clones
            .map{grp, vcf, idx -> grp}

        ch_mono_only_vcf = ch_grps1.concat(ch_grps2)
            .combine(ch_subpop_clonal_phased_imputed_vcf, by: 0)// tuple val(grp), path(vcf)

        ch_collated_vcf = ch_grouped_phased_vcf.predominant_clones
            .flatMap{grp, vcfs, idx->
                def a = vcfs.collate(200)
                def b = idx.collate(200)
                def n = a.size()
                def subgrps = (0..(n-1))
                return [subgrps, a, b].transpose().collect{subgrp, vcf_lst, idx_lst ->
                    def stdin_vcf_str = vcf_lst.collect{p->p.getName()}.join('\n') + '\n'
                    return [grp, subgrp, n, vcf_lst, idx_lst, stdin_vcf_str]
                }
            }

        // Merge level 1
        BCFTOOLS_MERGE_SUBSUBPOP(ch_collated_vcf)

        ch_merge_level1 = BCFTOOLS_MERGE_SUBSUBPOP.out.map {grp, subgrp, n, m_vcf, m_idx ->
                [ groupKey(grp, n), [subgrp, m_vcf, m_idx]]
            }
            .groupTuple()
            .map {grp, ll-> 
                ll.sort{a, b->a[0]<=>b[0]};                     // sort
                def sorted_vcfs = ll.collect{it -> it[1]}       // extract
                def sorted_vcf_idxs =ll.collect{it -> it[2]}    // extract
                def stdin_vcf_list_str = sorted_vcfs.join('\n') + '\n'
                return [grp, sorted_vcfs, sorted_vcf_idxs, stdin_vcf_list_str]
                }

        ch_merge_level1_with_monoclonals = ch_merge_level1
            .combine(ch_subpop_clonal_phased_imputed_vcf, by:0)
            .map{grp, poly_svcfs, poly_sidx, stdin_vcf_list_str, monoclonal_vcf->
                [grp, poly_svcfs + [monoclonal_vcf], poly_sidx, stdin_vcf_list_str + "${monoclonal_vcf.getName()}\n"]
            }

        // Merge level 2
        BCFTOOLS_MERGE_WITH_MONOCLONALS_FILETER_SUBPOP(ch_merge_level1_with_monoclonals)
        ch_merge_level2 = BCFTOOLS_MERGE_WITH_MONOCLONALS_FILETER_SUBPOP.out

        BEAGLE_IMPUTE(ch_merge_level2)

        PUBLISH_MONO_ONLY_VCF(ch_mono_only_vcf)

        ch_phased_imputed_polyclonal_monoclonal_vcf = BEAGLE_IMPUTE.out.concat(ch_mono_only_vcf)
        



    emit:
        ch_phased_imputed_polyclonal_monoclonal_vcf // Format: tuple val(grp), path("*_imputed.vcf.gz")
}

// --------------------- FOR TESTING -------------------------------
// will be ignored when included as subworkflow script
workflow {
    ch_fws = channel.fromPath("results/02_pop_fws/04_moimix_fws/*.fws.txt")
        .map {fws_p-> def grp=fws_p.getSimpleName(); return [grp, fws_p]}
    ch_subpop_vcf = channel.fromPath("results/02_pop_fws/02_subpop_vcf/*.vcf.gz")
        .map {subpop_vcf_p -> def grp=subpop_vcf_p.getSimpleName(); return [grp, subpop_vcf_p]}

    ch_subpop_deploidclonalvcf = channel.fromPath("results/03_phase_impute_monoclonal/04_beagle_imputed/*.vcf.gz")
        .map {subpop_deploidclonalvcf_p -> def grp=subpop_deploidclonalvcf_p.getSimpleName().replaceAll(/_imputed/, ''); return [grp, subpop_deploidclonalvcf_p]}

    WF_PHASE_IMPUTE_POLYCLONAL(ch_fws, ch_subpop_vcf, ch_subpop_deploidclonalvcf)
    
}

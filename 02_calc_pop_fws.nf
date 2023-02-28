nextflow.enable.dsl=2

def publish_dir = "${params.outdir ?: launchDir}/02_pop_fws"

// Comparator for sorting genome intervals
def cmp_interval(String i, String j){
    def split_i = i.split(':')
    def chr_i = split_i[0] // chromosome name of interval i
    def start_i = split_i.size() < 2 ? 0: split_i[1].split('-')[0].toInteger() // start of interval i
    def split_j = j.split(':') // chromosome name of inteval j
    def chr_j = split_j[0] // start of interval j
    def start_j = split_j.size() < 2 ? 0: split_j[1].split('-')[0].toInteger()
    return chr_i<=>chr_j ?: start_i <=> start_j
}

// This process serves two purposes:
// - Concatenate interval-VCF files across the genome.
// - Simplify annotations fields so file sizes are minimal
process BCFTOOLS_CONCAT_MERGE_REGIONS_SIMPLIFY_ANNOTATION {
    tag "merge_regions_samplify_annot"
    publishDir "${publish_dir}/01_merge_regions_simplify_annot", mode: 'symlink'
    input: tuple val(interval_lst), path(vcf_lst)
    output: tuple path("merged.bcf"), path("merged.bcf.*")
    script: """ bcftools concat ${vcf_lst} --threads ${task.cpus} | bcftools annotate -x INFO,^FORMAT/GT,^FORMAT/AD -Ob -o merged.bcf;
                bcftools index merged.bcf """
    // NOTE: moimix needs format/AD to be in the vcf files
}

// This process involves dividing the full dataset into population subsets or
// groups (grps) and preparing VCF files for Fws metric calculation.
//
// NOTE: use `--force-samples` option to ignore errors when some requested
// samples are not in the vcf file
process BCFTOOLS_VIEW_EXTRACT_SUBPOP_SAMPLES {
    tag "$grp"
    publishDir "${publish_dir}/02_subpop_vcf", mode: 'symlink'
    input: tuple val(grp), path(sample_list_file); tuple path(merged_vcf), path(merged_vcf_idx)
    output: tuple val(grp), path("${grp}.vcf.gz")
    script: """bcftools view -S ${sample_list_file} ${merged_vcf} -Oz -o ${grp}.vcf.gz --force-samples"""
}

// This process converts population VCF files into the GDS file format as
// required by the `moimix` package
//
// NOTE: seqVcf2Gds is very slow but can be parallelized by specifying the `parallel`
// argument.
process SEQARRAY_VCF_TO_GDS {
    tag "$grp"
    publishDir "${publish_dir}/03_seqarray_gds", mode: 'symlink'
    input: tuple val(grp), path(vcf)
    output: tuple val(grp), path("*.gds")
    script: """#! /usr/bin/env Rscript
        library(SeqArray); seqVCF2GDS("${vcf}", "${grp}.gds", parallel=${task.cpus})
        """
}

// Run moimix to generate Fws table for each samples of a population
//
// NOTE: Both seqVcf2Gds and GetFws need a lot memory shoud use a dyanmic meory disignation
process R_MOIMIX_FWS {
    tag "$grp"
    publishDir "${publish_dir}/04_moimix_fws", mode: 'symlink'
    input: tuple val(grp), path(gds)
    output: tuple val(grp), path("*.fws.txt")
    script: """#! /usr/bin/env Rscript
        library(SeqArray)
        library("moimix")
        isolates <- seqOpen("${gds}")
        sample.id <- seqGetData(isolates, "sample.id")
        fws_all <- getFws(isolates)
        df = data.frame(sample.id, fws_all)
        write.table(df, "${grp}.fws.txt", sep="\t", row.names=FALSE, quote=FALSE)
        """
}

workflow WF_CALC_POP_FWS {
    take: 
        ch_input_vcf // Format: tuple val(interval), path(vcf_gz)
        ch_input_subpop // Format: tumple val(group_name), path(sample_list_file)
    main:

    ch_input_vcf_sorted_list = ch_input_vcf.toSortedList{a, b-> cmp_interval(a[0], b[0])}
        .map {lst -> def interval_lst = lst.collect{it[0]}; def vcf_lst = lst.collect{it[1]}; return [interval_lst, vcf_lst]}

    BCFTOOLS_CONCAT_MERGE_REGIONS_SIMPLIFY_ANNOTATION(ch_input_vcf_sorted_list)

    ch_merged_vcf = BCFTOOLS_CONCAT_MERGE_REGIONS_SIMPLIFY_ANNOTATION.out
    BCFTOOLS_VIEW_EXTRACT_SUBPOP_SAMPLES(ch_input_subpop, ch_merged_vcf)

    ch_subpop_vcf = BCFTOOLS_VIEW_EXTRACT_SUBPOP_SAMPLES.out
    SEQARRAY_VCF_TO_GDS(ch_subpop_vcf)

    ch_gds = SEQARRAY_VCF_TO_GDS.out
    R_MOIMIX_FWS(ch_gds)

    ch_fws = R_MOIMIX_FWS.out

    emit:
        ch_fws //Format tuple val(grp), path(grp_fws_txt)
        ch_subpop_vcf //tuple val(grp), path("${grp}.vcf.gz") 
}

// --------------------- FOR TESTING -------------------------------
// will be ignored when included as subworkflow script
workflow{
    ch_input_vcf = channel.fromPath("results/01_filt_vcf/08_maf_bcf/*.vcf.gz", checkIfExists: true)
        .map{p-> def interval=p.getSimpleName().replaceAll(/_fmiss2_imiss2_maf/, ''); return [interval, p]}

    ch_input_subpop = channel.fromPath(params.input_subpop_fn)
        .splitCsv(sep: '\t')
        .map{
            assert it.size() == 2
            def grp = it[0]
            def sample_list_file = file(it[1], checkIfExists: true, type: 'file')
            return [grp, sample_list_file]
        }

    WF_CALC_POP_FWS(ch_input_vcf, ch_input_subpop)

}
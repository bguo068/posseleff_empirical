nextflow.enable.dsl=2

params.high_vmiss = 0.7
params.high_imiss = 0.7
params.low_vmiss = 0.3
params.low_imiss = 0.3
params.min_maf = 0.01
params.samples_included = '' // TODO: set a whitelist sample list as default

def publish_dir = "${params.outdir ?: launchDir}/01_filt_vcf"

// Filter VCF sites by the variant type and the number of alleles
// 
// NOTE: 
// The params.samples_included variable allows a  manual removal of samples that 
// - are duplicated  
// - or show low sequencing quality 
process BCFTOOLS_VIEW_PASS_BIALLELIC {
    tag "$interval"
    publishDir "${publish_dir}/01_pass_biallelic", mode: 'symlink'
    input: tuple val(interval), path(vcf) 
    output: tuple val(interval), path("*.bcf")
    script: 
    def S_args = params.samples_included ? "-S ${file(params.samples_included)}" : ""
    """ 
    if [ -f ${vcf}.csi ] ; then bcftools index ${vcf};  fi
    bcftools view ${S_args} ${vcf} -Ou | bcftools view -f PASS -Ou | bcftools view -v snps -Ou | bcftools view -m2 -M2 -Ob -o ${interval}.bcf
    """
}

// The process filters VCF sites by the site-missingness (what percentage of samples have 
// missing genotype for a given site) and then calculates per-sample stats. 
// Note that, as the input VCFs are only part of the genome (a range), only after full genome 
// information (per-sample-stats for all genome-interval VCFs) are gathered can
// we filter the VCF samples by missingness.
process BCFTOOLS_VIEW_F_MISSING1_WITH_PSC_STATS {
    tag "$interval"
    publishDir "${publish_dir}/02_fmissing1", mode: 'symlink'
    input: tuple val(interval), path(vcf)
    output: tuple val(interval), path("*_fmiss1.bcf"), path("*_fmiss1_stats.txt")
    script: """ bcftools view -i "F_MISSING<${params.high_vmiss}"  ${vcf} -Ob -o ${interval}_fmiss1.bcf
                bcftools stats -s- ${interval}_fmiss1.bcf > ${interval}_fmiss1_stats.txt
            """
} 

// Base on per-sample stats for all interval-VCF files across the whole genome, 
// the process generate a list of samples that have the percentage of genotype
// missiningness less than a given threshold.
process PANDAS_IMISS1_SAMPLES_TO_KEEP {
    tag "sample_to_keep_imiss1"
    publishDir "${publish_dir}/03_imiss1_samples", mode: 'symlink'
    input: path(vcf_stats)
    output: path("samples_to_keep_imiss1.txt")
    script: """ get_samples_with_max_imiss.py --bcftools_psc_stats_files *stats.txt --max_imiss ${params.high_imiss} --out samples_to_keep_imiss1.txt """
}

// This process does the actual sample filtering (based on missingness).
// After samples are filtering, a step added to ensure 
// the remaining sites are segregating.
process BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS1 {
    tag "$interval"
    publishDir "${publish_dir}/04_imiss1_bcf", mode: 'symlink'
    input: tuple val(interval), path(vmiss1_bcf); path(samples_to_keep_imiss1)
    output: tuple val(interval), path('*_fmiss1_imiss1.bcf')
    script: """ bcftools index ${vmiss1_bcf}; bcftools view -S ${samples_to_keep_imiss1} ${vmiss1_bcf} | bcftools view -i "MAC>0" -Ob -o ${interval}_fmiss1_imiss1.bcf"""
}

// This process is similar to `BCFTOOLS_VIEW_F_MISSING1_WITH_PSC_STATS` but uses
// a stricter threshold. 
process BCFTOOLS_VIEW_F_MISSING2_WITH_PSC_STATS {
    tag "$interval"
    publishDir "${publish_dir}/05_fmissing2", mode: 'symlink'
    input: tuple val(interval), path(vcf)
    output: tuple val(interval), path("*_fmiss2.bcf"), path("*_fmiss2_stats.txt")
    script: """ bcftools view -i "F_MISSING<${params.low_vmiss}"  ${vcf} -Ob -o ${interval}_fmiss2.bcf
                bcftools stats -s- ${interval}_fmiss2.bcf > ${interval}_fmiss2_stats.txt
            """
} 

// This process is similar to `PANDAS_IMISS1_SAMPLES_TO_KEEP` but uses
// a stricter threshold. 
process PANDAS_IMISS2_SAMPLES_TO_KEEP {
    tag "sample_to_keep_imiss2"
    publishDir "${publish_dir}/06_imiss2_samples", mode: 'symlink'
    input: path(vcf_stats)
    output: path("samples_to_keep_imiss2.txt")
    script: """ get_samples_with_max_imiss.py --bcftools_psc_stats_files *stats.txt --max_imiss ${params.low_imiss} --out samples_to_keep_imiss2.txt """
}

// This process is similar to `BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS1` but uses
// a stricter threshold. 
process BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS2 {
    tag "$interval"
    publishDir "${publish_dir}/07_imiss2_bcf", mode: 'symlink'
    input: tuple val(interval), path(vmiss2_bcf); path(samples_to_keep_imiss2)
    output: tuple val(interval), path('*_fmiss2_imiss2.bcf')
    script: """ bcftools index ${vmiss2_bcf}; bcftools view -S ${samples_to_keep_imiss2} ${vmiss2_bcf} | bcftools view -i "MAC>0" -Ob -o ${interval}_fmiss2_imiss2.bcf"""
}

// This process filter out rare variants
process BCFTOOLS_VIEW_FILTER_SITE_BY_MAF {
    tag "$interval"
    publishDir "${publish_dir}/08_maf_bcf", mode: 'symlink'
    input: tuple val(interval), path(imiss2_bcf)
    output: tuple val(interval), path("*_fmiss2_imiss2_maf.vcf.gz")
    script: """ bcftools view -q ${params.min_maf}:minor ${imiss2_bcf} -Oz -o ${interval}_fmiss2_imiss2_maf.vcf.gz"""
}

// --------------------------- SUBWORKFLOW --------------------------------------------
// This subworkflow filters the set of VCF files by variant type, number of
// alleles per site, and missingness per site and per sample.

workflow WF_FILTER_VCF {
    
    take: ch_input_vcf // FORMAT: tuple val(interval), path(vcf)

    main:
        // filter by filter=PASS and -m2 and -M2
        BCFTOOLS_VIEW_PASS_BIALLELIC(ch_input_vcf)

        // Filter by vmiss (round 1)
        ch_pass_biallelic_vcf = BCFTOOLS_VIEW_PASS_BIALLELIC.out
        BCFTOOLS_VIEW_F_MISSING1_WITH_PSC_STATS(ch_pass_biallelic_vcf)

        // Filter by imiss (round 1)
        ch_fmiss1_vcf_and_psc_stats = BCFTOOLS_VIEW_F_MISSING1_WITH_PSC_STATS.out
        ch_fmiss1_vcf = ch_fmiss1_vcf_and_psc_stats.map{interval, fmiss1_bcf, fmiss1_bcf_stat->[interval, fmiss1_bcf]}
        ch_fmiss1_stats_list = ch_fmiss1_vcf_and_psc_stats.map{interval, fmiss1_bcf, fmiss1_bcf_stat -> fmiss1_bcf_stat}.toList()

        PANDAS_IMISS1_SAMPLES_TO_KEEP(ch_fmiss1_stats_list)
        ch_imiss1_samples = PANDAS_IMISS1_SAMPLES_TO_KEEP.out

        BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS1(ch_fmiss1_vcf, ch_imiss1_samples)

        // Filter by vmiss (round 2)
        ch_fmiss1_imiss1_vcf =  BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS1.out

        // Filter by imiss (round 2)
        ch_fmiss2_vcf_and_psc_stats = BCFTOOLS_VIEW_F_MISSING2_WITH_PSC_STATS(ch_fmiss1_imiss1_vcf)
        ch_fmiss2_vcf = ch_fmiss2_vcf_and_psc_stats.map{interval, fmiss2_bcf, fmiss2_bcf_stat->[interval, fmiss2_bcf]}
        ch_fmiss2_stats_list = ch_fmiss2_vcf_and_psc_stats.map{interval, fmiss2_bcf, fmiss2_bcf_stat -> fmiss2_bcf_stat}.toList()

        PANDAS_IMISS2_SAMPLES_TO_KEEP(ch_fmiss2_stats_list)
        ch_imiss2_samples = PANDAS_IMISS2_SAMPLES_TO_KEEP.out

        BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS2(ch_fmiss2_vcf, ch_imiss2_samples)

        // Filter by maf
        ch_imiss2_vcf = BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS2.out
        BCFTOOLS_VIEW_FILTER_SITE_BY_MAF(ch_imiss2_vcf)

        ch_filterred_vcf = BCFTOOLS_VIEW_FILTER_SITE_BY_MAF.out

    emit: 
        ch_filterred_vcf // FORMAT: tuple val(interval), path(vcf)
}


// -------------------------- TESTING -----------------------------------------------
// The workflow is designed to efficiently process large datasets
// by breaking them down into smaller VCF files, each containing a specific
// genome interval or chunk for all samples. The workflow is designed to run in
// parallel, allowing for faster processing times, while also gathering related
// information as needed.

workflow {
    ch_input_vcf = channel
        .fromPath(params.input_vcf_fn, type: 'file', checkIfExists: true )
        .splitCsv(header: false, sep: '\t')
        .map{
            assert it.size() == 2
            def interval = it[0]
            def vcf = file(it[1], type: 'file', checkIfExists: true)
            return [interval, vcf]
        }.take(6)
    
    WF_FILTER_VCF(ch_input_vcf)
}

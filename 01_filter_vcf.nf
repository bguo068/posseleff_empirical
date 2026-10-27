nextflow.enable.types = true

// Parameters used when this script is run on its own (see the entry workflow
// at the bottom). When included by main.nf, main.nf's params block is used.
params {
    input_vcf_fn: Path = "${launchDir}/input_vcfs.tsv"
    // Optional whitelist of samples (see BCFTOOLS_VIEW_PASS_BIALLELIC)
    samples_included: Path?
    high_vmiss: Float = 0.7
    high_imiss: Float = 0.7
    low_vmiss: Float = 0.3
    low_imiss: Float = 0.3
    min_maf: Float = 0.01
}

record IntervalVcf {
    interval: String
    vcf: Path
}

record FilterVcfParams {
    samples_included: Path?
    high_vmiss: Float
    high_imiss: Float
    low_vmiss: Float
    low_imiss: Float
    min_maf: Float
}

// All files produced for one genome interval. `vcf` is the final filtered VCF.
record FilteredVcf {
    interval: String
    pass_bcf: Path
    fmiss1_bcf: Path
    fmiss1_stats: Path
    imiss1_bcf: Path
    fmiss2_bcf: Path
    fmiss2_stats: Path
    imiss2_bcf: Path
    vcf: Path
}

record FilterSampleLists {
    imiss1_samples: Path
    imiss2_samples: Path
}

// Filter VCF sites by the variant type and the number of alleles
//
// NOTE:
// The samples_included input allows a manual removal of samples that
// - are duplicated
// - or show low sequencing quality
process BCFTOOLS_VIEW_PASS_BIALLELIC {
    tag "${interval}"

    input:
    record(interval: String, vcf: Path)
    samples_included: Path?

    output:
    record(interval: interval, bcf: file("${interval}.bcf"))

    script:
    def S_args = samples_included ? "-S ${samples_included}" : ""
    """
    if [ -f ${vcf}.csi ] ; then bcftools index ${vcf};  fi
    bcftools view ${S_args} ${vcf} -Ou | bcftools view -f PASS -Ou | bcftools view -v snps -Ou | bcftools view -m2 -M2 -Ob -o ${interval}.bcf
    """

    stub:
    """
    touch ${interval}.bcf
    """
}

// The process filters VCF sites by the site-missingness (what percentage of samples have
// missing genotype for a given site) and then calculates per-sample stats.
// Note that, as the input VCFs are only part of the genome (a range), only after full genome
// information (per-sample-stats for all genome-interval VCFs) are gathered can
// we filter the VCF samples by missingness.
process BCFTOOLS_VIEW_F_MISSING1_WITH_PSC_STATS {
    tag "${interval}"

    input:
    record(interval: String, bcf: Path)
    max_vmiss: Float

    output:
    record(interval: interval, bcf: file("${interval}_fmiss1.bcf"), stats: file("${interval}_fmiss1_stats.txt"))

    script:
    """
    bcftools view -i "F_MISSING<${max_vmiss}"  ${bcf} -Ob -o ${interval}_fmiss1.bcf
    bcftools stats -s- ${interval}_fmiss1.bcf > ${interval}_fmiss1_stats.txt
    """

    stub:
    """
    touch ${interval}_fmiss1.bcf ${interval}_fmiss1_stats.txt
    """
}

// Base on per-sample stats for all interval-VCF files across the whole genome,
// the process generate a list of samples that have the percentage of genotype
// missiningness less than a given threshold.
process PANDAS_IMISS1_SAMPLES_TO_KEEP {
    tag "sample_to_keep_imiss1"

    input:
    vcf_stats: Bag<Path>
    max_imiss: Float

    output:
    file("samples_to_keep_imiss1.txt")

    script:
    """
    get_samples_with_max_imiss.py --bcftools_psc_stats_files *stats.txt --max_imiss ${max_imiss} --out samples_to_keep_imiss1.txt
    """

    stub:
    """
    touch samples_to_keep_imiss1.txt
    """
}

// This process does the actual sample filtering (based on missingness).
// After samples are filtering, a step added to ensure
// the remaining sites are segregating.
process BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS1 {
    tag "${interval}"

    input:
    record(interval: String, bcf: Path)
    samples_to_keep: Path

    output:
    record(interval: interval, bcf: file("${interval}_fmiss1_imiss1.bcf"))

    script:
    """
    bcftools index ${bcf}; bcftools view -S ${samples_to_keep} ${bcf} | bcftools view -i "MAC>0" -Ob -o ${interval}_fmiss1_imiss1.bcf
    """

    stub:
    """
    touch ${interval}_fmiss1_imiss1.bcf
    """
}

// This process is similar to `BCFTOOLS_VIEW_F_MISSING1_WITH_PSC_STATS` but uses
// a stricter threshold.
process BCFTOOLS_VIEW_F_MISSING2_WITH_PSC_STATS {
    tag "${interval}"

    input:
    record(interval: String, bcf: Path)
    max_vmiss: Float

    output:
    record(interval: interval, bcf: file("${interval}_fmiss2.bcf"), stats: file("${interval}_fmiss2_stats.txt"))

    script:
    """
    bcftools view -i "F_MISSING<${max_vmiss}"  ${bcf} -Ob -o ${interval}_fmiss2.bcf
    bcftools stats -s- ${interval}_fmiss2.bcf > ${interval}_fmiss2_stats.txt
    """

    stub:
    """
    touch ${interval}_fmiss2.bcf ${interval}_fmiss2_stats.txt
    """
}

// This process is similar to `PANDAS_IMISS1_SAMPLES_TO_KEEP` but uses
// a stricter threshold.
process PANDAS_IMISS2_SAMPLES_TO_KEEP {
    tag "sample_to_keep_imiss2"

    input:
    vcf_stats: Bag<Path>
    max_imiss: Float

    output:
    file("samples_to_keep_imiss2.txt")

    script:
    """
    get_samples_with_max_imiss.py --bcftools_psc_stats_files *stats.txt --max_imiss ${max_imiss} --out samples_to_keep_imiss2.txt
    """

    stub:
    """
    touch samples_to_keep_imiss2.txt
    """
}

// This process is similar to `BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS1` but uses
// a stricter threshold.
process BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS2 {
    tag "${interval}"

    input:
    record(interval: String, bcf: Path)
    samples_to_keep: Path

    output:
    record(interval: interval, bcf: file("${interval}_fmiss2_imiss2.bcf"))

    script:
    """
    bcftools index ${bcf}; bcftools view -S ${samples_to_keep} ${bcf} | bcftools view -i "MAC>0" -Ob -o ${interval}_fmiss2_imiss2.bcf
    """

    stub:
    """
    touch ${interval}_fmiss2_imiss2.bcf
    """
}

// This process filter out rare variants
process BCFTOOLS_VIEW_FILTER_SITE_BY_MAF {
    tag "${interval}"

    input:
    record(interval: String, bcf: Path)
    min_maf: Float

    output:
    record(interval: interval, vcf: file("${interval}_fmiss2_imiss2_maf.vcf.gz"))

    script:
    """
    bcftools view -q ${min_maf}:minor ${bcf} -Oz -o ${interval}_fmiss2_imiss2_maf.vcf.gz
    """

    stub:
    """
    touch ${interval}_fmiss2_imiss2_maf.vcf.gz
    """
}

// --------------------------- SUBWORKFLOW --------------------------------------------
// This subworkflow filters the set of VCF files by variant type, number of
// alleles per site, and missingness per site and per sample.

workflow WF_FILTER_VCF {
    take:
    input_vcfs: Channel<IntervalVcf>
    opts: FilterVcfParams

    main:
    // filter by filter=PASS and -m2 and -M2
    ch_pass_biallelic = BCFTOOLS_VIEW_PASS_BIALLELIC(input_vcfs, opts.samples_included)

    // Filter by vmiss (round 1)
    ch_fmiss1 = BCFTOOLS_VIEW_F_MISSING1_WITH_PSC_STATS(ch_pass_biallelic, opts.high_vmiss)

    // Filter by imiss (round 1)
    val_imiss1_samples = PANDAS_IMISS1_SAMPLES_TO_KEEP(ch_fmiss1.map { r -> r.stats }.collect(), opts.high_imiss)
    ch_imiss1 = BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS1(ch_fmiss1, val_imiss1_samples)

    // Filter by vmiss (round 2)
    ch_fmiss2 = BCFTOOLS_VIEW_F_MISSING2_WITH_PSC_STATS(ch_imiss1, opts.low_vmiss)

    // Filter by imiss (round 2)
    val_imiss2_samples = PANDAS_IMISS2_SAMPLES_TO_KEEP(ch_fmiss2.map { r -> r.stats }.collect(), opts.low_imiss)
    ch_imiss2 = BCFTOOLS_VIEW_FILTER_SAMPLE_BY_IMISS2(ch_fmiss2, val_imiss2_samples)

    // Filter by maf
    ch_maf = BCFTOOLS_VIEW_FILTER_SITE_BY_MAF(ch_imiss2, opts.min_maf)

    // Gather every file produced for an interval into a single record
    ch_filtered = ch_pass_biallelic.map { r -> record(interval: r.interval, pass_bcf: r.bcf) }
        .join(ch_fmiss1.map { r -> record(interval: r.interval, fmiss1_bcf: r.bcf, fmiss1_stats: r.stats) }, by: 'interval')
        .join(ch_imiss1.map { r -> record(interval: r.interval, imiss1_bcf: r.bcf) }, by: 'interval')
        .join(ch_fmiss2.map { r -> record(interval: r.interval, fmiss2_bcf: r.bcf, fmiss2_stats: r.stats) }, by: 'interval')
        .join(ch_imiss2.map { r -> record(interval: r.interval, imiss2_bcf: r.bcf) }, by: 'interval')
        .join(ch_maf, by: 'interval')

    val_sample_lists = val_imiss1_samples.combine(val_imiss2_samples)
        .map { imiss1, imiss2 -> record(imiss1_samples: imiss1, imiss2_samples: imiss2) }

    emit:
    filtered: Channel<FilteredVcf> = ch_filtered
    sample_lists: Value<FilterSampleLists> = val_sample_lists
}

// --------------------------- PUBLISHING ---------------------------------------------
// Where each output file goes, relative to the output directory. These are
// shared by the output blocks of this script and main.nf.

def filtered_vcf_publish_paths(r: FilteredVcf) -> Map {
    def d = '01_filt_vcf'
    return [
        (r.pass_bcf):     "${d}/01_pass_biallelic/",
        (r.fmiss1_bcf):   "${d}/02_fmissing1/",
        (r.fmiss1_stats): "${d}/02_fmissing1/",
        (r.imiss1_bcf):   "${d}/04_imiss1_bcf/",
        (r.fmiss2_bcf):   "${d}/05_fmissing2/",
        (r.fmiss2_stats): "${d}/05_fmissing2/",
        (r.imiss2_bcf):   "${d}/07_imiss2_bcf/",
        (r.vcf):          "${d}/08_maf_bcf/",
    ]
}

def filter_sample_lists_publish_paths(r: FilterSampleLists) -> Map {
    def d = '01_filt_vcf'
    return [
        (r.imiss1_samples): "${d}/03_imiss1_samples/",
        (r.imiss2_samples): "${d}/06_imiss2_samples/",
    ]
}

// -------------------------- TESTING -----------------------------------------------
// The workflow is designed to efficiently process large datasets
// by breaking them down into smaller VCF files, each containing a specific
// genome interval or chunk for all samples. The workflow is designed to run in
// parallel, allowing for faster processing times, while also gathering related
// information as needed.

workflow {
    main:
    def rows = params.input_vcf_fn.splitCsv(header: false, sep: '\t').take(6)
    ch_input_vcf = channel.fromList(rows).map { row ->
        assert row.size() == 2
        record(interval: row[0], vcf: file(row[1], type: 'file', checkIfExists: true))
    }

    def opts = record(
        samples_included: params.samples_included,
        high_vmiss: params.high_vmiss,
        high_imiss: params.high_imiss,
        low_vmiss: params.low_vmiss,
        low_imiss: params.low_imiss,
        min_maf: params.min_maf,
    )
    filt = WF_FILTER_VCF(ch_input_vcf, opts)

    publish:
    filt_vcf = filt.filtered
    filt_vcf_samples = filt.sample_lists
}

output {
    filt_vcf: Channel<FilteredVcf> {
        path { r -> filtered_vcf_publish_paths(r) }
    }
    filt_vcf_samples: FilterSampleLists {
        path { r -> filter_sample_lists_publish_paths(r) }
    }
}

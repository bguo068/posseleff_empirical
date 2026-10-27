nextflow.enable.types = true

// Parameters used when this script is run on its own (see the entry workflow
// at the bottom). When included by main.nf, main.nf's params block is used.
params {
    input_subpop_fn: Path = "${launchDir}/input_subpops.tsv"
}

record IntervalVcf {
    interval: String
    vcf: Path
}

record SubpopSamples {
    grp: String
    sample_list: Path
}

record MergedBcf {
    bcf: Path
    csi: Path
}

// One population (group): its VCF and the Fws table of its samples
record PopVcf {
    grp: String
    vcf: Path
    fws: Path
}

// Comparator for sorting genome intervals
def cmp_interval(i: String, j: String) -> Integer {
    def split_i = i.tokenize(':')
    def chr_i = split_i[0] // chromosome name of interval i
    def start_i = split_i.size() < 2 ? 0 : split_i[1].tokenize('-')[0].toInteger() // start of interval i
    def split_j = j.tokenize(':')
    def chr_j = split_j[0] // chromosome name of interval j
    def start_j = split_j.size() < 2 ? 0 : split_j[1].tokenize('-')[0].toInteger() // start of interval j
    return chr_i <=> chr_j ?: start_i <=> start_j
}

// This process serves two purposes:
// - Concatenate interval-VCF files across the genome.
// - Simplify annotations fields so file sizes are minimal
process BCFTOOLS_CONCAT_MERGE_REGIONS_SIMPLIFY_ANNOTATION {
    tag "merge_regions_samplify_annot"

    input:
    vcfs: List<Path> // sorted by genome interval

    output:
    record(bcf: file("merged.bcf"), csi: file("merged.bcf.csi"))

    script:
    // NOTE: moimix needs format/AD to be in the vcf files
    """
    bcftools concat ${vcfs.join(' ')} --threads ${task.cpus} | bcftools annotate -x INFO,^FORMAT/GT,^FORMAT/AD -Ob -o merged.bcf;
    bcftools index merged.bcf
    """

    stub:
    """
    touch merged.bcf merged.bcf.csi
    """
}

// This process involves dividing the full dataset into population subsets or
// groups (grps) and preparing VCF files for Fws metric calculation.
//
// NOTE: use `--force-samples` option to ignore errors when some requested
// samples are not in the vcf file
process BCFTOOLS_VIEW_EXTRACT_SUBPOP_SAMPLES {
    tag "${grp}"

    input:
    record(grp: String, sample_list: Path)
    merged: MergedBcf

    output:
    record(grp: grp, vcf: file("${grp}.vcf.gz"))

    script:
    """
    bcftools view -S ${sample_list} ${merged.bcf} -Oz -o ${grp}.vcf.gz --force-samples
    """

    stub:
    """
    touch ${grp}.vcf.gz
    """
}

// This process calculate moi based population VCF files using bcftk/calcfws
//
// Generate Fws table for each samples of a population
process BCFTK_MOIMIX_FWS {
    tag "${grp}"

    input:
    record(grp: String, vcf: Path)

    output:
    record(grp: grp, fws: file("${grp}.fws.txt"))

    script:
    """
    bcftk calcfws ${vcf} --output ${grp}.fws.txt
    """

    stub:
    """
    touch ${grp}.fws.txt
    """
}

workflow WF_CALC_POP_FWS {
    take:
    input_vcfs: Channel<IntervalVcf>
    subpops: Channel<SubpopSamples>

    main:
    val_sorted_vcfs = input_vcfs.collect()
        .map { rs -> rs.toSorted { a, b -> cmp_interval(a.interval, b.interval) }.collect { r -> r.vcf }.toList() }

    val_merged = BCFTOOLS_CONCAT_MERGE_REGIONS_SIMPLIFY_ANNOTATION(val_sorted_vcfs)

    ch_subpop_vcf = BCFTOOLS_VIEW_EXTRACT_SUBPOP_SAMPLES(subpops, val_merged)

    ch_fws = BCFTK_MOIMIX_FWS(ch_subpop_vcf)

    ch_pops = ch_subpop_vcf.join(ch_fws, by: 'grp')

    emit:
    merged: Value<MergedBcf> = val_merged
    pops: Channel<PopVcf> = ch_pops
}

// --------------------------- PUBLISHING ---------------------------------------------
// Where each output file goes, relative to the output directory. These are
// shared by the output blocks of this script and main.nf.

def merged_bcf_publish_paths(r: MergedBcf) -> Map {
    def d = '02_pop_fws'
    return [
        (r.bcf): "${d}/01_merge_regions_simplify_annot/",
        (r.csi): "${d}/01_merge_regions_simplify_annot/",
    ]
}

def pop_vcf_publish_paths(r: PopVcf) -> Map {
    def d = '02_pop_fws'
    return [
        (r.vcf): "${d}/02_subpop_vcf/",
        (r.fws): "${d}/04_moimix_fws/",
    ]
}

// --------------------- FOR TESTING -------------------------------
// will be ignored when included as subworkflow script
workflow {
    main:
    ch_input_vcf = channel.fromPath("results/01_filt_vcf/08_maf_bcf/*.vcf.gz", checkIfExists: true)
        .map { p -> record(interval: p.simpleName.replaceAll(/_fmiss2_imiss2_maf/, ''), vcf: p) }

    def rows = params.input_subpop_fn.splitCsv(sep: '\t')
    ch_input_subpop = channel.fromList(rows).map { row ->
        assert row.size() == 2
        record(grp: row[0], sample_list: file(row[1], checkIfExists: true, type: 'file'))
    }

    pop_fws = WF_CALC_POP_FWS(ch_input_vcf, ch_input_subpop)

    publish:
    pop_fws_merged = pop_fws.merged
    pop_fws = pop_fws.pops
}

output {
    pop_fws_merged: MergedBcf {
        path { r -> merged_bcf_publish_paths(r) }
    }
    pop_fws: Channel<PopVcf> {
        path { r -> pop_vcf_publish_paths(r) }
    }
}

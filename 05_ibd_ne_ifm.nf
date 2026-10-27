nextflow.enable.types = true

// Parameters used when this script is run on its own (see the entry workflow
// at the bottom). When included by main.nf, main.nf's params block is used.
params {
    // imputed vcf (glob string, work with pattern such as '*.vcf.gz')
    vcf: String = 'input/*.vcf.gz'  // test data can be found ./test_data/SAM_imputed.vcf.gz
    meta: Path?

    chrname_map: Path = "${projectDir}/assets/chrname_map.tsv"
    num_chrs: Integer = 14

    // hmmibd
    hmmibd_n: Integer = 100
    hmmibd_m: Integer = 5

    // vcf/ibd filtering parameters
    r: Float = 6.67e-7
    mincm: Float = 2.0

    // peak filtering methods
    peak_validate_meth: String = 'ihs' // 'xirs' or 'ihs'

    ibdne_no_diploid_convertion: String = "true"
    ibdne_mincm: String = "2" // "2,4"
    ibdne_minregion: String = "10" // "10,20,50"
    ibdne_flatmeth: String = 'none' // "none,merge,keep_hap_1_only"

    // infomap parameter grid
    ifm_transform_lst: String = "square" // "square,cube,none"
    ifm_mincm_lst: String = "4" // "2,4,6"
    ifm_mingwcm_lst: String = "5" //  "2,4,5,12"
    ifm_ntrials_lst: String = "1000"
    ifm_rmchr_lst: String = "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14" // 0 means not removing any chromosome when call infomap
    // otherwise iterate the list and each time remove a single chr before running infomap
}

record GroupVcf {
    grp: String
    vcf: Path
}

record IbdParams {
    chrname_map: Path
    num_chrs: Integer
    r: Float
    mincm: Float
    hmmibd_n: Integer
    hmmibd_m: Integer
    peak_validate_meth: String
    ibdne_no_diploid_convertion: String
    ibdne_mincm: String
    ibdne_minregion: String
    ibdne_flatmeth: String
    ifm_transform_lst: String
    ifm_mincm_lst: String
    ifm_mingwcm_lst: String
    ifm_ntrials_lst: String
    ifm_rmchr_lst: String
    meta: Path?
}

// The VCF of one chromosome of a group and the IBD segments called from it
record ChrIbd {
    grp: String
    chrno: Integer
    vcf: Path
    csi: Path
    ibd: Path
}

record NeInput {
    label: String
    orig_sh: Path
    orig_map: Path
    orig_ibd: Path
    rmpeaks_sh: Path
    rmpeaks_map: Path
    rmpeaks_ibd: Path
    ibddist_obj: Path
    cov_objs: Set<Path>
    ne_objs: Set<Path>
}

record NeOutput {
    label: String
    are_peaks_removed: Boolean
    ne: Path
}

record IfmInput {
    label: String
    orig: Path
    rmpeaks: Path
}

record IfmOutput {
    label: String
    are_peaks_removed: Boolean
    member: Path
}

def to_lst(value: String) -> List<String> {
    return value.tokenize(',')
}

// Rename chromosome names into integers and subset sites for a chromosome
process RENAME_CHR_FOR_VCF {
    tag "vcf_chr_${grp}"

    input:
    record(grp: String, vcf: Path, chrno: Integer)
    chr_name_map: Path

    output:
    record(grp: grp, chrno: chrno, vcf: file("${grp}_chr${chrno}.vcf.gz"), csi: file("${grp}_chr${chrno}.vcf.gz.csi"))

    script:
    """
    bcftools index -f ${vcf}
    bcftools annotate --rename-chr ${chr_name_map} ${vcf} -Oz -o ${grp}_renamed.vcf.gz
    bcftools index ${grp}_renamed.vcf.gz
    bcftools view -r ${chrno} ${grp}_renamed.vcf.gz -Oz -o ${grp}_chr${chrno}.vcf.gz
    bcftools index ${grp}_chr${chrno}.vcf.gz
    """

    stub:
    """
    touch ${grp}_chr${chrno}.vcf.gz{,.csi}
    """
}

// This process call IBD for a single chromosome for all samples within the input vcf files
process CALL_IBD {
    tag "hmmibd_${grp}_ch${chrno}"

    input:
    record(grp: String, chrno: Integer, vcf: Path, csi: Path)
    opts: Record

    output:
    record(grp: grp, chrno: chrno, vcf: vcf, csi: csi, ibd: file("*_hmmibd.ibd"))

    script:
    def args = [
        vcf: "${grp}_${chrno}.vcf.gz",
        chrno: chrno,
        label: grp,
        r: opts.r, mincm: opts.mincm, n: opts.hmmibd_n, m: opts.hmmibd_m,
    ]
    def args_local = args.keySet().collect { k -> "--${k} ${args[k]}" }.join(" ")
    """
    bcftools view -r ${chrno} ${vcf} -Oz -o ${grp}_${chrno}.vcf.gz
    call_ibd.py ${args_local}
    """

    stub:
    def prefix = "${grp}_${chrno}"
    """
    touch ${prefix}_hmmibd.ibd
    """
}

process PROC_DIST_NE {
    tag "${label}"

    input:
    record(
        label: String,
        ibd_files: List<Path>,
        vcf_files: List<Path>,
        ibdne_mincm: String,
        ibdne_minregion: String,
        ibdne_flatmeth: String
    )
    peak_validate_meth: String
    ibdne_no_diploid_conversion: String

    output:
    record(
        label: label,
        jar: file("ibdne.jar"),
        orig_sh: file("*_orig.sh"),
        orig_map: file("*_orig.map"),
        orig_ibd: file("*_orig.ibd.gz"),
        rmpeaks_sh: file("*_rmpeaks.sh"),
        rmpeaks_map: file("*_rmpeaks.map"),
        rmpeaks_ibd: file("*_rmpeaks.ibd.gz"),
        ibddist_obj: file("*_ibddist.ibdobj.gz"),
        cov_objs: files("*.ibdcov.ibdobj.gz"),
        ne_objs: files("*.ibdne.ibdobj.gz")
    )

    script:
    def args = [
        ibd_files: ibd_files.join(' '),
        vcf_files: vcf_files.join(' '),
        peak_validate_meth: peak_validate_meth,
        ibdne_no_diploid_conversion: ibdne_no_diploid_conversion,
        label: label,
        ibdne_mincm: ibdne_mincm,
        ibdne_minregion: ibdne_minregion,
        ibdne_flatmeth: ibdne_flatmeth,
    ]
    def args_local = args.keySet().collect { k -> "--${k} ${args[k]}" }.join(" ")
    """
    proc_dist_ne.py ${args_local}
    """

    stub:
    """
    touch ibdne.jar
    touch ${label}_{orig,rmpeaks}.{sh,map,ibd.gz}
    touch ${label}_ibddist.ibdobj.gz
    touch ${label}_{orig,rmpeaks}.ibdne.ibdobj.gz
    touch ${label}_orig_all.ibdcov.ibdobj.gz
    touch ${label}_orig_unrel.ibdcov.ibdobj.gz
    """
}

process PROC_INFOMAP {
    tag "${label}"

    input:
    record(label: String, ibd_files: List<Path>, vcf_files: List<Path>)
    peak_validate_meth: String
    ibdne_no_diploid_conversion: String

    output:
    record(label: label, orig: file("*_ifm_orig.ibdobj.gz"), rmpeaks: file("*_ifm_rmpeaks.ibdobj.gz"))

    script:
    def args = [
        ibd_files: ibd_files.join(' '),
        vcf_files: vcf_files.join(' '),
        peak_validate_meth: peak_validate_meth,
        ibdne_no_diploid_conversion: ibdne_no_diploid_conversion,
        label: label,
    ]
    def args_local = args.keySet().collect { k -> "--${k} ${args[k]}" }.join(" ")
    """
    proc_infomap.py ${args_local}
    """

    stub:
    """
    touch ${label}{_ifm_orig.ibdobj.gz,_ifm_rmpeaks.ibdobj.gz}
    """
}

process RUN_IBDNE {
    tag "${label}_${are_peaks_removed}"

    input:
    record(label: String, jar: Path, sh: Path, gmap: Path, ibd: Path, are_peaks_removed: Boolean)

    output:
    record(label: label, are_peaks_removed: are_peaks_removed, ne: file("*.ne"))

    script:
    """
    bash ${sh}
    """

    stub:
    def src = are_peaks_removed ? "rmpeaks" : "orig"
    """
    touch ${label}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag "${label}_${are_peaks_removed}"

    input:
    record(
        label: String,
        ibd_obj: Path,
        are_peaks_removed: Boolean,
        transform: String,
        ifm_mincm: String,
        ifm_mingwcm: String,
        ntrials: String,
        ifm_rmchr: String
    )
    meta: Path?

    output:
    record(label: label, are_peaks_removed: are_peaks_removed, member: file("*_member.pq"))

    script:
    def cut_mode = are_peaks_removed ? 'rmpeaks' : 'orig'
    def args = [
        transform: transform,
        ifm_mincm: ifm_mincm,
        ifm_mingwcm: ifm_mingwcm,
        ntrials: ntrials,
        ifm_rmchr: ifm_rmchr,
        ibd_obj: ibd_obj,
        meta: meta,
        label: label,
        cut_mode: cut_mode,
    ]
    def args_local = args.keySet().collect { k -> args[k] ? "--${k} ${args[k]}" : " " }.join(" ")
    """
    run_infomap.py ${args_local}
    """

    stub:
    def cut_mode = are_peaks_removed ? 'rmpeaks' : 'orig'
    def ifm_params_str = [transform, ifm_mincm, ifm_mingwcm, ntrials, ifm_rmchr].join("_")
    """
    touch ${label}_${cut_mode}_${ifm_params_str}_member.pq
    """
}



workflow WF_IBD_ANALYSES {
    take:
    grp_vcfs: Channel<GroupVcf>
    opts: IbdParams

    main:
    // one record per group and chromosome
    ch_grp_chr = grp_vcfs.flatMap { r ->
        (1..opts.num_chrs).collect { chrno -> record(grp: r.grp, vcf: r.vcf, chrno: chrno) }
    }

    ch_grp_chr_vcf = RENAME_CHR_FOR_VCF(ch_grp_chr, opts.chrname_map)

    ch_chr_ibd = CALL_IBD(ch_grp_chr_vcf, opts)

    // gather the chromosomes of each group, sorted by chromosome number
    ch_ibd_gw = ch_chr_ibd
        .map { r -> tuple(r.grp, opts.num_chrs, r) }
        .groupBy()
        .map { label, rs ->
            def sorted = rs.toSorted { a, b -> a.chrno <=> b.chrno }
            record(label: label, ibd_files: sorted.collect { r -> r.ibd }.toList(), vcf_files: sorted.collect { r -> r.vcf }.toList())
        }

    def proc_ne_params = to_lst(opts.ibdne_mincm).collectMany { ibdne_mincm ->
        to_lst(opts.ibdne_minregion).collectMany { ibdne_minregion ->
            to_lst(opts.ibdne_flatmeth).collect { ibdne_flatmeth ->
                record(ibdne_mincm: ibdne_mincm, ibdne_minregion: ibdne_minregion, ibdne_flatmeth: ibdne_flatmeth)
            }
        }
    }

    ch_in_proc_dist_ne = ch_ibd_gw.flatMap { r -> proc_ne_params.collect { ne_params -> r + ne_params } }

    ch_ne_input = PROC_DIST_NE(ch_in_proc_dist_ne, opts.peak_validate_meth, opts.ibdne_no_diploid_convertion)

    ch_ifm_input = PROC_INFOMAP(ch_ibd_gw, opts.peak_validate_meth, opts.ibdne_no_diploid_convertion)

    ch_ne_in = ch_ne_input.flatMap { r -> [
        record(label: r.label, jar: r.jar, sh: r.orig_sh, gmap: r.orig_map, ibd: r.orig_ibd, are_peaks_removed: false),
        record(label: r.label, jar: r.jar, sh: r.rmpeaks_sh, gmap: r.rmpeaks_map, ibd: r.rmpeaks_ibd, are_peaks_removed: true),
    ] }

    ch_ne_output = RUN_IBDNE(ch_ne_in)

    def ifm_params_grid = to_lst(opts.ifm_transform_lst).collectMany { trans ->
        to_lst(opts.ifm_mincm_lst).collectMany { mincm ->
            to_lst(opts.ifm_mingwcm_lst).collectMany { mingwcm ->
                to_lst(opts.ifm_ntrials_lst).collectMany { ntrials ->
                    to_lst(opts.ifm_rmchr_lst).collect { rmchr ->
                        record(transform: trans, ifm_mincm: mincm, ifm_mingwcm: mingwcm, ntrials: ntrials, ifm_rmchr: rmchr)
                    }
                }
            }
        }
    }.findAll { g -> g.ifm_mincm.toFloat() <= g.ifm_mingwcm.toFloat() }

    ch_ifm_in = ch_ifm_input.flatMap { r ->
        def ibd_objs = [
            record(label: r.label, ibd_obj: r.orig, are_peaks_removed: false),
            record(label: r.label, ibd_obj: r.rmpeaks, are_peaks_removed: true),
        ]
        ibd_objs.collectMany { ibd -> ifm_params_grid.collect { g -> ibd + g } }
    }

    ch_ifm_output = RUN_INFOMAP(ch_ifm_in, opts.meta)

    emit:
    chr_ibd: Channel<ChrIbd> = ch_chr_ibd
    ne_input: Channel<NeInput> = ch_ne_input
    ne_output: Channel<NeOutput> = ch_ne_output
    ifm_input: Channel<IfmInput> = ch_ifm_input
    ifm_output: Channel<IfmOutput> = ch_ifm_output
}

// --------------------------- PUBLISHING ---------------------------------------------
// Where each output file goes, relative to the output directory. These are
// shared by the output blocks of this script and main.nf.

def chr_ibd_publish_paths(r: ChrIbd) -> Map {
    def d = '05_ibdanalysis'
    return [
        (r.vcf): "${d}/01_vcf_chr/",
        (r.csi): "${d}/01_vcf_chr/",
        (r.ibd): "${d}/02_callibd/hmmibd/",
    ]
}

def ne_input_publish_paths(r: NeInput) -> Map {
    def d = "05_ibdanalysis/${r.label}"
    return [
        (r.orig_sh):     "${d}/ne_input/",
        (r.orig_map):    "${d}/ne_input/",
        (r.orig_ibd):    "${d}/ne_input/",
        (r.rmpeaks_sh):  "${d}/ne_input/",
        (r.rmpeaks_map): "${d}/ne_input/",
        (r.rmpeaks_ibd): "${d}/ne_input/",
        (r.ibddist_obj): "${d}/ibdobj/",
        (r.cov_objs):    "${d}/ibdobj/",
        (r.ne_objs):     "${d}/ibdobj/",
    ]
}

def ne_output_publish_paths(r: NeOutput) -> Map {
    return [(r.ne): "05_ibdanalysis/${r.label}/ne_output/"]
}

def ifm_input_publish_paths(r: IfmInput) -> Map {
    def d = "05_ibdanalysis/${r.label}"
    return [
        (r.orig):    "${d}/ifm_input/",
        (r.rmpeaks): "${d}/ifm_input/",
    ]
}

def ifm_output_publish_paths(r: IfmOutput) -> Map {
    return [(r.member): "05_ibdanalysis/${r.label}/ifm_output/"]
}

// Note this will not be called when included as a module
workflow {
    main:
    ch_vcf = channel.fromPath(params.vcf)
        .map { vcf -> record(grp: vcf.simpleName.replaceAll('_.*', ''), vcf: vcf) }

    def opts = record(
        chrname_map: params.chrname_map,
        num_chrs: params.num_chrs,
        r: params.r,
        mincm: params.mincm,
        hmmibd_n: params.hmmibd_n,
        hmmibd_m: params.hmmibd_m,
        peak_validate_meth: params.peak_validate_meth,
        ibdne_no_diploid_convertion: params.ibdne_no_diploid_convertion,
        ibdne_mincm: params.ibdne_mincm,
        ibdne_minregion: params.ibdne_minregion,
        ibdne_flatmeth: params.ibdne_flatmeth,
        ifm_transform_lst: params.ifm_transform_lst,
        ifm_mincm_lst: params.ifm_mincm_lst,
        ifm_mingwcm_lst: params.ifm_mingwcm_lst,
        ifm_ntrials_lst: params.ifm_ntrials_lst,
        ifm_rmchr_lst: params.ifm_rmchr_lst,
        meta: params.meta,
    )
    ibd = WF_IBD_ANALYSES(ch_vcf, opts)

    publish:
    ibd_params = params
    ibd_chr = ibd.chr_ibd
    ibd_ne_input = ibd.ne_input
    ibd_ne_output = ibd.ne_output
    ibd_ifm_input = ibd.ifm_input
    ibd_ifm_output = ibd.ifm_output
}

output {
    // the parameters of the run, saved as an index file
    ibd_params: Record {
        index {
            path '05_ibdanalysis/params.json'
        }
    }
    ibd_chr: Channel<ChrIbd> {
        path { r -> chr_ibd_publish_paths(r) }
    }
    ibd_ne_input: Channel<NeInput> {
        path { r -> ne_input_publish_paths(r) }
    }
    ibd_ne_output: Channel<NeOutput> {
        path { r -> ne_output_publish_paths(r) }
    }
    ibd_ifm_input: Channel<IfmInput> {
        path { r -> ifm_input_publish_paths(r) }
    }
    ibd_ifm_output: Channel<IfmOutput> {
        path { r -> ifm_output_publish_paths(r) }
    }
}

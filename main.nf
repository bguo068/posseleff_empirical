nextflow.enable.types = true

include { WF_FILTER_VCF ; filtered_vcf_publish_paths ; filter_sample_lists_publish_paths } from "./01_filter_vcf.nf"
include { WF_CALC_POP_FWS ; merged_bcf_publish_paths ; pop_vcf_publish_paths } from "./02_calc_pop_fws.nf"
include { WF_PHASE_IMPUTE_MONOCLONAL ; monoclonal_pop_publish_paths } from "./03_phase_impute_monoclonal.nf"
include {
    WF_PHASE_IMPUTE_POLYCLONAL ;
    polyclonal_split_publish_paths ;
    polyclonal_deploid_publish_paths ;
    polyclonal_merged_publish_paths ;
    monoclonal_only_publish_paths
} from "./04_phase_impute_polyclonal.nf"
include {
    WF_IBD_ANALYSES ;
    chr_ibd_publish_paths ;
    ne_input_publish_paths ;
    ne_output_publish_paths ;
    ifm_input_publish_paths ;
    ifm_output_publish_paths
} from "./05_ibd_ne_ifm.nf"

params {
    // ------------ Inputs -----------------
    // Tab-separated table without header: <interval> <vcf>
    input_vcf_fn: Path = "${launchDir}/input_vcfs.tsv"
    // Tab-separated table without header: <group name> <sample list file>
    input_subpop_fn: Path = "${launchDir}/input_subpops.tsv"
    // Relative paths in the two tables above are resolved against the project directory.

    // ------------ Submodule 1: filter VCF -----------------
    // Optional whitelist of samples, to manually remove duplicated or low quality samples
    samples_included: Path?
    high_vmiss: Float = 0.7
    high_imiss: Float = 0.7
    low_vmiss: Float = 0.3
    low_imiss: Float = 0.3
    min_maf: Float = 0.01

    // ------------ Submodules 3 and 4: phase and impute -----------------
    beagle_jar: Path = "${projectDir}/lib/beagle.08Feb22.fa4.jar"
    gmap: Path = "${projectDir}/assets/Pf3D7_all.txt"
    deploid_nopanel_options: String = "-rate 8 -burn 0.67 -k 4"
    beagle_options_monoclonal: String = ''
    deploid_panel_size: Integer = 16
    deploid_panel_options: String = "-nSample 250 -rate 8 -burn 0.67 -k 4 -ibd"
    beagle_options_polyclonal: String = 'overlap=3.0' // default is 2.0, cause error for some vcf files

    // ------------ Submodule 5: IBD analyses -----------------
    meta: Path?
    chrname_map: Path = "${projectDir}/assets/chrname_map.tsv"
    test: Boolean
    num_chrs: Integer? // default: 2 if test else 14

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

// Absolute paths are used as is; relative paths are relative to the project directory
def resolve_path(p: String) -> Path {
    return file(p.startsWith('/') || p.contains('://') ? p : "${projectDir}/${p}")
}

workflow {
    main:
    ch_input_vcf = channel.fromList(params.input_vcf_fn.splitCsv(header: false, sep: '\t'))
        .map { row ->
            assert row.size() == 2
            record(interval: row[0], vcf: resolve_path(row[1]))
        }

    ch_input_subpop = channel.fromList(params.input_subpop_fn.splitCsv(sep: '\t'))
        .map { row ->
            assert row.size() == 2
            record(grp: row[0], sample_list: resolve_path(row[1]))
        }

    // ------------ Run submodule 1 -----------------
    filt = WF_FILTER_VCF(ch_input_vcf, record(
        samples_included: params.samples_included,
        high_vmiss: params.high_vmiss,
        high_imiss: params.high_imiss,
        low_vmiss: params.low_vmiss,
        low_imiss: params.low_imiss,
        min_maf: params.min_maf,
    ))

    // ------------ Run submodule 2 -----------------
    pop_fws = WF_CALC_POP_FWS(filt.filtered, ch_input_subpop)

    // -----------  Run sumodule 3 -----------------
    mono = WF_PHASE_IMPUTE_MONOCLONAL(pop_fws.pops, record(
        deploid_nopanel_options: params.deploid_nopanel_options,
        beagle_jar: params.beagle_jar,
        gmap: params.gmap,
        beagle_options: params.beagle_options_monoclonal,
    ))

    // -----------  Run sumodule 4 -----------------
    poly = WF_PHASE_IMPUTE_POLYCLONAL(pop_fws.pops, mono, record(
        deploid_panel_size: params.deploid_panel_size,
        deploid_panel_options: params.deploid_panel_options,
        beagle_jar: params.beagle_jar,
        gmap: params.gmap,
        beagle_options: params.beagle_options_polyclonal,
    ))

    // -----------  Run sumodule 5 -----------------
    ibd = WF_IBD_ANALYSES(poly.phased, record(
        chrname_map: params.chrname_map,
        num_chrs: params.num_chrs ?: (params.test ? 2 : 14),
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
    ))

    publish:
    filt_vcf = filt.filtered
    filt_vcf_samples = filt.sample_lists
    pop_fws_merged = pop_fws.merged
    pop_fws = pop_fws.pops
    phase_impute_monoclonal = mono
    poly_split = poly.split
    poly_deploid = poly.deploid
    poly_merged = poly.merged
    poly_mono_only = poly.mono_only
    ibd_params = params
    ibd_chr = ibd.chr_ibd
    ibd_ne_input = ibd.ne_input
    ibd_ne_output = ibd.ne_output
    ibd_ifm_input = ibd.ifm_input
    ibd_ifm_output = ibd.ifm_output
}

// Each output is published under the output directory (`-output-dir` or the
// `outputDir` config setting, default: results). See the *_publish_paths
// functions of each submodule for the subdirectory of each file.
output {
    // 01_filt_vcf
    filt_vcf {
        path { r -> filtered_vcf_publish_paths(r) }
    }
    filt_vcf_samples {
        path { r -> filter_sample_lists_publish_paths(r) }
    }

    // 02_pop_fws
    pop_fws_merged {
        path { r -> merged_bcf_publish_paths(r) }
    }
    pop_fws {
        path { r -> pop_vcf_publish_paths(r) }
    }

    // 03_phase_impute_monoclonal
    phase_impute_monoclonal {
        path { r -> monoclonal_pop_publish_paths(r) }
    }

    // 04_phase_impute_polyclonal
    poly_split {
        path { r -> polyclonal_split_publish_paths(r) }
    }
    poly_deploid {
        path { r -> polyclonal_deploid_publish_paths(r) }
    }
    poly_merged {
        path { r -> polyclonal_merged_publish_paths(r) }
    }
    poly_mono_only {
        path { r -> monoclonal_only_publish_paths(r) }
    }

    // 05_ibdanalysis
    // the parameters of the run, saved as an index file
    ibd_params {
        index {
            path '05_ibdanalysis/params.json'
        }
    }
    ibd_chr {
        path { r -> chr_ibd_publish_paths(r) }
    }
    ibd_ne_input {
        path { r -> ne_input_publish_paths(r) }
    }
    ibd_ne_output {
        path { r -> ne_output_publish_paths(r) }
    }
    ibd_ifm_input {
        path { r -> ifm_input_publish_paths(r) }
    }
    ibd_ifm_output {
        path { r -> ifm_output_publish_paths(r) }
    }
}

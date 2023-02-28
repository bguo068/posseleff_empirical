nextflow.enable.dsl=2

include { WF_FILTER_VCF              } from "./01_filter_vcf.nf"
include { WF_CALC_POP_FWS            } from "./02_calc_pop_fws.nf"
include { WF_PHASE_IMPUTE_MONOCLONAL } from "./03_phase_impute_monoclonal.nf"
include { WF_PHASE_IMPUTE_POLYCLONAL } from "./04_phase_impute_polyclonal.nf"
include { WF_IBD_ANALYSES            } from "./05_ibd_ne_ifm.nf"

workflow {
  ch_input_vcf = channel
        .fromPath(params.input_vcf_fn, type: 'file', checkIfExists: true )
        .splitCsv(header: false, sep: '\t')
        .map{
            assert it.size() == 2
            def interval = it[0]
            def vcf = file(it[1], type: 'file', checkIfExists: params.test? false:true)
            return [interval, vcf]
        }

    ch_input_subpop = channel.fromPath(params.input_subpop_fn)
    .splitCsv(sep: '\t')
    .map{
        assert it.size() == 2
        def grp = it[0]
        def sample_list_file = file(it[1], checkIfExists: true, type: 'file')
        return [grp, sample_list_file]
    }
    
    // ------------ Run submodule 1 -----------------
    WF_FILTER_VCF(ch_input_vcf)


    // ------------ Run submodule 2 -----------------
    ch_filterred_vcf = WF_FILTER_VCF.out
    WF_CALC_POP_FWS(ch_filterred_vcf, ch_input_subpop)

/*
    // -----------  Run sumodule 3 -----------------
    ch_fws = WF_CALC_POP_FWS.out[0]
    ch_subpop_vcf = WF_CALC_POP_FWS.out[1]
    WF_PHASE_IMPUTE_MONOCLONAL(ch_fws, ch_subpop_vcf)

    // -----------  Run sumodule 4 -----------------
    ch_subpop_deploidclonalvcf = WF_PHASE_IMPUTE_MONOCLONAL.out
    WF_PHASE_IMPUTE_POLYCLONAL(ch_fws, ch_subpop_vcf, ch_subpop_deploidclonalvcf)

    // -----------  Run sumodule 5 -----------------
    ch_phased_imputed_grp_vcf = WF_PHASE_IMPUTE_POLYCLONAL.out
    // select one reps for high and low transmission settings
    ch_phased_imputed_grp_vcf.filter{grp, subpop_vcf-> grp in ["ESEA", "WAF"] }
    // Run
    WF_IBD_ANALYSES(ch_phased_imputed_grp_vcf)    
*/
}
import java.nio.file.*
import groovy.json.*

nextflow.enable.dsl=2

params.chrname_map = "$projectDir/assets/chrname_map.tsv"
params.num_chrs = params.test ? 2: 14

// hmmibd
params.hmmibd_n = 100
params.hmmibd_m = 5

// vcf/ibd filtering parameters
params.r = 6.67e-7
params.maf = params.test ? 0.00001: 0.01
params.mincm = 2.0

params.ibdne_mincm = 4
params.ibdne_minregion = 10

// infomap parameter grid
params.ifm_transform_lst = "square,cube,none"
params.ifm_mincm_lst = "2,4,6"
params.ifm_mingwcm_lst = "2,4,5,12"
params.ifm_ntrials_lst = "1000"
// params.ifm_rmchr_lst="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14"
params.ifm_rmchr_lst="0" // 0 means not removing any chromosome when call infomap
// otherwise iterate the list and each time remove a single chr before running infomap

// imputed vcf (glob string, work with pattern such as '*.vcf.gz')
params.vcf = 'input/*.vcf.gz'  // test data can be found ./test_data/SAM_imputed.vcf.gz
params.meta = ""

def publish_dir = "${params.outdir ?: launchDir}/05_ibdanalysis"

def to_lst(value){
    def lst = "${value}".split(",").collect{it};
    return lst
}

// logging params
process LOG_PARAMS {
    exec:
        def params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
        file("$publish_dir").mkdirs()
        file("$publish_dir/params.log").text = params_str
}


// Rename chromosome names into integers and subset sites for a chromosome
process RENAME_CHR_FOR_VCF {
    tag "vcf_chr_${grp}"
    publishDir "$publish_dir/01_vcf_chr"
    input: tuple val(grp), path(vcf), path(chr_name_map), val(chrno)
    output: tuple val(grp), val(chrno), \
        path('*_chr*.vcf.gz'), path("*_chr*.vcf.gz.csi")
    script: """
    bcftools index -f ${vcf}
    bcftools annotate --rename-chr ${chr_name_map} ${vcf} -Oz -o ${grp}_renamed.vcf.gz
    bcftools index ${grp}_renamed.vcf.gz
    bcftools view -r ${chrno} ${grp}_renamed.vcf.gz -Oz -o ${grp}_chr${chrno}.vcf.gz
    bcftools index ${grp}_chr${chrno}.vcf.gz
    """
    stub:
    """touch ${grp}_chr${chrno}.vcf.gz{,.csi}"""

}

// This process call IBD for a single chromosome for all samples within the input vcf files
process CALL_IBD {
    tag "hmmibd_${label}_ch${chrno}"
    publishDir "$publish_dir/02_callibd/hmmibd", pattern: "*_hmmibd.ibd", mode: 'symlink'

    input:
        tuple val(label), val(chrno), path(vcf), path(index)
    output:
        tuple val(label), val(chrno), path("*_hmmibd.ibd"), emit: hmmibd

    script:
    def args_local = [
        vcf: "${label}_${chrno}.vcf.gz",
        chrno: chrno,
        label: label,
        r: params.r, mincm: params.mincm, n:params.hmmibd_n, m:params.hmmibd_m,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    bcftools view -r ${chrno} ${vcf} -Oz -o ${label}_${chrno}.vcf.gz
    call_ibd.py ${args_local}
    """
    stub:
    def prefix="${label}_${chrno}"
    """
    touch ${prefix}_hmmibd.ibd
    """
}

process PROC_DIST_NE {
    tag "${label}"

    publishDir "${publish_dir}/${label}/ne_input/", pattern: "*.sh", mode: 'symlink'
    publishDir "${publish_dir}/${label}/ne_input/", pattern: "*.map", mode: 'symlink'
    publishDir "${publish_dir}/${label}/ne_input/", pattern: "*.ibd.gz", mode: 'symlink'
    publishDir "${publish_dir}/${label}/ibdobj/", pattern: "*.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label), path(ibd_lst), path(vcf_lst)
    output:
        tuple val(label), path("ibdne.jar"), path("*_orig.sh"), \
                path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
        tuple val(label), path("ibdne.jar"), path("*_rmpeaks.sh"),  \
                path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
        tuple val(label), path("*_ibddist.ibdobj.gz"), emit: ibddist_ibd_obj
        tuple val(label), path("*.ibdcov.ibdobj.gz"), emit: cov_ibd_obj
        tuple val(label), path("*.ibdne.ibdobj.gz"), emit: ne_ibd_obj
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        label: label,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
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

    publishDir "${publish_dir}/${label}/ifm_input/", pattern: "*.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label), path(ibd_lst), path(vcf_lst)

    output:
        tuple val(label), path("*_ifm_orig.ibdobj.gz"), emit: ifm_orig_ibd_obj
        tuple val(label), path("*_ifm_rmpeaks.ibdobj.gz"), emit: ifm_rmpeaks_ibd_obj
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        label: label,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
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

    publishDir "${publish_dir}/${label}/ne_output/",  mode: 'symlink'

    input:
        tuple val(label), path(ibdne_jar), path(ibdne_sh), path(gmap), path(ibd_gz), \
            val(are_peaks_removed)
    output:
        tuple val(label), val(are_peaks_removed), path("*.ne")
    script:
    """
    bash ${ibdne_sh}
    """
    stub:
    def src = are_peaks_removed ? "rmpeaks": "orig"
    """
    touch ${label}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag "${label}_${are_peaks_removed}"
    publishDir "${publish_dir}/${label}/ifm_output/",  mode: 'symlink'
    input:
        tuple val(label), path(ibd_obj), val(are_peaks_removed), val(ifm_params)
    output:
        tuple val(label), val(are_peaks_removed), path("*_member.pq")
    script:
    def meta = params.meta ? file(params.meta) : ''
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    def args_local = (ifm_params + [ // add infomap param here - ifm_params
        ibd_obj: ibd_obj,
        meta: meta,
        label: label,
        cut_mode: cut_mode,
    ]).collect{k, v-> v ? "--${k} ${v}": " "}.join(" ")
    """
    run_infomap.py ${args_local}
    """
    stub:
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    def ifm_params_str = [ifm_params.transform, ifm_params.ifm_mincm, 
        ifm_params.ifm_mingwcm, ifm_params.ntrials,ifm_params.ifm_rmchr].join("_")
    """
    touch ${label}_${cut_mode}_${ifm_params_str}_member.pq
    """
}



workflow WF_IBD_ANALYSES {
    take:
        ch_grp_vcf // tuple val(grp), path(vcf)
    main:

        LOG_PARAMS()

        RENAME_CHR_FOR_VCF(
            ch_grp_vcf.map{grp,vcf->[grp,vcf,file(params.chrname_map)]}
            .combine(Channel.fromList(1..params.num_chrs))
            ) // [grp, vcf, chrnamemap, chrno]

        ch_grp_chr_vcf = RENAME_CHR_FOR_VCF.out
             // [grp, chrno, vcf, idx]

        CALL_IBD(ch_grp_chr_vcf)

        ch_ibd_gw = CALL_IBD.out.hmmibd.combine(
                ch_grp_chr_vcf.map{grp,chrno,vcf,index->[grp,chrno,vcf]},
                by: [0, 1]
            )
            .map{label, chrno, ibd, vcf ->
                [ groupKey(label, params.num_chrs),  [chrno, ibd, vcf]  ]}
            .groupTuple(by: 0, sort: {a, b -> a[0]<=> b[0]} )
            .map{label, ll->
                [label, ll.collect{arr->arr[1]},ll.collect{arr->arr[2]} ] }

        PROC_DIST_NE(ch_ibd_gw)

        PROC_INFOMAP(ch_ibd_gw)

        ch_ne_in = (
            PROC_DIST_NE.out.ne_input_orig.map{label, ibdne_jar, script, gmap, ibd ->
            [label, ibdne_jar, script, gmap, ibd , false] }
         ) .concat(
            PROC_DIST_NE.out.ne_input_rmpeaks.map{label, ibdne_jar, script, gmap, ibd ->
            [label, ibdne_jar, script, gmap, ibd , true] }
         )

        RUN_IBDNE(ch_ne_in)

        ch_ifm_params_grid = Channel.fromList(to_lst(params.ifm_transform_lst))
            .combine( Channel.fromList(to_lst(params.ifm_mincm_lst)) )
            .combine( Channel.fromList(to_lst(params.ifm_mingwcm_lst)) )
            .combine( Channel.fromList(to_lst(params.ifm_ntrials_lst)) )
            .combine( Channel.fromList(to_lst(params.ifm_rmchr_lst)) )
            .filter{trans, mincm, mingwcm, ntrials, rmchr -> mincm<=mingwcm}
            .map {trans, mincm, mingwcm, ntrials, rmchr -> 
                [transform: trans, ifm_mincm: mincm, ifm_mingwcm: mingwcm, 
                    ntrials: ntrials, ifm_rmchr: rmchr]
            }//.view()

        ch_ifm_in = (
            PROC_INFOMAP.out.ifm_orig_ibd_obj.map{label, ibd->[label, ibd, false]}
        ).concat(
            PROC_INFOMAP.out.ifm_rmpeaks_ibd_obj.map{label, ibd->[label, ibd, true]}
        ).combine(ch_ifm_params_grid) // [label, ibd, t_or_f, ifm_params]

        RUN_INFOMAP(ch_ifm_in)

        RUN_INFOMAP.out//.view()
}

// Note this will not be called when included as a module
workflow {

    ch_vcf = Channel.fromPath(params.vcf)
        .map{vcf -> def grp = vcf.getSimpleName().replaceAll('_.*', ''); [grp, vcf]}

    WF_IBD_ANALYSES(ch_vcf)

}

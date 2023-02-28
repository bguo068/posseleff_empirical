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

// infomap parameters
params.ifm_transform = ["square", "cube", "none"][0]
params.ifm_ntrails = 1000

// imputed vcf (glob string, work with pattern such as '*.vcf.gz')
params.vcf = 'input/*.vcf.gz'  
params.meta = ""

def publish_dir = "${params.outdir ?: launchDir}/05_ibdanalysis"

// logging params
process LOG_PARAMS {
    exec:
        def params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
        file("$publish_dir").mkdirs()
        file("$publish_dir/params.log").text = params_str
}


// Rename chromosome names into integers
process RENAME_CHR_FOR_VCF {
    tag "rename_vcf_chr_${grp}"
    publishDir "$publish_dir/01_rename_vcf_chr"
    input: tuple val(grp), path(vcf), path(chr_name_map)
    output: tuple val(grp), path('*_chrs_renamed.vcf.gz'), path("*_chrs_renamed.vcf.gz.csi")
    script: """
    bcftools index -f ${vcf}
    bcftools annotate --rename-chr ${chr_name_map} ${vcf} -Oz -o ${grp}_chrs_renamed.vcf.gz 
    bcftools index ${grp}_chrs_renamed.vcf.gz 
    """
    stub: """touch ${grp}_chrs_renamed.vcf.gz{,.csi}"""

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
    touch ${prefix}{_hmmibd.ibd}
    """
}

process PROC_DIST_NE {
    tag "${label}"

    publishDir "${publish_dir}/${label}/ne_input/", pattern: "*.sh", mode: 'symlink'
    publishDir "${publish_dir}/${label}/ne_input/", pattern: "*.map", mode: 'symlink'
    publishDir "${publish_dir}/${label}/ne_input/", pattern: "*.ibd.gz", mode: 'symlink'
    publishDir "${publish_dir}/${label}/ibddist_ibd/", pattern: "*_ibddist_ibd.pq", mode: 'symlink'

    input:
        tuple val(label), path(ibd_lst)
    output:
        tuple val(label), path("ibdne.jar"), path("*_orig.sh"), \
                path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
        tuple val(label), path("ibdne.jar"), path("*_rmpeaks.sh"),  \
                path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
        tuple val(label), path("*_ibddist_ibd.pq"), emit: ibddist_ibd_pq
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        label: label,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_dist_ne.py ${args_local} 
    """
    stub:
    """
    touch ibdne.jar
    touch ${label}{_orig.sh,_orig.map,_orig.ibd.gz}
    touch ${label}{_rmpeaks.sh,_rmpeaks.map,_rmpeaks.ibd.gz}
    touch ${label}_ibddist_ibd.pq
    """
}

process PROC_INFOMAP {
    tag "${label}"

    publishDir "${publish_dir}/${label}/ifm_input/", pattern: "*_ibd.pq", mode: 'symlink'

    input:
        tuple val(label), path(ibd_lst)
    output:
        tuple val(label), path("*_ifm_orig_ibd.pq"), emit: ifm_orig_ibd_pq
        tuple val(label), path("*_ifm_rmpeaks_ibd.pq"), emit: ifm_rmpeaks_ibd_pq
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        label: label,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_infomap.py ${args_local}
    """
    stub:
    """
    touch ${label}{_ifm_orig_ibd.pq,_ifm_rmpeaks_ibd.pq}
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
    publishDir "${publish_dir}/${label}_${label}/ifm_output/",  mode: 'symlink'
    input:
        tuple val(label), path(ibd_pq), val(are_peaks_removed)
    output:
        tuple val(label), val(are_peaks_removed), path("*_member.pq")
    script:
    def meta = params.meta ? file(params.meta) : ''
    def args_local = [
        ibd_pq: ibd_pq,
        meta: meta,
        label: label,
        ntrails: params.ifm_ntrails,
        transform: params.ifm_transform,
    ].collect{k, v-> v ? "--${k} ${v}": " "}.join(" ")
    """
    run_infomap.py ${args_local}
    """
    stub:
    """
    touch ${label}_member.pq
    """
}



workflow WF_IBD_ANALYSES {
    take:
        ch_grp_vcf // tuple val(grp), path(vcf)
    main:

        LOG_PARAMS()

        RENAME_CHR_FOR_VCF( ch_grp_vcf.map{grp,vcf->[grp,vcf,file(params.chrname_map)]})
        
        ch_grp_chr_vcf = RENAME_CHR_FOR_VCF.out.combine(Channel.fromList(1..params.num_chrs))
            .map {grp, vcf, index, chrom-> [grp, chrom, vcf, index]}
        
        CALL_IBD(ch_grp_chr_vcf)

        ch_ibd_gw = CALL_IBD.out.hmmibd.map{label, chrno, ibd -> 
            [ groupKey(label, params.num_chrs),  [chrno, ibd]  ]}
            .groupTuple(by: 0, sort: {a, b -> a[0]<=> b[0]} )
            .map{label, ll-> [label, ll.collect{pair->pair[1]} ] }
        
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

        ch_ifm_in = (
            PROC_INFOMAP.out.ifm_orig_ibd_pq.map{label, ibd->[label, ibd, false]}
        ).concat(
            PROC_INFOMAP.out.ifm_rmpeaks_ibd_pq.map{label, ibd->[label, ibd, true]}
        )

        RUN_INFOMAP(ch_ifm_in)
}

// Testing only. It will not be called when included as a module
workflow {

    ch_vcf = Channel.fromPath(params.vcf)
        .map{vcf -> def grp = vcf.getSimpleName().replaceAll('_.*', ''); [grp, vcf]}

    WF_IBD_ANALYSES(ch_vcf)

}
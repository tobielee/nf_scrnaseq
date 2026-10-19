/*
You will need to define clusterres in nextflow config prior to running this pipeline. 
This pipeline will help define cell-annotation based on your specified cluster res.
*/  
project_dir = "${workflow.projectDir}"
scripts_dir = "${project_dir}/scripts"

// TODO need to maybe borrow method from previous pipeline for integration using seurat object as input
// problem is that there are multiple .h5ad but single .rds
// maybe I save .pkl to correspond to .rds for backprop
process integrate_singlecell {
    debug true
    publishDir "${params.intoutdir}/${datlabel}", mode: 'copy'

    input:
    val indir
    val datlabel
    val mode
    val annofield1
    val annofield2

    output:
    path "*_merged.h5ad", emit: integrated_files

    script:
    """
    conda run -n scseq python ${scripts_dir}/nf_scanvi_integration-9.py \
        --datlabel "${datlabel}" \
        --indir "${indir}" \
        --outdir "." \
        --mode "${mode}" \
        --annofield1 "${annofield1}" \
        --annofield2 "${annofield2}"
    """
}

workflow {
    integrate_singlecell(
        params.indir,
        params.datlabel,
        params.mode,
        params.annofield1,
        params.annofield2
    )
    integrate_singlecell.out.integrated_files.view()
}

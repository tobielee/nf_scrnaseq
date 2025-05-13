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
    publishDir params.intoutdir, mode: 'copy', saveAs: { file -> "${file}" }

    input:
    val indir

    output:
    path "scanvi_integration/*", emit: intplot_files
    path '*_merged.h5ad', emit: integrated_files
    //  // output plots

    script:
    """
    conda run -n scseq python $scripts_dir/nf_scanvi_integration-9.py \
        --datlabel '${params.datlabel}' \
        --indir '${indir}' \
        --outdir 'scanvi_integration' \
    """
}

workflow {
    test = integrate_singlecell(params.intoutdir)
    test.integrated_files.view()
}

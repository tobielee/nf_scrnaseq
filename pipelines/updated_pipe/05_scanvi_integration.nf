nextflow.enable.dsl = 2

/*
 * Integrates annotated h5ad or Seurat RDS inputs from intoutdir using ScanVI.
 * mode defaults to h5ad. annofield1 is the reviewed source label; annofield2
 * names the derived label column used to train ScanVI after unknown labels
 * have been normalized.
 */


// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process integrate_singlecell {

    debug true
    tag "${params.datlabel}"
    publishDir params.intoutdir, mode: 'copy'

    input:
    val indir
    path integration_script

    output:
    path "scanvi_integration", emit: integration_results

    script:
    """
    echo "Integration input directory : ${indir}"
    echo "Integration script          : ${integration_script}"

    python ${integration_script} \
        --datlabel '${params.datlabel}' \
        --indir '${indir}' \
        --outdir 'scanvi_integration' \
        --mode '${params.mode}' \
        --annofield1 '${params.annofield1}' \
        --annofield2 '${params.annofield2}'
    """
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    // ------------------------------------------------------------------
    // Directories & scripts
    // ------------------------------------------------------------------

    def project_dir = workflow.projectDir
    def scripts_dir = "${project_dir}/pipelines/scripts"

    def integration_script = file(
        "${scripts_dir}/nf_scanvi_integration-9.py"
    )


    // ------------------------------------------------------------------
    // Logging
    // ------------------------------------------------------------------

    log.info """
    SCANVI INTEGRATION PIPELINE
    =======================================
    dataset_name : ${params.datlabel}
    input        : ${params.intoutdir}
    output       : ${params.intoutdir}
    mode         : ${params.mode}
    labels       : ${params.annofield1} -> ${params.annofield2}
    """


    // ------------------------------------------------------------------
    // Run integration
    // ------------------------------------------------------------------

    integrated = integrate_singlecell(
        params.intoutdir,
        integration_script
    )


    // ------------------------------------------------------------------
    // Display integrated files
    // ------------------------------------------------------------------

    integrated.integration_results.view()
}



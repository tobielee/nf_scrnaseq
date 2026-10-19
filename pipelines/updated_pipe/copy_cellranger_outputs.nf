nextflow.enable.dsl = 2

/*
 * Stage filtered CellRanger count matrices for one or more CellRanger runs.
 *
 * Required parameters:
 *   --cellrangerroot   Parent directory containing one subdirectory per run ID
 *   --datlabel         Run ID to copy (or shared config default)
 *   --cellrangermulti  1 for `cellranger multi`, 0 for `cellranger count`
 *   --outdir           Destination for the staged matrix files
 *
 * For multiple runs, pass --datlabels 'run1,run2'. Each label is also used as
 * the run directory name under cellrangerroot and as a filename prefix.
 */

params.cellrangerroot = null
params.cellrangermulti = null
params.datlabels = null

process copy_cellranger_outputs {
    tag "${datlabel}"
    publishDir params.outdir, mode: 'copy', saveAs: { file -> file.name }

    input:
    tuple val(datlabel), val(cellrangermulti), val(cellrangeroutdir)
    path copy_script

    output:
    path "${datlabel}/*", emit: copied_matrices

    script:
    """
    python3 ${copy_script} \
        --datlabel '${datlabel}' \
        --cellrangermulti '${cellrangermulti}' \
        --cellrangeroutpath '${cellrangeroutdir}'
    """
}

workflow {
    if (params.cellrangerroot == null || params.cellrangerroot.toString().isBlank()) {
        error "Set --cellrangerroot to the parent directory containing the CellRanger run directories."
    }
    if (params.cellrangermulti == null || !(params.cellrangermulti.toString() in ['0', '1'])) {
        error "Set --cellrangermulti to 0 for cellranger count or 1 for cellranger multi."
    }
    if (params.outdir == null || params.outdir.toString().isBlank() || params.outdir == '/path/to/analysis/outs') {
        error "Set --outdir to the directory where the copied count matrices should be published."
    }

    def datlabels = params.datlabels
        ? params.datlabels.toString().split(',').collect { it.trim() }.findAll { it }
        : [params.datlabel?.toString()?.trim()].findAll { it }

    if (!datlabels || datlabels.any { !it }) {
        error "Provide at least one run ID with --datlabel or --datlabels."
    }

    def scriptsDir = "${workflow.projectDir}/pipelines/scripts"
    def copyScript = file("${scriptsDir}/nf_cellrangerouts_copy.py")
    def runInputs = datlabels.collect { datlabel ->
        tuple(datlabel, params.cellrangermulti.toString(), "${params.cellrangerroot}/${datlabel}")
    }

    copy_cellranger_outputs(
        Channel.fromList(runInputs),
        file(copyScript)
    )
}

nextflow.enable.dsl = 2

/*
    You will need to define clusterres in nextflow config prior to running
    this pipeline.

    This pipeline will help define cell annotation based on your specified
    cluster resolution.
*/


// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process annotate_h5ad {

    debug true
    tag "${params.datlabel}"
    publishDir params.intoutdir, mode: 'copy'

    input:
    tuple path(anno_csv), path(anndata)
    path annotation_script

    output:
    path '*_annotated1.h5ad', emit: annotated_files

    script:
    """
    echo "Annotation CSV : ${anno_csv}"
    echo "AnnData file   : ${anndata}"
    echo "Python script  : ${annotation_script}"

    python ${annotation_script} \
        --annocsv '${anno_csv}' \
        --adatafile '${anndata}' \
        --annofield1 '${params.annofield1}'
    """
}


process backpropagate_annotations {

    debug true
    tag "${params.datlabel}"
    publishDir params.intoutdir, mode: 'copy'

    input:
    path anndataAnno
    path seuratAnno
    path backprop_script
    val  scripts_dir

    output:
    path '*_annotated1.rds', emit: backpropagated_files

    script:
    """
    echo "Seurat RDS     : ${seuratAnno}"
    echo "Annotated h5ad : ${anndataAnno}"
    echo "Scripts dir    : ${scripts_dir}"
    echo "R script       : ${backprop_script}"

    Rscript ${backprop_script} \
        '${scripts_dir}' \
        --seuratcountfile '${seuratAnno}' \
        --annfile '${anndataAnno}' \
        --integrate '${params.integrate}' \
        --sample_id '${params.sample_id}' \
        --annofield1 '${params.annofield1}'
    """
}


process integrate_singlecell {

    debug true
    tag "${params.datlabel}"
    publishDir params.intoutdir, mode: 'copy'

    input:
    val indir
    path integration_script

    output:
    path "scanvi_integration/*", emit: intplot_files
    path '*_merged.h5ad', emit: integrated_files

    script:
    """
    echo "Input directory : ${indir}"
    echo "Python script   : ${integration_script}"

    python ${integration_script} \
        --datlabel '${params.datlabel}' \
        --indir '${indir}' \
        --outdir 'scanvi_integration'
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

    def annotation_script = file(
        "${scripts_dir}/nf_annotate-8.py"
    )

    def backprop_script = file(
        "${scripts_dir}/nf_annotations_h5ad_to_rds-8.R"
    )

    def integration_script = file(
        "${scripts_dir}/nf_scanvi_integration-9.py"
    )

    def norm_filtered_data = params.removedoub == 1
        ? "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}_dubremoved.rds"
        : "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}.rds"
    def clustername = params.subcluster
        ? "${params.subclustername}_clustered"
        : '_clustered'
    def annotate1 = norm_filtered_data.replace(
        '_unlabeled', clustername + params.clusterres + '_sctype'
    )
    def expected_h5ad = norm_filtered_data
        .replace(
            '_unlabeled', clustername + params.clusterres + '_sctype_scanvi_celltypist'
        )
        .replace('.rds', '.h5ad')

    if (params.integrate != 'NONE') {
        error "Manual annotation currently requires integrate='NONE' so each sample CSV maps to its own h5ad."
    }


    // ------------------------------------------------------------------
    // Locate annotation inputs
    // ------------------------------------------------------------------

    def h5ad_files = file(
        "${params.outdir}/${params.subdirectory_name}/*.h5ad"
    ).collect { it }

    def csv_files = file(
        "${params.outdir}/*_annotations.csv"
    ).collect { it }

    if (h5ad_files.isEmpty()) {
        error "No annotated h5ad files found in ${params.outdir}/${params.subdirectory_name}"
    }
    if (csv_files.isEmpty()) {
        error "No manual annotation CSV files found in ${params.outdir}"
    }


    // ------------------------------------------------------------------
    // Pair each annotation CSV with its corresponding h5ad
    // ------------------------------------------------------------------

    def file_pairs = csv_files.collect { csv_file ->

        def sample_id = csv_file.getName().replaceFirst(/_annotations\.csv$/, '')
        def expected_name = "${sample_id}_${expected_h5ad}"
        def matching_h5ad = h5ad_files.findAll { it.getName() == expected_name }
        if (matching_h5ad.size() != 1) {
            error "Expected exactly one h5ad for '${csv_file.getName()}', found ${matching_h5ad*.getName()}"
        }
        [csv_file, matching_h5ad[0]]

    }.findAll { it != null }

    def paired_h5ad_names = file_pairs.collect { it[1].getName() }.toSet()
    def unpaired_h5ads = h5ad_files.findAll { !paired_h5ad_names.contains(it.getName()) }
    if (!unpaired_h5ads.isEmpty()) {
        error "Annotated h5ad files have no matching manual annotation CSV: ${unpaired_h5ads*.getName()}"
    }


    log.info """
    SCANVI / CELLTYPIST ANNOTATION PIPELINE
    =======================================
    dataset_name : ${params.datlabel}
    outdir       : ${params.outdir}
    intoutdir    : ${params.intoutdir}
    h5ad files   : ${h5ad_files.size()}
    CSV files    : ${csv_files.size()}
    file pairs   : ${file_pairs.size()}
    """


    // ------------------------------------------------------------------
    // Step 1: Annotate h5ad files
    // ------------------------------------------------------------------

    ch_files = Channel.fromList(file_pairs)

    annotated_h5ad_ch = annotate_h5ad(
        ch_files,
        annotation_script
    )


    // ------------------------------------------------------------------
    // Pair each annotated h5ad with the single ScType Seurat RDS produced
    // by the previous stage. That RDS contains the sample objects by name.
    // ------------------------------------------------------------------

    def seurat_annotation_rds = file(
        "${params.outdir}/${annotate1}",
        checkIfExists: true
    )


    // ------------------------------------------------------------------
    // Step 3: Backpropagate annotations to Seurat objects
    // ------------------------------------------------------------------

    backpropagated = backpropagate_annotations(
        annotated_h5ad_ch,
        Channel.value(seurat_annotation_rds),
        backprop_script,
        scripts_dir
    )


    // ------------------------------------------------------------------
    // Optional Step 4: Integrate single-cell data
    // ------------------------------------------------------------------

    // Uncomment if/when the integration step is needed:
    //
    // integrated = integrate_singlecell(
    //     params.intoutdir,
    //     integration_script
    // )

}


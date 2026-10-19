nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process clusterres_selection {

    debug true
    tag "${params.datlabel}"
    publishDir params.outdir, mode: 'copy'

    input:
    path clustered_data
    path clusterres_script
    val  scripts_dir
    val  clusterres_data

    output:
    path "${clusterres_data}", emit: clusterres_data

    script:
    """
    echo "Using R script       : ${clusterres_script}"
    echo "Scripts directory    : ${scripts_dir}"
    echo "Clustered data       : ${clustered_data}"
    echo "Clusterres output    : ${clusterres_data}"

    Rscript ${clusterres_script} \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${clustered_data}' \
        --outfile '${clusterres_data}' \
        --integrate '${params.integrate}' \
        --clusterreslist '${params.clusterreslist}' \
        --clusterres '${params.clusterres}' \
        --subcluster '${params.subcluster}'
    """
}


process annotate_sctype {

    debug true
    tag "${params.datlabel}"
    publishDir params.outdir, mode: 'copy'

    input:
    path clusterres_data
    path annotate_script
    val  scripts_dir
    val  annotate1
    val  annotate_outprefix

    output:
    path "${annotate1}", emit: annotate1
    path "sctype_anno/*"
    path "marker_diffgenes/*"

    script:
    """
    echo "Using R script       : ${annotate_script}"
    echo "Scripts directory    : ${scripts_dir}"
    echo "Clusterres data      : ${clusterres_data}"
    echo "Annotation output    : ${annotate1}"

    Rscript ${annotate_script} \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${clusterres_data}' \
        --norm '${params.norm}' \
        --integrate '${params.integrate}' \
        --outfile '${annotate1}' \
        --outfilepre '${annotate_outprefix}' \
        --subcluster '${params.subcluster}' \
        --subclustername '${params.subclustername}' \
        --sample_id '${params.sample_id}' \
        --clusterres '${params.clusterres}' \
        --markerdir '${params.markerdir}'
    """
}


process annotate_scanvi_celltypist {

    debug true
    tag "${params.datlabel}"
    publishDir params.outdir, mode: 'copy'

    input:
    path annotate1
    path annotation_script
    val  annotate2
    val  queryOutPath

    output:
    path "${params.subdirectory_name}/*", emit: annotated_data

    script:
    """
    echo "Using Python script  : ${annotation_script}"
    echo "Input annotation     : ${annotate1}"
    echo "Output                : ${queryOutPath}"

    python3 ${annotation_script} \
        --datlabel '${params.datlabel}' \
        --infile '${annotate1}' \
        --outfile '${queryOutPath}' \
        --integrate '${params.integrate}' \
        --reflab '${params.reflab}' \
        --scanvi_model_path '${params.scanvi_model_path}' \
        --scanvi_model_batchkey '${params.scanvi_model_batchkey}' \
        --scanvi_model_unknownlab '${params.scanvi_model_unknownlab}' \
        --celltypist_model_path '${params.celltypist_model_path}' \
        --sample_id '${params.sample_id}' \
        --subdirectory_name '${params.subdirectory_name}'
    """
}


process plot_annotations {

    debug true
    tag "${params.datlabel}"
    publishDir params.outdir, mode: 'copy'

    input:
    path plot_input
    path plot_script

    output:
    path "plots/*"
    path "*.csv"

    script:
    if (params.integrate == "NONE") {
        """
        echo "Looking for .pkl file in:"
        echo "${plot_input}"

        PKL_FILE=""

        for file in ${plot_input}; do
            if [[ \$file == *.pkl ]]; then
                PKL_FILE=\$file
                break
            fi
        done

        if [[ -z \$PKL_FILE ]]; then
            echo "ERROR: No .pkl file found in ${plot_input}"
            exit 1
        fi

        echo "Using .pkl file: \$PKL_FILE"

        python3 ${plot_script} \
            --datlabel '${params.datlabel}' \
            --infile "\$PKL_FILE" \
            --sample_id '${params.sample_id}' \
            --reflab '${params.reflab}' \
            --integrate '${params.integrate}'
        """
    } else {
        """
        echo "Looking for .h5ad file in:"
        echo "${plot_input}"

        H5AD_FILE=""

        for file in ${plot_input}; do
            if [[ \$file == *.h5ad ]]; then
                H5AD_FILE=\$file
                break
            fi
        done

        if [[ -z \$H5AD_FILE ]]; then
            echo "ERROR: No .h5ad file found in ${plot_input}"
            exit 1
        fi

        echo "Using .h5ad file: \$H5AD_FILE"

        python3 ${plot_script} \
            --datlabel '${params.datlabel}' \
            --infile "\$H5AD_FILE" \
            --sample_id '${params.sample_id}' \
            --reflab '${params.reflab}' \
            --integrate '${params.integrate}'
        """
    }
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

    def clusterres_script = file("${scripts_dir}/nf_selectclusterres-4.R")
    def annotate_script   = file("${scripts_dir}/nf_annomarkerbased-5.R")
    def annotation_script = file("${scripts_dir}/nf_cellannotation_scanvi_celltypist_refbased_gpu-6.py")
    def plot_script       = file("${scripts_dir}/nf_plot_cellannotations_eda-7.py")


    // ------------------------------------------------------------------
    // Expected filenames
    // ------------------------------------------------------------------

    def norm_filtered_data = params.removedoub == 1
        ? "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}_dubremoved.rds"
        : "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}.rds"

    def clustername = params.subcluster
        ? "${params.subclustername}_clustered"
        : "_clustered"

    def clustered_data = norm_filtered_data.replace(
        '_unlabeled',
        clustername
    )

    def clusterres_data = norm_filtered_data.replace(
        '_unlabeled',
        clustername + params.clusterres
    )

    def annotate1 = norm_filtered_data.replace(
        '_unlabeled',
        clustername + params.clusterres + '_sctype'
    )

    def annotate_outprefix = annotate1.substring(
        0,
        annotate1.indexOf("_sctype")
    )

    def annotate2 = norm_filtered_data
        .replace(
            '_unlabeled',
            clustername + params.clusterres + '_sctype_scanvi_celltypist'
        )
        .replace('.rds', '.h5ad')

    def queryOutPath = "${params.subdirectory_name}/${annotate2}"


    // ------------------------------------------------------------------
    // Logging
    // ------------------------------------------------------------------

    log.info """
    R scRNA-seq PIPELINE
    =======================================
    dataset_name    : ${params.datlabel}
    species         : ${params.species}
    outdir          : ${params.outdir}

    clustered       : ${clustered_data}
    clusterres      : ${clusterres_data}
    annotate1       : ${annotate1}
    annotate2       : ${annotate2}
    """


    // ------------------------------------------------------------------
    // Input from previous pipeline step
    //
    // clusterres_selection starts from the clustered RDS produced by
    // the previous pipeline/process.
    // ------------------------------------------------------------------

    def clustered_input = Channel.fromPath(
        "${params.outdir}/${clustered_data}",
        checkIfExists: true
    )


    // ------------------------------------------------------------------
    // Step 1: Select clustering resolution
    // ------------------------------------------------------------------

    clusterres = clusterres_selection(
        clustered_input,
        clusterres_script,
        scripts_dir,
        clusterres_data
    )


    // ------------------------------------------------------------------
    // Step 2: Cell-type annotation with marker-based ScType
    // ------------------------------------------------------------------

    annotated_sctype = annotate_sctype(
        clusterres,
        annotate_script,
        scripts_dir,
        annotate1,
        annotate_outprefix
    )


    // ------------------------------------------------------------------
    // Step 3: ScanVI + CellTypist annotation
    // ------------------------------------------------------------------

    annotated_data = annotate_scanvi_celltypist(
        annotated_sctype.annotate1,
        annotation_script,
        annotate2,
        queryOutPath
    )


    // ------------------------------------------------------------------
    // Step 4: Plot annotations
    // ------------------------------------------------------------------

    plot_annotations(
        annotated_data,
        plot_script
    )
}


nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process qcFilter {

    debug true
    tag "${params.datlabel}"
    publishDir params.outdir, mode: 'copy'

    input:

    path qc_script
    val  scripts_dir
    val  filtered_rds

    output:

    path "${filtered_rds}", emit: filtered_data

    script:

    """
    echo "Using R script     : ${qc_script}"
    echo "Scripts directory  : ${scripts_dir}"
    echo "Expected output    : ${filtered_rds}"

    Rscript ${qc_script} \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --datadir '${params.datadir}' \
        --outdir '${params.outdir}' \
        --fcellmin ${params.fcellmin} \
        --ffeatmin ${params.ffeatmin} \
        --fquantmin ${params.fquantmin} \
        --fquantmax ${params.fquantmax} \
        --fmito ${params.fmito} \
        --species '${params.species}' \
        --tissue '${params.tissue}' \
        --detectdoub ${params.detectdoub} \
        --removedoub ${params.removedoub} \
        --inputseurat '${params.inputseurat}' \
        --sample_id ${params.sample_id} \
        --parseregex '${params.parseregex}' \
        --parsesplit '${params.parsesplit}'
    """
}


process normIntegrate {

    debug true
    tag "${params.datlabel}"
    publishDir params.outdir, mode: 'copy'

    input:

    path filtered_data
    path norm_script
    val  scripts_dir
    val  norm_rds

    output:

    path "${norm_rds}", emit: norm_filtered_data

    script:

    """
    echo "Using R script     : ${norm_script}"
    echo "Scripts directory  : ${scripts_dir}"
    echo "Filtered data      : ${filtered_data}"
    echo "Expected output    : ${norm_rds}"

    Rscript ${norm_script} \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${filtered_data}' \
        --outfile '${norm_rds}' \
        --norm ${params.norm} \
        --integrate ${params.integrate} \
        --detectdoub ${params.detectdoub} \
        --removedoub ${params.removedoub}
    """
}


process clustering {

    debug true
    tag "${params.datlabel}"
    publishDir params.outdir, mode: 'copy'

    input:

    path norm_filtered_data
    path cluster_script
    val  scripts_dir
    val  clustered_rds

    output:

    path "${clustered_rds}", emit: clustered_data
    path "cluster_deg/*", emit: cluster_deg

    script:

    """
    echo "Using R script     : ${cluster_script}"
    echo "Scripts directory  : ${scripts_dir}"
    echo "Normalized data    : ${norm_filtered_data}"
    echo "Expected output    : ${clustered_rds}"

    Rscript ${cluster_script} \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${norm_filtered_data}' \
        --outfile '${clustered_rds}' \
        --integrate ${params.integrate} \
        --npca ${params.npca} \
        --clusterreslist ${params.clusterreslist} \
        --topnfeat ${params.topn_feat} \
        --reduction ${params.reduction} \
        --sample_id ${params.sample_id} \
        --subcluster ${params.subcluster} \
        --subclustername ${params.subclustername}
    """
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    // ------------------------------------------------------------------
    // Directories & scripts
    // ------------------------------------------------------------------

    project_dir = workflow.projectDir
    scripts_dir = "${project_dir}/pipelines/scripts"

    qc_script = file("${scripts_dir}/nf_preprocessing_doubletdetect-1.R")
    norm_script = file("${scripts_dir}/nf_normintegrate-2.R")
    cluster_script = file("${scripts_dir}/nf_clustering-3.R")


    // ------------------------------------------------------------------
    // Pre-compute expected output filenames
    // ------------------------------------------------------------------

    filtered_rds = params.removedoub == 1 \
        ? "${params.datlabel}_removedDoublets_filtered_list.rds" \
        : params.detectdoub == 1 \
            ? "${params.datlabel}_detectedDoublets_filtered_list.rds" \
            : "${params.datlabel}_filtered_list.rds"

    norm_rds = params.removedoub == 1 \
        ? "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}_dubremoved.rds" \
        : "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}.rds"

    cluster_suffix = params.subcluster \
        ? "${params.subclustername}_clustered" \
        : "_clustered"

    clustered_rds = norm_rds.replace('_unlabeled', cluster_suffix)


    // ------------------------------------------------------------------
    // Logging
    // ------------------------------------------------------------------

    log.info """
    R scRNA-seq PIPELINE
    =======================================
    dataset_name : ${params.datlabel}
    species      : ${params.species}
    filtered     : ${filtered_rds}
    normalized   : ${norm_rds}
    clustered    : ${clustered_rds}
    """


    // ------------------------------------------------------------------
    // Run pipeline
    // ------------------------------------------------------------------

    filtered = qcFilter(
        qc_script,
        scripts_dir,
        filtered_rds
    )

    normalized = normIntegrate(
        filtered,
        norm_script,
        scripts_dir,
        norm_rds
    )

    clustering(
        normalized,
        cluster_script,
        scripts_dir,
        clustered_rds
    )
}

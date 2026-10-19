// Define the qcFilter process
process qcFilter {
    debug true // to output script info to console 
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    input:
    each sample
    path qc_script
    val scripts_dir

    output:
    path "${params.datlabel}_${sample}_${params.fcellmin}-${params.ffeatmin}-${params.fquantmin}-${params.fquantmax}-${params.fmito}_qcviolin.png", emit: plots

    script:
    """
    echo "Using R script: ${qc_script}"
    echo "Scripts directory: ${scripts_dir}"    
    Rscript ${qc_script} \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --datadir '${params.datadir}' \
        --outdir '${params.outdir}' \
        --inputseurat '${params.inputseurat}' \
        --fcellmin ${params.fcellmin} \
        --ffeatmin ${params.ffeatmin} \
        --fquantmin ${params.fquantmin} \
        --fquantmax ${params.fquantmax} \
        --fmito ${params.fmito} \
        --species '${params.species}' \
	--parseregex '${params.parseregex}' \
        --sample ${sample}
    """
}

workflow {
    // Set paths relative to the Nextflow project directory
    project_dir = workflow.projectDir
    scripts_dir = "${project_dir}/pipelines/scripts"
    qc_script = file("${scripts_dir}/nf_qcfiltering_check-0.R")
    
    log.info """
    R scRNA-seq PIPELINE
    =======================================
    dataset_name : ${params.datlabel}
    species      : ${params.species}
    datadir      : ${params.datadir}
    outdir       : ${params.outdir}
    """    

    if (params.inputseurat != null && !params.inputseurat.isBlank()) {
        // TODO need to process data differently 
        // data needs to be split and renormalized but I can just plot it as one dataset for now
        println "Using input seurat object"
	qcFilter(
            Channel.of(params.datlabel),
            qc_script,
	    scripts_dir
        )        
    } else {
        // Sample parsing
        def allFiles = new File(params.datadir).list() // List all files in the specified directory
        def relevantFiles = allFiles.findAll { it =~ /(\.mtx|\.mtx\.gz|\.tsv\.gz|\.tsv)$/ } // Filter relevant files
        println relevantFiles // Output the list of unique sample names
        def regex = ~params.parseregex
        def parsedNames = relevantFiles.collect { fname ->
            def m = fname =~ regex
            m.find() ? m.group(1) : null
        }.findAll { it != null }.unique()

        println "Parsed sample names: $parsedNames"
        qcFilter(
            Channel.fromList(parsedNames),
            qc_script,
	    scripts_dir
        )        
    }
}

project_dir = "${workflow.projectDir}"

log.info """
    R scRNA-seq  P I P E L I N E
    =======================================
    dataset_name    : ${params.datlabel}
    species       : ${params.species}
    outdir          : ${params.outdir}
    """

// processGetSampleNames {
//     publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
//     tag "${params.datlabel}"

//     input:
//     file datadir
//     file parseregex
//     file inputseurat

//     output:
//     file 'sampleList.txt', emit: sampleList

//     script:
//     """
//     echo $project_dir
//     Rscript $project_dir/nf_get_sample_names.R \
//         --datadir '${params.datadir}' \
//         --parseregex '${params.parseregex}' \
//         --inputseurat '${params.inputseurat}' \
//         --outdir '${params.outdir}'
//     """
// }

// Define the qcFilter process
process qcFilter {
    debug true // to output script info to console 
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    input:
    each sample

    output:
    path "${params.datlabel}_${sample}_${params.fcellmin}-${params.ffeatmin}-${params.fquantmin}-${params.fquantmax}-${params.fmito}_qcviolin.png", emit: plots

    script:
    """
    echo $project_dir
    conda run -n scseq Rscript $project_dir/nf_qcfiltering_check-0.R \
        '${project_dir}' \
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
        --sample ${sample}
    """
}

workflow {
    if (params.inputseurat != null && !params.inputseurat.isBlank()) {
        // TODO need to process data differently 
        // data needs to be split and renormalized but I can just plot it as one dataset for now
        println "Using input seurat object"
        qcFilter(params.datlabel)
    } else {
        // Sample parsing
       def allFiles = new File(params.datadir).list() // List all files in the specified directory
        def relevantFiles = allFiles.findAll { it =~ /(\.mtx|\.mtx\.gz|\.tsv\.gz|\.tsv)$/ } // Filter relevant files
        println relevantFiles // Output the list of unique sample names

        def parsedNames = relevantFiles.collect { it.replaceAll(params.parseregex, '$1') } // Parse file names
        def sampleList = parsedNames.unique() // Get unique sample names
        // typechan = Channel.from(sampleList)

        println sampleList // Output the list of unique sample names
        qcFilter(sampleList)
    }
}

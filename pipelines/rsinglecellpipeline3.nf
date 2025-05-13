/*
You will need to define clusterres in nextflow config prior to running this pipeline. 
This pipeline will help define cell-annotation based on your specified cluster res.
*/  
project_dir = "${workflow.projectDir}"
scripts_dir = "${project_dir}/scripts"

// -- TODO these filenameing definitions probably can be elsewhere
def norm_filtered_data
if (params.removedoub == 1) {
    norm_filtered_data = "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}_dubremoved.rds"
} else {
    norm_filtered_data = "${params.datlabel}_unlabeled_${params.integrate}_${params.norm}.rds"
}
def clustered_data
if (params.subcluster) {
    clustername = "${params.subclustername}_clustered"
} else {
    clustername = "_clustered"
}
clustered_data = "${norm_filtered_data.replace('_unlabeled', clustername)}"

def clusterres_data = "${norm_filtered_data.replace('_unlabeled', clustername + params.clusterres)}"

def annotate1 = "${norm_filtered_data.replace('_unlabeled', clustername + params.clusterres + '_sctype')}"


// annotate1 is the one to start with for backpropagating annotations to seurat

def parseFileName(seuratAnno) {
    File infile = new File(seuratAnno)
    String baseFilename = infile.getName()
    String parentDirName = infile.getParent()
    String prefixedFilename = parentDirName + "/scanvi-celltypist_outs/" + baseFilename
    return prefixedFilename.replace('sctype', 'sctype_scanvi_celltypist').replace('.rds', '_annotated1.h5ad')    
}
def parseFileName2(anndataAnno) {
    File infile = new File(anndataAnno)
    String baseFilename = infile.getName().replace('.h5ad', '.rds')
    return baseFilename    
}


process map_seurat {
    input:
    val anndataAnno

    output:
    tuple val(new File(params.outdir, annotate1).toString()), val(anndataAnno), emit: paired_files

    exec:
    println "Processing annotation file: ${anndataAnno}"
}

process annotate_h5ad {
    debug true 

    publishDir params.intoutdir, mode: 'copy', saveAs: { file -> "${file}" }

    input:
    tuple path(anno_csv), path(anndata)
    
    output:
    path '*_annotated1.h5ad', emit: annotated_files

    script:
    """
    echo "Processing annotation file: ${anno_csv} with data file: ${anndata}"
    conda run -n scseq python $scripts_dir/nf_annotate-8.py \
        --annocsv '${anno_csv}' \
        --adatafile '${anndata}' \
        --annofield1 '${params.annofield1}' \
    """
}

process backpropagate_annotations {
    debug true 
    publishDir params.intoutdir, mode: 'copy', saveAs: { file -> "${file}" }

    input:
    tuple path(seuratAnno), path(anndataAnno)
    
    output:
    path '*_annotated1.rds'

    script:
    """
    conda run -n scseq Rscript $scripts_dir/nf_annotations_h5ad_to_rds-8.R \
        '${scripts_dir}' \
        --seuratcountfile '${seuratAnno}' \
        --annfile '${anndataAnno}' \
        --integrate '${params.integrate}' \
        --sample_id '${params.sample_id}' \
        --annofield1 '${params.annofield1}' \
    """
}



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
    // Step 1: Collect h5ad and CSV files for annotation
    def h5ad_files = file("${params.outdir}/scanvi-celltypist_outs/*.h5ad").collect { it }
    def csv_files = file("${params.outdir}/scanvi-celltypist_outs/*.csv").collect { it }

    // Step 2: Generate file pairs based on sample ID
    def file_pairs = csv_files.collect { csv_file ->
        def sample_id = (csv_file.getName() =~ /(.*)_annotations\.csv/)[0][1]
        def matching_h5ad = h5ad_files.find { it.getName().startsWith(sample_id) }
        matching_h5ad ? [csv_file, matching_h5ad] : null
    }.findAll { it != null }

    ch_files = Channel.fromList(file_pairs)

    // Step 3: Annotate h5ad files
    annotated_files = annotate_h5ad(ch_files)
    annotated_files.view()

    // Step 4: Map annotated files to Seurat paths
    tuples = map_seurat(annotated_files)

    // Step 5: Backpropagate annotations
    test = backpropagate_annotations(tuples)
    test.view()
}

// workflow {
//     test = integrate_singlecell(params.intoutdir)
//     test.integrated_files.view()
// }

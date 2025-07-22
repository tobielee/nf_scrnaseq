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
// paste0(OUTSDIR, OUTFILE_PRE, "_list_cellmarker2_panglao_toppcell_sctype_annotated.rds")
def annotate_outprefix = annotate1.substring(0, annotate1.indexOf("_sctype"))

//TODO fix this
def annotate2 = "${norm_filtered_data.replace('_unlabeled', clustername + params.clusterres + '_sctype_scanvi_celltypist').replace('.rds', '.h5ad')}"
def queryOutPath = "${params.subdirectory_name}/${annotate2}"
log.info """
    R scRNA-seq  P I P E L I N E
    =======================================
    dataset_name    : ${params.datlabel}
    species       : ${params.species}
    outdir          : ${params.outdir}
    annotate1      : ${annotate1}
    """
process clusterres_selection {
    debug true // to output script info to console 
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"
    def clustered_data_path = new File(params.outdir, clustered_data).toString() // not sure why I need two of these

    input:
    path(clustered_data_path)

    output:
    path "${clusterres_data}", emit: clusterres_data

    script: // TODO still have dependency of filtered_data annotation and rscript (so need to set both if adjusting)
    """
    conda run -n scseq Rscript $scripts_dir/nf_selectclusterres-4.R \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${clustered_data_path}' \
        --outfile '${clusterres_data}' \
        --integrate ${params.integrate} \
        --clusterreslist ${params.clusterreslist} \
        --clusterres ${params.clusterres} \
        --subcluster ${params.subcluster} \
    """
}

process annotate_sctype {
    debug true // to output script info to console 

    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    input:
    path(clusterres_data)

    output:
    path "${annotate1}", emit: annotate1
    path "sctype_anno/*" // output plots
    path "marker_diffgenes/*" // output csv and rds for degs 

    script:
    """
    conda run -n scseq Rscript $scripts_dir/nf_annomarkerbased-5.R \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${clusterres_data}' \
        --norm ${params.norm} \
        --integrate ${params.integrate} \
        --outfile '${annotate1}' \
        --outfilepre '${annotate_outprefix}' \
        --sample_id '${params.sample_id}' \
        --subcluster ${params.subcluster} \
        --subclustername ${params.subclustername} \
        --clusterres ${params.clusterres} \
        --markerdir ${params.markerdir} 
    """
}

process annotate_scanvi_celltypist {
    debug true // to output script info to console 
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    input:
    path(annotate1)

    output:
    path "${params.subdirectory_name}/*", emit: annotated_data // sample names are prepended right after directory if unintegrated

    script:
    """
    echo $scripts_dir
    echo $norm_filtered_data
    conda run -n scseq python3 $scripts_dir/nf_cellannotation_scanvi_celltypist_refbased_gpu-6.py \
        --datlabel '${params.datlabel}' \
        --infile '${annotate1}' \
        --outfile '${queryOutPath}' \
        --integrate '${params.integrate}' \
        --reflab '${params.reflab}' \
        --scanvi_model_path '${params.scanvi_model_path}' \
        --celltypist_model_path '${params.celltypist_model_path}' \
        --sample_id '${params.sample_id}' \
        --subdirectory_name '${params.subdirectory_name}'
    """
}

process plot_annotations {
    debug true // to output script info to console 
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    input:
    path(plot_input)

    output:
    path "plots/*"
    path "*.csv"

    script:
    if (params.integrate == "NONE") {
        """
        echo $scripts_dir
        PKL_FILE=""
        for file in \$(echo ${plot_input}); do
            if [[ \$file == *.pkl ]]; then
                PKL_FILE=\$file
                break
            fi
        done
        if [[ -z \$PKL_FILE ]]; then
            echo "No .pkl file found in ${plot_input}"
            exit 1  # Exit with an error code to indicate failure
        else
            echo "Using .pkl file \$PKL_FILE from ${plot_input}"
            conda run -n scseq python3 $scripts_dir/nf_plot_cellannotations_eda-7.py \
                --datlabel '${params.datlabel}' \
                --infile "\$PKL_FILE" \
                --sample_id '${params.sample_id}' \
                --reflab '${params.reflab}' \
                --integrate '${params.integrate}'
        fi
        """
    } else {
        """
        INT_FILE=""
        for file in \$(echo ${plot_input}); do
            if [[ \$file == *.h5ad ]]; then
                INT_FILE=\$file
                echo "Using .h5ad file \$INT_FILE from ${plot_input}"
                conda run -n scseq python3 $scripts_dir/nf_plot_cellannotations_eda-7.py \
                    --datlabel '${params.datlabel}' \
                    --infile "\$INT_FILE" \
                    --sample_id '${params.sample_id}' \
                    --reflab '${params.reflab}' \
                    --integrate '${params.integrate}'
            fi
        done
        """
    }
}



workflow {
    def clustered_data_path = new File(params.outdir, clustered_data).toString()
    
    // Step 1: Run clusterres_selection
    test = clusterres_selection(clustered_data_path)
    
    // Step 2: Run annotate_sctype
    test2 = annotate_sctype(test)
    
    // Step 3: Run annotate_scanvi_celltypist
    test3 = annotate_scanvi_celltypist(test2.annotate1)
    
    // Step 4: Determine plot input based on integration status
    plot_annotations(test3)
    // def dic_path = new File("/${queryOutPath}".replace('.h5ad', '.pkl')).toString()
    // def plot_input = params.integrate == "NONE" ? dic_path : test3.annotated_data
    // plot_annotations(plot_input)

    // test3 = "${params.annotated_datapath}"
    // plot_annotations(test3)
}

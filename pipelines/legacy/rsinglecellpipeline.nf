project_dir = "${workflow.projectDir}"
scripts_dir = "${project_dir}/scripts"

process qcFilter {
    debug true // to output script info to console 
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    // Define filtered_data outputfile
    def filtered_data
    if (params.removedoub == 1) {
        filtered_data = "${params.datlabel}_removedDoublets_filtered_list.rds"
    } else if (params.detectdoub == 1) {
        filtered_data = "${params.datlabel}_detectedDoublets_filtered_list.rds"
    } else {
        filtered_data = "${params.datlabel}_filtered_list.rds"
    }

    output:
    path "${filtered_data}", emit: filtered_data
     
    script: // TODO still have dependency of filtered_data annotation and rscript (so need to set both if adjusting)
    """
    echo $scripts_dir
    conda run -n scseq Rscript $scripts_dir/nf_preprocessing_doubletdetect-1.R \
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
        --sample_id ${params.sample_id} \
        --detectdoub ${params.detectdoub} \
        --removedoub ${params.removedoub} \
        --inputseurat '${params.inputseurat}' \
        --parseregex '${params.parseregex}' \
        --parsesplit '${params.parsesplit}' 
    """
}

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

process normIntegrate {
    debug true // to output script info to console 

    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    input:
    path(filtered_data)

    output:
    path "${norm_filtered_data}", emit: norm_filtered_data

    script:
    """
    echo $scripts_dir
    echo $filtered_data
    echo $norm_filtered_data
    conda run -n scseq Rscript $scripts_dir/nf_normintegrate-2.R \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${filtered_data}' \
        --outfile '${norm_filtered_data}' \
        --norm ${params.norm} \
        --integrate ${params.integrate} \
        --detectdoub ${params.detectdoub} \
        --removedoub ${params.removedoub}
    """
}

process clustering {
    debug true // to output script info to console 
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }
    tag "${params.datlabel}"

    input:
    path(norm_filtered_data)

    output:
    path "${clustered_data}", emit: clustered_data
    path "cluster_deg/*" // output plots

    script:
    """
    echo $scripts_dir
    echo $norm_filtered_data
    conda run -n scseq Rscript $scripts_dir/nf_clustering-3.R \
        '${scripts_dir}' \
        --datlabel '${params.datlabel}' \
        --infile '${norm_filtered_data}' \
        --outfile '${clustered_data}' \
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

workflow {
    test = qcFilter() // Assign the output of qcFilter to a variable
    test2 = normIntegrate(test)
    clustering(test2)
    // def norm_filtered_data_path = new File(params.outdir, norm_filtered_data).toString()
    // log.info("Starting workflow with ${norm_filtered_data_path}")
    // clustering(norm_filtered_data_path)
}

// Set script_dir to the current directory
script_dir = "${workflow.projectDir}"

params.cellrangermulti = 1 // 1 true or 0 false depending on your case
params.cellrangeroutpath = '/path/to/cellrangerouts'
params.outdir = '/path/to/select/cellrangerouts'

process copyOuts {
    publishDir params.outdir, mode: 'copy', saveAs: { file -> "${file}" }

    input:
    val datlabel

    script:
    """
    python3 $script_dir/nf_cellrangerouts_copy.py \
        --datlabel '${datlabel}' \
        --cellrangermulti '${params.cellrangermulti}' \
        --cellrangeroutpath '${params.cellrangeroutpath}/${datlabel}'
    """
}

workflow {
    // Create a channel to emit each datlabel
    datlabel_channel = Channel.of(
        'XZ71233_XZ71234'
        //'XZ62858-XZ62859',
        //'XZ68423-XZ68424',
        //'XZ63068-XZ63069',
        //'XZ65028_XZ65032',
        //'XZ69317-XZ69318',
        //'XZ64126-XZ64127',
        //'XZ69320-XZ69321',
        //'XZ64255-XZ64256',
        //'XZ65029_XZ65033',
        //'XZ69803-XZ69804',
        //'XZ65026_XZ65030',
        //'XZ67501-XZ67502',
        //'XZ71045-XZ71046',
        //'XZ65027_XZ65031'
    ) 
    // Start the copyOuts process for each datlabel
    copyOuts(datlabel_channel)
}

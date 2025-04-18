#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Check if help flag was passed (my hacky solution)
help1 = "${params.help}" == "nohelp" ? "nohelp" : "help"
help2 = "${params.h}" == "nohelp" ? "nohelp" : "help"

// Print help message and quit
def printHelp() {
    println """
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FISHY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Global default params:

  --out            Set name for output folder
  --tmp_dir        Manually specify a TMP directory for pybedtools output
  --help/--h       Display this help menu
  """

    System.exit(0)
}

if (help1 == "help") {
    printHelp()
} else if(help2 =="help"){
    printHelp()
}

workflow{
    // input_data = fetchData()
    //input_query_fasta = input_data.query_fasta
    ("${params.fasta}" != "" ? getQueryAssemblies("${params.fasta}") : Channel.empty()).set{query_fasta}
    ("${params.reference_fasta}" != "" ? getReferenceAssemblies("${params.reference_fasta}") : Channel.empty()).set{reference_fasta}
}


workflow getQueryAssemblies{

    take:
    fasta_loc

    emit:
    fasta_data
    
    main:

    if(fasta_loc == ""){
        error "No assembly data provided via --fasta/--ref_fasta"
    } else{

        fasta_dir = file(fasta_loc)

        // If --fasta is a directory...
        if(fasta_dir.isDirectory()){
            ch_fasta = Channel.fromPath(["${fasta_dir}/*.fa","${fasta_dir}/*.fasta","${fasta_dir}/*.fna"])
        } 
        // If --fasta is a file...
        else if(fasta_dir.isFile()){
            
            // Check if it is a single fasta file...
            if(fasta_dir.getExtension() == "fa" || fasta_dir.getExtension() == "fna" || fasta_dir.getExtension() == "fasta"){
                ch_fasta = Channel.from(fasta_dir).map{it-> file(it)}
            } 
            // Otherwise, assume a file with paths to FASTAs
            else{
                ch_fasta = Channel.from(fasta_dir.readLines()).filter{ file -> file =~ /\.(fa|fasta|fna)$/}.map{it-> file(it)}
            }
        } else{
            error "$fasta_dir is not a valid directory or file..."
        }

        fasta_data = ch_fasta
        .map { filePath ->
            if (!file(filePath).exists()) { error "$filePath is not a valid directory or file..." }
            return filePath }
        .map { filePath ->
            def fileName = file(filePath).getBaseName()
            tuple(fileName, filePath)
        }
    }
}

workflow getReferenceAssemblies{

    take:
    fasta_loc

    emit:
    fasta_data
    
    main:

    if(fasta_loc == ""){
        error "No assembly data provided via --fasta/--ref_fasta"
    } else{

        fasta_dir = file(fasta_loc)

        // If --fasta is a directory...
        if(fasta_dir.isDirectory()){
            ch_fasta = Channel.fromPath(["${fasta_dir}/*.fa","${fasta_dir}/*.fasta","${fasta_dir}/*.fna"])
        } 
        // If --fasta is a file...
        else if(fasta_dir.isFile()){
            
            // Check if it is a single fasta file...
            if(fasta_dir.getExtension() == "fa" || fasta_dir.getExtension() == "fna" || fasta_dir.getExtension() == "fasta"){
                ch_fasta = Channel.from(fasta_dir).map{it-> file(it)}
            } 
            // Otherwise, assume a file with paths to FASTAs
            else{
                ch_fasta = Channel.from(fasta_dir.readLines()).filter{ file -> file =~ /\.(fa|fasta|fna)$/}.map{it-> file(it)}
            }
        } else{
            error "$fasta_dir is not a valid directory or file..."
        }

        fasta_data = ch_fasta
        .map { filePath ->
            if (!file(filePath).exists()) { error "$filePath is not a valid directory or file..." }
            return filePath }
        .map { filePath ->
            def fileName = file(filePath).getBaseName()
            tuple(fileName, filePath)
        }
    }
}

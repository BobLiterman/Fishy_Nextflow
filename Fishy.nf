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
    ("${params.ab1}" != "" ? processAB1Files("${params.ab1}") : Channel.empty()).set{query_ab1}
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

workflow processAB1Files{

    take:
    ab1_loc

    emit:
    ab1_files
    
    main:

    if(ab1_loc == ""){
        error "No data provided to --ab1"
    } else{
        ab1_dir = file(ab1_loc)

        // If --ab1 is a single directory, get all ab1 from that directory
        if(ab1_dir.isDirectory()){
            ab1_files = Channel.fromPath(["${ab1_dir}/*.ab1"])
        }         
        // If --ab1 is a file including paths to files, process all ab1s
        else if(ab1_dir.isFile()){
            if(ab1_dir.getExtension() == "ab1"){
                ab1_files = Channel.from(ab1_dir).map{it-> file(it)}
            } else{
                ab1_files = Channel.from(ab1_dir.readLines()).filter{ file -> file =~ /\.(ab1)$/}.map{it-> file(it)}
            }
        }
        // Error if --ab1 doesn't point to a valid file or directory
        else{
            error "$ab1_dir is neither a valid file or directory..."
        }

        view(ab1_files)
        ab1_data = ab1_files.map { filePath ->
            if (!file(filePath).exists()) { error "$filePath is not a valid directory or file..." }
            return filePath 
            }      


        paired_ab1_data = ab1_data |
        branch{it ->
            forward: "${it}".endsWith("${params.ab1_forward}")
                def fileName = file("${it}").getName()
                def sampleName = fileName.replaceAll("${params.ab1_forward}", "")
                return(tuple([sampleName,"${it}"]))
            reverse: "${it}".endsWith("${params.ab1_reverse}")
                def fileName = file("${it}").getName()
                def sampleName = fileName.replaceAll("${params.ab1_reverse}", "")
                return(tuple([sampleName,"${it}"]))
            single: true
                def fileName = file("${it}").getName()
                return(tuple(["${it}".replaceAll(".ab1", ""),"${it}"]))}

        paired_ab1_data.forward.subscribe{println("Forward: ${it}")}
        paired_ab1_data.reverse.subscribe{println("Reverse: ${it}")}
        paired_ab1_data.single.subscribe{println("Single: ${it}")}
    }
}

process fetchAB1{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val dir // Directory containing read files
    val read_ext // Extention for read files (e.g., fastq.gz or fq)
    val forward_suffix // Identifier for forward reads (e.g., _1.fastq or _R1_001.fq.gz)
    val reverse_suffix // Identifier for reverse reads (e.g., _2.fastq or _R2_001.fq.gz)

    output:
    stdout

    script:

    if(!file(dir).isDirectory()){
        error "$dir is not a valid directory..."
    } else{
    """
    $params.load_python_module
    python ${findReads} --read_dir ${dir} --read_filetype ${read_ext} --forward_suffix ${forward_suffix} --reverse_suffix ${reverse_suffix} --trim_name ${params.trim_name}
    """
    }
}

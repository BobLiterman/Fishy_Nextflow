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
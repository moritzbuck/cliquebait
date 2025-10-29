#!/usr/bin/env nextflow

include { download } from "./processes.nf"
include { fastani } from "./processes.nf"
include { cliquebait } from "./processes.nf"

workflow {
   

    taxa = Channel.fromPath(params.taxa)
    taxa 
	| splitCsv
	| map { x -> x[0] }
    | branch { taxon ->
        done : file("${params.genomes_out}/${taxon}", checkIfExists: false).exists()
        todo : true
        } 
    | set { genome_sets }

    genome_sets.todo
         | download
         | set { downloaded_genomes }
    
    genome_sets.done.map { taxon -> file("${params.genomes_out}/${taxon}") }
        | concat(downloaded_genomes)
        | branch { taxon ->
            done : file("${taxon}_fastani.tsv".replace("${params.genomes_out}", "${params.fastani_out}"), checkIfExists: false).exists()
            todo : true
            } 
        | set{ fastani_sets}

//    fastani_sets.todo
//        | fastani
//        | set { fastanied }

    fastani_sets.done.map { taxon -> file("${taxon}_fastani.tsv".replace("genomes/", "fastanis/")) }
//        | concat(fastanied)
        | cliquebait
    

}


// End of file
            

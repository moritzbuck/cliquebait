#!/usr/bin/env nextflow


include { download } from "./processes.nf"
include { fastani } from "./processes.nf"
include { cliquebait } from "./processes.nf"
include { motupan } from "./processes.nf"
include { genecall  } from "./processes.nf"
include { clusterAAs } from "./processes.nf"

include { parseCliqueBaitJson } from "./utils.nf"

workflow {
   

    taxa = channel.fromPath(params.taxa)
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
        | map { json_file -> parseCliqueBaitJson(json_file) }
        | flatMap
        | filter{ t -> t[1].size() > params.min_clique_size  }
        | set { kept_motus }

    kept_motus
        | flatMap{t -> t[1].collect{ tt -> Tuple.tuple(t[0], tt) } }
        | genecall
        | groupTuple
        | motupan

//    kept_motus
//        | motupan
}


// End of file
            

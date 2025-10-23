process download {
    cpus  1
    publishDir "${params.genomes_out}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 20
    input :
    val(taxon)
    output:
    file("${taxon}")
    script:
    """

    grep "${taxon};" ${params.gtdb_md} > gtdb_md.tsv
    cut -f1 -d\$'\\t' gtdb_md.tsv  | cut -f2- -d"_" > to_dls
    datasets  download genome accession --inputfile to_dls
    unzip ncbi_dataset.zip
    rm ncbi_dataset.zip
    mkdir ${taxon}
    find ncbi_dataset -name "*.fna" -exec mv {} ${taxon}/ \\; 
   """
}

process fastani{
    cpus  20
    publishDir "${params.fastani_out}", mode: 'copy'
    conda '0079'
    input:  
    file(taxon)
    output:
    file("${taxon}_fastani.tsv")
    script:
    """
    find -L ${taxon}  -name "*_genomic.fna" > tmp
    paste  <(cut -f9 -d"/"  tmp)  tmp > genome_list
#    fastANI -t ${task.cpus} --ql tmp  --rl tmp  -o ${taxon}_fastani.tsv    
    mOTUlize.py --threads ${task.cpus} --keep-simi-file ${taxon}_fastani.tsv -T --fna genome_list --output /dev/null --force
    """

}

// End of file
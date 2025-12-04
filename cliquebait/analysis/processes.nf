process download {
    cpus  1
    publishDir "${params.genomes_out}", mode: 'copy'
    errorStrategy 'ignore'
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

process cliquebait{
    publishDir "${params.cliquebait_out}", mode: 'copy'
    input:  
    file(fastani_file)
    output:
    file("*.json")
    script:
    """
    cwd=`pwd`
    sed 's#.*/\\(GC[AF]_[0-9.]*\\).*\\(GC[AF]_[0-9.]*\\).*_genomic.fna\\(.*\\)#\\1\\t\\2\\t\\3# ' ${fastani_file} > fixed_fastani.tsv
    cd ${params.cliquebait_path}
    outfile=`echo ${fastani_file} | sed 's/_fastani.tsv/_baits.json/'`
    python -m cliquebait.bin.clique_bait --similarities \${cwd}/fixed_fastani.tsv -o \${cwd}/\${outfile} --min_size ${params.min_clique_size} --gap_size  ${params.gap_size}  --checkm "${params.checkm_path};${params.checkm_completeness};${params.checkm_contamination}"
    cd \${cwd}/
    """

}


process genecall{
    cpus 1
    publishDir "${params.genecall_out}", mode: 'copy'
    input:  
    val(genome)    
    output:
    tuple val("${genome[0]}"), file("${genome[1]}.faa")
    script:
    """
    fna=`find ${params.genomes_out} -name "${genome[1]}*_genomic.fna"`
    if [ -f ${params.genecall_out}/${genome[1]}.faa  ];
    then
        cp ${params.genecall_out}/${genome[1]}.faa .
        cp ${params.genecall_out}/${genome[1]}.gbk .
    else 
       prodigal -a ${genome[1]}.faa -c  -i \${fna} -o ${genome[1]}.gbk 
    fi
    """
}

process clusterAAs{
    cpus 20
    publishDir "${params.base_out}", mode: 'copy'
    input:  
    file(faa_files)
    output:
    file("mmsseqs_cluster.tsv")
    script:
    """
    for f in `ls | grep .faa`
    do
        sed 's/>.*_\\([0-9]*\\) #.*/>'\${f}'_\\1/ ; s/.faa//' \${f} >> all_genes.faa
    done
    mmseqs easy-cluster --min-seq-id ${params.seqid} --cov-mode ${params.covmode} -c ${params.cov} --threads ${task.cpus} all_genes.faa  mmsseqs mmseqs_temp    
    """
}


process motupan{
    publishDir "${params.motupan_out}", mode: 'copy'
    errorStrategy 'ignore'
    cpus 8
    input:  
    val(input)
    output:
    tuple(file("${input[0]}_motupan.csv"), file("${input[0]}_gcs.json") )
    script:
    """
    mOTUpan.py --genome2gene_clusters_only  --threads ${task.cpus} --faas ${input[1].join(" ")} --output ${input[0]}_gcs.json --name ${input[0]} 
    mOTUpan.py --gene_clusters_file  ${input[0]}_gcs.json --checkm ${params.checkm_path}  --threads ${task.cpus} --output ${input[0]}_motupan.csv --name ${input[0]} 
    """

}


// End of file
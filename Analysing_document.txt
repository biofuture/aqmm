# aqmm
AQMM analyzing document

Many analyses were performed with MetaP (https://github.com/biofuture/MetaP), a package for both 16S amplicion, metagenomics and metatranscriptmics analyses

1. quality control of DNA and RNA sequences with trimmomatic using MetaP script .   
    perl illumina_pe_qc.pl <fq.list> <output dir> <number of threads>

2. Using SortMerRNA to remove the rRNA from metatranscriptome data .   

    time ./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db  --reads /MDHome/FOAMINGMETATRAN/FCHCF5JBBXX-WHENVfbwEAAARAAPEI-201_L2_1.fq --sam --num_alignments 1 --fastx --aligned /MDHome/FOAMINGMETATRAN/FCHCF5JBBXX-WHENVfbwEAAARAAPEI-201_L2_1.rRNA --other /MDHome/FOAMINGMETATRAN/FCHCF5JBBXX-WHENVfbwEAAARAAPEI-201_L2_1.non_rRNA --log -v

    perform this step for every pair-end metatranscriptomics data set.  

3. Using AQMM to generate the absolute quantification index for different samples   
    perl ../GitProject/AQMM/aqmm/dist/aqmm -a <Input_metagenome_dir> -b <Input_metatranscriptome_dir> -m <Experimental_meta_data> -n <num_threads> -o <ouput_dir> .   
    The Input_metagenome_dir contains all the metagenomic data after quality control for all samples .   
    The Input_metatranscriptome_dir contains all the metatranscriptomic data after quality control and removing rRNA for all samples .        
    The Experimental_meta_dta contains all the meta data information for the extraction of DNA and RNA .   
    
4. Co-assembly of metagenome data sets and generate the contigs with CLC Bio .   
   Using CLC Bio to co-assemble all the metagneome data set and generate the contig set "Contigs_all.mixassemble.fa" .   
    
5. gene prediction from the contigs and obtain the catalog of genes .   
    using prodigal to predict genes from the assembled contigs "Contigs_all.mixassemble.fa" to get the gene file "all_nucleotide_nu.fna", and then unique the gene set to get the unique gene set of this gene pool with usearch.    
    
    prodigal -i Contigs_all.mixassemble.fa -o allout.prodigal.txt -a all.pro.aa.fna -d all.gene.nu.fna -f gbk -p meta .      
    usearch -fastx_uniques all.gene.nu.fna -fastaout uniques.fasta -sizeout -relabel Uniq .   
    
6. mapping metagenome and metatranscriptome to the gene catalog with CLC Bio .   
    Using CLC Bio to map metagenomic and metatranscriptomic reads back to the unique contig/gene set to obtain the coverage of each contig and gene .   
    
7. taxonomy assignment of genes and contigs using MEGAN6 .    
    Using diamond to align all the unique gene to NR database (version 201608) and then using MEGAN6 to process the alignment result to get the taxonomy assignment of each hit gene.        
    
8. ARGs identification with ARGs-OAP .   
   Alignment all the gene to SARGG database to identify the type and subtype of ARGs .   
    
9. ARGs matrix and genus/phylum matrix of the samples .       
   Using the absolute quantification of samples obtained by AQMM, coverage information of genes, ARGs and taxonomy annotation information of genes, the absolute quantificaiton of ARGs and genus were obtained and merge into matrix for further statistical analysis.     
    
10. differential expression genes/genes identification with the absolute quantification .   
   For all the absolute quantificaiton matrix, using T-test to test the significance between foaming and non-foaming .  

11. visulization with R .  
   visulize the top ARGs/genus with boxplot using R .  

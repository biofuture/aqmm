# aqmm

Document for AQMM (Absolute quantification of metagenome and metatranscriptome) .  

The AQMM was developed to perform absolute quantification of mulitipile metagenome and it's paralell metatranscriptome. In order to use this algorithm, the experiment should be designed with both metagenome and metatranscriptome data. From initial stage, the molecular experimental data should be recorded to help on estimating the overall DNA or RNA of a unit (ml/gram) of sample. The AQMM was demonstrated to obtain better results of differential experssion genes identification in comparative metatranscriptomic studies.       

## clone the source code to local server .   

git clone https://github.com/biofuture/aqmm.git .   

## prepare the meta_data.txt which contains all the molecular experimental data .   

The meta_data.txt should .

SID	|  SNAME  | DNA_weight    |  DNA_volume   |   DNA_extr_eff  |  RNA_weight  |    RNA_volume   |   RNA_extr_eff |   mRNA_ratio 
--------|---------|---------------|---------------|-----------------|--------------|-----------------|----------------|-
1       | AS1     | 24570         | 0.5           | 0.282           | 6135         |     0.5         |    0.075       |  0.03 
2       | AS2     |27300          |0.5            | 0.282           |  6355        |    0.5          |    0.075       |  0.03   

Each line contains information of one sample .

SNAME is the name of the sample, which stored under the input directory .   
DNA_weight is the weight of DNA for the extraction (ng) .   
DNA_volume is the total volume of sample used for the DNA extraction .   
DNA_extr_eff is the extracting efficiency for the sample, which is an emperical value for a DNA extraction KIT for a sample .   
RNA_weight is the weight of RNA for the extraction .   
RNA_volume is the total volume of sample used for the RNA extraction .   
RNA_extr_eff is the extracting efficiency for the sample, which is an emperical value for a RNA extraction KIT for a type of sample .   
mRNA ratio is the ratio of mRNA to total RNA, which is a emperical value for different environment          

Notice: Before using this algorithm, users need to get the above meta data information as accurately as possible, although the algorithm show robusty to the variance to the extraction efficiency among different samples, users should keep all the condition as consistent as possible in order to get more accurate/reliable quantification and comparativety. 

## Running the aqmm just by one command .  

By typing ./aqmm -h the help information will show here .    
    
	perl ./aqmm -a <Input_metagenome_dir> -b <Input_metatranscriptome_dir> -m <Experimental_meta_data> -n <num_threads> -o <ouput_dir>
	Author Jiang Xiao-Tao 2017-08
	Email  biofuture.jiang@gmail.com

where:    

       -a <Input_metagenome_dir> including all the fq files of your metagenomic samples after quality filtering
       
       -b <Input_metatranscriptome_dir> including all the paralell metatranscriptomic fq files with metagenomic samples after . remvoing the rRNA with tools like SortMeRNA.
       
       -m meta_data.txt including all the experimental extraction information and the extraction efficiency information of your sample type
       
       -o <ouput_dir> will include all the output results for the quantification, including the esimated cell numbers of metagenomic data of each sample and the estimated sequenced cell numbers your metatranscriptomic data included. With the quantification results obtained, the abundance of gene/species could be futher normalized to per cell/volume/gram etc. al.   

## Example to run

There is a small demo to show how to use AQMM to generate absolute quantification of your data under the example directory. 

Enter the example directory and run the following command 

	../dist/aqmm  -a AS_metagenome -b AS_metatranscriptome -m meta_data.txt -n 12 -o testout

After runnning, check the results 

	cat dna_cell_info.txt .   

SID	| SNAME	| library_size	| lib_cell	| cell_per_ml	| A_ratio	| T_ratio	| C_ratio	| G_ratio    
--------|-------|---------------|---------------|---------------|---------------|---------------|---------------|----------
1	|AS2	|50000000	|8.56742513314301	| 3.23275242169008e+19	| 0.22371768	| 0.22220922	| 0.27649788	| 0.27757522   
2	|AS1	|50000000	|9.93418747046929	| 3.37319232480681e+19	| 0.21272236	| 0.21128754	| 0.28672574	| 0.28926436   

	cat rna_cell_info.txt .

SID	|SNAME	|library_size	|lib_cell	|cell_per_ml|	A_ratio	|T_ratio|	C_ratio	|G_ratio    
--------|-------|---------------|---------------|-----------|-----------|-------|---------------|-------
1	|AS2	|50000000	|8.56742513314301 |	3.23275242169008e+19|	0.22371768|	0.22220922	|0.27649788	|0.27757522    
2	|AS1	|50000000	|9.93418747046929 |	3.37319232480681e+19|	0.21272236|	0.21128754	|0.28672574	|0.28926436    

These numbers could be used to normalize your genes/species to per cell/volume level in the end!       


### Compare with RQ methods

A script was developed to process the results for both AQMM and RQ methods like RPKM/TPM/edgeR to identify the differential expression genes between groups

	perl normalization_RNA_sequencing.pl <RNA_Depth> <gene_length> <RNA_Table> <Oprefix>

where 

	RNA_Depth	is the meta_data file contains the the library size of each sample 
	gene_length	is the file contains all the gene length information
	RNA_Table	is the coverage matrix of each gene in different samples generated by mapping process
	Oprefix	        is the output prefix for all the normalization methods for the RNA_Table 


I will show the our foaming vs nonfoaming activated sludge as the example

cat RNA_DEPTH.txt

	NAME	lib.size	Group
	RNA-201.gene.txt	148949591	Foaming
	RNA-202.gene.txt	181507132	Foaming
	RNA-203.gene.txt	135988557	Foaming
	RNA-205.gene.txt	152020602	NonFoaming
	RNA-206.gene.txt	155510444	Foaming
	RNA-207.gene.txt	161566079	NonFoaming
	RNA-208.gene.txt	148860009	NonFoaming
	RNA-209.gene.txt	132982810	Foaming
	RNA-210.gene.txt	137007411	NonFoaming

head -10 len_gene.txt

	contig_1_1	1071
	contig_1_2	552
	contig_1_3	1284
	contig_1_4	279
	contig_1_5	699
	contig_1_6	867
	

### Supporting for time series metatranscriptomics studies 

As many metatranscriptomics studies only investiagte the activity of a system without change of the DNA part, hence this is a very important application scence. To perform absolute quantification of this condition, a optional parameters were integrated into the aqmm tool.  Users process this type of data could quantify the transcript in this way.  Inputing the relative abundance informaition, per species quantification could be achieved. 



Copyright: LG209, Environmental biotechnology laborotory HKU.    

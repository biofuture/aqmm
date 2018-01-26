
#!/usr/bin/env perl

use strict;

#This program is to normalize integrated metagenome and metatranscriptome (IMM) generated DNA table and RNA expression mother table to get differential expressed genes (DEGs) in control and experimental groups. 
#This script integrate differnt methods to perform normalizaton including RPKM, FPKM, TPM, TMM and also negative bionormal distribution based normalization methods such as DESeq2 and edgeR.   
##Before running these methods the R package DESeq2 and edgeR should be installed on your unix/linux platform 
##edgrR DEseq2 should be install in your local R in dadvance to use this script
##normalize RNA script with DNA samples gene relative abundacne (DNA abundance were normalizaed with the gene coverage against averagely coverage of universal signle copy marker genes in each sample, so that the abundance are comparable to each other), as different Groups of samples should have many species with high different number of cells, the relative abundance of gene in DNA samples are divided to remove this factor, in order to get a comparable quantification among groups, DNA samples are fistly normalize to the coverage of USCMGs.
##Genes in one contigs should have the same DNA abundance in one sample, so we use the aveage abundance of one contigs to all the genes in that contig, we using the RPKM to divide the gene relative abundance in paralelly DNA sample

##Parse input  parameters 
##Prelimilary requirements for this script to run correctly
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);

our (@dirset,$aqmmdir);
BEGIN {
    @dirset = split(/\//,$Bin);
    $aqmmdir = join("/", @dirset);
    unshift @INC, "$aqmmdir/bin";
}

our ($opt_h, $opt_a, $opt_b, $opt_c, $opt_d, $opt_o, $opt_l, $opt_h, $opt_s) = "";
my  $usage = <<USE;
        Author: Xiao-Tao JIANG
        Date: 
        Modified : 19-01-2019
        Email: biofuture.jiang\@gmail.com
	 "perl $0 -a <meta_data_rna.txt> -b <unique_gene_set> -c <rna_reads_count.matrix> -o <Oprefix> -h \n" 
        
	-a input meta data file for all the metatranscriptomic data <required>
	-b the gene set fasta file <required>
	-c matrix of reads number mapped in RNA samples for each gene  <required>
	-d matrix of transcript per copy gene  run AQMM required.
	-o output prefix, all output result files name are prefix with this string <required>  
	-l the length of reads default strategies PE 150 sequencing, default 150 bps 
	-s this option is to select the normalization methods, RPKM, TPM, EDGER, AQMM, ALL. default AQMM 
	-h print this help information 
USE

#Get all the input parameters
getopts('a:b:c:d:o:l:s:h');
if($opt_h  ||  ($opt_a eq "") ){
	die "$usage\n";
}

#// represent unit test code 
##Variables -------------------------------------------------------------
my %len;  #store the lenght of each gene in the unique gene set  global varaibles 
my %scalingfac; ##store the RNA sample sequencing depth
my $oprefix = $opt_o;
my %scalingsmapped; ##store the mapped number of reads for each sample 
my %samplegroup; #store groups information for samples 

#default value 
$opt_l ||= 150;
$opt_s ||= "ALL";

#Defined sub funtions 
# s_total_mapped len_gene parse_meta  dge_iden rpkm tpm edger percellnorm perunitnorm

##------------------------------------------------------------------------------------
#Main flow starting 
#1. parse the gene set to get lens of each genes 

	len_gene($opt_b, \%len);
#	print "len_gene finish\n";
#//for my $test (sort {$len{$b} <=> $len{$a}} keys %len){
#//	if($len{$test} > 5000){
#//		print "$test\t$len{$test}\n";
#//	}
#//}
#//die;

#2. parse the meta data rna file to input the meta_data information into %hash

	parse_meta($opt_a, 0, 1, \%scalingfac); #hash address as parametes
	parse_meta($opt_a, 0, 4, \%samplegroup);
						#0 name 1 depth of the RNA library 
	#print "parse meta finish\n";
#//for my $key (keys %scalingfac){
#//	print "$key\t$scalingfac{$key}\n";
#//}
#//die;
#3 normalization the RNA matrix with different methods 
	if($opt_s eq "RPKM"){		

		s_total_mapped($opt_c, \%scalingsmapped); ##store total number of reads mapped to each sample for this gene set 
		##define the output RPKM matrix file, and use rpkm function to normalize
		my $orpkm = "$oprefix.RPKM.txt";
  		rpkm($opt_c, \%scalingsmapped,\%len, $orpkm);
		my $odgerpkm = "$oprefix.odge.rpkm.txt";
		dge_iden($orpkm, \%samplegroup, $odgerpkm);	
		##differnetial gene expression identification using RPKM matrix with statistical test 

	}elsif($opt_s eq "TPM"){
		my $otpm = "$oprefix.TPM.txt";
		tpm($opt_c, \%len, $otpm);	
		my $odgetpm = "$oprefix.odge.tpm.txt";
		dge_iden($otpm, \%samplegroup, $odgetpm);

	}elsif($opt_s eq "DESEQ2"){

	}elsif($opt_s eq "EDGER"){
		
		##do all the edgeR analyses on the reads count matrix to identity differential expression geens 
		edger($opt_c, $opt_a);
				
	}elsif($opt_s eq "AQMM"){

		#easily normalize the mapped genes coverage with total cell number 	
		my %seqcells; #store the sequenced RNA cells for each sample
		my %perunitcells; #store the per unit RNA cells for each sample
		die "Without AQMM result: \n" unless (-e $opt_a);
		parse_meta($opt_a, 0, 2, \%seqcells);	#initial %seqcells	
		parse_meta($opt_a, 0, 3, \%perunitcells); #initial %perunitcells	
		
		#norlaize mRNA abundance to global absolute scaling factors 
		#$opt_l is also needed to calculate the avearge coveage  
		my $opercell = "$oprefix.global.percell.txt";
		my $operunit = "$oprefix.global.perunit.txt";
		percell_norm($opt_c, $opt_l, \%seqcells, $opercell); #normalize the mRNA to total cell numbers in the library
		perunit_norm($opt_c, $opt_l, \%perunitcells, \%seqcells, $operunit);	#normalize the mRNA to absolute unit like ml/gram 
		
		##differential gene expression idenntification using percell or perunit 
		my $odgepercell = "$oprefix.dge.percell.txt";
		my $odgeperunit = "$oprefix.dge.perunit.txt";
	
		##per copy transcript indicate that for each copy of this gene, how many transcript are there for it. As the RNA is sequenced much deeper than the DNA, sometimes we detected DNA zero abundance, however with RNA expression, these RNA acturally could help to get better assemble results for the DNA data

		dge_iden($opercell, \%samplegroup, $odgepercell);		
		dge_iden($operunit, \%samplegroup, $odgeperunit);

	}elsif($opt_s eq "TPCG"){
		
		my $odgepertpcg = "$oprefix.dge.tpcg.txt";
		dge_iden($opt_d, \%samplegroup, $odgepertpcg);
	}elsif($opt_s eq "ALL"){

		## perform all the analyses 
		s_total_mapped($opt_c, \%scalingsmapped); ##store total number of reads mapped to each sample for this gene set 
		##define the output RPKM matrix file, and use rpkm function to normalize
		my $orpkm = "$oprefix.RPKM.txt";
  		rpkm($opt_c, \%scalingsmapped,\%len, $orpkm);
		my $odgerpkm = "$oprefix.odge.rpkm.txt";
		dge_iden($orpkm, \%samplegroup, $odgerpkm);	
		##differnetial gene expression identification using RPKM matrix with statistical test 
		my $otpm = "$oprefix.TPM.txt";
		tpm($opt_c, \%len, $otpm);	
		my $odgetpm = "$oprefix.odge.tpm.txt";
		dge_iden($otpm, \%samplegroup, $odgetpm);
		##do all the edgeR analyses on the reads count matrix to identity differential expression geens 
		edger($opt_c, $opt_a);
		
		#easily normalize the mapped genes coverage with total cell number 	
		my %seqcells; #store the sequenced RNA cells for each sample
		my %perunitcells; #store the per unit RNA cells for each sample
		die "Without AQMM result: \n" unless (-e $opt_a);
		parse_meta($opt_a, 0, 2, \%seqcells);	#initial %seqcells	
		parse_meta($opt_a, 0, 3, \%perunitcells); #initial %perunitcells	
		
		#norlaize mRNA abundance to global absolute scaling factors 
		#$opt_l is also needed to calculate the avearge coveage  
		my $opercell = "$oprefix.global.percell.txt";
		my $operunit = "$oprefix.global.perunit.txt";
		percell_norm($opt_c, $opt_l, \%seqcells, $opercell); #normalize the mRNA to total cell numbers in the library
		perunit_norm($opt_c, $opt_l, \%perunitcells, \%seqcells, $operunit);	#normalize the mRNA to absolute unit like ml/gram 
		
		##differential gene expression idenntification using percell or perunit 
		my $odgepercell = "$oprefix.dge.percell.txt";
		my $odgeperunit = "$oprefix.dge.perunit.txt";
	
		my $odgepertpcg = "$oprefix.dge.tpcg.txt";
		dge_iden($opercell, \%samplegroup, $odgepercell);		
		dge_iden($operunit, \%samplegroup, $odgeperunit);
		dge_iden($opt_d, \%samplegroup, $odgepertpcg);
	}else{
		die "Module selection error:\n";
	}

##--------------------------------------END of Main Flow Function-------------------###

##sub function
##---------------------------------------------------------------------------------------
sub len_gene {
	##Inport the gene length of your gene expressed and detected in the metatranscriptomics data
	my ($fapath, $hashlen) = @_; ##This first parameter is the path of fa, the second is the address of hash to store the gene length info
	die "$!\n" unless open(GL, "$fapath"); ##Parse the fasta gene file to get length info
	my $id = '';
	while(<GL>){
		chomp;
		if($_ =~ /^>(\S+)/){
			$id = $1;
		}else{
			#${$hashlen}{$id} .= $_;
			${$hashlen}{$id} += length($_);
		}
	}
	close GL;
} ##len_gene

sub parse_meta {
	
	my ($mfile, $in1, $in2, $hashref) = @_;

	die "$! $mfile\n" unless open(RNADEPTH, "$mfile");
	<RNADEPTH>;
	while(<RNADEPTH>){
		chomp;
		my @tem = split /\t/;
		${$hashref}{$tem[$in1]} = $tem[$in2]; ##keep the total reads of sample 
	}
	close RNADEPTH;
}##parse meta

sub s_total_mapped {

	my ($rnamatrix, $hashref) = @_;
	
	die "$!\n" unless open(RNAREAD, "$rnamatrix");
	my $head = <RNAREAD>;
	chomp($head);
	my @namearray = split("\t", $head);
	#initial %{$hashref} with zero
	for(my $i = 1; $i<= $#namearray; $i++){
		${$hashref}{$namearray[$i]} = 0;
	}

	while(<RNAREAD>){
		chomp;
		my @tmp = split /\t/;
		for(my $i=1; $i <= $#tmp; $i++){
			${$hashref}{$namearray[$i]} += $tmp[$i]; 
		}	
	}
	close RNAREAD;
}##statistic total mapped reads number for each sample 

##-------------------------------------------RPKM only on RNA table-----------------------------------##
##reads per kilo bps per million mapped reads
##The calculation of RPKM does not consider the abundance of genes in the DNA sample
##It is quite easy to calculate this index by supplying the length of gene 
#These three metrics attempt to normalize for sequencing depth and gene length. Here show you do it for RPKM:
#1. Count up the total mapped reads in a sample and divide that number by 1,000,000 this is our per million scaling factor.
#2. Divide the read counts by the per million scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
#3. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.

sub rpkm {
	
	my ($rmatrix, $hashref, $lenhash, $rpkmo) = @_;
	#1. RNA reads mapping matrix
	#2. hash ref for the each smaple mapped reads 
	#3. output rpkm matrix for this RNA reads count matrix 
	die "$!\n" unless open(RNA, "$rmatrix");
	#my $rpkmo = "$oprefix.RPKM.txt";
	die "$!\n" unless open(ORPKM, ">$rpkmo");
	my $headrna = <RNA>;chomp($headrna);
	my @names = split("\t", $headrna);
#print "$names[1]\n";
	print ORPKM "$headrna\n";
	while(<RNA>){
		chomp;
		my @line = split /\t/;
##calculate RPKM
		print ORPKM "$line[0]";
		for(my $i =1; $i <= $#line; $i++){
			die "$i\t$names[$i]" unless(exists ${$hashref}{$names[$i]});
			die "$i\t$names[$i]" unless(exists ${$lenhash}{$line[0]});
			my $rpkmnorm = $line[$i] * 1000 * 1000000 / (${$lenhash}{$line[0]} * ${$hashref}{$names[$i]});
			print ORPKM "\t$rpkmnorm";
		}
		print ORPKM "\n";
	}
	close RNA;
	close ORPKM;
}#rpkm
##----------------------------------------------------close RPKM-------------------------------------##

##--------------------------------------------------------TPM----------------------------------------##
#TPM is very similar to RPKM, The only difference is the order of operations. Here show you calculate TPM:
#1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your per million scaling factor.
#3. Divide the RPK values by the per millions scaling factor. This gives you TPM.

sub tpm {

	my ($infile, $lenhref, $otpm) = @_;
	die "$!\n" unless open(RNA, "$infile");
	die "$!\n" unless open(OTPM, ">$otpm");

	my $headrna = <RNA>;chomp($headrna);
	my @names = split("\t", $headrna);
	my %scalerpk; ##store the total rpk for the sample 
#initial scalrpk with zero
		for(my $i = 1; $i <=$#names; $i++){
			$scalerpk{$names[$i]} = 0;
		}

	print OTPM "$headrna\n";
	while(<RNA>){
		chomp;
		my @line = split /\t/;
		if(exists ${$lenhref}{$line[0]}){
			for(my $i =1; $i <= $#line; $i++){
				my $rpk = $line[$i] * 1000  / (${$lenhref}{$line[0]});
				$scalerpk{$names[$i]} += $rpk; ##get the total rpk for one sample 
			}
		}
	}
	close RNA;

##output each TPM
##while calculating the TPM, the total number of reads is not used, the sumup reads number of the sample is used 
	die "$!\n" unless open(RNA, "$infile");
	<RNA>;
	while(<RNA>){
		chomp;
		my @line = split /\t/;
		if(exists ${$lenhref}{$line[0]}){
##calculate RPKM
			print OTPM "$line[0]";
			for(my $i =1; $i <= $#line; $i++){
				my $rpk = $line[$i] * 1000  / (${$lenhref}{$line[0]});

				my $tpmv =   $rpk / ( $scalerpk{$names[$i]} / 1000000);  ##get the total rpk for one sample 
					print OTPM "\t$tpmv";
			}

			print OTPM "\n";
		}
	}
	close RNA;
	close OTPM;
}#tpm function 

##------------------------------------------------------END TPM------------------------------------##

##Using edgeR to normalization of RNA_seq data to detect differential expression genes in system
##Generate R script for edgeR to do the Differnt expression genes analysis
sub edger {
	
my ($readsmatrix, $metadata) = @_;
my $rfile="$oprefix.edgeR.R";
die "$!" unless open(RS, ">$rfile");
 
my $erscript =<<EDGER;
library(edgeR)
table <- read.delim("$readsmatrix", sep="\\t", header=TRUE,row.name=1)
fac <- read.table("$metadata", sep="\\t",header=TRUE)
y<- DGEList(counts=table, group=fac\$Group, lib.size=fac\$lib.size) 
##filtering low abundant genes 
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y<- calcNormFactors(y) ##using trimmed mean of M values (TMM) to adjust lib.size 
pdf("$oprefix.edgeR.MDS.pdf")
plotMDS(y)
dev.off()
design <- model.matrix(~fac\$Group)
y <- estimateDisp(y,design) #estimate the dispersion 
##generate BCV plot to show the biological coefficient variance of genes  (BCV biological coefficient vaiation)
pdf("$oprefix.edgeR.BCV.pdf")
plotBCV(y)
dev.off()
##generate differential expression gene plots 
fit <- glmFit(y,design) ##Generalize linear model (GLM) likelihood ratio test  
lrt <- glmLRT(fit)	##LogFC, FC means fold change logFC eq 1 means 2 fold change 
de <- decideTestsDGE(lrt)
detags <- rownames(y)[as.logical(de)]
pdf("$oprefix.edgeR.Smear.pdf")
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()
rownames(de) <- rownames(y)
out <- data.frame(lrt\$table, de)
write.csv(file="$oprefix.edgeR.DEGs.list", out)  ##CPM copy per million The output  LogCPM is the log2 of CPM after TMM adjustion of the total lib.size 
write.csv(file="$oprefix.edgeR.GLM.fitted.txt",lrt\$fitted.values)
#qlf <- glmQLFTest(fit)   ##quasi-likelihood QL F-test 
#fit <- glmFit(y,design)

EDGER
print RS $erscript;
close RS;

`R CMD BATCH $rfile`;

}#edger
##------------------------------------END of edgeR-------------------------------#

sub percell_norm {  ##This function generate the tpc matrix 

	my ($rnafile, $readlen, $seqcellref, $outputpercel) = @_;
	
	die "$rnafile\n$!" unless open(RNAF, "$rnafile");
	die "$outputpercel\n" unless open(OFCELL, ">$outputpercel");
	
	my $headrna = <RNAF>;chomp($headrna);
	my @names = split("\t", $headrna);
#print "$names[1]\n";
	print OFCELL "$headrna\n";
	
	while(<RNAF>){
		chomp;
		my @tem = split /\t/;
		print OFCELL "$tem[0]";
		if(exists $len{$tem[0]}){
			for(my $i=1; $i <= $#tem; $i++){
##normalize the reads count to coveage of gene per cell
				die "$tem[0]" if($len{$tem[0]} ==  0);
				my $geneavcov = ($tem[$i] * $readlen) / $len{$tem[0]};
##coverage / cell numbers of this sample  to get the TPC 
				my $tpc = $geneavcov / ${$seqcellref}{$names[$i]};  
				print OFCELL "\t$tpc";
			}
			print OFCELL "\n"; #END of this gene/mRNA/transcript 

		}
	}	
	close OFCELL;
	close RNAF;
}#per cell norm

sub perunit_norm {

	#($opt_c, $opt_l, \%perunitcells, \%seqcells, $operunit);
	my ($rnafile, $readlen, $perunitcells, $seqcellref, $outputperunit) = @_;
	
	die "$rnafile\n$!" unless open(RNAF, "$rnafile");
	die "$outputperunit\n" unless open(OFUNIT, ">$outputperunit");
	
	my $headrna = <RNAF>;chomp($headrna);
	my @names = split("\t", $headrna);
#print "$names[1]\n";
	print OFUNIT "$headrna\n";
	
	while(<RNAF>){
		chomp;
		my @tem = split /\t/;
		print OFUNIT "$tem[0]";
		if(exists $len{$tem[0]}){
			
			if(exists $len{$tem[0]}){
			for(my $i=1; $i <= $#tem; $i++){
				die "$tem[0]\n" unless (exists $len{$tem[0]});
##normalize the reads count to coveage of gene per cell
				my $geneavcov = $tem[$i] * $readlen / $len{$tem[0]};
##coverage / cell numbers of this sample  to get the TPC 
				my $tpunit = ${$perunitcells}{$names[$i]} *  $geneavcov / ${$seqcellref}{$names[$i]};  
				print OFUNIT "\t$tpunit";
			}
			print OFUNIT "\n"; #END of this gene/mRNA/transcript 
			}
		}
	}
	
	close OFUNIT;
	close RNAF;
}#per unit norm

sub dge_iden {
	
	#differential gene expression for a matrix
	my ($matrix, $hashgroup, $odgefile) = @_;
	#wilcox.test(mpg ~ am, data=mtcars) 

	##generate the R script and use the wilxoc.test to identify the DGEs
	##correct the P value with B-H method
	##Kruskal-Wallis test as the non parametric methods for one way ANOVA
	my $temprscript = "$oprefix.temp.dge.R";
	my $head = `head -1 $matrix`; ##put the first line names into a vector
	chomp($head);
	my @headarray = split(/\t/, $head);
	my @grouparray = ();
	
	for(my $i=1; $i <=$#headarray; $i++){
		push @grouparray, ${$hashgroup}{$headarray[$i]};
	}	
	my $garray = join("\t", @grouparray);
	
	die "$temprscript\t$!\n" unless open(RT, ">$temprscript");
my $rscripts =<<RSP;
mat <- read.table(file="$matrix", sep="\\t", header=TRUE, row.names=1)
#divide the matrix into two by groups information
n <- dim(mat)[1] 
m <- dim(mat)[2]
name <- colnames(mat)
groupstring <- \"$garray\"
garray <- strsplit(groupstring, "\\\\s+")
vga <- matrix(unlist(garray), byrow=T)[,1]
#levels(garray)

##vector store all dge 
dgevector <- matrix(, nrow=n)

j <- 1;
for(i in 1:n){

	ratio <- table(mat[i,] != 0)[1]  / m
	if(ratio >= 0.5){		
		#data.frame(t(mat[25,]), vga)
		kt <- kruskal.test(t(mat[i,]) ~ as.factor(vga))
		dgevector[j,1] <- kt\$p.value
		#dgevector[j,2] <- lgfold
		rownames(dgevector) <- rownames(mat)
		j <- j+1
	}
}

write.table(file="$odgefile", dgevector,sep="\\t")

RSP
	print RT $rscripts;
	close RT;

	`R CMD BATCH $temprscript`;
}#dge_iden


1;
__END__

#!/usr/bin/perl -w
use strict;

die "perl $0 <DNA_PerCell> <RNA_PerCell> <OUT_TPCG>
	This script transform the DNA matrix normalized to per cell and RNA matrix normalized to per cell into the TPCG, after running aqmm, the number of sequenced cells in DNA sample and the number of eqivalent expressed cells are obtained. Users can easly normalized the read count matrix of genes with these absolute quantifications \n" unless (@ARGV == 3);

die "$!\n" unless open(D, "$ARGV[0]");
die "$!\n" unless open(R, "$ARGV[1]");
die "$!\n" unless open(O, ">$ARGV[2]");

## DNA G/Cell            
## RNA Exoressuib/Cell  whether the total cell is the same? 

my %dna;
<D>;
while(<D>){
	chomp;
	my @m = split(/\t/, $_, 2);
	$dna{$m[0]} = $m[1];	
}
close D;

my %excep; 
<R>;
while(<R>){
	chomp;
	my @m = split(/\t/, $_, 2);
	my @tm = split(/_/, $m[0]);
	pop @tm;
	my $name = join("_", @tm);
	if(exists $dna{$name}){
		my @dm = split(/\t/, $dna{$name});			
		my @rm = split(/\t/, $m[1]);
		
		my @temp = ();
		my $flag = 0;
		for(my $i=0;$i<=$#dm;$i++){
			if($dm[$i] != 0){
				my $gpc = $rm[$i] / $dm[$i];
				push @temp, $gpc;
			}else{
				if($rm[$i] == 0){
					push @temp, 0;
				}else{
					push @temp, 0;
					$flag = 1; 
				}
			}
		}

		if($flag == 1){
			$excep{$m[0]} = $dna{$name};
		}

		my $o = join("\t", @temp);
		print O "$m[0]\t$o\n";
	}else{
		die "$m[0]\n";
	}

}
close R;
close O;

for (keys %excep){
	print "$_\t$excep{$_}\n";
}

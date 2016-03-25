#!/usr/bin/env perl
#########################################################################
# File Name: CIRCpseudo.pl
# Author: Rui
# mail: dongruipicb@gmail.com
# Created Time: Thu 03 Mar 2016 02:55:30 PM CST
#########################################################################

use strict;
use warnings;
use Getopt::Long;

sub usage {
	print <<"END_USAGE";

CIRCpseudo.pl 1.0 -- circRNA-derived pseudogene analysis tool

Usage: perl $0 [options]

Required:
	--circ 		CircRNA file (CIRCexplorer format)
	--ref 		Reference file (refFlat format)
	--genome 	Genome file (fa format)
	--bwaidx 	bwa index of reference genome
	--output 	output file
Optional:
	--mismatch	max mismaches between fusion sequences and genome, defalt 4
	--fusionlen	fusion lenth of back-splice exon-exon junctions defalt 40
END_USAGE
	exit;
}

my ($circ,$ref,$genome,$output,$bwaidx,$mismatch,$fusionlen);
GetOptions (
	'circ=s'=>\$circ,
	'ref=s'=>\$ref,
	'genome=s'=>\$genome,
	'bwaidx=s'=>\$bwaidx,
	'output=s'=>\$output,
	'mismatch=i'=>\$mismatch,
	'fusionlen=i'=>\$fusionlen,
	'output=s'=>\$output,
) or usage();
usage() if (!$circ or !$ref or !$genome or !$output or !$bwaidx);

$mismatch=4 if (!$mismatch);
$fusionlen=40 if (!$fusionlen);
my $half_len=$fusionlen/2;
open my $circ_in,"$circ";
mkdir "temp",0755;
my (%hash_judge_circ,%hash_judge_linear_splice,%hash_judge_exon_intron);
open my $temp_loc,">temp/temp_loc.bed";
print "Generate circRNA back-splice exon-exon junctions\n";
while (<$circ_in>){
	chomp;
	my @F=split;
	next if $F[13]=~/ciRNA|yes/;
	my $name=join "-",@F[13,14,0..2];
	my $dis=$F[2]-$half_len-$F[1];
	print $temp_loc "$F[0]\t$F[1]\t$F[2]\t$name\t0\t$F[5]\t$F[6]\t$F[7]\t0,0,0\t2\t$half_len,$half_len\t0,$dis\n" if (!exists $hash_judge_circ{"$F[0]_$F[1]_$F[2]"});
	$hash_judge_circ{"$F[0]_$F[1]_$F[2]"}=1;
}
close $circ_in;
print "Done\n";

open my $ref_in,"$ref";
open my $temp_linear_splice_loc,">temp/temp_linear_splice.bed";
open my $temp_exon_intron_loc,">temp/temp_exon_intron.bed";
print "Generate linear exon-exon junctions and exon-intron junctions\n";
while (<$ref_in>){
	chomp;
	my @F=split;
	my @start=split(/,/,$F[9]);
	my @end=split(/,/,$F[10]);
	for my $num1(0..$#start){
		my $left_exon_intron_start=$start[$num1]-$half_len;
		my $left_exon_intron_end=$start[$num1]+$half_len;
		my $left_name=join "-",@F[0..2],$left_exon_intron_start,$left_exon_intron_end;
		my $right_exon_intron_start=$end[$num1]-$half_len;
		my $right_exon_intron_end=$end[$num1]+$half_len;
		my $right_name=join "-",@F[0..2],$right_exon_intron_start,$right_exon_intron_end;
		print $temp_exon_intron_loc "$F[2]\t$left_exon_intron_start\t$left_exon_intron_end\t$left_name\t0\t$F[3]\n" if ($left_exon_intron_start>0 and !exists $hash_judge_exon_intron{"$F[0]_${left_exon_intron_start}_${left_exon_intron_end}"});
		print $temp_exon_intron_loc "$F[2]\t$right_exon_intron_start\t$right_exon_intron_end\t$right_name\t0\t$F[3]\n" if ($right_exon_intron_start>0 and !exists $hash_judge_exon_intron{"$F[0]_${right_exon_intron_start}_${right_exon_intron_end}"});
		$hash_judge_exon_intron{"$F[0]_${left_exon_intron_start}_${left_exon_intron_end}"}=1;
		$hash_judge_exon_intron{"$F[0]_${right_exon_intron_start}_${right_exon_intron_end}"}=1;
		next if $num1==$#start;
		for my $num2($num1+1..$#start){
			my $linear_start=$end[$num1]-$half_len;
			my $linear_end=$start[$num2]+$half_len;
			my $dis=$start[$num2]+$half_len-$end[$num1];
			my $name=join "-",@F[0..2],$start[$num1],$end[$num2];
			print $temp_linear_splice_loc "$F[2]\t$linear_start\t$linear_end\t$name\t0\t$F[3]\t$linear_start\t$linear_start\t0,0,0\t2\t$half_len,$half_len\t0,$dis\n" if ($linear_start>0 and !exists $hash_judge_linear_splice{"$F[2]_$start[$num1]_$end[$num2]"});
			$hash_judge_linear_splice{"$F[2]_$start[$num1]_$end[$num2]"}=1;
		}
	}
}
close $ref_in;
print "Done\n";

print "Generate fa files\n";
`bedtools getfasta -fi $genome -bed temp/temp_exon_intron.bed -fo temp/temp_exon_intron.fa -name -s`;
`bedtools getfasta -fi $genome -bed temp/temp_linear_splice.bed -fo temp/temp_linear_splice.fa -name -split -s`;
`bedtools getfasta -fi $genome -bed temp/temp_loc.bed -fo temp/temp_loc.fa -name -split -s`;

open my $back_junc_in,"temp/temp_loc.fa";
open my $back_splice,">temp/temp_back_splice.fa";
my $back_splice_name;
while (<$back_junc_in>){
	chomp;
	my @F=split;
	if ($_=~/^>/){
		$back_splice_name=$_;
		next;
	}
	my $left_seq=substr($_,0,$half_len);
	my $right_seq=substr($_,$half_len,$half_len);
	print $back_splice "$back_splice_name\n$right_seq$left_seq\n";
}
close $back_junc_in;
close $back_splice;
print "Done\n";

print "Generate Bowtie index\n";
`bowtie-build temp/temp_linear_splice.fa temp/temp_linear_splice.fa`;
`bowtie-build temp/temp_exon_intron.fa temp/temp_exon_intron.fa`;
print "Done\n";

print "Bowtie and BWA mapping\n";
`bowtie temp/temp_linear_splice.fa -f temp/temp_back_splice.fa -p 5 -m 1 -v 2 -S temp/map2linear.sam --un temp/unmap2linear.fa &>temp/temp_log`;
`bowtie temp/temp_exon_intron.fa -f temp/unmap2linear.fa -p 5 -m 1 -v 2 -S temp/map2ei.sam --un temp/unmap2ei.fa &>>temp/temp_log`;
`bwa aln -n $mismatch -t 25 -f temp/map2genome.sai $bwaidx temp/unmap2ei.fa &>>temp/temp_log`;
`bwa samse $bwaidx temp/map2genome.sai temp/unmap2ei.fa 1>temp/map2genome_bwa.sam 2>temp/temp_log`;
print "Done\n";

print "Generate putative circRNA-derived pseudogene sites\n";
open my $bwa_map_in,"temp/map2genome_bwa.sam";
open my $output_file,">$output";
print $output_file "Host_location\tHost_gene\tFusion_sequence\tPseudo_location\tMismatches\n";
while (<$bwa_map_in>){
	chomp;
	my @F=split;
	next if $_=~/^\@/ or $F[1]==4;
	$F[-1]=~s/MD:Z://;
	my $count=$F[-1]=~s/A|T|C|G//g;
	my $end_loc=$F[3]+40;
	my $pseudoloc=join "",$F[2],":",$F[3],"-",$end_loc;
	my @host=split(/-/,$F[0]);
	my $hostloc=join "",$host[2],":",$host[3],"-",$host[4];
	print $output_file "$hostloc\t$host[1]\t$F[9]\t$pseudoloc\t$count\n";
}
close $bwa_map_in;
close $output_file;
print "Done\n";

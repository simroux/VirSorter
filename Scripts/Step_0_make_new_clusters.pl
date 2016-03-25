#!/usr/bin/env perl

use strict;
use autodie;
use File::Spec::Functions;
use File::Path 'mkpath';
use File::Which 'which';

# Script to generate a new db with putative new clusters
# Argument 0 : revision directory
# Argument 1 : Fasta file of the predicted proteins
# Argument 2 : Fasta file of the unclustered from previous Runs
# Argument 3 : Liste of prots to try to cluster
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[3])))
{
	print "# Script to generate a new db with putative new clusters
# Argument 0 : revision directory
# Argument 1 : Fasta file of the predicted proteins
# Argument 2 : Fasta file of the unclustered from previous Runs
# Argument 3 : Liste of prots to try to cluster\n";
	die "\n";
}

my $path_to_blastp      = which("blastp")      or die "No blastp\n";
my $MCX_LOAD            = which("mcxload")     or die "No mcxload\n";
my $MCL                 = which("mcl")         or die "No mcl\n";
my $path_to_muscle      = which("muscle")      or die "No muscle\n";
my $path_to_hmmbuild    = which("hmmbuild")    or die "No hmmbuild\n";
my $path_to_hmmpress    = which("hmmpress")    or die "No hmmpress\n";
my $path_to_makeblastdb = which("makeblastdb") or die "No makeblastdb\n";

my $r_dir=$ARGV[0];
$r_dir=~/(r_\d*)\/?$/;
my $r_number=$1;
print "Revision $r_number\n";
my $fasta_prot_contigs=$ARGV[1];
my $fasta_prot_unclustered=$ARGV[2];
my $blast_unclustered=$fasta_prot_unclustered;
$blast_unclustered=~s/Pool_unclustered.faa/Blast_unclustered.tab/;
my $liste=$ARGV[3];
my $liste=$ARGV[3];
my $min_seq_in_a_cluster=3;

my %check;
open(LI,"<$liste") || die ("pblm opening liste $liste\n");
while (<LI>){
	chomp($_);
	my @tab=split(",",$_);
	foreach(@tab){$check{$_}=1;}
}
close LI;

my $pool_new= catfile($r_dir, "pool_new_proteins.fasta");
open(S1,">$pool_new") || die ("pblm opening file $pool_new");
open(FA,"<$fasta_prot_contigs") || die ("pblm opening file $fasta_prot_contigs");
my $tag=0;
while (<FA>){
	chomp($_);
	if ($_=~/^>(.*)/){
		my $seq=$1;
		$tag=0;
		if ($check{$seq}==1){
			print S1 "$_\n";
			$tag=1;
		}
	}
	elsif($tag==1){
		print S1 "$_\n";
	}
}
close FA;
close S1;

my $db= catfile($r_dir, "pool_new_proteins");
my $cmd_format="$path_to_makeblastdb -dbtype prot -in $pool_new -out $db";
print "$cmd_format\n";
my $out=`$cmd_format`;
print "makeblastdb : $out\n";
my $cmd_cat="cat $fasta_prot_unclustered >> $pool_new";
print "$cmd_cat\n";
$out=`$cmd_cat`;
print "Cat : $cmd_cat\n";
# BLAST des unclustered et des new contre les new
my $out_blast=catfile($r_dir, "pool_unclustered-and-new-proteins-vs-new-proteins.tab");
my $cmd_blast="$path_to_blastp -query $pool_new -db $db -out $out_blast -outfmt 6 -num_threads 10 -evalue 0.00001"; # On 10 cores to keep a few alive for the rest of the scripts
print "$cmd_blast\n";
$out=`$cmd_blast`;
print "Blast : $out\n";
$cmd_cat="cat $blast_unclustered >> $out_blast";
print "$cmd_cat\n";
$out=`$cmd_cat`;
print "Cat : $out\n";
print "Generating abc file\n";
my $out_abc= catfile($r_dir, "new_clusters.abc");
my $th_score=50;
my $th_evalue=0.00001;
my $max=200; # Max on sig
open(S1,">$out_abc") || die ("pblm opening file $out_abc\n");
open(BL,"<$out_blast") || die ("pblm opening file $out_blast\n");
while(<BL>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[11]>$th_score && $tab[10]<$th_evalue && $tab[0] ne $tab[1]){
		my $evalue=$tab[10];
# 		$evalue=-log10($evalue);
# 		if ($evalue>$max){$evalue=$max;}
		print S1 "$tab[0]\t$tab[1]\t$evalue\n";
	}
}
close BL;
close S1;
my $out_mci=catfile($r_dir, "new_clusters.mci");
my $out_tab=catfile($r_dir, "new_clusters.tab");
my $cmd_mcxload="$MCX_LOAD -abc $out_abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $out_mci -write-tab $out_tab";
print "$cmd_mcxload\n";
$out=`$cmd_mcxload`;
print "Mxc Load : $out\n";
my $dump_file=catfile($r_dir, "new_clusters.csv");
my $cmd_mcl="$MCL $out_mci -I 2  -use-tab $out_tab -o $dump_file";
print "$cmd_mcl\n";
$out=`$cmd_mcl`;
print "Mcl : $out\n";

my %unclustered;
my %clusters;
my %check_cluster;
my $last_cluster_id=0;
# toutes les séquences clusterisées dans des groupes de plus de 2 (3 et plus) -> on prend / Toutes les autres on les garde en tant qu'unclustered
open(DUMP,"<$dump_file") || die "pblm ouverture fichier $dump_file\n";
while(<DUMP>){
	chomp($_);
	my @tab=split("\t",$_);
	my $n_s_c=$#tab+1;
	if ($n_s_c>=$min_seq_in_a_cluster){
		# on a trouvé un cluster de plus de deux
		my $cluster_id=$last_cluster_id+1;
		$cluster_id="Phage_cluster_".$cluster_id."-".$r_number;
		print "We found a cluster with $n_s_c sequences => Cluster $cluster_id\n";
		$last_cluster_id++;
		foreach(@tab){
			$clusters{$cluster_id}{$_}=1;
			$check_cluster{$_}=1;
		}
	}
	else{
		foreach(@tab){
			$unclustered{$_}=1;
			$check_cluster{$_}=1;
		}
	}
}
close DUMP;
my %seq_temp;
my $id_c="";
open(FA,"<$pool_new") || die "pblm ouverture fichier $pool_new\n";
while(<FA>){
	chomp($_);
	if ($_=~/^>(\S*)/){
		$id_c=$1;
		if (!defined($check_cluster{$id_c})){$unclustered{$id_c}=1;$check_cluster{$id_c}=1;}
	}
	else{$seq_temp{$id_c}.=$_;}
}
close FA;

`mkdir $r_dir/clusts`;
foreach(keys %clusters){
	my $cluster_id=$_;
	my $out_file=catfile($r_dir, "clusts", $cluster_id . ".faa");
	open(S1,">$out_file") || die "pblm ouverture fichier $out_file\n";
	foreach(keys %{$clusters{$cluster_id}}){
		print S1 ">$_\n$seq_temp{$_}\n";
	}
	close S1;
}

mkpath(catdir($r_dir, 'db'));
my $pool_unclustered=catfile($r_dir, "db", "Pool_unclustered.faa");
my $blast_unclustered=catfile($r_dir, "db", "Blast_unclustered.tab");
my $pool_new_unclustered=catfile($r_dir, "db", "Pool_new_unclustered.faa");
my $blastable_new_unclustered=$pool_new_unclustered;
$blastable_new_unclustered=~s/\.faa//;

open(S2,">$pool_new_unclustered") || die "pblm ouverture fichier $pool_new_unclustered\n";
open(S1,">$pool_unclustered") || die "pblm ouverture fichier $pool_unclustered\n";
foreach(keys %unclustered){
	print S1 ">$_\n$seq_temp{$_}\n";
	if ($check{$_}==1){
		print S2 ">$_\n$seq_temp{$_}\n";
	}
}
close S1;
close S2;
print "making a blastable db from the new unclustered\n";
$out=`$path_to_makeblastdb -dbtype prot -in $pool_new_unclustered -out $blastable_new_unclustered`;
# on réduit aussi le fichier blast qu'on ajoute au blast des unclustered
open(BL,"<$out_blast") || die "pblm ouverture fichier $out_blast\n";
open(S1,">$blast_unclustered") || die "pblm ouverture fichier $blast_unclustered\n";
while(<BL>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($unclustered{$tab[0]}==1 && $unclustered{$tab[1]}==1){
		print S1 "$_\n";
	}
}
close BL;
close S1;

my $tag=0;
foreach(sort keys %clusters){
	$tag=1;
	my $ali_id=$_;
	my $path_to_file= catfile($r_dir, "clusts", $ali_id);
	my $path_to_fasta=catfile($r_dir, "clusts", $ali_id . ".faa");
	my $path_to_ali=catfile($r_dir, "clusts", $ali_id . ".ali_faa");
	my $path_to_hmm=catfile($r_dir, "clusts", $ali_id . "_ali.hmm");
	if (-e $path_to_ali){
		`rm $path_to_ali $path_to_hmm`;
	}
	my $muscle_out= catfile($r_dir, "log_out_muscle");
	my $muscle_err= catfile($r_dir, "log_err_muscle");
	`$path_to_muscle -in $path_to_fasta -out $path_to_ali > $muscle_out 2> $muscle_err`;
	my $out_stokcholm=$path_to_ali.".stockholm";
	open(S1,">$out_stokcholm") || die "pblm opening $out_stokcholm\n";
	print S1 "# STOCKHOLM 1.0\n";
	open(FA,"<$path_to_ali") || die "pblm ouverture $path_to_ali\n";
	while(<FA>){
		chomp($_);
		if ($_=~/^>(.*)/){
			my $id=$1;
			$id=~s/\s/_/g;
			print S1 "\n$id  ";
			
		}
		else{print S1 "$_";}
	}
	close FA;
	print S1 "\n//\n";
	`$path_to_hmmbuild --amino $path_to_hmm $out_stokcholm`;
}

my @tab_hmm=<$r_dir/clusts/*.hmm>;
if ($#tab_hmm>=0){
	# we gather all hmm and fasta (if any)
	$out=`cat $r_dir/clusts/*.hmm >> $r_dir/db/Pool_clusters.hmm`;
	print "cat new hmm : $out\n";
	# and create a new database for hmmsearch
	$out=`$path_to_hmmpress $r_dir/db/Pool_clusters.hmm`;
	print "hmm press :$out\n";
}
else{
	$out=`touch $r_dir/db/Pool_clusters.hmm`;
}

#!/usr/bin/env perl

use strict;
use autodie;

# Script to make a summary of the predictions to add to previous predictions
# Argument 0 : summary file of the phage fragments
# Argument 1 : global summary to be completed
# Argument 2 : Out file for new prot list
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[3])))
{
	print "# Script to make a summary of the predictions to add to previous predictions
# Argument 0 : affiliation file of the contigs
# Argument 1 : summary file of the phage fragments
# Argument 2 : global summary to be completed
# Argument 3 : Out file for new prot list\n";
	die "\n";
}

my $affi_contigs   = $ARGV[0];
my $new_summary    = $ARGV[1];
my $global_summary = $ARGV[2];
my $new_prot_list  = $ARGV[3];

my %infos;
my $tag=0;
my %check_prot_old;
my %check_contig_old;
if (-e $global_summary){
	 # Get info from global_summary
	open SUM, '<', $global_summary;
	while (<SUM>){
		chomp($_);
		if ($_=~/^## (\d+)/){
			$tag=$1;
		}
		elsif($_=~/^##/ || $_ eq ""){}
		elsif($tag<=3){
# 			print "we had $_ -> tag $tag\n";
			my @tab=split(",",$_);
			$infos{$tag}{$tab[2]}{"nb_gene"}=$tab[1];
			$infos{$tag}{$tab[2]}{"category"}=$tab[4];
			$infos{$tag}{$tab[2]}{"hallmark"}=$tab[5];
			$infos{$tag}{$tab[2]}{"phage"}=$tab[6];
			$infos{$tag}{$tab[2]}{"noncaudo"}=$tab[7];
			$infos{$tag}{$tab[2]}{"pfam"}=$tab[8];
			$infos{$tag}{$tab[2]}{"unch"}=$tab[9];
			$infos{$tag}{$tab[2]}{"switch"}=$tab[10];
			$infos{$tag}{$tab[2]}{"size"}=$tab[11];
			$check_contig_old{$tab[2]}=1;
			print "checking old $tab[2]\n";
		}
		else{
			my @tab=split(",",$_);
			print "we had $tab[2] -> tag $tag\n";
			$infos{$tag}{$tab[0]}{$tab[2]}{"nb_gene_contig"}=$tab[1];
			$infos{$tag}{$tab[0]}{$tab[2]}{"nb_gene"}=$tab[3];
			$infos{$tag}{$tab[0]}{$tab[2]}{"category"}=$tab[4];
			$infos{$tag}{$tab[0]}{$tab[2]}{"hallmark"}=$tab[5];
			$infos{$tag}{$tab[0]}{$tab[2]}{"phage"}=$tab[6];
			$infos{$tag}{$tab[0]}{$tab[2]}{"noncaudo"}=$tab[7];
			$infos{$tag}{$tab[0]}{$tab[2]}{"pfam"}=$tab[8];
			$infos{$tag}{$tab[0]}{$tab[2]}{"unch"}=$tab[9];
			$infos{$tag}{$tab[0]}{$tab[2]}{"switch"}=$tab[10];
			$infos{$tag}{$tab[0]}{$tab[2]}{"size"}=$tab[11];
			if($infos{$tag}{$tab[0]}{$tab[2]}{"category"}==1){ # If the category is 4, we check all the prot from this fragment
				$tab[2]=~/.*-gene_(\d*)-gene_(\d*)/;				
				for (my $i=$1;$i<=$2;$i++){
					my $prot_id=$tab[0]."-gene_".$i;
					$check_prot_old{$prot_id}=1;
				}
			}
		}
	}
	close SUM;
}
else{
	print "This is the first global summary that we'll do\n";
}

my %check_prot_new;
my %check_contig_new;
open SUM, '<', $new_summary;
while (<SUM>){
	chomp($_);
	$_=~s/,/;/g;
	my @tab=split("\t",$_);
	if ($tab[4] eq "complete_phage"){
#    0   /       1        /    2       /  3  /     4          /    5      /      6        /     7        /     8       /       9        /     10   /      11     /    
# Contig / Total Nb Genes /  Fragment / Size / Type detection / Category /  Enrich Phage / Enrich Pfam / Enrich Unch / Enrich Switch / Avg_g_size / Nb Hallmark
		# Determine order in which this contig will be displayed
		my $class=3;
		if ($tab[5]==1){# If the category is 1, we check all the prot from this contig
			$class=1;
			$check_contig_new{$tab[0]}=1;
		}
		elsif ($tab[5]==2){$class=2;}
 		for(my $i=5;$i<=$#tab;$i++){
 			if ($tab[$i]=~/(.*);$/){
 				$tab[$i]=$1;
 			}
 		}
# 		print "$_ => tag $class\n";
		$infos{$class}{$tab[0]}{"nb_gene"}=$tab[1];
		$infos{$class}{$tab[0]}{"category"}=$tab[5];
		$infos{$class}{$tab[0]}{"phage"}=$tab[6];
		$infos{$class}{$tab[0]}{"noncaudo"}=$tab[7];
		$infos{$class}{$tab[0]}{"pfam"}=$tab[8];
		$infos{$class}{$tab[0]}{"unch"}=$tab[9];
		$infos{$class}{$tab[0]}{"switch"}=$tab[10];
		$infos{$class}{$tab[0]}{"size"}=$tab[11];
		$infos{$class}{$tab[0]}{"hallmark"}=$tab[12];
	}
	else{
		my $class=6;
		if ($tab[5]==1){
			$class=4;
			# If the category is 1, we check all the prot from this fragment
			$tab[2]=~/.*-gene_(\d*)-gene_(\d*)/;
			for (my $i=$1;$i<=$2;$i++){
				my $prot_id=$tab[0]."-gene_".$i;
				$check_prot_new{$prot_id}=1;
				print "we check new $prot_id\n";
			}
		}
		elsif($tab[5]==2){$class=5;}
		# Remove all former prophages (if any) is there is an overlap
		for (my $i=4;$i<=6;$i++){
			if (defined $infos{$i}{$tab[0]}){
				foreach (keys %{$infos{$i}{$tab[0]}}){
					if (overlap($tab[2],$_)==1){
						print "Overlap between $tab[1] and $_, we remove $_ ($tab[0] - 3)\n";
						delete($infos{$i}{$tab[0]}{$_});
					}
				}
			}
		}
 		for(my $i=5;$i<=$#tab;$i++){
 			if ($tab[$i]=~/(.*);$/){
 				$tab[$i]=$1;
 			}
 		}
# 		print "Prophage $class / $tab[0] - $tab[2]\n";
		$infos{$class}{$tab[0]}{$tab[2]}{"nb_gene_contig"}=$tab[1];
		$infos{$class}{$tab[0]}{$tab[2]}{"nb_gene"}=$tab[3];
		$infos{$class}{$tab[0]}{$tab[2]}{"category"}=$tab[5];
		$infos{$class}{$tab[0]}{$tab[2]}{"phage"}=$tab[6];
		$infos{$class}{$tab[0]}{$tab[2]}{"noncaudo"}=$tab[7];
		$infos{$class}{$tab[0]}{$tab[2]}{"pfam"}=$tab[8];
		$infos{$class}{$tab[0]}{$tab[2]}{"unch"}=$tab[9];
		$infos{$class}{$tab[0]}{$tab[2]}{"switch"}=$tab[10];
		$infos{$class}{$tab[0]}{$tab[2]}{"size"}=$tab[11];
		$infos{$class}{$tab[0]}{$tab[2]}{"hallmark"}=$tab[12];
	}
}

# Remove redundancy 
foreach(sort {$a <=> $b } keys %infos){
	my $class=$_;
	my @liste_contigs=keys %{$infos{$class}};
	if ($class<=3){ ## For complete phages, remove all predictions with higher categories
		for (my $i=$class+1;$i<=6;$i++){
			foreach(@liste_contigs){
				if (defined($infos{$i}{$_})){
					print "$_ defined in $class, so we remove its info for $i\n";
					delete($infos{$i}{$_});
				}
			}
		}
	}
	else{
		foreach(@liste_contigs){ ## For prophages, remove the prediction of the same prophages with higher categories
			my $contig=$_;
			foreach(keys %{$infos{$class}{$contig}}){
				for (my $i=$class+1;$i<=6;$i++){
					if (defined($infos{$i}{$contig}{$_})){
						print "$_ defined in $class, so we remove its info for $i\n";
						delete($infos{$i}{$contig}{$_});
					}
				}
			}
			
			
		}
	}
}


open S1, '>', $global_summary;
for (my $class=1;$class<=6;$class++){
	if ($class==1){
		print S1 "## 1 - Complete phage contigs - category 1 (sure)\n";
		print S1 "## Contig_id,Nb genes contigs,Fragment,Nb genes,Category,Nb phage hallmark genes,Phage gene enrichment sig,Non-Caudovirales phage gene enrichment sig,Pfam depletion sig,Uncharacterized enrichment sig,Strand switch depletion sig,Short genes enrichment sig\n";
	}
	if ($class==2){
		print S1 "## 2 - Complete phage contigs - category 2 (somewhat sure)\n";
		print S1 "## Contig_id,Nb genes contigs,Fragment,Nb genes,Category,Nb phage hallmark genes,Phage gene enrichment sig,Non-Caudovirales phage gene enrichment sig,Pfam depletion sig,Uncharacterized enrichment sig,Strand switch depletion sig,Short genes enrichment sig\n";
	}
	if ($class==3){
		print S1 "## 3 - Complete phage contigs - category 3 (not so sure)\n";
		print S1 "## Contig_id,Nb genes contigs,Fragment,Nb genes,Category,Nb phage hallmark genes,Phage gene enrichment sig,Non-Caudovirales phage gene enrichment sig,Pfam depletion sig,Uncharacterized enrichment sig,Strand switch depletion sig,Short genes enrichment sig\n";
	}
	if ($class==4){
		print S1 "## 4 - Prophages - category 1 (sure)\n";
		print S1 "## Contig_id,Nb genes contigs,Fragment,Nb genes,Category,Nb phage hallmark genes,Phage gene enrichment sig,Non-Caudovirales phage gene enrichment sig,Pfam depletion sig,Uncharacterized enrichment sig,Strand switch depletion sig,Short genes enrichment sig\n";
	}
	if ($class==5){
		print S1 "## 5 - Prophages - category 2 (somewhat sure)\n";
		print S1 "## Contig_id,Nb genes contigs,Fragment,Nb genes,Category,Nb phage hallmark genes,Phage gene enrichment sig,Non-Caudovirales phage gene enrichment sig,Pfam depletion sig,Uncharacterized enrichment sig,Strand switch depletion sig,Short genes enrichment sig\n";	
	}
	if ($class==6){
		print S1 "## 6 - Prophages - category 3 (not so sure)\n";
		print S1 "## Contig_id,Nb genes contigs,Fragment,Nb genes,Category,Nb phage hallmark genes,Phage gene enrichment sig,Non-Caudovirales phage gene enrichment sig,Pfam depletion sig,Uncharacterized enrichment sig,Strand switch depletion sig,Short genes enrichment sig\n";	
	}
	foreach(sort keys %{$infos{$class}}){
		my $contig=$_;
		if ($class<=3){
			print S1 "$_,$infos{$class}{$contig}{nb_gene},$_,$infos{$class}{$contig}{nb_gene},$infos{$class}{$contig}{category},$infos{$class}{$contig}{hallmark},$infos{$class}{$contig}{phage},$infos{$class}{$contig}{noncaudo},$infos{$class}{$contig}{pfam},$infos{$class}{$contig}{unch},$infos{$class}{$contig}{switch},$infos{$class}{$contig}{size}\n";
		}
		else{
			foreach (sort keys %{$infos{$class}{$contig}}) {
				print S1 "$contig,$infos{$class}{$contig}{$_}{nb_gene_contig},$_,$infos{$class}{$contig}{$_}{nb_gene},$infos{$class}{$contig}{$_}{category},$infos{$class}{$contig}{$_}{hallmark},$infos{$class}{$contig}{$_}{phage},$infos{$class}{$contig}{$_}{noncaudo},$infos{$class}{$contig}{$_}{pfam},$infos{$class}{$contig}{$_}{unch},$infos{$class}{$contig}{$_}{switch},$infos{$class}{$contig}{$_}{size}\n";
			}
		}
	}
}
close S1;

# Check if they could be new clusters among the new proteins
my @liste_to_add=();
my $th_evalue=0.0000000001; # Big threshold, to prevent too much false positive
open AFI, '<', $affi_contigs;
my $contig_c="";
while (<AFI>){
	chomp($_);
	if ($_=~/>(.*)/){
		my @tab=split(/\|/,$1);
		$contig_c=$tab[0];
	}
	else{
		#     0  | 1   | 2  |  3   |  4   |    5     |  6  |   7  |   8     |   9    | 10   | 11
		# gene_id|start|stop|length|strand|affi_phage|score|evalue|category|affi_pfam|score|evalue|
		my @tab=split(/\|/,$_);
		my $gene=$tab[0];
		if (($check_prot_new{$gene}==1 && !defined($check_prot_old{$gene})) || (($check_contig_new{$contig_c}==1) && !defined($check_contig_old{$contig_c}))){
# 			print "Ah, a new prot, putatively a new cluster\n";
			if ($tab[5] eq "-"){ 
# 				print "\t oh yep, no phage cluster, so we take this\n";
				push(@liste_to_add,$gene);
			}
			elsif($tab[7]<=$th_evalue && $tab[5]=~/^Phage_cluster_\d+/){
# 				print "\t nope, evalue of $tab[7] on a cluster, so will likely cluster with an existing PC\n";
			}
			else{
# 				print "\t oh yep, we take this, because it's either a bigger evalue than the th or a non-clustered phage protein -> $tab[5]\n";
				push(@liste_to_add,$gene);
			}
		}
	}
}

if ($#liste_to_add>=0){
	print "Listing the new prots to add\n";
	my $l=join(",",@liste_to_add);
	open S2, '>', $new_prot_list;
	print S2 "$l\n";
	close S2;
}
else{
	print "No new prots, no list\n";
}


sub overlap { # To check if a prediction is not within another of the same type, in which case we don't really care
	my $pred_wide=$_[0];
	my $pred_short=$_[1];
	$pred_wide=~/.*-gene_(\d*)-gene_(\d*)/;
	my $start_pred_wide=$1; 
	my $end_pred_wide=$2;
	$pred_short=~/.*-gene_(\d*)-gene_(\d*)/;
	my $start_pred_short=$1; 
	my $end_pred_short=$2;
	my $o=0;
	if ($start_pred_short>=$start_pred_wide && $end_pred_short<=$end_pred_wide && (!($start_pred_short==$start_pred_wide && $end_pred_short==$end_pred_wide))){
		$o=1;
	}
# 	print "$pred_short / $pred_wide $start_pred_short>=$start_pred_wide && $end_pred_short<=$end_pred_wide => $o\n";
	return $o;
}

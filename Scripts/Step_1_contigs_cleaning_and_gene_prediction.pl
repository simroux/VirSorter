#!/usr/bin/env perl

use strict;
use autodie;
use Bio::Seq;
use File::Spec::Functions;
use File::Which 'which';

# Script to detect circular contigs, nett sequences, and predict genes with mga
# Argument 0 : Fasta file of contigs
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[3])))
{
	print "# Script to detect circular contigs, nett sequences, and predict genes with mga
# Argument 0 : Id du dataset
# Argument 1 : Working dir
# Argument 2 : Fasta file of contigs
# Argument 3 : Threshold on the number of genes \n";
	die "\n";
}

my $id                = $ARGV[0];
my $tmp_dir           = $ARGV[1];
my $fasta_contigs     = $ARGV[2];
my $th_nb_genes       = $ARGV[3];
my $path_to_mga       = which('mga_linux_ia64') or die "Cannot find mga_linux_ia64\n";
my $in_file           = catfile($tmp_dir, $id . "_nett.fasta");
my $circu_file        = catfile($tmp_dir, $id . "_circu.list");
my $out_special_circu = catfile($tmp_dir, $id . "_contigs_circu_temp.fasta");

# Reading fasta file of the contigs
open my $fa, '<', $fasta_contigs;
my %seq_base;
my $id_seq="";
while(<$fa>){
	$_=~s/\r\n/\n/g; #Cas d'un fichier windows ##AJOUT
	chomp($_);
	if ($_=~/^>(\S*)/){$id_seq=$1;}
	else{$seq_base{$id_seq}.=$_;}
}
close $fa;

## DETECTION OF CIRCULAR CONTIG AND CLEANING OF THESE CIRCULAR (REMOVE THE MATCHING ENDS)
my $minimum_size=1500;
my %order_contig;
my %length;
my $n1=0;

open my $s1, '>', $in_file;
open my $s2, '>', $circu_file;
for my $id_contig (
    sort {length($seq_base{$b}) <=> length($seq_base{$a})} keys %seq_base){
	$order_contig{$id_contig}=$n1;
	$n1++;
	my $s=$seq_base{$id_contig};
	$length{$id_contig}=length($seq_base{$id_contig});
	# Test its circularity
	my $prefix=substr($seq_base{$id_contig},0,10);
	if ($s=~/(.+)($prefix.*?)$/){
# 		print "on a retrouvé prefix ($prefix) plus loin dans la séquence de $_\n";
		my $sequence=$1;
		my $suffixe=$2;
		my $test=substr($seq_base{$id_contig},0,length($suffixe));
# 		print "$suffixe\n$test\n";
		if ($suffixe eq $test){
# 			print " et il match bien $suffixe, donc c'est un contig circulaire\n";
			my $l=$length{$id_contig};
			$id_contig=$id_contig."-circular";
			$length{$id_contig}=$l;
			print $s2 "$id_contig\t$length{$id_contig}\n";
			$seq_base{$id_contig}=$sequence;
		}
	}
	# Update the length of the contig
	$length{$id_contig}=length($seq_base{$id_contig});
	print $s1 ">$id_contig\n$seq_base{$id_contig}\n";
}
close $s1;
close $s2;

# Gene prediction for all contigs
my $out_file= $tmp_dir."/".$id."_mga.predict";
print "mga ($path_to_mga) $in_file -m > $out_file\n";
my $mga=`$path_to_mga $in_file -m > $out_file`;

# Special prediction for circular contigs if we have some
my $out_file_circu="";
my %circu;
if (-e $circu_file){
	open my $tsv, '<', $circu_file;
	while(<$tsv>){
		chomp($_);
		my @tab=split("\t",$_);
		my $id_c=$tab[0];
		$circu{$id_c}=1;
	}
	close $tsv;
	open my $s3, '>', $out_special_circu;
	my $long=1000; # we cp the 1000 first bases to the end of the contig
	my $seuil_long=1000;
	my $n_circu=0;
	foreach(sort {$order_contig{$a} <=> $order_contig{$b} } keys %circu){
		my $id_c=$_;
		my $s=$seq_base{$id_c}.substr($seq_base{$id_c},0,$long);
		print $s3 ">$id_c\n$s\n";
		$n_circu++;
	}
	close $s3;
	$out_file_circu= $tmp_dir."/".$id."_special_circus_mga.predict";
	if ($n_circu>0){
		my $mga=`$path_to_mga $out_special_circu -m > $out_file_circu`;
	}
	else{
		`touch $out_file_circu`;
	}
}

# Mix 'n match of the two results of gene prediction
my %order_gene;
my $n2=0;
open my $fts, '<', $out_file;
my %predict;
my %type;
my $id_c="";
while(<$fts>){
	chomp($_);
	if($_=~/^# gc/){}
	elsif($_=~/^# self: (.*)/){$type{$id_c}=$1;}
	elsif ($_=~/^# (.*)/){
		$id_c=$1;
		$n2=0;
	}
	else{
		my @tab=split("\t",$_);
		$predict{$id_c}{$tab[0]}=$_;
		if (!defined($order_gene{$id_c}{$tab[0]})){$order_gene{$id_c}{$tab[0]}=$n2;$n2++;}
	}
}
close $fts;
if (-e $circu_file){
	open my $fts_c, '<', $out_file_circu;
	my $tag=0;
	while(<$fts_c>){
		chomp($_);
		if($_=~/^# gc/){}
		elsif($_=~/^# self: (.*)/){$type{$id_c}=$1;}
		elsif ($_=~/^# (.*)/){
			if($tag==1){
				my %to_start;
				# Some ORFs were modified, we clean up
				foreach(sort {$order_gene{$a} <=> $order_gene{$b} } keys %{$predict{$id_c}}){
					my @tab=split("\t",$predict{$id_c}{$_});
					if ($tab[5]!=11){
						# $tab[0] miss start and/or stop codon
						if(($tab[1]<3) || ($tab[2]>($length{$id_c}-3))){
							# And it spans the origin, so we can remove it
							if ($tab[1]<3){
								$to_start{$tab[0]}{"start"}=$tab[1];
								$to_start{$tab[0]}{"stop"}=$tab[2];
							}
							delete($predict{$id_c}{$tab[0]});
						}
						elsif(($tab[2]>997) && ($tab[2]<1001)){ # if we are around the zone of ~ 1000
							foreach(keys %to_start){
								my $total=($length{$id_c}-$tab[1]+1)+($to_start{$_}{"stop"}); 
								if ($total % 3 == 0){
									$tab[2]=$to_start{$_}{"stop"};
									$tab[5]=11;
									my $new_line=join("\t",@tab);
									$predict{$id_c}{$tab[0]}=$new_line;
								}
							}
						}
					}
				}
			}
			$id_c=$1;
			$tag=0;
		}
		else{
			my @tab=split("\t",$_);
			if (defined($predict{$id_c}{$tab[0]})){
				my @tab2=split("\t",$predict{$id_c}{$tab[0]});
				if (($tab2[1]==$tab[1]) && ($tab2[2]==$tab[2])){}# same prediction, we don't change anything
				else{
					if (($tab[1]<$length{$id_c}) && ($tab[2]>$length{$id_c})){
						# we span the origin, we replace the prediction
						$tag=1;
						my $stop=$tab[2]-$length{$id_c};
						$tab[2]=$stop;
						my $new_line=join("\t",@tab);
						$predict{$id_c}{$tab[0]}=$new_line;
					}
				}
			}
			else{
				# we predict a new gene, we keep only if at the start / end of the contig
				if (($tab[1]<$length{$id_c}) && ($tab[2]>$length{$id_c})){
					$tag=1;
					my $stop=$tab[2]-$length{$id_c};
					$tab[2]=$stop;
					my $new_line=join("\t",@tab);
					$predict{$id_c}{$tab[0]}=$new_line;
					$tag=1;
				}
			}
		}
	}
	if($tag==1){
		my %to_start;
		# we changed some things, we clean up
		foreach(sort {$order_gene{$a} <=> $order_gene{$b} } keys %{$predict{$id_c}}){
			my @tab=split("\t",$predict{$id_c}{$_});
			if ($tab[5]!=11){
				if(($tab[1]<3) || ($tab[2]>($length{$id_c}-3))){
					if ($tab[1]<3){
						$to_start{$tab[0]}{"start"}=$tab[1];
						$to_start{$tab[0]}{"stop"}=$tab[2];
					}
					delete($predict{$id_c}{$tab[0]});
				}
				elsif(($tab[2]>997) && ($tab[2]<1001)){
					foreach(keys %to_start){
						my $total=($length{$id_c}-$tab[1]+1)+($to_start{$_}{"stop"}); 
						if ($total % 3 == 0){
							$tab[2]=$to_start{$_}{"stop"};
							$tab[5]=11;
							my $new_line=join("\t",@tab);
							$predict{$id_c}{$tab[0]}=$new_line;
						}
					}
				}
			}
		}
	}
	close $fts_c;
}

## Generation of the final files
## One with all sequences nett and filtered (based on number of genes) - Fasta
## One of the associated gene prediction - MGA-like
## One of the predicted protein sequences - Fasta
my $final_file=$tmp_dir."/".$id."_nett_filtered.fasta";
my $out_final=$tmp_dir."/".$id."_mga_final.predict";
my $prot_file=$tmp_dir."/".$id."_prots.fasta";

open my $fa_s,  '>', $final_file;
open my $out_s, '>', $out_final;
open my $prot_s,'>', $prot_file;
my $n=0;
foreach(sort {$order_contig{$a} <=> $order_contig{$b} } keys %predict){
	$n++;
	if ($n % 10000 == 0){print "$n-ieme contig\n";}
	my $id=$_;
	my @tab_genes=sort {$order_gene{$id}{$a} <=> $order_gene{$id}{$b} } keys %{$predict{$id}};
	my $n_complete_genes=0;
	for (my $i=0;$i<=$#tab_genes;$i++){
		my @tab=split("\t",$predict{$id}{$tab_genes[$i]});
		if ($tab[5]!=11){}
		else{$n_complete_genes++;}
	}
	if ($n_complete_genes<$th_nb_genes){
# 		print "$id is excluded because too short ($n_complete_genes) \n";
	}
	else{
		## We check the first gene and modify it if needed
		my @tab_first=split("\t",$predict{$id}{$tab_genes[0]});
		my @tab_second=split("\t",$predict{$id}{$tab_genes[1]});
		$tab_first[0]=~/gene_(\d*)/;
		my $n_1=$1;
		$tab_second[0]=~/gene_(\d*)/;
		my $n_2=$1;
		if ($n_1>$n_2){
			print "We probably have a circular contig ($id) as the first gene $tab_first[0] is beyond the second gene $tab_second[0] ($n_1>$n_2), so we switch $tab_first[0] ";
			$tab_first[0]="gene_0";
			print "to $tab_first[0]\n";
			$predict{$id}{$tab_genes[0]}=join("\t",@tab_first);
		}
# 		if ($n_complete_genes<$th_nb_genes){print "$id is saved because of its circularity\n";}
# 		else{print "We keep $id = $#tab_genes +1 genes\n";}
		print $out_s ">$id\t$length{$id}\n";
		print $fa_s ">$id\n";
		my $seq_c=$seq_base{$id};
		print $fa_s "$seq_c\n";
		foreach(@tab_genes){
			my @tab=split("\t",$predict{$id}{$_});
			if ($tab[5]!=11){
				# soit on est au début de séquence, soit en toute fin (théoriquement)
				if ($tab[4]!=0){
					if ($tab[3] eq "-"){
						$tab[2]-=$tab[4];
					}
					elsif($tab[3] eq "+"){
						$tab[1]+=$tab[4];
					}
					else{
						print "%%%%%% pblm on a pas de sens pour $id : @tab\n";
					}
				}
				my $new_line=join("\t",@tab);
				$predict{$id}{$_}=$new_line;
			}
			print $out_s "$predict{$id}{$_}\n";
			@tab=split("\t",$predict{$id}{$_});
	                my $name  = $tab[0];
	                my $start = $tab[1];
	                my $stop  = $tab[2];
	                my $sens  = $tab[3];
	                my $frame = $tab[4];
	                my $frag  = "";
			# Regular case (not spanning the origin)
			if ($start<$stop){
				my $long=$stop-$start+1;
				$frag=substr($seq_c,$start-1,$long);
			}
			# Exceptional case, we span the origin
			else{
				my $l1=length($seq_c)-$start+1;
				$frag=substr($seq_c,$start-1,$l1);
				$frag.=substr($seq_c,0,$stop);
			}
			## WE GET THE PREDICTED PROTEIN SEQUENCE
			if ($frag eq ""){
				print "!!!! FRAG IS $frag\n";
			}
			my $seq_bio = Bio::Seq->new(-id => "dummy_id" , -seq =>$frag, -alphabet => 'dna', -verbose => -1);
			# Test to catch the Bio SeqUtils warning
			my @seqs;
			eval{
				@seqs = Bio::SeqUtils->translate_6frames($seq_bio, -verbose => -1);
			};
			if ( $@ ){
				print "We got the error $@\n";
			}
			#my @seqs = Bio::SeqUtils->translate_6frames($seq_bio, -verbose => -1);
			# End of test
			my $cadre=0;
			if ($sens eq "-"){$cadre=3;}
			my $prot=$seqs[$cadre];
			my $prot_sequence=$prot->seq;
			if ($prot_sequence=~/\*$/){
				# we remove the stop codon
				chop($prot_sequence);
			}
			my $id_out=$id."-".$name;
			if (($prot_sequence=~/X{50,}/) || ($prot_sequence=~/F{50,}/) || ($prot_sequence=~/A{50,}/) || ($prot_sequence=~/K{50,}/) || ($prot_sequence=~/P{50,}/)){
				print "we exclude $id_out because there is a pblm with the sequence -> too many succesive X, F, A, K or P\n";
			}
			else{
				print $prot_s ">$id_out\n$prot_sequence\n";
			}
		}
	}
}
close $fa_s;
close $out_s;
close $prot_s;

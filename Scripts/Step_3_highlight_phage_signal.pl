#!/usr/bin/env perl

use strict;
use autodie;
use File::Spec::Functions;
use FindBin '$Bin';

# Script to measure metrics on the sliding window
# Argument 0 : csv file of the contigs
# Argument 1 : summary file of the phage fragments
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[1])))
{
	print "# Script to measure metrics on the sliding window
# Argument 0 : csv file of the contigs
# Argument 1 : summary file of the phage fragments
# Argument 2 (optional) : a file with the refs values that we could use instead of estimating them \n";
	die "\n";
}
$| = 1;
my $csv_file = $ARGV[0];
my $out_file = $ARGV[1];
if ( -e $out_file ) { `rm $out_file`; }
my $ref_file = $ARGV[0];
$ref_file =~ s/\.csv/.refs/g;
my $do_ref_estimation = 0;
if (defined($ARGV[2])){
#	$ref_file=$ARGV[2];
	`cp $ARGV[2] $ref_file`; # That way, the ref file is in the result directory if a use wants to check it
	$do_ref_estimation=1;
}

## ABSOLUTE THRESHOLDS ##
my $th_viral_hallmark=1;
my $th_sig=2;
my $th_sig_2=4;
my $th_nb_genes_covered=0.80;
my $th_nb_genes_noncaudo=1;
## END OF ABSOLUTE THRESHOLDS ##
my $script_dir= catfile($Bin);
my $path_to_c_script= catfile($script_dir, "Sliding_windows_3");

print "## Taking information from the contig info file ($csv_file)\n";
open F1, '<', $csv_file;
my $n=0;
my $id_c=$_;
my %infos;
my @liste_contigs;
my %nb_genes;
while(<F1>){
	chomp($_);
	if ($_=~/>(.*)/){
		my @tab=split(/\|/,$1);
		$id_c=$tab[0];
		push(@liste_contigs,$id_c);
		$nb_genes{$id_c}=$tab[1];
		$n++;
	}
	else{
		#     0  | 1   | 2  |  3   |  4   |    5     |  6  |   7  |   8     |   9    | 10   | 11
		# gene_id|start|stop|length|strand|affi_phage|score|evalue|category|affi_pfam|score|evalue|
		my @tab=split(/\|/,$_);
		my $gene=$tab[0];
		$gene=~/.*-(gene_\d*)/;
		$gene=$1;
		$infos{$id_c}{$gene}{"order"}=$n;
		$infos{$id_c}{$gene}{"length"}=$tab[3];
		$infos{$id_c}{$gene}{"strand"}=$tab[4];
		$infos{$id_c}{$gene}{"category"}=-1;
		if ($tab[5] eq "-"){ ## no Phage Cluster affiliation
			if ($tab[9] eq "-"){ ## no PFAM either, ok.. 
				$infos{$id_c}{$gene}{"best_domain_hit"}="-";
			}
			else{
				$infos{$id_c}{$gene}{"best_domain_hit"}="PFAM-".$tab[9];
			}
		}
		else{
			if ($tab[9] eq "-" || $tab[6]>$tab[10]){ ## no PFAM or Phage Cluster better than PFAM (score comparison)
				$infos{$id_c}{$gene}{"best_domain_hit"}="PC-".$tab[5];
				if ($tab[9] ne "-"){$infos{$id_c}{$gene}{"hit_PFAM"}="PFAM-".$tab[5];}
				if ($tab[8] eq "-"){$infos{$id_c}{$gene}{"category"}=-1;}
				else{$infos{$id_c}{$gene}{"category"}=$tab[8];}
			}
			else{ ## So we have a PFAM, which is clearly better than Phage Cluster, so we keep it
				$infos{$id_c}{$gene}{"best_domain_hit"}="PFAM-".$tab[9];
				$infos{$id_c}{$gene}{"hit_PC"}="PC-".$tab[5];
			}
		}
		$n++;
	}
}
close F1;

my $th_gene_size=0;
# WE HAVE A REF FILE, WE DONT ESTIMATE
if ($do_ref_estimation==1){
	print "## We have a ref file : $ref_file , so will use it\n";
	open F1, '<', $ref_file;
	while (<F1>){
		chomp($_);
		my @tab=split("\t",$_);
		$th_gene_size=$tab[4];
	}
	close F1;
}
### ELSE, IF WE ESTIMATE THE PARAMETERS FROM THE DATASET
else{
	my %total;
	my @store_avg_g_size;
	# look at all contigs at once for the global metrics
	print "## First look at everything to get the totals\n";
	foreach(@liste_contigs){
		my $contig=$_;
		my @tab_genes=sort {$infos{$contig}{$a}{"order"} <=> $infos{$contig}{$b}{"order"}} keys %{$infos{$contig}};
		my $total_nb_genes=$#tab_genes+1;
		my $n_f=0;
		# First, taking all the metrics for the totals
		my $last_strand=$infos{$contig}{$tab_genes[0]}{"strand"};
		for (my $i=0;$i<$total_nb_genes;$i++){
			if (defined($tab_genes[$i])){
				$total{"n_obs"}++;
				if ($infos{$contig}{$tab_genes[$i]}{"best_domain_hit"}=~/^PC-/){ # look at best domain hit on a phage
					$total{"phage"}++;
					if (defined($infos{$contig}{$tab_genes[$i]}{"hit_PFAM"})){$total{"pfam"}++;}
					if ($infos{$contig}{$tab_genes[$i]}{"category"}>=3){
						$total{"noncaudo"}++;
					}
				}
				elsif($infos{$contig}{$tab_genes[$i]}{"best_domain_hit"}=~/^PFAM-/){
					$total{"pfam"}++;
					if (defined($infos{$contig}{$tab_genes[$i]}{"hit_PC"})){$total{"phage"}++;}
				}
				elsif($infos{$contig}{$tab_genes[$i]}{"best_domain_hit"} eq "-"){
					$total{"unch"}++;
				}
				if ($infos{$contig}{$tab_genes[$i]}{"strand"} ne $last_strand){
					$total{"switch"}++;
					$last_strand=$infos{$contig}{$tab_genes[$i]}{"strand"};
				}
				push(@store_avg_g_size,$infos{$contig}{$tab_genes[$i]}{"length"});
			}
		}
	}

	print "## Transform it into probability and gene size decile\n";
	# Transform it into probability / ratios and sort the gene size table
	$total{"phage"}/=$total{"n_obs"};
	$total{"noncaudo"}/=$total{"n_obs"};
	$total{"pfam"}/=$total{"n_obs"};
	$total{"unch"}/=$total{"n_obs"};
	$total{"switch"}/=$total{"n_obs"};
	# Determine d1 (first decile) of the gene size distribution, so we divide the distribution in 10 parts
	$th_gene_size=get_th_gene_size(\@store_avg_g_size,10);
	open S2, '>', $ref_file;
	print S2 $total{"phage"}."\t".$total{"pfam"}."\t".$total{"unch"}."\t".$total{"switch"}."\t".$th_gene_size."\t".$total{"noncaudo"};
	close S2;
}

my $nb_gene_th=2;
# Now the sliding windows
print "## Then look at each contig and each sliding window\n";
open S1, '>', $out_file;
close S1;
my $i=0;
foreach(@liste_contigs){
	my $contig_c=$_;
	my @tab_genes=sort {$infos{$contig_c}{$a}{"order"} <=> $infos{$contig_c}{$b}{"order"}} keys %{$infos{$contig_c}};
	my $total_nb_genes=$#tab_genes+1;
	### Preparing data for C program
	my $out_file_c=$ref_file;
	$out_file_c=~s/\.refs/.tmp_$i/g;
	my $out_file_c2=$ref_file;
	$out_file_c2=~s/\.refs/.out_$i/g;
	my $out_file_c3=$ref_file;
	$out_file_c3=~s/\.refs/.out_$i-sorted/g;
# 	print "we have $out_file_c $out_file_c2 $out_file_c3\n";
	open MAP_C, '>', $out_file_c;
	print MAP_C "$nb_genes{$contig_c}\n";
	my $last_strand="0";
	my $total_hallmark=0;
	my $total_noncaudo=0;
	foreach(@tab_genes){
		my $gene=$_;
		my $tag="";
		# Line : PC / PFAM / UNCH / SIZE / STRAND / HALLMARK
		if($infos{$contig_c}{$gene}{"best_domain_hit"}=~/^PC/){
			if ($infos{$contig_c}{$gene}{"category"}>=3){$tag="1\t1\t0\t0\t";$total_noncaudo++;}
			else{$tag="1\t0\t0\t0\t";}
		}
		elsif($infos{$contig_c}{$gene}{"best_domain_hit"}=~/^PFAM/){$tag="0\t0\t1\t0\t";}
		else{$tag="0\t0\t0\t1\t";}
		if ($infos{$contig_c}{$gene}{"length"}<$th_gene_size){$tag.="1\t";}
		else{$tag.="0\t";}
		if (($last_strand eq "0") || ($infos{$contig_c}{$gene}{"strand"} eq $last_strand)){$tag.="0\t";}
		else{$tag.="1\t";}
		$last_strand=$infos{$contig_c}{$gene}{"strand"};
		if (($infos{$contig_c}{$gene}{"category"}==0) || ($infos{$contig_c}{$gene}{"category"}==3)){
			$tag.="1\t";$total_hallmark++;
			print "Gene $contig_c / $gene -> category $infos{$contig_c}{$gene}{category} -> putative hallmark\n";
		} # look at putative hallmarklmark
		else{$tag.="0\t";}
		print MAP_C "$tag\n";
	}
	close MAP_C;
	### Now go execute the C program
	my $c_cmd="$path_to_c_script $ref_file $out_file_c $out_file_c2";
#        print "Step 1 - $c_cmd\n";
	my $out=`$c_cmd`;
# 	print "$out\n";
	$c_cmd="sort -r -n -k 4 $out_file_c2 > $out_file_c3";
#        print "Step 2 - $c_cmd\n";	
	$out=`$c_cmd`;
# 	print "$out\n";
	### reading the c program output to fill the match hash table / and removing overlap
	my %match;
	my %check;
	my @check_gene;
	open OUT_C, '<', $out_file_c3;
	while(<OUT_C>){
		chomp($_);
		my @tab=split("\t",$_);
		my $start=$tab[0];
		my $last=$tab[0]+$tab[1]-1;
		my $fragment_id=$contig_c."-".$tab_genes[$start]."-".$tab_genes[$last];
		my $tag=0;
		# Code : 0 phage / 1 pfam / 2 unch / 3 size / 4 strand switch
		if ($tab[2]==0){
			if (overlap($fragment_id,$check{"phage"})==0){
				$match{$fragment_id}{"proof"}{"phage"}=$tab[3];
				$check{"phage"}{$fragment_id}=1;
				$tag=1;
				for (my $i=$start;$i<=$last;$i++){$check_gene[$i]++;}
			}
		}
		if ($tab[2]==1){
			if (overlap($fragment_id,$check{"pfam"})==0){
				$match{$fragment_id}{"proof"}{"pfam"}=$tab[3];
				$check{"pfam"}{$fragment_id}=1;
				$tag=1;
				for (my $i=$start;$i<=$last;$i++){$check_gene[$i]++;}
			}
		}
		if ($tab[2]==2){
			if (overlap($fragment_id,$check{"unch"})==0){
				$match{$fragment_id}{"proof"}{"unch"}=$tab[3];
				$check{"unch"}{$fragment_id}=1;
				$tag=1;
				for (my $i=$start;$i<=$last;$i++){$check_gene[$i]++;}
			}
		}
		if ($tab[2]==3){
			if (overlap($fragment_id,$check{"avg_g_size"})==0){
				$match{$fragment_id}{"proof"}{"avg_g_size"}=$tab[3];
				$check{"avg_g_size"}{$fragment_id}=1;
				$tag=1;
			}
		}
		if ($tab[2]==4){
			if (overlap($fragment_id,$check{"switch"})==0){
				$match{$fragment_id}{"proof"}{"switch"}=$tab[3];
				$check{"switch"}{$fragment_id}=1;
				$tag=1;
			}
		}
		if ($tab[2]==5){
			if (overlap($fragment_id,$check{"noncaudo"})==0){
				$match{$fragment_id}{"proof"}{"noncaudo"}=$tab[3];
				$check{"noncaudo"}{$fragment_id}=1;
				$tag=1;
				for (my $i=$start;$i<=$last;$i++){$check_gene[$i]++;}
			}
		}
		if ($tag==1){
			# If a match, we also take the nb of hallmark genes, and the size
			if ($tab[4]>0){$match{$fragment_id}{"hallmark"}=$tab[4];}
			$match{$fragment_id}{"size"}=$tab[1];
		}
	}
	close OUT_C;
	### Ok, we read the C output, no we try (neatly) to merge all predictions for this sequence
	my $n=0;
	my %merged_match;
	my $th_contig_size=$th_nb_genes_covered*$total_nb_genes;
	my @tab_matches=sort { $match{$b}{"size"} <=> $match{$a}{"size"} } keys %match;
	if (!defined($match{$tab_matches[0]}{"size"})){} # Not even an interesting region, skip to the next sequence
	else{
		my $tag_complete=0;
		my $i=0;
		while ($match{$tab_matches[$i]}{"size"}>$th_contig_size && $tag_complete==0){
			if ($match{$tab_matches[$i]}{"size"}>$th_contig_size && (defined($match{$tab_matches[$i]}{"proof"}{"pfam"}) || defined($match{$tab_matches[$i]}{"proof"}{"phage"}) || defined($match{$tab_matches[$i]}{"proof"}{"unch"}) || defined($match{$tab_matches[$i]}{"proof"}{"noncaudo"}))){ # SEEMS LIKE WE HAVE A COMPLETE PHAGE SEQUENCE 
				$tag_complete=1;
				my $fragment_id=$contig_c."-".$tab_genes[0]."-".$tab_genes[$#tab_genes];
				if (defined($match{$fragment_id})){
					$merged_match{$fragment_id}=$match{$fragment_id}; # If we indeed have complete metrics, we take themn
				}
				else{
					$merged_match{$fragment_id}{"size"}=$total_nb_genes;# Otherwise we store just the size
					$merged_match{$fragment_id}{"hallmark"}=$total_hallmark;# And the total number of hallmark genes on this fragment
				}
				$merged_match{$fragment_id}{"type"}="complete_phage";# And we store the type of fragment
				foreach(@tab_matches){
					my $fragment_id=$_;
					if ($match{$fragment_id}{"size"}<$total_nb_genes){
						my $r=get_overlap($fragment_id,\%merged_match);
						if ($r eq "no"){ # if no overlap
							$merged_match{$fragment_id}=$match{$fragment_id}; # NO OVERLAP WITH THE COMPLETE 
							print "!!!!!!!!!!!!!!!!!!! THIS SHOULD NOT BE POSSIBLE\n";
						}
						else{
							# Overlap, we propagate the proof and note it "partial"
							foreach(keys %{$match{$fragment_id}{"proof"}}){
								if (defined($merged_match{$r}{"proof"}{$_})){
									if ($merged_match{$r}{"proof"}{$_}=~/:/){
										$fragment_id=~/.*-(gene_\d*-gene_\d*)/;
										$merged_match{$r}{"proof"}{$_}.=$1.":".$match{$fragment_id}{"proof"}{$_}.",";
									}
									else{} # already a score for the entire match, no pblm
								} 
								else {
									$fragment_id=~/.*-(gene_\d*-gene_\d*)/;
									$merged_match{$r}{"proof"}{$_}=$1.":".$match{$fragment_id}{"proof"}{$_}.",";
								}
							}
						}
					}
				}
			}
			$i++;
		}
		if($tag_complete==0){ # No complete phage, putatively one or several prophages
			# First get all the phage region
			# We look for interesting regions   my $fragment_id=$contig_c."-".$tab_genes[0]."-".$tab_genes[$#tab_genes];
			my $tag=-1;
			my $tag_h=0;
			for (my $i=0;$i<$total_nb_genes;$i++){
				if ($tag>=0 && (!defined($check_gene[$i]) || $check_gene[$i]<1)){ # end of an interesting region
					my $fragment_id.=$contig_c."-".$tab_genes[$tag]."-".$tab_genes[$i-1];
					if ($merged_match{$fragment_id}{"size"}>$th_contig_size){ # Complete phage
						$fragment_id=$contig_c."-".$tab_genes[0]."-".$tab_genes[$#tab_genes];
						$merged_match{$fragment_id}{"type"}="complete_phage";
						$merged_match{$fragment_id}{"size"}=$total_nb_genes;
						$merged_match{$fragment_id}{"hallmark"}=$tag_h;
					} 
					else{ # Prophage
						$merged_match{$fragment_id}{"size"}=$i-$tag;
						$merged_match{$fragment_id}{"type"}="prophage";
						$merged_match{$fragment_id}{"hallmark"}=$tag_h;
					}
					$tag=-1;
					$tag_h=0;
				}
				elsif ($tag==-1 && $check_gene[$i]>=1){
					$tag=$i;
					$tag_h=0;
				}
				if ($infos{$contig_c}{$tab_genes[$i]}{"category"}==0 || $infos{$contig_c}{$tab_genes[$i]}{"category"}==3){$tag_h++;} # look at putative hallmark
			}
			if ($tag>=0){
				my $fragment_id.=$contig_c."-".$tab_genes[$tag]."-".$tab_genes[$#tab_genes];
				print "Region is $fragment_id ..";
				if ($merged_match{$fragment_id}{"size"}>$th_contig_size){ # Complete phage
					print "which is a complete phage\n";
					$fragment_id=$contig_c."-".$tab_genes[0]."-".$tab_genes[$#tab_genes];
					$merged_match{$fragment_id}{"type"}="complete_phage";
					$merged_match{$fragment_id}{"size"}=$total_nb_genes;
					$merged_match{$fragment_id}{"hallmark"}=$tag_h;
				} 
				else{ # Prophage
					print "which is a prophage\n";
					$merged_match{$fragment_id}{"size"}=$total_nb_genes-$tag;
					$merged_match{$fragment_id}{"type"}="prophage";
					$merged_match{$fragment_id}{"hallmark"}=$tag_h;
				}
			}
			# Now we merge the annotation in these regions
			foreach(@tab_matches){
				my $fragment_id=$_;
				# Check if overlap
				my $r=get_overlap($fragment_id,\%merged_match);
				if ($r eq "no"){ } # if no overlap # not in an interesting region 
				else{
					# Overlap, we propagate the proof and note it "partial"
					foreach(keys %{$match{$fragment_id}{"proof"}}){
						if (defined($merged_match{$r}{"proof"}{$_})){
							if ($merged_match{$r}{"proof"}{$_}=~/:/){
								$fragment_id=~/.*-(gene_\d*-gene_\d*)/;
								$merged_match{$r}{"proof"}{$_}.=$1.":".$match{$fragment_id}{"proof"}{$_}.",";
							}
							else{} # already a score for the entire match, no pblm
						} 
						else {
							$fragment_id=~/.*-(gene_\d*-gene_\d*)/;
							$merged_match{$r}{"proof"}{$_}=$1.":".$match{$fragment_id}{"proof"}{$_}.",";
						}
					}
					delete($match{$fragment_id});
				}
			}
			## New addition that should help to get the prophage coordinates correctly !
			# And now check if one of the prophage map to the whole sequence
			foreach(keys %merged_match){
				print "This is a prophage\n";
				my $fragment_id=$_;
				if ($merged_match{$fragment_id}{"size"}>$th_contig_size){
					$tag_complete=1;
					my $new_fragment_id=$contig_c."-".$tab_genes[0]."-".$tab_genes[$#tab_genes];
					print "We have a complete prophage -- we add it $new_fragment_id !\n";
					# <STDIN>;
					foreach(keys %{$merged_match{$fragment_id}}){
						$merged_match{$new_fragment_id}{$_}=$merged_match{$fragment_id}{$_};
					}
					$merged_match{$new_fragment_id}{"type"}="complete_phage";# And we store the type of fragment
				}
			}
			if ($tag_complete==1){
				# We can remove all the prophages
				my @tab_temp=keys %merged_match;
				foreach(@tab_temp){
					if ($merged_match{$_}{"type"} eq "complete_phage"){}
					else{
						delete($merged_match{$_});
					}
				}
			}
			## END OF THE NEW ADDITION
		}
		open S1, '>>', $out_file;
		foreach(sort { $merged_match{$b}{"size"} <=> $merged_match{$a}{"size"} } keys %merged_match){ ## IMPORTANT, HAVE TO BE SIZE ORDERED
			my $fragment_id=$_;
			$fragment_id=~/.*-(gene_\d+-gene_\d+)/;
			my $zone=$1;
			my $type_detection=$merged_match{$fragment_id}{"type"};
			print "$fragment_id\t$merged_match{$fragment_id}{size}\t$merged_match{$fragment_id}{hallmark}\t$merged_match{$fragment_id}{proof}{phage}\t$merged_match{$fragment_id}{proof}{pfam}\t$merged_match{$fragment_id}{proof}{unch}\t$merged_match{$fragment_id}{proof}{switch}\t$merged_match{$fragment_id}{proof}{avg_g_size}\n";
			my $category=3;
			if ($merged_match{$fragment_id}{"hallmark"}==0){delete($merged_match{$fragment_id}{"hallmark"});}
			# Determine the category. To check this, we want several good indicators - And also remove prediction based on one single indicator, unless it's a strong one (sig >2)
			# New categories : 
			# Cat 1 - hallmark + gene phage enrichment
			# Cat 2 - gene phage or hallmark without gene phage
			# Cat 3 - no hallmark or gene phage, but other signal
			my @tab_proof=keys %{$merged_match{$fragment_id}{"proof"}};
			if ($merged_match{$fragment_id}{"hallmark"}>0){
				if (defined($merged_match{$fragment_id}{"proof"}{"noncaudo"}) || defined($merged_match{$fragment_id}{"proof"}{"phage"})){
					if ($merged_match{$fragment_id}{"proof"}{"noncaudo"}=~/(gene_\d+-gene_\d+):(\d+)/){
						my $match_region=$1;
						my $score=$2;
						if ($match_region eq $zone && $score>=$th_sig){$category=1;} # Phage metric on the whole region
					}
					elsif ($merged_match{$fragment_id}{"proof"}{"noncaudo"}>=$th_sig){$category=1;} # if we have hallmark or gene_size + a phage metric on the whole fragment -> should be quite sure $category=1; ## THRESHOLD TO REMOVE THE NONCAUDO ON THE SMALL SMALL CONTIGS
					if ($merged_match{$fragment_id}{"proof"}{"phage"}=~/(gene_\d+-gene_\d+):(\d+)/){
						my $match_region=$1;
						my $score=$2;
						if ($match_region eq $zone && $score>=$th_sig){$category=1;} # Phage metric on the whole region
					}
					elsif ($merged_match{$fragment_id}{"proof"}{"phage"}>=$th_sig){$category=1;} # if we have hallmark or gene_size + a phage metric on the whole fragment -> should be quite sure $category=1;
					if ($category==3){ # no match complete, so category 2
						$category=2;
					}
				}
				else{
					foreach(@tab_proof){
						if ($merged_match{$fragment_id}{"proof"}{$_}=~/(gene_\d+-gene_\d+):(\d+)/){
							my $match_region=$1;
							my $score=$2;
							print "Hallmark but no phage or noncaudo, but other proof $_ -> $match_region / $score ($merged_match{$fragment_id}{proof}{$_})\n";
							if ($match_region eq $zone && $score>=$th_sig){$category=2;} # other metric on the whole region
							elsif($score>=$th_sig_2){$category=2;} # metric partial only but strong enough so we keep it
						}
						elsif ($merged_match{$fragment_id}{"proof"}{$_}>=$th_sig){
							if ($_ eq "pfam" || $_ eq "unch"){
								$category=2; # if we have hallmark or gene_size + a metric pfam or unch on the whole fragment -> should be quite sure
							}
						}
					}
				}
			}
			elsif (defined($merged_match{$fragment_id}{"proof"}{"phage"}) || defined($merged_match{$fragment_id}{"proof"}{"noncaudo"})){# If we have some phage signal, 
				if ($merged_match{$fragment_id}{"proof"}{"phage"}=~/:(\d*)/){
					if ($1>=$th_sig){
						$category=2; # Good, phage signal significant -> should be quite sure
					}
				} 
				elsif($merged_match{$fragment_id}{"proof"}{"phage"}>=$th_sig){
					$category=2; # Good, phage signal significant -> should be quite sure
				}
				if ($merged_match{$fragment_id}{"proof"}{"noncaudo"}=~/:(\d*)/){ ## THRESHOLD TO AVOID SHORT CONTIGS BIAS
					if ($1>=$th_sig && $total_noncaudo>$th_nb_genes_noncaudo){
						$category=2; # Good, phage signal significant -> should be quite sure
					}
				} 
				elsif($merged_match{$fragment_id}{"proof"}{"noncaudo"}>=$th_sig && $total_noncaudo>$th_nb_genes_noncaudo){ ## THRESHOLD TO AVOID SHORT CONTIGS BIAS
					$category=2; # Good, phage signal significant -> should be quite sure
				}
			}
			if ($category==3){ # If the category is still 3, meaning that the phage signal (if there was any) was not that strong ..
				if ($#tab_proof==0){
					$category=0; # No phage signal nor hallmark gene, and only one metric, we remove
				}
				else{
					my $tag1=0;
					foreach(@tab_proof){
						if ($merged_match{$fragment_id}{"proof"}{$_}=~/:(\d*)/){
							if ($1>=$th_sig_2){
								$tag1=1; # Good, one signal very significant 
							}
						} 
						elsif ($merged_match{$fragment_id}{"proof"}{$_}>=$th_sig_2){
							$tag1=1; # Good, one signal very significant 
						}
					}
					if ($tag1==0){ # If none of the metrics is really strong ...
						$category=0; # .. we remove the detection
					}
				}
			}
			# Columns index :  0  /     1        /     2    /   3   /      4        /    5     /     6        /      7            /    8         /      9       /     10     /    11     /       12
			# Columns : Contig / Total Nb Genes / Fragment / Size / Type detection / Category /  Enrich Phage / Enrich Non Caudo / Enrich Pfam / Enrich Unch / Enrich Switch / Avg_g_size / Nb Hallmark
			if ($category>0){
				print S1 "$contig_c\t$total_nb_genes\t$fragment_id\t$merged_match{$fragment_id}{size}\t$type_detection\t$category\t$merged_match{$fragment_id}{proof}{phage}\t$merged_match{$fragment_id}{proof}{noncaudo}\t$merged_match{$fragment_id}{proof}{pfam}\t$merged_match{$fragment_id}{proof}{unch}\t$merged_match{$fragment_id}{proof}{switch}\t$merged_match{$fragment_id}{proof}{avg_g_size}\t$merged_match{$fragment_id}{hallmark}\n";
			}
		}
		close S1;
	}
	$i++;
	`rm $out_file_c $out_file_c2 $out_file_c3`;
}

sub factorial { # factorial $n
	my $n = shift;
	my $f = 1;
	$f *= $n-- while $n > 0;    # Multiply, then decrement
	return $f;
}

sub combine { # combination of $k elements in $n ensemble
	my $k=$_[0];
	my $n=$_[1];
	my $f=factorial($n)/(factorial($k) * factorial($n-$k));
	return $f;
}

sub proba { # probability of x=$i knowing nb_obs $n and p $p
	my $i=$_[0];
	my $n=$_[1];
	my $p=$_[2];
	my $f=combine($i,$n)*($p**$i)*((1-$p)**($n-$i));
	return $f;
}

sub proba_more_than { # probability of x>=$s knowing nb_obs $n and p $p
	my $s=$_[0];
	my $n=$_[1];
	my $p=$_[2];
	my $f=0;
	for (my $i=$s;$i<=$n;$i++){
		$f+=proba($i,$n,$p);
	}
	return $f;
}

sub proba_less_than { # probability of x<=$s knowing nb_obs $n and p $p
	my $s=$_[0];
	my $n=$_[1];
	my $p=$_[2];
	my $f=0;
	for (my $i=0;$i<=$s;$i++){
		$f+=proba($i,$n,$p);
	}
	return $f;
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

sub overlap { # To check if a prediction is not within another of the same type, in which case we don't really care
	my $pred=$_[0];
	$pred=~/.*-gene_(\d*)-gene_(\d*)/;
	my $start_pred=$1; my $end_pred=$2;
	my $p_hash=$_[1];
	my $o=0;
	foreach(keys %{$p_hash}){
		$_=~/.*-gene_(\d*)-gene_(\d*)/;
		if (($start_pred<=$1 && $1<$end_pred) || ($start_pred<$2 && $2<=$end_pred) || ($1<=$start_pred && $2>=$end_pred)){
			$o=1;
		}
	}
	return $o;
}


sub get_overlap { # To get the overlapping if any
	my $pred=$_[0];
	$pred=~/.*-gene_(\d*)-gene_(\d*)/;
	my $start_pred=$1; my $end_pred=$2;
	my $p_hash=$_[1];
	my $o="no";
	foreach(keys %{$p_hash}){
		$_=~/.*-gene_(\d*)-gene_(\d*)/;
		if ($start_pred>=$1 && $end_pred<=$2){
			$o=$_;
		}
	}
	return $o;
}

sub is_local_max {
	my $p_metrics=$_[0];
	my $s=$_[1];
	my $w=$_[2];
	my $c=$_[3];
	my $v=$$p_metrics{$s}{$w}{$c};
	my $f=1;
	my $how_much_to_look=5;
	for (my $i=-$how_much_to_look;$i<=$how_much_to_look;$i++){
		for (my $j=-$how_much_to_look;$j<=$how_much_to_look;$j++){
			if ($i==0 && $j==0){}
			elsif(defined($$p_metrics{$s+$i}{$w+$j})){
				if(defined($$p_metrics{$s+$i}{$w+$j}{$c})){
					if ($$p_metrics{$s+$i}{$w+$j}{$c}>$v){
						$f=0; # we found a neightbor with a greater value, not a local maxima
					}
				}
			}
		}
	}
	return $f;
}

sub get_position{
	my $value=$_[0];
	my $p_tab=$_[1];
	my @tab=sort {$a <=> $b} @$p_tab;
	print "looking for $value in the gene size table\n";
	my $index = 0;
	while($tab[$index]<$value && $index<$#tab){$index++;}
	print "found at index $index - $tab[$index] - total : $#tab\n";
	my $ratio=$index/$#tab;
	print "which gives a ratio of $ratio\n";
	return $ratio;
}

sub get_th_gene_size{
	my @tab=sort {$a <=> $b} (@{$_[0]});
	my $div=$_[1];
	my $m=0;
	if ($#tab % $div == 0){return ($tab[$#tab/$div]);}
	else{return (($tab[($#tab-1)/$div]+$tab[($#tab+1)/$div])/2);}
}

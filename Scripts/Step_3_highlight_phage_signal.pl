#!/usr/bin/env perl

use strict;
use autodie;
use File::Spec::Functions;
use FindBin '$Bin';
use Parallel::ForkManager;
use Getopt::Long 'GetOptions';
# use List::MoreUtils qw(natatime);
# Script to measure metrics on the sliding window
my $csv_file="";
my $out_file="";
my $n_cpus=1;
my $ref_file="";
my $external_ref="";
my $no_c=0;

GetOptions(
   'csv=s'     => \$csv_file, # Csv files of the contigs
   'out=s'     => \$out_file, # Summary file of virus fragments (output)
   'n_cpu=i'   => \$n_cpus, # Number of threads
   'ref=s'     => \$external_ref, # Reference file (optional)
   'no_c=i'    => \$no_c # Tag to use the perl version rather than the C version
);

if ($csv_file eq "" || $out_file eq ""){
	print "# Script to measure metrics on the sliding window
# -csv : csv file of the contigs
# -out : summary file of the phage fragments
# -n_cpu : number of CPUs to use in parallel processing of sliding window analysis
# -ref : a file with the refs values that we could use instead of estimating them 
# -no_c : 0 for using the C program, 1 for using perl\n";
	die "\n";
}
$| = 1;

if ( -e $out_file ) { `rm $out_file`; }

# Deal with the ref file depending if we are in a virome decontamination or not
$ref_file=$csv_file;
$ref_file =~ s/\.csv/.refs/g;
my $do_ref_estimation = 0;
if ($external_ref ne ""){ # If we specify a ref file, that's because we are in a virome decontamination mode
	`cp $external_ref $ref_file`;
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
# 		print "## First look at $contig\n";
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

my $num_contigs = scalar(@liste_contigs);
my $chunk_len = int($num_contigs/$n_cpus)+1;
#print "$chunk_len\n";
my @grouped = @liste_contigs;
my $it = &natatime($chunk_len, @grouped);

my $pm = Parallel::ForkManager->new($n_cpus); #Starts the parent process for parallelizing the next foreach loop, sets max number of parallel processes
#$pm->set_waitpid_blocking_sleep(0);
while (my @vals = $it->()){
	$pm->start and next; #do the fork
	foreach(@vals){
		my $contig_c=$_;
# 		print "## Looking for the second time at $contig_c\n";
		my @tab_genes=sort {$infos{$contig_c}{$a}{"order"} <=> $infos{$contig_c}{$b}{"order"}} keys %{$infos{$contig_c}};
		my $total_nb_genes=$#tab_genes+1;
		### Preparing data for C program
		my $out_file_c=$ref_file;
		$out_file_c=~s/\.refs/.tmp_$$/g; #The variable $$ is the PID for the fork
		my $out_file_c2=$ref_file;
		$out_file_c2=~s/\.refs/.out_$$/g;
		my $out_file_c3=$ref_file;
		$out_file_c3=~s/\.refs/.out_$$-sorted/g;
		my $line_input="$ref_file\n$nb_genes{$contig_c}\n"; # We start a "line" instead, so put everything in a string ($line_input)
		my $last_strand="0";
		my $total_hallmark=0;
		my $total_noncaudo=0;
		foreach(@tab_genes){
			my $gene=$_;
			my $tag="";
			# Line : PC / noncaudo / PFAM / UNCH / SIZE / STRAND / HALLMARK
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
# 				print "Gene $contig_c / $gene -> category $infos{$contig_c}{$gene}{category} -> putative hallmark\n";
			} # look at putative hallmarklmark
			else{$tag.="0\t";}
			$line_input.=$tag."\n";
		}
		
		my @tab_lines;
		if ($no_c == 1){
			@tab_lines=&get_metrics($line_input);
		}
		else{
			### reading the c program output to fill the match hash table / and removing overlap
			## Now go execute the C program
			print $line_input."\n";
			my $path_to_c_script= catfile($script_dir, "Sliding_windows_3");
			my $c_cmd="printf \"$line_input\" | $path_to_c_script | sort -r -n -k 4 ";
			print $c_cmd."\n";
			my $out=`$c_cmd`;
# 			print "### RESULT FROM C: $out\n";
			@tab_lines=split("\n",$out);
		}
		
		
		my %match;
		my %check;
		my @check_gene;
		foreach my $line (@tab_lines){
			print $line."\n";
			my @tab=split("\t",$line);
			if ($tab[3]<2){next;} # Just in case
			my $start=$tab[0];
			my $last=$tab[0]+$tab[1]-1;
			my $fragment_id=$contig_c."-".$tab_genes[$start]."-".$tab_genes[$last];
			my $tag=0;
			# Code : 0 phage / 1 pfam / 2 unch / 3 size / 4 strand switch
			if ($tab[2]==0){
				if (overlap($fragment_id,$check{"phage"})==0 || $tab[1]==100){ # Exception - we check if the window size is as large as the max window, i.e. 100 genes, regardless of overlap
					$match{$fragment_id}{"proof"}{"phage"}=$tab[3];
					$check{"phage"}{$fragment_id}=1;
					$tag=1;
					for (my $i=$start;$i<=$last;$i++){$check_gene[$i]++;}
				}
			}
			if ($tab[2]==1 || $tab[1]==100){
				if (overlap($fragment_id,$check{"pfam"})==0){
					$match{$fragment_id}{"proof"}{"pfam"}=$tab[3];
					$check{"pfam"}{$fragment_id}=1;
					$tag=1;
					for (my $i=$start;$i<=$last;$i++){$check_gene[$i]++;}
				}
			}
			if ($tab[2]==2 || $tab[1]==100){
				if (overlap($fragment_id,$check{"unch"})==0){
					$match{$fragment_id}{"proof"}{"unch"}=$tab[3];
					$check{"unch"}{$fragment_id}=1;
					$tag=1;
					for (my $i=$start;$i<=$last;$i++){$check_gene[$i]++;}
				}
			}
			if ($tab[2]==3 || $tab[1]==100){
				if (overlap($fragment_id,$check{"avg_g_size"})==0){
					$match{$fragment_id}{"proof"}{"avg_g_size"}=$tab[3];
					$check{"avg_g_size"}{$fragment_id}=1;
					$tag=1;
				}
			}
			if ($tab[2]==4 || $tab[1]==100){
				if (overlap($fragment_id,$check{"switch"})==0){
					$match{$fragment_id}{"proof"}{"switch"}=$tab[3];
					$check{"switch"}{$fragment_id}=1;
					$tag=1;
				}
			}
			if ($tab[2]==5 || $tab[1]==100){
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
	#	close OUT_C;
		### Ok, we read the C/Perl output, now we try (neatly) to merge all predictions for this sequence
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
								print "Overlap -- we propagate the proof and not it partial\n";
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
						print "We have a complete genome -- we add it $new_fragment_id !\n";
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
	#	$i++;
		#`rm $out_file_c $out_file_c2 $out_file_c3`;
	}
	$pm->finish(0); # do the exit in the child process
}
$pm->wait_all_children; # wait until everything in the above foreach loop is done before moving on

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
	if ($n<0){print "########## PBLM WE ARE TRYING TO TAKE A LOG OF A NEGATIVE \n";}
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
# 	print "looking for $value in the gene size table\n";
	my $index = 0;
	while($tab[$index]<$value && $index<$#tab){$index++;}
# 	print "found at index $index - $tab[$index] - total : $#tab\n";
	my $ratio=$index/$#tab;
# 	print "which gives a ratio of $ratio\n";
	return $ratio;
}

sub get_th_gene_size{
	my @tab=sort {$a <=> $b} (@{$_[0]});
	my $div=$_[1];
	my $m=0;
	if ($#tab % $div == 0){return ($tab[$#tab/$div]);}
	else{return (($tab[($#tab-1)/$div]+$tab[($#tab+1)/$div])/2);}
}

sub natatime ($@)
{
	my $n    = shift;
	my @list = @_;
	return sub {
		return splice @list, 0, $n;
	}
}


sub get_th{
	my $sw=$_[0];
	my $corr=$_[1];
	my $proba=$_[2];
	my $th_nb_gene=$sw+1;
# 	print "Sliding window $sw - Th nb genes starting at ".$th_nb_gene." with proba ".$proba." - we aim for $corr\n";
	my $p_t=0;
	while($p_t<=$corr && $th_nb_gene>0){
		$th_nb_gene--;
		$p_t+=&proba($th_nb_gene,$sw,$proba);
# 		print "\tp(x>=$th_nb_gene) = $p_t\n";
	}
# 	print "Final is -> $th_nb_gene, which has proba $p_t\n";
	return($th_nb_gene);
}

sub get_th_less{
	my $sw=$_[0];
	my $corr=$_[1];
	my $proba=$_[2];
	my $th_nb_gene=0;
# 	print "Sliding window $sw - Th nb genes starting at ".$th_nb_gene." with proba ".$proba." - we aim for $corr\n";
	my $p_t=0;
	while($p_t<=$corr && $th_nb_gene<=$sw){
		$th_nb_gene++;
		$p_t+=&proba($th_nb_gene,$sw,$proba);
# 		print "\tp(x>=$th_nb_gene) = $p_t\n";
	}
# 	print "Final is -> $th_nb_gene, which has proba $p_t\n";
	return($th_nb_gene);
}


sub is_local_maximum(int start,int size,int type, int p_nb_genes, int p_max,double ***store){
	my $start=$_[0];
	my $size=$_[1];
	my $cat=$_[2];
	my $nb_genes=$_[3];
	my $p_max=$_[4];
	my $store=$_[5];
	my $result=1;
	my $hood=5;
	for (my $i=$start-$hood;$i<=$start+$hood;$i++){
		for (my $j=$size-$hood;$j<=$size+$hood;$j++){
# 			print "-- Looking at $i - $j\n";
			if ($i>=0 && $j>=0 && $i<$nb_genes && $j<=$p_max){
# 				print "-- Really Looking at $i - $j\n";
				if ($$store{$i}{$j}{$cat}>$$store{$start}{$size}{$cat}){
					$result=0;
					$i=$start+$hood+1;
					$j=$size+$hood+1;
				}
			}
		}
	}
	return $result;
}



sub get_metrics{
	my @lines=split("\n",$_[0]);
	my $file_ref=shift(@lines);
	my $nb_genes=shift(@lines);
	# Load ref probabilities
# 	print "Load ref probabilities\n";
	my @tab_categories=("phage","pfam","unch","size","strand","noncaudo","hallmark");
	my %p;
	open my $tsv,"<",$file_ref;
	my @ref_line=split(" ",<$tsv>);
	$p{"phage"}=$ref_line[0];
	$p{"pfam"}=$ref_line[1];
	$p{"unch"}=$ref_line[2];
	$p{"strand"}=$ref_line[3];
	$p{"size"}=0.1; ## This is the true probability, the number in the file is for reference, the size of genes considered as "short", which corresponds to the first decile, i.e. 10% chance to be that size or shorter
	$p{"noncaudo"}=$ref_line[5];
	# Load sequence information
# 	print "Load sequence information\n";
	my %store;
	for(my $i=0;$i<=$#lines;$i++){
		my @t=split("\t",$lines[$i]);
		$store{$i}{"phage"}=$t[0];
		$store{$i}{"noncaudo"}=$t[1];
		$store{$i}{"pfam"}=$t[2];
		$store{$i}{"unch"}=$t[3];
		$store{$i}{"size"}=$t[4];
		$store{$i}{"strand"}=$t[5];
		$store{$i}{"hallmark"}=$t[6];
	}
	# Counting number of sliding windows for correction
# 	print "Counting sliding windows\n";
	my $nb_sw;
	my $min=10;
	my $max=100;
	if ($min>$nb_genes){$min=$nb_genes;}
	if ($max>$nb_genes){$max=$nb_genes;}
	for(my $k=$max;$k>=$min;$k--){$nb_sw+=$nb_genes-$k+1;}
	my $corr=0.01/$nb_sw;
	# Now do sliding windows
# 	print "Going through sliding windows\n";
	my %store_result;
	my %th;
	my %nb;
	my $tag=0;
	for(my $k=$max;$k>=$min;$k--){
# 		print "\twindow of size $k\n";
		%th=();
		$th{"phage"}=&get_th($k,$corr,$p{"phage"});
		$th{"pfam"}=&get_th_less($k,$corr,$p{"pfam"});
		$th{"unch"}=&get_th($k,$corr,$p{"unch"});
		$th{"size"}=&get_th($k,$corr,0.1);
		$th{"strand"}=&get_th_less($k,$corr,$p{"strand"});
		$th{"noncaudo"}=&get_th($k,$corr,$p{"noncaudo"});
		%nb=();
		for (my $i=0;$i<($nb_genes-$k+1);$i++){
# 			print "\t\tstarting at position $i\n";
			$tag=0;
			if ($i==0){
				# initialize first window
				for (my $j=0;$j<$k;$j++){
					foreach my $cat (@tab_categories){
						$nb{$cat}+=$store{$j}{$cat};
					}
				}
			}
			else{ # shift window
				my $dropped=$i-1;
				my $last=$i+$k-1;
				foreach my $cat (@tab_categories){
					$nb{$cat}-=$store{$dropped}{$cat};
					$nb{$cat}+=$store{$last}{$cat};
				}
			}
			if ($nb{"phage"}>$th{"phage"}){
# 				print "Testing phage - $nb{phage} == $th{phage} - $k / $p{phage} / $nb_sw\n";
				$store_result{$i}{$k}{"phage"}= -1 * &log10(&proba_more_than($nb{"phage"},$k,$p{"phage"})*$nb_sw); $tag=1;
			}
			if ($nb{"pfam"}<$th{"pfam"}){
# 				print "Testing pfam - $nb{pfam} == $th{pfam} - $k / $p{pfam} / $nb_sw\n";
				$store_result{$i}{$k}{"pfam"}= -1 * &log10(&proba_less_than($nb{"pfam"},$k,$p{"pfam"})*$nb_sw); $tag=1;
			}
			if ($nb{"unch"}>$th{"unch"}){
# 				print "Testing unch - $nb{unch} == $th{unch} - $k / $p{unch} / $nb_sw\n";
				$store_result{$i}{$k}{"unch"}= -1 * &log10(&proba_more_than($nb{"unch"},$k,$p{"unch"})*$nb_sw); $tag=1;
			}
			if ($nb{"size"}>$th{"size"}){
# 				print "Testing size - $nb{size} == $th{size} - $k / $p{size} / $nb_sw\n";
				$store_result{$i}{$k}{"size"}= -1 * &log10(&proba_more_than($nb{"size"},$k,$p{"size"})*$nb_sw); $tag=1;
			}
			if ($nb{"strand"}<$th{"strand"}){
# 				print "Testing strand - $nb{strand} == $th{strand} - $k / $p{strand} / $nb_sw\n";
				$store_result{$i}{$k}{"strand"}= -1 * &log10(&proba_less_than($nb{"strand"},$k,$p{"strand"})*$nb_sw); $tag=1;
			}
			if ($nb{"noncaudo"}>$th{"noncaudo"}){
# 				print "Testing noncaudo - $nb{noncaudo} == $th{noncaudo} - $k / $p{noncaudo} / $nb_sw\n";
				$store_result{$i}{$k}{"noncaudo"}= -1 * &log10(&proba_more_than($nb{"noncaudo"},$k,$p{"noncaudo"})*$nb_sw); $tag=1;
			}
			if ($tag==1){
				$store_result{$i}{$k}{"hallmark"}=$nb{"hallmark"};
			}
		}
	}
	my %store_output;
	my $n=0;
	for (my $k=$max;$k>=$min;$k--){
		for ($i=0;$i<($nb_genes-$k+1);$i++){
			for (my $j=0;$j<6;$j++){
				my $cat=$tab_categories[$j];
				if ($store_result{$i}{$k}{$cat} >= 2){ # the stored value is worth looking at
#  					print "potential local maximum ".$i." ".$k." ".$cat." ".$store{$i}{$k}{$cat}." ".$store_result{$i}{$k}{"hallmark"}."\n";
					if (&is_local_maximum($i,$k,$cat,$nb_genes-1,$max,\%store_result)==1){ # and is a local maximum
						# so we print it, with the nb_hallmark (start / window size / type / sig / nb_hallmark)
#  						print "local maximum ! $i $k $j ($cat) $store_result{$i}{$k}{$cat}\n";
 						$n++;
 						$store_output{$n}{"score"}=$store_result{$i}{$k}{$cat};
 						$store_output{$n}{"start"}=$i;
 						$store_output{$n}{"line"}=$i."\t".$k."\t".$j."\t".sprintf("%.14lf",$store_result{$i}{$k}{$cat})."\t".$store_result{$i}{$k}{"hallmark"};
					}
				}
			}
		}
	}
	my @tab_out;
	foreach my $key (sort {$store_output{$b}{"score"} <=> $store_output{$a}{"score"} || $store_output{$b}{"start"} <=> $store_output{$a}{"start"}} keys %store_output){
		push(@tab_out,$store_output{$key}{"line"});
	}
	return (@tab_out);
}

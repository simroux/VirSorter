#!/usr/bin/perl
use strict;
use lib '/usr/local/bin/Virsorter/Tools/BioPerl-1.6.1';
use Bio::Seq;
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


my $id=$ARGV[0];
my $tmp_dir=$ARGV[1];
my $fasta_contigs=$ARGV[2];
my $th_nb_genes=$ARGV[3];

my $path_to_mga="/usr/local/bin/Virsorter/Tools/Metagene_annotator/mga_linux_ia64";

my $in_file=$tmp_dir."/".$id."_nett.fasta";
my $circu_file=$tmp_dir."/".$id."_circu.list";
my $out_special_circu=$tmp_dir."/".$id."_contigs_circu_temp.fasta";

open(F1,"<$fasta_contigs") || die "pblm ouverture fichier $fasta_contigs\n";
my %seq_base;
my $id_seq="";
while(<F1>){
	$_=~s/\r\n/\n/g; #Cas d'un fichier windows ##AJOUT
	chomp($_);
	if ($_=~/^>(\S*)/){$id_seq=$1;}
	else{$seq_base{$id_seq}.=$_;}
}
close F1;

## DETECTION OF CIRCULAR CONTIG AND CLEANING OF THESE CIRCULAR (REMOVE THE MATCHING ENDS)
my $minimum_size=1500;
my %order_contig;
my %length;
my $n1=0;
open(S1,">$in_file") || die "pblm ouverture fichier $in_file\n";
open(S2,">$circu_file") || die "pblm ouverture fichier $circu_file\n";
foreach(sort {length($seq_base{$b}) <=> length($seq_base{$a})} keys %seq_base){
	my $id_contig=$_;
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
			print S2 "$id_contig\t$length{$id_contig}\n";
			$seq_base{$id_contig}=$sequence;
		}
	}
	# Update the length of the contig
	$length{$id_contig}=length($seq_base{$id_contig});
	print S1 ">$id_contig\n$seq_base{$id_contig}\n";
}
close S1;
close S2;

# Gene prediction for all contigs
my $out_file= $tmp_dir."/".$id."_mga.predict";
print "$path_to_mga $in_file -m > $out_file\n";
my $mga=`$path_to_mga $in_file -m > $out_file`;

# Special prediction for circular contigs
my $out_file_circu="";
my %circu;
if (-e $circu_file){
	open(CI,"<$circu_file") || die "pblm ouverture fichier $circu_file\n";
	while(<CI>){
		chomp($_);
		my @tab=split("\t",$_);
		my $id_c=$tab[0];
		$circu{$id_c}=1;
	}
	close CI;
	open(S2,">$out_special_circu") || die "pblm ouverture fichier $out_special_circu\n";
	my $long=1000; # we cp the 1000 first bases to the end of the contig
	my $seuil_long=1000;
	my $n_circu=0;
	foreach(sort {$order_contig{$a} <=> $order_contig{$b} } keys %circu){
		my $id_c=$_;
		my $s=$seq_base{$id_c}.substr($seq_base{$id_c},0,$long);
		print S2 ">$id_c\n$s\n";
		$n_circu++;
	}
	close S2;
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
open(RESU,"<$out_file")  || die "pblm ouverture fichier $out_file\n";
my %predict;
my %type;
my $id_c="";
while(<RESU>){
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
close RESU;
if (-e $circu_file){
	open(CIR,"<$out_file_circu") || die "pblm ouverture fichier $out_file_circu\n";
	my $tag=0;
	while(<CIR>){
		chomp($_);
		if($_=~/^# gc/){}
		elsif($_=~/^# self: (.*)/){$type{$id_c}=$1;}
		elsif ($_=~/^# (.*)/){
			if($tag==1){
				my %to_start;
# 				print "on a remplacé certains orfs pour $id_c, on fait du ménage en cherchant des gènes qui seraient partiellement prédits, soit au début, soit à la fin du contig\n";
				foreach(sort {$order_gene{$a} <=> $order_gene{$b} } keys %{$predict{$id_c}}){
					my @tab=split("\t",$predict{$id_c}{$_});
					if ($tab[5]!=11){
# 						print "\t$tab[0] n'a pas le(s) codon(s) start et / ou stop\n";
						if(($tab[1]<3) || ($tab[2]>($length{$id_c}-3))){
# 							print "\tet il est autour de l'origine, on enlève\n";
							if ($tab[1]<3){
								$to_start{$tab[0]}{"start"}=$tab[1];
								$to_start{$tab[0]}{"stop"}=$tab[2];
							}
							delete($predict{$id_c}{$tab[0]});
						}
						elsif(($tab[2]>997) && ($tab[2]<1001)){ # si on est dans la zone des 1000
							foreach(keys %to_start){
								my $total=($length{$id_c}-$tab[1]+1)+($to_start{$_}{"stop"}); # longueur du fragment terminal + longueur du fragment initial potentiel - pour ce deuxième, on omet le -1 +1 puisque de toute facon ce fragment commencera forcèment en +1
								if ($total % 3 == 0){
		#							print "\t et on peut l'étendre grâce à une autre prédiction : -> $to_start{$_}{stop}\n";
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
# 			print "## on passe à $id_c\n";
			$tag=0;
		}
		else{
			my @tab=split("\t",$_);
			if (defined($predict{$id_c}{$tab[0]})){
				my @tab2=split("\t",$predict{$id_c}{$tab[0]});
				if (($tab2[1]==$tab[1]) && ($tab2[2]==$tab[2])){}#print "même prédiction pour $tab[0], on change rien\n";}
				else{
# 					print "prédiction différente pour $tab[0] de $id_c, on regarde si on est dans la zone qui nous interesse $length{$id_c} == $tab[1] $tab[2] vs $tab2[1]  $tab2[2]\n";
					# si on chevauche l'origine
					if (($tab[1]<$length{$id_c}) && ($tab[2]>$length{$id_c})){
# 						print "on chevauche l'origine ! on remplace la prédiction\n";
						$tag=1;
						my $stop=$tab[2]-$length{$id_c};
						$tab[2]=$stop;
						my $new_line=join("\t",@tab);
						$predict{$id_c}{$tab[0]}=$new_line;
					}
				}
			}
			else{
				# print "on prédit un nouveau gène, qu'on ne garde que si il débute avant la fin du vrai contig\n";
				if (($tab[1]<$length{$id_c}) && ($tab[2]>$length{$id_c})){
					$tag=1;
					my $stop=$tab[2]-$length{$id_c};
					$tab[2]=$stop;
					my $new_line=join("\t",@tab);
					$predict{$id_c}{$tab[0]}=$new_line;
# 					print "on chevauche l'origine avec un nouveau gene ($tab[0]) ! on rajoute la prédiction\n";
					$tag=1;
				}
			}
		}
	}
	if($tag==1){
		my %to_start;
# 		print "on a remplacé certains orfs pour $id_c, on fait du ménage en cherchant des gènes qui seraient partiellement prédits, soit au début, soit à la fin du contig\n";
		foreach(sort {$order_gene{$a} <=> $order_gene{$b} } keys %{$predict{$id_c}}){
			my @tab=split("\t",$predict{$id_c}{$_});
			if ($tab[5]!=11){
# 				print "\t$tab[0] n'a pas le(s) codon(s) start et / ou stop\n";
				if(($tab[1]<3) || ($tab[2]>($length{$id_c}-3))){
# 					print "\tet il est autour de l'origine, on enlève\n";
					if ($tab[1]<3){
						$to_start{$tab[0]}{"start"}=$tab[1];
						$to_start{$tab[0]}{"stop"}=$tab[2];
					}
					delete($predict{$id_c}{$tab[0]});
				}
				elsif(($tab[2]>997) && ($tab[2]<1001)){ # si on est dans la zone des 1000
					foreach(keys %to_start){
						my $total=($length{$id_c}-$tab[1]+1)+($to_start{$_}{"stop"}); # longueur du fragment terminal + longueur du fragment initial potentiel - pour ce deuxième, on omet le -1 +1 puisque de toute facon ce fragment commencera forcèment en +1
						if ($total % 3 == 0){
#							print "\t et on peut l'étendre grâce à une autre prédiction : -> $to_start{$_}{stop}\n";
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
	close CIR;
}

## Generation of the final files
## One with all sequences nett and filtered (based on number of genes) - Fasta
## One of the associated gene prediction - MGA-like
## One of the predicted protein sequences - Fasta
my $final_file=$tmp_dir."/".$id."_nett_filtered.fasta";
my $out_final=$tmp_dir."/".$id."_mga_final.predict";
my $prot_file=$tmp_dir."/".$id."_prots.fasta";
open(FNA,">$final_file") || die "pblm opening file $final_file\n";
open(SF,">$out_final") || die "pblm ouverture fichier $out_final\n";
open(PROT,">$prot_file") || die "pblm ouverture fichier $prot_file\n";
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
		print SF ">$id\t$length{$id}\n";
		print FNA ">$id\n";
		my $seq_c=$seq_base{$id};
		print FNA "$seq_c\n";
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
			print SF "$predict{$id}{$_}\n";
			my @tab=split("\t",$predict{$id}{$_});
			my $name=$tab[0];
			my $start=$tab[1];
			my $stop=$tab[2];
			my $sens=$tab[3];
			my $frame=$tab[4];
			my $frag="";
			# Cas "normal" (on chevauche pas l'origine)
			if ($start<$stop){
				my $long=$stop-$start+1;
				$frag=substr($seq_c,$start-1,$long);
			}
			# Cas exceptionnel, on chevauche l'origine du contig
			else{
				my $l1=length($seq_c)-$start+1;
				$frag=substr($seq_c,$start-1,$l1);
				$frag.=substr($seq_c,0,$stop);
			}
			## POUR RECUPERER LA SEQ PROT
			my $seq_bio = Bio::Seq->new(-seq =>$frag,-alphabet => 'dna' );
			my @seqs = Bio::SeqUtils->translate_6frames($seq_bio);
			my $cadre=0;
			if ($sens eq "-"){$cadre=3;}
			my $prot=$seqs[$cadre];
			my $prot_sequence=$prot->seq;
			if ($prot_sequence=~/\*$/){
	# 				print "on enlève le codon stop final pour muscle\n";
				chop($prot_sequence);
			}
			my $id_out=$id."-".$name;
			if (($prot_sequence=~/X{50,}/) || ($prot_sequence=~/F{50,}/) || ($prot_sequence=~/A{50,}/) || ($prot_sequence=~/K{50,}/) || ($prot_sequence=~/P{50,}/)){
				print "we exclude $id_out because there is a pblm with the sequence -> too many succesive X, F, A, K or P\n";
			}
			else{
				print PROT ">$id_out\n$prot_sequence\n";
			}
		}
	}
}
close FNA;
close SF;
close PROT;

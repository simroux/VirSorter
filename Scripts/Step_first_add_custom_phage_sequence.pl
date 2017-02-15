#!/usr/bin/env perl

use strict;
use autodie;
use FindBin '$Bin';
use Bio::Seq;
use File::Spec::Functions;
use File::Path 'mkpath';
use File::Which 'which';
# Script to generate a new db with putative new clusters
# Argument 0 : Fasta file of the new phages

if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[2])))
{
	print "# Script to generate a new db with putative new clusters
# Argument 0 : fasta of custom phages
# Argument 1 : db-in directory
# Argument 2 : db-out directory
\n";
	die "\n";
}

my $path_to_makeblastdb = which("makeblastdb")    or die "No makeblastdb\n";
my $path_to_blastp      = which("blastp")         or die "No blastp\n";
my $path_to_muscle      = which("muscle")         or die "No muscle\n";
my $path_to_hmmbuild    = which("hmmbuild")       or die "No hmmbuild\n";
my $path_to_hmmpress    = which("hmmpress")       or die "No hmmpress\n";
my $path_hmmsearch      = which("hmmsearch")      or die "No hmmsearch\n";
my $path_to_mga         = which("mga_linux_ia64") or die "No mga\n";
my $MCX_LOAD            = which("mcxload")        or die "No mcxloafd\n";
my $MCL                 = which("mcl")            or die "No mcl\n";

my $min_seq_in_a_cluster=3;
my $n_cpus=8;

my $fasta_contigs=$ARGV[0];
my $db_in=$ARGV[1];
my $db_out=$ARGV[2];

my $tmp_dir=$db_out."/initial_db";
`mkdir $tmp_dir`;
`cp $db_in/* $tmp_dir/`;
$tmp_dir.="/";
my $log_out=$tmp_dir."log_out_step_custom_phage";
my $log_err=$tmp_dir."log_err_step_custom_phage";

my $db_phage              = $tmp_dir . "Pool_clusters.hmm";
my $blastable_unclustered = $tmp_dir . "Pool_new_unclustered";
my $fasta_unclustered     = $tmp_dir . "Pool_new_unclustered.faa";
my $ref_phage_clusters    = $tmp_dir . "Phage_Clusters_current.tab";
my $blast_unclustered     = $tmp_dir . "Blast_unclustered.tab";

open(F1,"<$fasta_contigs") || die "pblm ouverture fichier $fasta_contigs\n";
my %seq_base;
my $id_seq="";
my $i=0;
my %order_contig;
while(<F1>){
	$_=~s/\r\n/\n/g; #Cas d'un fichier windows ##AJOUT
	chomp($_);
	if ($_=~/^>(\S*)/){$id_seq=$1;$order_contig{$id_seq}=$i;$i++}
	else{$seq_base{$id_seq}.=$_;}
}
close F1;


# Predict genes on the new phages
my $out_file= $db_out."/Custom_phages_mga.predict";
print "$path_to_mga $fasta_contigs -m > $out_file\n";
my $mga=`$path_to_mga $fasta_contigs -m > $out_file`;
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
	elsif ($_=~/^# (\S*)/){
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
my %check_prot;
my $prot_file=$tmp_dir."/Custom_phages_mga_prots.fasta";
open(PROT,">$prot_file") || die "pblm ouverture fichier $prot_file\n";
my $n=0;
foreach(sort {$order_contig{$a} <=> $order_contig{$b} } keys %predict){
	$n++;
	my $id=$_;
	my @tab_genes=sort {$order_gene{$id}{$a} <=> $order_gene{$id}{$b} } keys %{$predict{$id}};
	## We check the first gene and modify it if needed
	my $seq_c=$seq_base{$id};
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

		@tab=split("\t",$predict{$id}{$_});
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
		print PROT ">$id_out\n$prot_sequence\n";
		$check_prot{$id_out}=1;
	}
}
close PROT;
# Clustering the proteins
# - 1 - Hmmsearch vs the original db
my $out_hmmsearch=$tmp_dir."New_prots_vs_Phagedb.tab";
my $out_hmmsearch_bis=$tmp_dir."New_prots_vs_Phagedb.out";
my $cmd_hmm_phage="$path_hmmsearch --tblout $out_hmmsearch --cpu $n_cpus -o $out_hmmsearch_bis --noali $db_phage $prot_file >> $log_out 2>> $log_err";
print "Step 0.9 : $cmd_hmm_phage\n";
`echo $cmd_hmm_phage >> $log_out 2>> $log_err`;
my $out=`$cmd_hmm_phage`;
print "$out\n";
open(HMM,"<$out_hmmsearch") || die ("pblm opening file $out_hmmsearch\n");
my $score_th=200;
my $evalue_th=0.0000000001;
while(<HMM>){
	chomp($_);
	if ($_=~m/^#/){
		next;
	}
	else{
		my @splign=split(m/\s+/,$_);
		my $seq=$splign[0];
		my $match=$splign[2];
		my $evalue=$splign[4];
		my $score=$splign[5];
		if ($score>=$score_th && $evalue<=$evalue_th){
			$check_prot{$seq}=0;
		}
	}
}
close HMM;
# - 2 - All which does not match a known -> get it
my $prot_file_to_cluster=$tmp_dir."/Custom_phages_mga_prots-to-cluster.fasta";
my $tag=0;
my %seq_temp;
open(PROT,"<$prot_file") || die ("pblm opening file $prot_file\n");
open(NEWPROT,">$prot_file_to_cluster") || die ("pblm opening file $prot_file_to_cluster\n");
my $id_c="";
while (<PROT>){
	chomp($_);
	if ($_=~/^>(.*)/){
		my $id=$1;
		$id_c=$id;
		$tag=0;
		if ($check_prot{$id}==1){
			print NEWPROT "$_\n";
			$tag=1;
		}
	}
	elsif($tag==1){
		print NEWPROT "$_\n";
		$seq_temp{$id_c}.=$_;
	}
}
close PROT;
close NEWPROT;
# - 3 - and make new clusters
my $db=$tmp_dir."Custom_phages_mga_prots-to-cluster";
my $cmd_format="$path_to_makeblastdb -in $prot_file_to_cluster -out $db -dbtype prot";
print "$cmd_format\n";
my $out=`$cmd_format`;
print "Formatdb : $out\n";
my $cmd_cat="cat $fasta_unclustered >> $prot_file_to_cluster";
print "$cmd_cat\n";
$out=`$cmd_cat`;
print "Cat : $cmd_cat\n";
#     - blast vs themselves and the unclustered
my $out_blast=$tmp_dir."pool_unclustered-and-custom-phages-vs-custom-phages.tab";
my $cmd_blast="$path_to_blastp -query $prot_file_to_cluster -db $db -out $out_blast -outfmt 6 -num_threads 10 -evalue 0.00001"; # On 10 cores to keep a few alive for the rest of the scripts
print "$cmd_blast\n";
$out=`$cmd_blast`;
print "Blast : $out\n";
$cmd_cat="cat $blast_unclustered >> $out_blast";
print "$cmd_cat\n";
$out=`$cmd_cat`;
print "Cat : $out\n";
print "Generating abc file\n";
#     - mcl 
my $out_abc=$tmp_dir."new_clusters.abc";
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
my $out_mci=$tmp_dir."new_clusters.mci";
my $out_tab=$tmp_dir."new_clusters.tab";
my $cmd_mcxload="$MCX_LOAD -abc $out_abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $out_mci -write-tab $out_tab";
print "$cmd_mcxload\n";
$out=`$cmd_mcxload`;
print "Mxc Load : $out\n";
my $dump_file=$tmp_dir."new_clusters.csv";
my $cmd_mcl="$MCL $out_mci -I 2  -use-tab $out_tab -o $dump_file";
print "$cmd_mcl\n";
$out=`$cmd_mcl`;
print "Mcl : $out\n";
#     - make new cluster
my %unclustered;
my %clusters;
my %check_cluster;
my $last_cluster_id=0;
# All predicted proteins clustered in PCs of 3 and more members are added to the db, the others are added to the pool of unclustered
open(DUMP,"<$dump_file") || die "pblm ouverture fichier $dump_file\n";
while(<DUMP>){
	chomp($_);
	my @tab=split("\t",$_);
	my $n_s_c=$#tab+1;
	if ($n_s_c>=$min_seq_in_a_cluster){
		# on a trouvé un cluster de plus de deux
		my $cluster_id=$last_cluster_id+1;
		$cluster_id="Phage_cluster_".$cluster_id."-c";
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
open(FA,"<$prot_file_to_cluster") || die "pblm ouverture fichier $prot_file_to_cluster\n";
while(<FA>){
	chomp($_);
	if ($_=~/^>(\S*)/){
		$id_c=$1;
		if (!defined($check_cluster{$id_c})){$unclustered{$id_c}=1;$check_cluster{$id_c}=1;}
	}
	else{$seq_temp{$id_c}.=$_;}
}
close FA;
`mkdir $tmp_dir/clusts`;
foreach(keys %clusters){
	my $cluster_id=$_;
	my $out_file=$tmp_dir."clusts/".$cluster_id.".faa";
	open(S1,">$out_file") || die "pblm ouverture fichier $out_file\n";
	foreach(keys %{$clusters{$cluster_id}}){
		print S1 ">$_\n$seq_temp{$_}\n";
	}
	close S1;
}
# - 4 - Plus add the unclustered to the unclustered database
my $final_pool_unclustered=$db_out."/Pool_new_unclustered.faa";
my $final_blastable_unclustered=$final_pool_unclustered;
$final_blastable_unclustered=~s/\.faa//;
my $final_blast_unclustered=$db_out."/Blast_unclustered.tab";
open(S1,">$final_pool_unclustered") || die "pblm ouverture fichier $final_pool_unclustered\n";
foreach(keys %unclustered){
	print S1 ">$_\n$seq_temp{$_}\n";
}
close S1;
print "making a blastable db from the new unclustered\n";
$out=`$path_to_makeblastdb -in $final_pool_unclustered -out $final_blastable_unclustered -dbtype prot`;
# We subset the BLAST result to only unclustered proteins, and add it to the previous unclustered blast result
open(BL,"<$out_blast") || die "pblm ouverture fichier $out_blast\n";
open(S1,">$final_blast_unclustered") || die "pblm ouverture fichier $final_blast_unclustered\n";
while(<BL>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($unclustered{$tab[0]}==1 && $unclustered{$tab[1]}==1){
		print S1 "$_\n";
	}
}
close BL;
close S1;
# Generating the new database
my $tag=0;
foreach(sort keys %clusters){
	$tag=1;
	my $ali_id=$_;
	my $path_to_file=$tmp_dir."clusts/".$ali_id;
	my $path_to_fasta=$tmp_dir."clusts/".$ali_id.".faa";
	my $path_to_ali=$tmp_dir."clusts/".$ali_id.".ali_faa";
	my $path_to_hmm=$tmp_dir."clusts/".$ali_id."_ali.hmm";
	if (-e $path_to_ali){
		`rm $path_to_ali $path_to_hmm`;
	}
	my $muscle_out=$tmp_dir."log_out_muscle";
	my $muscle_err=$tmp_dir."log_err_muscle";
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
# We pool all hmm / fasta from all PCs
$out=`cat $db_phage > $db_out/Pool_clusters.hmm`;
print "cat previous hmm : $out\n";
$out=`cat $tmp_dir/clusts/*.hmm >> $db_out/Pool_clusters.hmm`;
print "cat new hmm : $out\n";
# We create a db for hmmscan
$out=`$path_to_hmmpress $db_out/Pool_clusters.hmm`;
print "hmm press :$out\n";
# And update the phage clusters catalog
my $final_catalog=$db_out."/Phage_Clusters_current.tab";
$out=`cat $ref_phage_clusters > $final_catalog`;
print "Cat old catalog : $out\n";
open(CA,">>$final_catalog") || die ("pblm opening file $final_catalog\n");
foreach(keys %clusters){
	my $liste=join(" ",keys %{$clusters{$_}});
	print CA "$_|2||$liste\n";
}
close CA;

#!/usr/bin/env perl

use strict;
use autodie;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Location::Simple;
use Bio::Location::Split;
use Cwd 'cwd';
use File::Spec::Functions;
use File::Path 'mkpath';

# Script to get fasta file from VirSorter results
# Argument 0 : code of the run
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[0])))
{
	print "# Script to get fasta file from VirSorter results
# Argument 0 : code of the run\n";
	die "\n";
}


my $code    = $ARGV[0] or die 'No code';
my $wdir    = $ARGV[1] || cwd();
my $dir_out = catdir($wdir, "Predicted_viral_sequences");

unless (-d $dir_out) {
    mkpath($dir_out);
}

# We decal each zone by 50 nt before and beyond
my $decal=50;
print "Code $code\n";
my $out_file_1  = catfile( $dir_out, $code . '_cat-1.fasta' );
my $out_file_2  = catfile( $dir_out, $code . '_cat-2.fasta' );
my $out_file_3  = catfile( $dir_out, $code . '_cat-3.fasta' );
my $out_file_p1 = catfile( $dir_out, $code . '_prophages_cat-4.fasta' );
my $out_file_p2 = catfile( $dir_out, $code . '_prophages_cat-5.fasta' );
my $out_file_p3 = catfile( $dir_out, $code . '_prophages_cat-6.fasta' );
my $gb_file_1   = catfile( $dir_out, $code . '_cat-1.gb' );
my $gb_file_2   = catfile( $dir_out, $code . '_cat-2.gb' );
my $gb_file_3   = catfile( $dir_out, $code . '_cat-3.gb' );
my $gb_file_p1  = catfile( $dir_out, $code . '_prophages_cat-4.gb' );
my $gb_file_p2  = catfile( $dir_out, $code . '_prophages_cat-5.gb' );
my $gb_file_p3  = catfile( $dir_out, $code . '_prophages_cat-6.gb' );
print join("\n", "The sequences will be put in:",
    ( map { " - $_" } 
        $out_file_1,
        $out_file_2,
        $out_file_3,
        $out_file_p1,
        $out_file_p2,
        $out_file_p3,
    ),
    ''
);

my $summary       = catfile($wdir, $code . '_global-phage-signal.csv');
my $last_affi     = catfile($wdir, $code . '_phage-signal.csv');
my $affi_contigs  = catfile($wdir, $code . '_affi-contigs.csv');
my $fasta_contigs = catfile($wdir, 'fasta', $code . '_nett_filtered.fasta');
my $fasta_prot    = catfile($wdir, 'fasta', $code . '_prots.fasta');

print "Checking '$last_affi'\n";

if (-e $last_affi){
	my %compte=(1=>0,2=>0,3=>0,4=>0,5=>0,6=>0);
	my %check;
	my $current_c="";
	open(SUM, '<', $summary);

	while (<SUM>){
		chomp($_);
		if ($_=~/## (\d)/){$current_c=$1;}
		elsif ($_=~/##/){}
		else{
			my @tab=split(",",$_);
			$tab[0]=~s/\(/_/g;
			$tab[0]=~s/\)/_/g;
			$tab[2]=~s/\(/_/g;
			$tab[2]=~s/\)/_/g;
			$tab[0]=~s/\[/_/g;
			$tab[0]=~s/\]/_/g;
			$tab[2]=~s/\[/_/g;
			$tab[2]=~s/\]/_/g;
			if($tab[0]=~/(.*)-circular/){
# 				$tab[2]=~s/-circular//g;
				$check{$tab[0]}{$tab[2]}{"circular"}=1;
			}
			if($tab[0] eq ""){
				print "!!!!! void\n";
			}
			else{
				$check{$tab[0]}{$tab[2]}{"prophage"}=$tab[2];
				$check{$tab[0]}{$tab[2]}{"line"}=$_;
				$compte{$current_c}++;
				$check{$tab[0]}{$tab[2]}{"category"}=$current_c;
			}
		}
	}
	close SUM;
	print "$code\t$compte{1}\t$compte{2}\t$compte{3}\t$compte{4}\t$compte{5}\t$compte{6}\n";
	# Get the sequence annotation
	my $id_c="";
	my %infos;
	my $i=0;
	open(ANOT,"<$affi_contigs") || die ("pblm opening file $affi_contigs\n");
	while(<ANOT>){
		chomp($_);
		if ($_=~/>(.*)/){
			my @tab=split(/\|/,$1);
			$tab[0]=~s/\(/_/g;
			$tab[0]=~s/\)/_/g;
			$tab[0]=~s/\[/_/g;
			$tab[0]=~s/\]/_/g;
			$id_c=$tab[0];
		}
		else{
			#     0  | 1   | 2  |  3   |  4   |    5     |  6  |   7  |   8     |   9    | 10   | 11
			# gene_id|start|stop|length|strand|affi_phage|score|evalue|category|affi_pfam|score|evalue|
			my @tab=split(/\|/,$_);
			my $gene=$tab[0];
			$gene=~/.*-(gene_\d*)/;
			$gene=$1;
			$infos{$id_c}{$gene}{"start"}=$tab[1];
			$infos{$id_c}{$gene}{"stop"}=$tab[2];
			$infos{$id_c}{$gene}{"length"}=$tab[3];
			$infos{$id_c}{$gene}{"strand"}=$tab[4];
			$infos{$id_c}{$gene}{"category"}=-1;
			$infos{$id_c}{$gene}{"order"}=$i;
			$i++;
			if ($tab[5] eq "-"){ ## no Phage Cluster affiliation
				if ($tab[9] eq "-"){ ## no PFAM either, ok.. 
					$infos{$id_c}{$gene}{"affi"}="hypothetical protein";
				}
				else{
					$infos{$id_c}{$gene}{"affi"}="PFAM-".$tab[9];
				}
			}
			else{
				if ($tab[9] eq "-"){ ## no PFAM or Phage Cluster better than PFAM (score comparison)
					$infos{$id_c}{$gene}{"affi"}=$tab[5];
				}
				else{ ## So we have a PFAM, which is clearly better than Phage Cluster, so we keep it
					$infos{$id_c}{$gene}{"affi"}=$tab[5]."_"."PFAM-".$tab[9];
				}
			}
		}
	}
	close ANOT;
	open(FA,"<$fasta_prot") || die ("pblm opening file $fasta_prot");
	my $gene_c="";
	my $tag=0;
	while (<FA>){
		chomp($_);
		if ($_=~/^>(\S*)/){
			$tag=0;
			my $gene_temp=$1;
			$gene_temp=~/(.*)-(gene_\d*)/;
			$id_c=$1;
			$gene_c=$2;
			if(defined($infos{$id_c}{$gene_c})){$tag=1;}
		}
		elsif($tag==1){
			$infos{$id_c}{$gene_c}{"seq"}.=$_;
		}
	}
	close FA;
	# Now get all the fasta cut of the contigs
	open(SP1,">$out_file_p1") || die ("pblm opening file $out_file_p1\n");
	open(SP2,">$out_file_p2") || die ("pblm opening file $out_file_p2\n");
	open(SP3,">$out_file_p3") || die ("pblm opening file $out_file_p3\n");
	open(S1,">$out_file_1") || die ("pblm opening file $out_file_1\n");
	open(S2,">$out_file_2") || die ("pblm opening file $out_file_2\n");
	open(S3,">$out_file_3") || die ("pblm opening file $out_file_3\n");
	my $output_1 = Bio::SeqIO->new(-file => ">$gb_file_1",-format => 'GenBank');
	my $output_2 = Bio::SeqIO->new(-file => ">$gb_file_2",-format => 'GenBank');
	my $output_3 = Bio::SeqIO->new(-file => ">$gb_file_3",-format => 'GenBank');
	my $output_p1 = Bio::SeqIO->new(-file => ">$gb_file_p1",-format => 'GenBank');
	my $output_p2 = Bio::SeqIO->new(-file => ">$gb_file_p2",-format => 'GenBank');
	my $output_p3 = Bio::SeqIO->new(-file => ">$gb_file_p3",-format => 'GenBank');
	my $sequence=0;
	open(FASTA,"<$fasta_contigs") || die ("pblm opening file $fasta_contigs\n");
	$id_c="";
	my $seq_c="";
	while (<FASTA>){
		chomp($_);
		if ($_=~/^>(.*)/){
			my $id_c_temp=$1;
			$id_c_temp=~s/\(/_/g;
			$id_c_temp=~s/\)/_/g;
			$id_c_temp=~s/\[/_/g;
			$id_c_temp=~s/\]/_/g;
			if (defined($check{$id_c})){
				my $id_red=$id_c;
				print "We had checked $id_c -> $id_red\n";
				foreach(keys %{$check{$id_c}}){
					$id_red=$id_c;
					my @tab=split(",",$check{$id_c}{$_}{"line"});
					$tab[0]=~s/\(/_/g;
					$tab[0]=~s/\)/_/g;
					$tab[2]=~s/\(/_/g;
					$tab[2]=~s/\)/_/g;
					$tab[0]=~s/\[/_/g;
					$tab[0]=~s/\]/_/g;
					$tab[2]=~s/\[/_/g;
					$tab[2]=~s/\]/_/g;
					# $tab[2]=~s/-circular//g;
					my $desc="Putative phage sequence (category $check{$id_c}{$_}{category}), predicted by PhageSorter";
					my $iscirc=0;
					if ($check{$id_c}{$tab[2]}{"circular"}==1){
						# $id_red.="-circ";
						$iscirc=1;
					}
					if ($check{$id_c}{$_}{"category"}<=3){
						print ".. predicted to be a complete phage..\n";
						$sequence = Bio::Seq::RichSeq->new(-display_id => "$id_red", -accession_number => "$id_red", -desc => $desc ,-seq =>"$seq_c",-is_circular =>$iscirc,-division => "ENV",-alphabet => "dna");
						$sequence->add_date(`date +%D`);
						my $featsource = Bio::SeqFeature::Generic->new(-start => 1,-end => length($seq_c),-primary => "source",-tag => {'organism' => "$desc"});
						$sequence->add_SeqFeature($featsource);
						foreach(sort { $infos{$id_c}{$a}{"order"} <=> $infos{$id_c}{$b}{"order"} } keys %{$infos{$id_c}}){
							my $gene=$_;
							my $splitlocation = Bio::Location::Split->new();
							my $strand=0;
							if ($infos{$id_c}{$gene}{"strand"} eq "-"){$strand=-1;}
							# si on est sur un join, etc..
							if ($infos{$id_c}{$gene}{"stop"} < $infos{$id_c}{$gene}{"start"}){
								$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>$infos{$id_c}{$gene}{"start"},-end=>length($seq_c),-strand=>$strand));
								$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>1,-end=>$infos{$id_c}{$gene}{"stop"},-strand=>$strand));
							}
							else{
								$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>$infos{$id_c}{$gene}{"start"},-end=>$infos{$id_c}{$gene}{"stop"},-strand=>$strand));
							}
							my $featgene = Bio::SeqFeature::Generic->new(-location => $splitlocation,-primary => "gene",-tag => {'gene' => "$gene",'locus_tag' => $id_c."_".$gene});
							$sequence->add_SeqFeature($featgene);
							my $product=$infos{$id_c}{$gene}{"affi"};
							my $note="Predicted by MGA";
							my $featcds = Bio::SeqFeature::Generic->new(-location=>$splitlocation,-primary => "CDS",-tag => {'product' => "$product",'note' => "$note",'locus_tag' => $id_c."_".$gene,'codon_start' => "1",'gene' => "$gene",'transl_table' => "11"});
							$sequence->add_SeqFeature($featcds);
							$featcds->add_tag_value('translation',$infos{$id_c}{$gene}{"seq"});
						}
						if ($check{$id_c}{$_}{"category"}==1){
							print S1 ">$id_red-cat_$check{$id_c}{$_}{category}\n$seq_c\n";
							$output_1->write_seq($sequence);
						}
						elsif ($check{$id_c}{$_}{"category"}==2){
							print S2 ">$id_red-cat_$check{$id_c}{$_}{category}\n$seq_c\n";
							$output_2->write_seq($sequence);
						}
						elsif ($check{$id_c}{$_}{"category"}==3){
							print S3 ">$id_red-cat_$check{$id_c}{$_}{category}\n$seq_c\n";
							$output_3->write_seq($sequence);
						}
					}
					else{
						print ".. predicted to be a prophage..\n";
						if ($tab[2]=~/^$id_c-(gene_\d+)-(gene_\d+)/){
							my $gene_start=$1;
							my $start=$infos{$id_c}{$gene_start}{"start"}-$decal;
							if ($start<0){$start=0;}
							my $gene_stop=$2;
							my $stop=$infos{$id_c}{$gene_stop}{"stop"}+$decal;
							my $length=$stop-$start;
							print "  from $1 to $2 .. from $start to $stop ($length)\n";
							my $substr=substr($seq_c,$start,$length);
							$iscirc=0; # An integrated prophage cannot be circular, so set this to linear 
							my $display_id=$id_red."_".$gene_start."_".$gene_stop."-".$start."-".$stop."-cat_".$check{$id_c}{$_}{"category"};
							$sequence = Bio::Seq::RichSeq->new(-display_id => "$display_id", -accession_number => "$display_id", -desc => $desc ,-seq =>"$substr",-is_circular =>$iscirc,-division => "ENV",-alphabet => "dna");
							$sequence->add_date(`date +%D`);
							my $featsource = Bio::SeqFeature::Generic->new(-start => 1,-end => length($substr),-primary => "source",-tag => {'organism' => "$desc"});
							$sequence->add_SeqFeature($featsource);
							foreach(sort { $infos{$id_c}{$a}{"order"} <=> $infos{$id_c}{$b}{"order"} } keys %{$infos{$id_c}}){
								my $gene=$_;
								# Check if the gene is in the fragment entirely
								if (($infos{$id_c}{$gene}{"start"}>=$start) && ($infos{$id_c}{$gene}{"start"}<=$stop) && ($infos{$id_c}{$gene}{"stop"}>=$start) && ($infos{$id_c}{$gene}{"stop"}<=$stop)){
									my $splitlocation = Bio::Location::Split->new();
									my $strand=0;
									if ($infos{$id_c}{$gene}{"strand"} eq "-"){$strand=-1;}
									$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>$infos{$id_c}{$gene}{"start"}-$start,-end=>$infos{$id_c}{$gene}{"stop"}-$start,-strand=>$strand));
									my $featgene = Bio::SeqFeature::Generic->new(-location => $splitlocation,-primary => "gene",-tag => {'gene' => "$gene",'locus_tag' => $id_c."_".$gene});
									$sequence->add_SeqFeature($featgene);
									my $product=$infos{$id_c}{$gene}{"affi"};
									my $note="Predicted by MGA";
									my $featcds = Bio::SeqFeature::Generic->new(-location=>$splitlocation,-primary => "CDS",-tag => {'product' => "$product",'note' => "$note",'locus_tag' => $id_c."_".$gene,'codon_start' => "1",'gene' => "$gene",'transl_table' => "11"});
									$sequence->add_SeqFeature($featcds);
									$featcds->add_tag_value('translation',$infos{$id_c}{$gene}{"seq"});
								}
							}
							if ($check{$id_c}{$_}{"category"}==4){
								print SP1 ">".$display_id."\n".$substr."\n";
								$output_p1->write_seq($sequence);
							}
							elsif ($check{$id_c}{$_}{"category"}==5){
								print SP2 ">".$display_id."\n".$substr."\n";
								$output_p2->write_seq($sequence);
							}
							elsif($check{$id_c}{$_}{"category"}==6){
								print SP3 ">".$display_id."\n".$substr."\n";
								$output_p3->write_seq($sequence);
							}
						}
						else{
							print "Pblm with $tab[2] - tab 2\n";
						}
					}
				}
			}
			$id_c=$id_c_temp;
			$seq_c="";
		}
		else{$seq_c.=$_;}
	}
	close FASTA;
	# We do not forget the last one
	if (defined($check{$id_c})){
		my $id_red=$id_c;
		print "We had checked $id_c -> $id_red\n";
		foreach(keys %{$check{$id_c}}){
			$id_red=$id_c;
			my @tab=split(",",$check{$id_c}{$_}{"line"});
			$tab[0]=~s/\(/_/g;
			$tab[0]=~s/\)/_/g;
			$tab[2]=~s/\(/_/g;
			$tab[2]=~s/\)/_/g;
			$tab[0]=~s/\[/_/g;
			$tab[0]=~s/\]/_/g;
			$tab[2]=~s/\[/_/g;
			$tab[2]=~s/\]/_/g;
#			$tab[2]=~s/-circular//g;
			my $desc="Putative phage sequence (category $check{$id_c}{$_}{category}), predicted by PhageSorter";
			my $iscirc=0;
			if ($check{$id_c}{$tab[2]}{"circular"}==1){
#				$id_red.="-circ";
				$iscirc=1;
			}
			if ($check{$id_c}{$_}{"category"}<=3){
				print ".. predicted to be a complete phage..\n";
				$sequence = Bio::Seq::RichSeq->new(-display_id => "$id_red", -accession_number => "$id_red", -desc => $desc ,-seq =>"$seq_c",-is_circular =>$iscirc,-division => "ENV",-alphabet => "dna");
				$sequence->add_date(`date +%D`);
				my $featsource = Bio::SeqFeature::Generic->new(-start => 1,-end => length($seq_c),-primary => "source",-tag => {'organism' => "$desc"});
				$sequence->add_SeqFeature($featsource);
				foreach(sort { $infos{$id_c}{$a}{"order"} <=> $infos{$id_c}{$b}{"order"} } keys %{$infos{$id_c}}){
					my $gene=$_;
					my $splitlocation = Bio::Location::Split->new();
					my $strand=0;
					if ($infos{$id_c}{$gene}{"strand"} eq "-"){$strand=-1;}
					# si on est sur un join, etc..
					if ($infos{$id_c}{$gene}{"stop"} < $infos{$id_c}{$gene}{"start"}){
						$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>$infos{$id_c}{$gene}{"start"},-end=>length($seq_c),-strand=>$strand));
						$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>1,-end=>$infos{$id_c}{$gene}{"stop"},-strand=>$strand));
					}
					else{
						$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>$infos{$id_c}{$gene}{"start"},-end=>$infos{$id_c}{$gene}{"stop"},-strand=>$strand));
					}
					my $featgene = Bio::SeqFeature::Generic->new(-location => $splitlocation,-primary => "gene",-tag => {'gene' => "$gene",'locus_tag' => $id_c."_".$gene});
					$sequence->add_SeqFeature($featgene);
					my $product=$infos{$id_c}{$gene}{"affi"};
					my $note="Predicted by MGA";
					my $featcds = Bio::SeqFeature::Generic->new(-location=>$splitlocation,-primary => "CDS",-tag => {'product' => "$product",'note' => "$note",'locus_tag' => $id_c."_".$gene,'codon_start' => "1",'gene' => "$gene",'transl_table' => "11"});
					$sequence->add_SeqFeature($featcds);
					$featcds->add_tag_value('translation',$infos{$id_c}{$gene}{"seq"});
				}
				if ($check{$id_c}{$_}{"category"}==1){
					print S1 ">$id_red-cat_$check{$id_c}{$_}{category}\n$seq_c\n";
					$output_1->write_seq($sequence);
				}
				if ($check{$id_c}{$_}{"category"}==2){
					print S2 ">$id_red-cat_$check{$id_c}{$_}{category}\n$seq_c\n";
					$output_2->write_seq($sequence);
				}
				elsif ($check{$id_c}{$_}{"category"}==3){
					print S3 ">$id_red-cat_$check{$id_c}{$_}{category}\n$seq_c\n";
					$output_3->write_seq($sequence);
				}
			}
			else{
				print ".. predicted to be a prophage..\n";
				if ($tab[2]=~/^$id_c-(gene_\d+)-(gene_\d+)/){
					my $gene_start=$1;
					my $start=$infos{$id_c}{$gene_start}{"start"}-$decal;
					if ($start<0){$start=0;}
					my $gene_stop=$2;
					my $stop=$infos{$id_c}{$gene_stop}{"stop"}+$decal;
					my $length=$stop-$start;
					print "  from $1 to $2 .. from $start to $stop ($length)\n";
                                        $iscirc=0; # An integrated prophage cannot be circular, so set this to linear 
                                        my $display_id=$id_red."_".$gene_start."_".$gene_stop."-".$start."-".$stop."-cat_".$check{$id_c}{$_}{"category"};
					my $substr=substr($seq_c,$start,$length);
					$sequence = Bio::Seq::RichSeq->new(-display_id => "$display_id", -accession_number => "$display_id", -desc => $desc ,-seq =>"$substr",-is_circular =>$iscirc,-division => "ENV",-alphabet => "dna");
					$sequence->add_date(`date +%D`);
					my $featsource = Bio::SeqFeature::Generic->new(-start => 1,-end => length($substr),-primary => "source",-tag => {'organism' => "$desc"});
					$sequence->add_SeqFeature($featsource);
					foreach(sort { $infos{$id_c}{$a}{"order"} <=> $infos{$id_c}{$b}{"order"} } keys %{$infos{$id_c}}){
						my $gene=$_;
						if ((($infos{$id_c}{$gene}{"start"}-$start)>0) && (($infos{$id_c}{$gene}{"start"}-$start)<=$stop) && (($infos{$id_c}{$gene}{"stop"}-$start)>0) && (($infos{$id_c}{$gene}{"stop"}-$start)<=$stop)){
							my $splitlocation = Bio::Location::Split->new();
							my $strand=0;
							if ($infos{$id_c}{$gene}{"strand"} eq "-"){$strand=-1;}
							$splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>$infos{$id_c}{$gene}{"start"}-$start,-end=>$infos{$id_c}{$gene}{"stop"}-$start,-strand=>$strand));
							my $featgene = Bio::SeqFeature::Generic->new(-location => $splitlocation,-primary => "gene",-tag => {'gene' => "$gene",'locus_tag' => $id_c."_".$gene});
							$sequence->add_SeqFeature($featgene);
							my $product=$infos{$id_c}{$gene}{"affi"};
							my $note="Predicted by MGA";
							my $featcds = Bio::SeqFeature::Generic->new(-location=>$splitlocation,-primary => "CDS",-tag => {'product' => "$product",'note' => "$note",'locus_tag' => $id_c."_".$gene,'codon_start' => "1",'gene' => "$gene",'transl_table' => "11"});
							$sequence->add_SeqFeature($featcds);
							$featcds->add_tag_value('translation',$infos{$id_c}{$gene}{"seq"});
						}
					}
					if ($check{$id_c}{$_}{"category"}==4){
						print SP1 ">".$display_id."\n".$substr."\n";
						$output_p1->write_seq($sequence);
					}
					elsif ($check{$id_c}{$_}{"category"}==5){
						print SP2 ">".$display_id."\n".$substr."\n";
						$output_p2->write_seq($sequence);
					}
					elsif($check{$id_c}{$_}{"category"}==6){
						print SP3 ">".$display_id."\n".$substr."\n";
						$output_p3->write_seq($sequence);
					}
				}
				else{
					print "Pblm with $tab[2] - tab 2\n";
				}
			}
		}
	}
	close S1;
	close S2;
	close S3;
	close SP1;
	close SP2;
	close SP3;
	$output_1->close();
	$output_2->close();
	$output_3->close();
	$output_p1->close();
	$output_p2->close();
	$output_p3->close();
}
else{
	print "$code\tin progress\n";
}

#!/usr/bin/env perl

=head1 USAGE 

  ./wrapper_phage_contigs_sorter_iPlant.pl -d Code_dataset --fna Fasta file of contigs --db 1 --wdir /path/to/working_directory
  Database codes : 1 for RefseqABVir only, 2 for RefseqABVir + Viromes
  An additional set of reference sequences can be added to the database as a fasta file with the argument cp (--cp /path/to/fasta_file)
  
=cut

use strict;
use warnings;
use autodie;
use File::Spec::Functions;
use File::Path 'mkpath';
use Getopt::Long 'GetOptions';
use Pod::Usage;
use Cwd 'cwd';

my $help              = '';
my $code_dataset      = $0;
my $original_fna_file = '';
my $choice_database   = '';
my $custom_phage      = '';
my $wdir              = cwd();

GetOptions(
   'd|dataset=s' => \$code_dataset,
   'fna=s'       => \$original_fna_file,
   'db=i'        => \$choice_database,
   'wdir=s'      => \$wdir,
   'cp=s'        => \$custom_phage,
   'h|help'      => \$help,
);

if ($help || $original_fna_file eq '' || $choice_database eq '') {
   pod2usage();
}

print "Dataset : $code_dataset, Fna file : $original_fna_file, Db : $choice_database, Wdir : $wdir, Custom phages : $custom_phage\n";
# We check if the custom phage is an actual fasta file, or if it's the working dir (which means -> no custom phage).
if ($custom_phage=~/.*\.f.*/){}
else{
	if($wdir=~/.*\/$custom_phage$/){
		print "The custom phage is actually the wdir id, so we remove it\n";
		$custom_phage="";
	}
	else{
		die("we do not understand this custom phage : $custom_phage");
	}
}

# Need 2 databases
# PCs from Refseq (phages) or PCs from Refseq+Viromes
# PFAM (26.0)

my $n_cpus = 8;
my $virsorter_dir="/usr/local/bin/Virsorter/";

# my $code_dataset      = $ARGV[0];
# my $original_fna_file = $ARGV[1];
# my $choice_database   = $ARGV[2];
# my $wdir              = $ARGV[3];
# my $custom_phage      = "";
# if ( defined( $ARGV[4] ) ) { $custom_phage = $ARGV[4]; }
print "#%#%#%#%#%# Processing $code_dataset....\n";
my $microbial_base_needed = 0;
## replace this directory with the iPlant dir
#my $wdir=$wdir."/".$code_dataset."/";
my $log_out = "log_out";
my $log_err = "log_err";

my $path_to_mga = catfile($virsorter_dir,"Tools/Metagene_annotator/mga_linux_ia64");
my $path_hmmsearch = catfile($virsorter_dir,"Tools/hmmer-3.0-linux-intel-x86_64/binaries/hmmsearch");
my $path_blastall    = "/usr/bin/blastall";
my $path_to_formatdb = "/usr/bin/formatdb";
my $script_dir       = catdir($virsorter_dir,"Scripts/");
my $dir_Phage_genes  = catdir($virsorter_dir,"Database/Phage_gene_catalog/");
my $ref_phage_clusters = catfile($virsorter_dir,"Database/Phage_gene_catalog/Phage_Clusters_current.tab");
my $readme_file = catfile($virsorter_dir,"VirSorter_Readme.txt");

if ( $choice_database == 2 ) {
    $dir_Phage_genes = catdir($virsorter_dir,"Database/Phage_gene_catalog_plus_viromes/");
    $ref_phage_clusters = catfile($virsorter_dir,"Database/Phage_gene_catalog_plus_viromes/Phage_Clusters_current.tab");
}
my $db_PFAM_a = catfile($virsorter_dir,"Database/PFAM_27/Pfam-A.hmm");
my $db_PFAM_b = catfile($virsorter_dir,"Database/PFAM_27/Pfam-B.hmm");

my $out = "";

## SETTING UP THE WORKING DIRECTORY
if (!-d 'log') {
    mkpath('log');
}
$out = `rm -r log/* *.csv`;
print "rm -r log* *.csv => $out\n";

# cp fasta file in the wdir
my $fastadir = "fasta";
if ( -d $fastadir ) {
    # Should not be used on iPlant, but you never know
    $out = "Already a fasta directory .. skipping (the great guru)\n";
}
else {
    mkpath($fastadir);
    my $fna_file = catfile( $fastadir, "input_sequences.fna" );
    open my $fa, '<', $original_fna_file;
    open my $s1, '>', $fna_file;
    while (<$fa>) {
        chomp($_);
        if ( $_ =~ /^>(.*)/ ) {
            my $id = $1;
            $id =~ s/[\/\.,\|\s?!\*%]/_/g;
            my $new_id = $code_dataset . "_" . $id;
            print $s1 ">$new_id\n";
        }
        else {
            print $s1 "$_\n";
        }

    }
    close $fa;
    close $s1;
    # detect circular, predict genes on contigs and extract proteins, as well
    # as filtering on size (nb genes) and/or circular
    my $nb_gene_th = 2; # At least two complete genes on the contig
    my $path_script_step_1 = catfile($script_dir,"Step_1_contigs_cleaning_and_gene_prediction.pl");
    my $cmd_step_1 = "$path_script_step_1 $code_dataset $fastadir $fna_file $nb_gene_th >> $log_out 2>> $log_err";
    print "Step 0.5 : $cmd_step_1\n";
    `echo $cmd_step_1 >> $log_out 2>> $log_err`;
    { $out = `$cmd_step_1`; }
}

print "\t$out\n";
my $fasta_contigs_nett = catfile( $fastadir, $code_dataset . "_nett_filtered.fasta" );
my $fasta_file_prots = catfile( $fastadir, $code_dataset . "_prots.fasta" );

# Match against PFAM, once for all
# compare to PFAM a then b (hmmsearch)
my $out_hmmsearch_pfama     = "Contigs_prots_vs_PFAMa.tab";
my $out_hmmsearch_pfama_bis = "Contigs_prots_vs_PFAMa.out";
my $cmd_hmm_pfama =
"$path_hmmsearch --tblout $out_hmmsearch_pfama --cpu $n_cpus -o $out_hmmsearch_pfama_bis --noali $db_PFAM_a $fasta_file_prots >> $log_out 2>> $log_err";
print "Step 0.8 : $cmd_hmm_pfama\n";

`echo $cmd_hmm_pfama >> $log_out 2>> $log_err`;

if ( !( -e $out_hmmsearch_pfama ) ) {
    $out = `$cmd_hmm_pfama`;
    print "\t$out\n";
}
else { 
    $out = "Already a results for PFAM A .. skipping (the great guru)\n"; 
}
my $out_hmmsearch_pfamb     = "Contigs_prots_vs_PFAMb.tab";
my $out_hmmsearch_pfamb_bis = "Contigs_prots_vs_PFAMb.out";
my $cmd_hmm_pfamb =
"$path_hmmsearch --tblout $out_hmmsearch_pfamb --cpu $n_cpus -o $out_hmmsearch_pfamb_bis --noali $db_PFAM_b $fasta_file_prots >> $log_out 2>> $log_err";
print "Step 0.9 : $cmd_hmm_pfamb\n";
`echo $cmd_hmm_pfamb >> $log_out 2>> $log_err`;

if ( !( -e $out_hmmsearch_pfamb ) ) {
    $out = `$cmd_hmm_pfamb`;
    print "\t$out\n";
}
else { 
    $out = "Already a results for PFAM B .. skipping (the great guru)\n"; 
}

# Now work on the phage gene catalog

# Files that will stay along the computations
my $predict_file = catfile( $fastadir, $code_dataset . "_mga_final.predict" );
my $out_hmmsearch = "Contigs_prots_vs_Phage_Gene_Catalog.tab";
my $out_hmmsearch_bis        = "Contigs_prots_vs_Phage_Gene_Catalog.out";
my $out_blast_unclustered    = "Contigs_prots_vs_Phage_Gene_unclustered.tab";
my $out_file_affi            = $code_dataset . "_affi-contigs.csv";
my $out_file_phage_fragments = $code_dataset . "_phage-signal.csv";
my $global_out_file          = $code_dataset . "_global-phage-signal.csv";
my $new_prots_to_cluster     = $code_dataset . "_new_prot_list.csv";

# Constant scripts
my $script_merge_annot = catfile($script_dir,"Step_2_merge_contigs_annotation.pl");
my $cmd_merge =
"$script_merge_annot $predict_file $out_hmmsearch $out_blast_unclustered $out_hmmsearch_pfama $out_hmmsearch_pfamb $ref_phage_clusters $out_file_affi >> $log_out 2>> $log_err";
# my $cmd_merge = "$script_merge_annot -m $predict_file -hmm_pc $out_hmmsearch -blast_pc $out_blast_unclustered -hmm_pfa $out_hmmsearch_pfama -hmm_pfb $out_hmmsearch_pfamb -ref_pc $ref_phage_clusters -out_f $out_file_affi >> $log_out 2>> $log_err";

my $script_detect = catfile($script_dir,"Step_3_highlight_phage_signal.pl");
my $cmd_detect =
"$script_detect $out_file_affi $out_file_phage_fragments >> $log_out 2>> $log_err";

my $script_summary =catfile($script_dir,"Step_4_summarize_phage_signal.pl");
my $cmd_summary =
"$script_summary $out_file_affi $out_file_phage_fragments $global_out_file $new_prots_to_cluster >> $log_out 2>> $log_err";

# # Get the final result file ready
`touch $global_out_file`;
my $r_n = -1;
# Si on a des nouvelles prots a clusteriser ou si on est dans la premiere
# revision
while ( (-e $new_prots_to_cluster || $r_n == -1) && ($r_n<=10) ) {
    $r_n++;    # New revision of the prediction
    my $dir_revision = "r_" . $r_n;
    print "### Revision $r_n\n";
    if ( !-d $dir_revision ) {
        ## mkdir de la db de cette revision
        #print "mkdir $dir_revision >> $log_out 2>> $log_err\n";
        #$out=`mkdir $dir_revision >> $log_out 2>> $log_err`;
        mkpath($dir_revision);
        print "Out : $out\n";
        ## Clustering of the new prots with the unclustered
        my $script_new_cluster = catfile($script_dir,"Step_0_make_new_clusters.pl");
        # First revision, we just import the Refseq database
        if ( $r_n == 0 ) {
            #`mkdir $dir_revision/db`;
            mkpath( catdir( $dir_revision, 'db' ) );

            ## Adding custom sequences to the database if required by the user
            if ( $custom_phage ne "" ) {
                my $script_custom_phage = catfile($script_dir,"Step_first_add_custom_phage_sequence.pl");
                $out =`$script_custom_phage $custom_phage $dir_Phage_genes/ $dir_revision/db >> $log_out 2>> $log_err`;
                print "Adding custom phage to the database : $out\n";
            }
            # should replace Pool_cluster / Pool_unclustered and
            # Pool_new_unclustered else , we just import the Refseq database
            else { `cp $dir_Phage_genes/* $dir_revision/db/`; }
        }
        else {
            my $previous_r = $r_n - 1;
            my $previous_fasta_unclustered =
              catfile( "r_" . $previous_r, "db", "Pool_unclustered.faa" );
            my $cmd_new_clusters = join(' ',
                "$script_new_cluster $dir_revision $fasta_file_prots",
                "$previous_fasta_unclustered",
                "$new_prots_to_cluster >> $log_out 2>> $log_err"
            );

            print "$cmd_new_clusters\n";
            $out = `$cmd_new_clusters`;
            print "Step 1.1 new clusters and new database : $out\n";
            # Rm the list of prots to be clustered now that they should be
            # clustered
            $out = `rm $new_prots_to_cluster`;
            print "rm $new_prots_to_cluster -> $out\n";
        }

        # Check if there are some data in these new clusters, or if all the new
        # proteins are unclustered
        my $new_db_profil = catfile( $dir_revision, "db", "Pool_clusters.hmm" );
        my $check = 0;
        open my $DB, '<', $new_db_profil;
        while (<$DB>) {
	    chomp($_);
            if ( $_ =~ /^NAME/ ) { $check++; print "there is a cluster $_ in the database, so we're good\n"; }
        }
        close $DB;
        if ( $check == 0 ) {
            print "There is no clusters in the database, so we skip the hmmsearch\n";
        }
        else {
            my $out_hmmsearch_new =
              catfile( $dir_revision, "Contigs_prots_vs_New_clusters.tab" );
            my $out_hmmsearch_bis_new =
              catfile( $dir_revision, "Contigs_prots_vs_New_clusters.out" );
            my $cmd_hmm_cluster = join(' ',
                "$path_hmmsearch --tblout $out_hmmsearch_new --cpu $n_cpus",
                "-o $out_hmmsearch_bis_new --noali $new_db_profil",
                "$fasta_file_prots >> $log_out 2>> $log_err"
            );

            print "Step 1.2 : $cmd_hmm_cluster\n";

            `echo $cmd_hmm_cluster >> $log_out 2>> $log_err`;

            $out = `$cmd_hmm_cluster`;
            print "\t$out\n";

            $out = `cat $out_hmmsearch_new >> $out_hmmsearch`;
            print "\t$out\n";
        }

        my $out_blast_new_unclustered =
          catfile( $dir_revision, "Contigs_prots_vs_New_unclustered.tab" );
        my $blastable_unclustered =
          catfile( $dir_revision, 'db', 'Pool_new_unclustered' );
        my $cmd_blast_unclustered = join(' ',
            "$path_blastall -p blastp -i $fasta_file_prots -d",
            "$blastable_unclustered -o $out_blast_new_unclustered -a $n_cpus", 
            "-m 8 -e 0.001 >> $log_out 2>> $log_err"
        );

        print "\nStep 1.3 : $cmd_blast_unclustered\n";
        `echo $cmd_blast_unclustered >> $log_out 2>> $log_err`;
        $out = `$cmd_blast_unclustered`;
        print "\t$out\n";
        $out = `cat $out_blast_new_unclustered >> $out_blast_unclustered`;
        print "\t$out\n";
        ## Make backup of the previous files to have trace of the different steps
        my $backup_affi = catfile( $dir_revision, "affi_backup.csv" );
        my $backup_phage_signal =
          catfile( $dir_revision, "phage_signal_backup.csv" );
        my $backup_global_signal =
          catfile( $dir_revision, "global_signal_backup.csv" );
        if ( -e $out_file_affi ) { `cp $out_file_affi $backup_affi`; }
        if ( -e $out_file_phage_fragments ) {
            `cp $out_file_phage_fragments $backup_phage_signal`;
        }
        if ( -e $global_out_file ) {
            `cp $global_out_file $backup_global_signal`;
        }
    }

    ## Complete the affi
    print "Step 2 : $cmd_merge\n";
    `echo $cmd_merge >> $log_out 2>> $log_err`;
    $out = `$cmd_merge`; 
    ## This generate a csv table including the map of each contig, with PFAM
    #and Viral PCs annotations, as well as strand and length of genes

    print "\t$out\n";
    ## Complete the summary
    print "Step 3 : $cmd_detect\n";
    `echo $cmd_detect >> $log_out 2>> $log_err`;
    $out = `$cmd_detect`;
    print "\t$out\n";

    # Decide which contigs are entirely viral and which are prophages, and
    # which of both of these categories are phage enough to be added to the
    # databases
    print "Setting up the final result file\n";
    print "Step 4 : $cmd_summary\n";
    `echo $cmd_summary >> $log_out 2>> $log_err`;
    $out = `$cmd_summary`;
    print "\t$out\n";
}

# Last step -> extract all sequences as fasta files and gb
my $script_generate_output = catfile($script_dir,"Step_5_get_phage_fasta-gb.pl");
my $cmd_step_5 = "$script_generate_output $code_dataset >> $log_out 2>> $log_err";
print "\nStep 5 : $cmd_step_5\n";

`echo $cmd_step_5 >> $log_out 2>> $log_err`;

$out = `$cmd_step_5`;
print "\t$out\n";


# Plus clean the output directory
print "Cleaning the output directory\n";
# We rm the first db to not overload user disk space
my $db_revision_0="r_0/db";
$out=`rm -r $db_revision_0`;
print "rm -r $db_revision_0 : $out\n";
`mv fasta/ Fasta_files/`;
# We put all results from Hmmsearch and BLAST files in a separate directory
my $store_database_comparison="Tab_files";
mkpath($store_database_comparison);
`mv $out_hmmsearch $store_database_comparison/`;
# `mv $out_hmmsearch_bis $store_database_comparison/`;
`mv $out_blast_unclustered $store_database_comparison/`;
`mv $out_hmmsearch_pfama $store_database_comparison/`;
`mv $out_hmmsearch_pfama_bis $store_database_comparison/`;
`mv $out_hmmsearch_pfamb $store_database_comparison/`;
`mv $out_hmmsearch_pfamb_bis $store_database_comparison/`;
`mv error.log logs/`;
`mv formatdb.log logs/`;
`mv log_err logs/Virsorter_error_log`;
`mv log_out logs/Virsorter_stdout_log`;
# We put all the files linked to the metric computation in a new directory
my $store_metric_files="Metric_files";
mkpath($store_metric_files);
`mv $out_file_affi $store_metric_files/VIRSorter_affi-contigs.tab`;
my $out_file_affi_ref  = $code_dataset . "_affi-contigs.refs";
`mv $out_file_affi_ref $store_metric_files/`;
`mv $out_file_phage_fragments $store_metric_files/VIRSorter_phage_signal.tab`;
`mv $new_prots_to_cluster $store_metric_files/`;
# And we customize and add the readme file in the output directory
my $datestring=localtime();
my $local_readme_file="Readme.txt";
open my $s1,'>',$local_readme_file;
print $s1 "VirSorter parameters used :\n\n";
print $s1 "--> Fasta file mined for viral sequences : $original_fna_file\n";
print $s1 "--> Viral database used : ";
if ($choice_database==2){print $s1 "Viromes : all bacterial and archaeal virus genomes in Refseq, as of January 2014, plus non-redundant predicted genes from viral metagenomes (including seawater, freshwater, and human-related samples)\n"}
else{print $s1 "RefseqABVir (all bacterial and archaeal virus genomes in Refseq, as of January 2014\n";}
if ($custom_phage eq ""){print $s1 "--> No custom reference sequence was added to the database\n";}
else{print $s1 "--> Custom reference sequences from fasta file $custom_phage were added to the database\n";}
print $s1 "\nThis VirSorter computation finished on $datestring\n";
close $s1;
`cat $readme_file >> $local_readme_file`;

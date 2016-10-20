#!/usr/bin/env perl

=head1 SYNOPSIS

  wrapper_phage_contigs_sorter_iPlant.pl --fasta sequences.fa

Required Arguments:

  -f|--fna       Fasta file of contigs

Options: 

  -d|--dataset   Code dataset (DEFAULT "VIRSorter")
  --cp           Custom phage sequence 
  --db           Either "1" (DEFAULT Refseqdb) or "2" (Viromedb)
  --wdir         Working directory (DEFAULT cwd)
  --ncpu         Number of CPUs
  --virome       Virome decontamination mode, for datasets mostly viral, force the use of generic metrics instead of calculated from the whole dataset

  --help         Show help and exit

=head1 DESCRIPTION

Wrapper for detection of viral contigs

=cut

use strict;
use warnings;
use feature 'say';
use autodie;
use FindBin '$Bin';
use File::Spec::Functions;
use File::Path 'mkpath';
use File::Which 'which';
use Getopt::Long 'GetOptions';
use Pod::Usage;
use Cwd 'cwd';

my $help            = '';
my $code_dataset    = 'VIRSorter';
my $input_file      = '';
my $choice_database = 1;
my $tag_virome      = 0;
my $custom_phage    = '';
my $data_dir        = '/data';
my $n_cpus          = 16;
my $wdir            = cwd();

GetOptions(
   'f|fna=s'     => \$input_file,
   'd|dataset:s' => \$code_dataset,
   'db:i'        => \$choice_database,
   'virome:i'    => \$tag_virome,
   'wdir:s'      => \$wdir,
   'cp:s'        => \$custom_phage,
   'data-dir:s'  => \$data_dir,
   'ncpu:i'      => \$n_cpus,
   'h|help'      => \$help,
);

if ($help) {
    pod2usage();
}

unless ($input_file) {
    pod2usage('Missing FASTA file');
}

if ($choice_database < 1 || $choice_database > 3) {
    pod2usage('choice_database must be 1, 2, or 3');
}

say map { sprintf "%-15s: %s\n", @$_ } (
    ['Bin',           $Bin],
    ['Dataset',       $code_dataset],
    ['Input file',    $input_file],
    ['Db',            $choice_database],
    ['Working dir',   $wdir],
    ['Custom phages', $custom_phage],
);

if ($tag_virome == 1) {
    say "WARNING: THIS WILL BE A VIROME DECONTAMINATION RUN";
}

# Need 2 databases
# PCs from Refseq (phages) or PCs from Refseq+Viromes
# PFAM (26.0)

my $path_hmmsearch     = which('hmmsearch') or die "Missing hmmsearch\n";
my $path_blastp        = which('blastp')    or die "Missing blastp\n";
my $script_dir         = catdir($Bin, 'Scripts');
my $dir_Phage_genes    = catdir($data_dir,'Phage_gene_catalog');
my $readme_file        = catfile($data_dir, 'VirSorter_Readme.txt');
my $ref_phage_clusters = catfile($data_dir,
                         'Phage_gene_catalog', 'Phage_Clusters_current.tab');

if ($tag_virome == 1) {
    $readme_file = catfile($data_dir, 'VirSorter_Readme_viromes.txt');
}

my $generic_ref_file = catfile($data_dir,'Generic_ref_file.refs');

if ($choice_database == 2) {
    $dir_Phage_genes    = catdir($data_dir, 'Phage_gene_catalog_plus_viromes');
    $ref_phage_clusters = catfile($data_dir,
        'Phage_gene_catalog_plus_viromes', 'Phage_Clusters_current.tab');
}
elsif ($choice_database == 3) {
    $dir_Phage_genes    = catdir($data_dir, 'euk-virus');
    # ??? what goes here?  I don't have this file
    $ref_phage_clusters = catfile($data_dir,
        'euk-virus', 'Phage_Clusters_current.tab');
}

my $db_PFAM_a = catfile($data_dir, 'PFAM_27', 'Pfam-A.hmm');
my $db_PFAM_b = catfile($data_dir, 'PFAM_27', 'Pfam-B.hmm');

my $out = "";

## SETTING UP THE WORKING DIRECTORY
my $log_dir = catdir($wdir, 'logs');
if (-d $log_dir) {
## Commented on iPlant, but can be useful when running VirSorter on a directory already processed 
## (to avoid recomputing the gene prediction and comparison to PFAM especially)
#    $out = `rm -r $log_dir/* *.csv`; 
#    print "rm -r log* *.csv => $out\n";
} 
else {
    mkpath($log_dir);
}
my $log_out = catfile($log_dir, 'out');
my $log_err = catfile($log_dir, 'err');

# cp fasta file in the wdir
my $fastadir = catdir($wdir, 'fasta');
if ( !-d $fastadir ) {
    mkpath($fastadir);
    my $fna_file = catfile($fastadir, 'input_sequences.fna');
    open my $fa, '<', $input_file;
    open my $s1, '>', $fna_file;

    while (<$fa>) {
        chomp($_);
        if ( $_ =~ /^>(.*)/ ) {
            my $id = $1;
            $id =~ s/[\/\.,\|\s?!\*%]/_/g;
            my $new_id = $code_dataset . "_" . $id;
            say $s1 ">$new_id";
        }
        else {
            say $s1 $_;
        }

    }
    close $fa;
    close $s1;

    # detect circular, predict genes on contigs and extract proteins, as well
    # as filtering on size (nb genes) and/or circular
    my $nb_gene_th = 2; # At least two complete genes on the contig
    my $path_script_step_1 
        = catfile($script_dir,"Step_1_contigs_cleaning_and_gene_prediction.pl");
    my $cmd_step_1 
        = "$path_script_step_1 $code_dataset $fastadir $fna_file $nb_gene_th "
        . ">> $log_out 2>> $log_err";
    say "Step 0.5 : $cmd_step_1";
    `echo $cmd_step_1 >> $log_out 2>> $log_err`;
    $out = `$cmd_step_1`;
}

say "\t$out";

my $fasta_contigs_nett 
    = catfile($fastadir, $code_dataset . "_nett_filtered.fasta");
my $fasta_file_prots = catfile($fastadir, $code_dataset . "_prots.fasta");

# Match against PFAM, once for all
# compare to PFAM a then b (hmmsearch)
my $out_hmmsearch_pfama     = catfile($wdir, 'Contigs_prots_vs_PFAMa.tab');
my $out_hmmsearch_pfama_bis = catfile($wdir, 'Contigs_prots_vs_PFAMa.out');
my $cmd_hmm_pfama
    = "$path_hmmsearch --tblout $out_hmmsearch_pfama --cpu $n_cpus "
    . "-o $out_hmmsearch_pfama_bis --noali $db_PFAM_a $fasta_file_prots "
    . ">> $log_out 2>> $log_err";

say "Step 0.8 : $cmd_hmm_pfama";

`echo $cmd_hmm_pfama >> $log_out 2>> $log_err`;

if (!(-e $out_hmmsearch_pfama)) {
    $out = `$cmd_hmm_pfama`;
    say "\t$out";
}

my $out_hmmsearch_pfamb     = catfile($wdir, 'Contigs_prots_vs_PFAMb.tab');
my $out_hmmsearch_pfamb_bis = catfile($wdir, 'Contigs_prots_vs_PFAMb.out');
my $cmd_hmm_pfamb 
    = "$path_hmmsearch --tblout $out_hmmsearch_pfamb --cpu $n_cpus "
    . "-o $out_hmmsearch_pfamb_bis --noali $db_PFAM_b $fasta_file_prots "
    . ">> $log_out 2>> $log_err";
say "Step 0.9 : $cmd_hmm_pfamb";
`echo $cmd_hmm_pfamb >> $log_out 2>> $log_err`;

if (!(-e $out_hmmsearch_pfamb)) {
    $out = `$cmd_hmm_pfamb`;
    say "\t$out";
}

# Now work on the phage gene catalog

# Files that will stay along the computations
my $predict_file  = catfile($fastadir, $code_dataset . '_mga_final.predict' );
my $out_hmmsearch = catfile($wdir, 'Contigs_prots_vs_Phage_Gene_Catalog.tab');
my $out_hmmsearch_bis        
    = catfile($wdir, 'Contigs_prots_vs_Phage_Gene_Catalog.out');
my $out_blast_unclustered    
    = catfile($wdir, 'Contigs_prots_vs_Phage_Gene_unclustered.tab');
my $out_file_affi  
    = catfile($wdir, $code_dataset . '_affi-contigs.csv');
my $out_file_phage_fragments 
    = catfile($wdir, $code_dataset . '_phage-signal.csv');
my $global_out_file          
    = catfile($wdir, $code_dataset . '_global-phage-signal.csv');
my $new_prots_to_cluster     
    = catfile($wdir, $code_dataset . '_new_prot_list.csv');

# Constant scripts
my $script_merge_annot 
    = catfile($script_dir,"Step_2_merge_contigs_annotation.pl");
my $cmd_merge
    = "$script_merge_annot $predict_file $out_hmmsearch $out_blast_unclustered "
    . "$out_hmmsearch_pfama $out_hmmsearch_pfamb $ref_phage_clusters "
    . "$out_file_affi >> $log_out 2>> $log_err";

my $script_detect = catfile($script_dir, "Step_3_highlight_phage_signal.pl");
my $cmd_detect 
    = "$script_detect $out_file_affi $out_file_phage_fragments "
    . ">> $log_out 2>> $log_err";

if ($tag_virome == 1) {
    $cmd_detect 
        = "$script_detect $out_file_affi $out_file_phage_fragments "
        . "$generic_ref_file >> $log_out 2>> $log_err";
}

my $script_summary = catfile($script_dir, "Step_4_summarize_phage_signal.pl");
my $cmd_summary 
    = "$script_summary $out_file_affi $out_file_phage_fragments "
    . "$global_out_file $new_prots_to_cluster >> $log_out 2>> $log_err";

# # Get the final result file ready
`touch $global_out_file`;
my $r_n = -1;
# Si on a des nouvelles prots a clusteriser ou si on est dans la premiere
# revision
#
while ( (-e $new_prots_to_cluster || $r_n == -1) && ($r_n<=10) ) {
    $r_n++;    # New revision of the prediction
    my $dir_revision = catdir($wdir, 'r_' . $r_n);
    say "### Revision $r_n";

    if (!-d $dir_revision) {
        ## mkdir for this revision
        mkpath($dir_revision);
        say "Out : $out";

        ## Clustering of the new prots with the unclustered
        my $script_new_cluster 
            = catfile($script_dir, "Step_0_make_new_clusters.pl");
        
        # First revision, we just import the Refseq database
        if ( $r_n == 0 ) {
            #`mkdir $dir_revision/db`;
            mkpath(catdir($dir_revision, 'db'));

            ## Adding custom sequences to the database if required by the user
            if ($custom_phage ne '') {
                my $script_custom_phage = catfile(
                    $script_dir, "Step_first_add_custom_phage_sequence.pl"
                );
                $out = `$script_custom_phage $custom_phage $dir_Phage_genes/ $dir_revision/db >> $log_out 2>> $log_err`;

                say "Adding custom phage to the database : $out";
            }
            # should replace Pool_cluster / Pool_unclustered and
            # Pool_new_unclustered else , we just import the Refseq database
            else { 
                `cp $dir_Phage_genes/* $dir_revision/db/`; 
            }
        }
        else {
            my $previous_r = $r_n - 1;
            my $previous_fasta_unclustered =
              catfile($wdir, 'r_'. $previous_r, 'db', 'Pool_unclustered.faa');

            my $cmd_new_clusters = join(' ',
                "$script_new_cluster $dir_revision $fasta_file_prots",
                "$previous_fasta_unclustered",
                "$new_prots_to_cluster >> $log_out 2>> $log_err"
            );

            say $cmd_new_clusters;
            $out = `$cmd_new_clusters`;

            say "Step 1.1 new clusters and new database : $out";
            # Rm the list of prots to be clustered now that they should be
            # clustered
            $out = `rm $new_prots_to_cluster`;
            #print "rm $new_prots_to_cluster -> $out\n";
        }

        # Check if there are some data in these new clusters, or if all the new
        # proteins are unclustered
        my $new_db_profil = catfile($dir_revision, 'db', 'Pool_clusters.hmm');
        my $check = 0;

        if (-s $new_db_profil) {
            open my $DB, '<', $new_db_profil;

            while (<$DB>) {
                chomp($_);
                if ( $_ =~ /^NAME/ ) { 
                    $check++; 
                }
            }
            close $DB;
        }

        if ($check == 0) {
            say "There are no clusters in the database, so skip the hmmsearch";
        }
        else {
            my $out_hmmsearch_new =
              catfile($dir_revision, 'Contigs_prots_vs_New_clusters.tab');

            my $out_hmmsearch_bis_new =
              catfile($dir_revision, 'Contigs_prots_vs_New_clusters.out');

            my $cmd_hmm_cluster = join(' ',
                "$path_hmmsearch --tblout $out_hmmsearch_new --cpu $n_cpus",
                "-o $out_hmmsearch_bis_new --noali $new_db_profil",
                "$fasta_file_prots >> $log_out 2>> $log_err"
            );

            say "Step 1.2 : $cmd_hmm_cluster";

            `echo $cmd_hmm_cluster >> $log_out 2>> $log_err`;

            $out = `$cmd_hmm_cluster`;
            say "\t$out";

            $out = `cat $out_hmmsearch_new >> $out_hmmsearch`;
            say "\t$out";
        }

        my $out_blast_new_unclustered =
          catfile($dir_revision, 'Contigs_prots_vs_New_unclustered.tab');

        my $blastable_unclustered =
          catfile( $dir_revision, 'db', 'Pool_new_unclustered' );

        my $cmd_blast_unclustered = join(' ',
            $path_blastp,
            "-query $fasta_file_prots",
            "-db $blastable_unclustered",
            "-out $out_blast_new_unclustered",
            "-num_threads $n_cpus", 
            "-outfmt 6",
            "-evalue 0.001 >> $log_out 2>> $log_err"
        );

        say "\nStep 1.3 : $cmd_blast_unclustered";

        `echo $cmd_blast_unclustered >> $log_out 2>> $log_err`;
        $out = `$cmd_blast_unclustered`;

        say "\t$out";
        $out = `cat $out_blast_new_unclustered >> $out_blast_unclustered`;

        say "\t$out";

        # Make backup of the previous files to have 
        # trace of the different steps
        my $backup_affi = catfile($dir_revision, 'affi_backup.csv');
        my $backup_phage_signal =
          catfile($dir_revision, 'phage_signal_backup.csv');
        my $backup_global_signal =
          catfile($dir_revision, 'global_signal_backup.csv');

        if (-e $out_file_affi) {
            `cp $out_file_affi $backup_affi`;
        }

        if (-e $out_file_phage_fragments) {
            `cp $out_file_phage_fragments $backup_phage_signal`;
        }

        if (-e $global_out_file) {
            `cp $global_out_file $backup_global_signal`;
        }
    }

    ## Complete the affi
    say "Step 2 : $cmd_merge";
    `echo $cmd_merge >> $log_out 2>> $log_err`;
    $out = `$cmd_merge`; 
    ## This generate a csv table including the map of each contig, with PFAM
    #and Viral PCs annotations, as well as strand and length of genes

    say "\t$out";
    ## Complete the summary
    say "Step 3 : $cmd_detect";
    `echo $cmd_detect >> $log_out 2>> $log_err`;
    $out = `$cmd_detect`;
    say "\t$out";

    # Decide which contigs are entirely viral and which are prophages, and
    # which of both of these categories are phage enough to be added to the
    # databases
    say "Setting up the final result file";
    say "Step 4 : $cmd_summary";
    `echo $cmd_summary >> $log_out 2>> $log_err`;
    $out = `$cmd_summary`;
    say "\t$out";
}

# Last step -> extract all sequences as fasta files and gb
my $script_generate_output 
    = catfile($script_dir, 'Step_5_get_phage_fasta-gb.pl');

my $cmd_step_5 
    = "$script_generate_output $code_dataset $wdir >> $log_out 2>> $log_err";

say "\nStep 5 : $cmd_step_5";

`echo $cmd_step_5 >> $log_out 2>> $log_err`;

$out = `$cmd_step_5`;
say "\t$out";

# Plus clean the output directory
say "Cleaning the output directory";

# We rm the first db to not overload user disk space
my $db_revision_0 = catdir($wdir, 'r_0', 'db');
if (-d $db_revision_0) {
    $out = `rm -r $db_revision_0`;
    say "rm -r $db_revision_0 : $out";
}

#`mv $fastadir $wdir/Fasta_files`;

# We put all results from Hmmsearch and BLAST files in a separate directory
my $store_database_comparison = catdir($wdir, "Tab_files");
mkpath($store_database_comparison);
`mv $out_hmmsearch $store_database_comparison`;

# `mv $out_hmmsearch_bis $store_database_comparison/`;
`mv $out_blast_unclustered $store_database_comparison`;
`mv $out_hmmsearch_pfama $store_database_comparison`;
`mv $out_hmmsearch_pfama_bis $store_database_comparison`;
`mv $out_hmmsearch_pfamb $store_database_comparison`;
`mv $out_hmmsearch_pfamb_bis $store_database_comparison`;

#`mv error.log $log_dir`;
#`mv formatdb.log $log_dir`;
#
#my $final_error_log = catfile($log_dir, 'Virsorter_stderr_log');
#`mv log_err $final_error_log`;
#my $final_out_log = catfile($log_dir, 'Virsorter_stdout_log');
#`mv log_out $final_out_log`;

# We put all the files linked to the metric computation in a new directory
my $store_metric_files = catdir($wdir, 'Metric_files');

if (!-d $store_metric_files) {
    mkpath($store_metric_files);
}

`mv $out_file_affi $store_metric_files/VIRSorter_affi-contigs.tab`;
my $out_file_affi_ref = $code_dataset . "_affi-contigs.refs";
`mv $out_file_affi_ref $store_metric_files`;
`mv $out_file_phage_fragments $store_metric_files/VIRSorter_phage_signal.tab`;

if (-e $new_prots_to_cluster) {
    `mv $new_prots_to_cluster $store_metric_files/`;
}

# And we customize and add the readme file in the output directory
my $datestring        = localtime();
my $local_readme_file = catfile($wdir, 'Readme.txt');

open my $s1, '>', $local_readme_file;
say $s1 "VirSorter parameters used :\n";
say $s1 "--> Fasta file mined for viral sequences : $input_file";
say $s1 "--> Viral database used : ";

if ($choice_database == 2) {
    say $s1 join(' ', 
        "Viromes : all bacterial and archaeal virus genomes in Refseq,",
        "as of January 2014, plus non-redundant predicted genes from viral",
        "metagenomes (including seawater, freshwater, and human-related",
        "samples)"
    );
}
elsif ($choice_database == 3) {
    print $s1 "Eukaryotic";
}
else {
    say $s1 "RefseqABVir (all bacterial and archaeal virus genomes " .
        "in Refseq, as of January 2014)";
}

if ($custom_phage eq "") {
    say $s1 "--> No custom reference sequence was added to the database";
}
else {
    say $s1 "--> Custom reference sequences from fasta file $custom_phage " .
        "were added to the database";
}

if ($tag_virome == 1) {
    say $s1 join(' ',
        "VirSorter was run with the in the 'Virome Decontamination' mode:",
        "overall metrics for microbial sequences were not evaluated from the",
        "complete dataset, but instead pre-computed values based on bacterial",
        "and archaeal genomes from Refseq were used."
    );
}

say $s1 "This VirSorter computation finished on $datestring";
close $s1;

`cat $readme_file >> $local_readme_file`;

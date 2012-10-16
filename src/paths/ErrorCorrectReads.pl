#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# ErrorCorrectReads.pl
#
# Script to run the combination of RemoveDodgyReads, FindErrors and
# CleanCorrectedReads that make up the error correction phase of ALLPATHS-LG.
# Takes and returns a fastq file.
# Generates a read ID map (read names are not preserved)
# Optionally generates kmer spectra.

use strict;
use FindBin;
use File::Basename;
use File::Path;

# ---- Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
    (READS_IN              => { value => undef,
				help  => "Fastq file to error correct." },
     READS_OUT             => { value => undef,
				help  => "Error corrected fastq file." },
     PHRED_ENCODING        => { value => undef,
                                help  => "Fastq Phred quals are encoded using an offset of either 33 or 64" },
     THREADS               => { value => "max",
                                help  => "Number of threads to use in error correction" },
     MAX_MEMORY_GB         => { value => "0",
                                help  => "Restrict memory usage - default is to use it all" },
     PLOIDY                => { value => "2",
                                help  => "The ploidy. Used for evaluation purposes only." },
     FORCE_PHRED           => { value => "False",
                                help  => "True: accepts specified PHRED encoding. False: Aborts is incorrect PHRED detected." },
     REMOVE_DODGY_READS    => { value => "True",
                                help  => "Run RemoveDodgyReads as part of the ALLPATHS-LG error correction pipeline." },
     KEEP_KMER_SPECTRA     => { value => "False",
                                help  => "Save the kmer spectra." },
     KEEP_INTERMEDIATES    => { value => "False",
                                help  => "Whether or not to keep intermediate files." });


# ---- Validate arguments

my $phred_encoding = $args{PHRED_ENCODING};
if ($phred_encoding != 33 && $phred_encoding != 64) {
  die "Invalid value for PHRED_ENCODING. Must be either 33 or 64\n";
}

my $ploidy = $args{PLOIDY};
if ($ploidy != 1 && $ploidy != 2) {
  die "Invalid value for PLOIDY. Must be either 1 or 2\n";
}

# ---- Determine filenames and directories to use

my($fastq_head_out, $output_dir, $suffix) = fileparse($args{READS_OUT}, ".fastq");
my $fastq_in = $args{READS_IN};
my $fastq_out = "$output_dir$fastq_head_out.fastq";
my $lg_ref_name = "$fastq_head_out.allpaths-lg";
my $lg_ref_dir = $output_dir . $lg_ref_name;
my $lg_data_dir = "$lg_ref_dir/data";
my $lg_run_dir = "$lg_data_dir/run";

# ---- Make sure input FASTQ file exists

abort("Fastq file '$args{READS_IN}' does not exist.") unless (-e $args{READS_IN});


# ---- Create allpaths-lg pipeline directory

if (-e $lg_ref_dir) {
  die "Found existing allpaths-lg pipeline directory: " .
  $lg_ref_dir . "\nPlease erase to continue.\n";
}

mkdir $lg_ref_dir;
mkdir $lg_data_dir;

my $cmd;

# ---- Setup Ploidy

run_or_die("echo $args{PLOIDY} > $lg_data_dir/ploidy");

# ---- Create dummy jump files to fool ALLPATHS-LG

run_or_die("touch $lg_data_dir/jump_reads_orig.fastb");
run_or_die("touch $lg_data_dir/jump_reads_orig.qualb");
run_or_die("touch $lg_data_dir/jump_reads_orig.pairs");


# ---- Run FastqToFastbQualb

print "\nRunning FastqToFastbQualb.\n";
$cmd = ("$FindBin::Bin/FastqToFastbQualb" .
	" FASTQ=$fastq_in" .
	" OUT_HEAD=$lg_data_dir/frag_reads_orig" .
	" PHRED_64=" . ($args{PHRED_ENCODING} == 64 ? "True" : "False") .
	" FORCE_PHRED=" . ($args{FORCE_PHRED} == 0 ? "False" : "True") );
run_or_die($cmd);

# ---- Run PairsFake

print "\nRunning PairsFake.\n";
$cmd = ("$FindBin::Bin/PairsFake" .
	" HEAD=$lg_data_dir/frag_reads_orig" );
run_or_die($cmd);


# ---- Run ALLPATHS-LG

print "\nRunning ALLPATHS-LG.\n";
$cmd = ("$FindBin::Bin/RunAllPathsLG" .
	" PRE=$output_dir" .
	" REFERENCE_NAME=$lg_ref_name" .
	" DATA_SUBDIR=data" .
	" RUN=run" .
	" THREADS=$args{THREADS}" .
	" MAX_MEMORY_GB=$args{MAX_MEMORY_GB}" .
	" VALIDATE_INPUTS=False" .
	" TARGETS=error_correction" .
	" REMOVE_DODGY_READS_FRAG=" . ($args{REMOVE_DODGY_READS} != 0 ? "True" : "False")  );
run_or_die($cmd);

# ---- Run FastbQualbToFastq

print "\nRunning FastbQualbToFastq.\n";
$cmd = ("$FindBin::Bin/FastbQualbToFastq" .
	" HEAD_IN=$lg_run_dir/frag_reads_corr" .
	" HEAD_OUT=$output_dir$fastq_head_out" .
	" NAMING_PREFIX=read" .
	" PAIRED=false" .
	" PHRED_OFFSET=$phred_encoding" );
run_or_die($cmd);

# ---- Run ReadTrack

print "\nRunning ReadTrack. (silent)\n";
$cmd = ("$FindBin::Bin/ReadTrack" .
	" NH=True" .
	" READS=$lg_run_dir/frag_reads_corr" .
	" TRACK_BACK=True" .
	" > $lg_run_dir/frag_reads_corr.tracked");
run_or_die($cmd);

# ---- Extract read ID tracking information

print "\nParsing read tracking information. (silent)\n";
parse_readtrack("$lg_run_dir/frag_reads_corr.tracked", "$fastq_out.ids");

# ---- Keep kmer spectra

if ($args{KEEP_KMER_SPECTRA} != 0) {
  mkdir "$fastq_out.kspec";
  run_or_die("cp $lg_run_dir/*.kspec $fastq_out.kspec");
}

# ---- Clean up

if ($args{KEEP_INTERMEDIATES} == 0) {
  rmtree($lg_ref_dir);
}

print "\nDone.\n";




sub parse_readtrack {

  my ($file_in, $file_out) = @_;

  my $input;
  open ($input, "<", $file_in) or die "Unable to open: $file_in\n";

  my $output;
  open ($output, ">", $file_out) or die "Unable to write to: $file_out\n";

  print $output "#New_ID,  Original_ID\n";

  my $line;
  while ($line = <$input>) {
    chomp $line;
    if ($line =~ /^frag_reads_corr:(\d+)\s+.*:(\d+)\s*$/ ) {
      print $output "$1,  $2\n";
    }
  }

  close $output;
  close $input;
  
}

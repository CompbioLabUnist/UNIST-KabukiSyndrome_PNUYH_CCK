#!/usr/bin/perl

use strict;
use File::Basename;
use Cwd;

die " Usage : " .__FILE__. " < sorted bam dir > < output dir > \n" if @ARGV != 2;

my $bam_dir = $ARGV[0];
my $MethylKit_dir = $ARGV[1];

my @bam_files = glob("$bam_dir/sample_ID/*.st.bam");

foreach my $bam_file (@bam_files) {
	my $sample_id = basename(dirname($bam_file));

	&check_dir("$MethylKit_dir/$sample_id");
	my $out_dir = "$MethylKit_dir/$sample_id";
	&check_dir("$out_dir/sh");

	open (SH,">$out_dir/sh/methylkit_bismarkAln_$sample_id.sh");
	print SH "#!/bin/sh\n";
        print SH "#SBATCH --ntasks=1 --cpus-per-task=2 --mem-per-cpu=5g\n";
        print SH "#SBATCH -o $out_dir/sh/methylkit_bismarkAln_$sample_id.sh.%j\n";

	print SH "Rscript /path/to/your/bismarkAln_script/2_2_bismarkAln.R --input $bam_file --sample_id $sample_id --output $out_dir\n";
	close SH;

	system ("sbatch $out_dir/sh/methylkit_bismarkAln_$sample_id.sh");
}

exit;

sub check_dir   {
        my $temp_dir = shift;

        unless (-d $temp_dir)   {
                system ("mkdir $temp_dir");
        }
}


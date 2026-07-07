#!/usr/bin/perl

use strict;
use File::Basename;
use Config::Simple;

die "Usage : " .__FILE__. " < config file > < bam file dir > < output dir > \n" if @ARGV != 3;


#####	set parameter
my ($cfg_file, $bam_dir, $out_dir) = @ARGV;


my %config = ();
Config::Simple->import_from($cfg_file, \%config);

my $java = $config{"TOOLS.JAVA"};
my $gatk = $config{"TOOLS.GATK_4.3.0"};
my $ref_file = $config{"REF.REF_FA_hg38"};
my $dbsnp = $config{"REF.DBSNP_hg38"};

$out_dir = &check_dir("$out_dir");
$out_dir = &check_dir("$out_dir/GVCF");
$out_dir = &check_dir("$out_dir/WGS");


#####	Run  
my @file_array = <$bam_dir/*\/*.rc.bam>;

foreach my $file (@file_array)  {
	my ($file_name, $file_dir) = fileparse($file);
	my $id = basename($file_dir);
  print "$id\n";
	&check_dir("$out_dir/$id");
	&check_dir("$out_dir/$id/sh");
	&check_dir("$out_dir/$id/tmp");
  
  
  for(my $chr_num = 1; $chr_num < 26; $chr_num++) {
                my $chr = $chr_num;
                if($chr_num == 23) {$chr = "X";}
                if($chr_num == 24) {$chr = "Y";}
		if($chr_num == 25) {$chr = "M";}
		
  my $sh_file = "$out_dir/$id/sh/HaplotypeCaller_$id.chr$chr.sh";
  
	#$thread--;
	open (SH,">$sh_file");
	print SH "#!/bin/sh\n";
	print SH "#SBATCH --ntasks=1 --cpus-per-task=1 --mem-per-cpu=10g\n";
	print SH "$java -Xmx10g ";
	print SH "-XX:ParallelGCThreads=1 ";
	print SH "-jar $gatk ";
	print SH "HaplotypeCaller \\\n";
	print SH "-R $ref_file \\\n";
	print SH "-I $file \\\n";
	print SH "-L chr$chr ";
	print SH "-O $out_dir/$id/$id.chr$chr.gvcf.gz ";
	print SH "-ERC GVCF \\\n";
	print SH "--dbsnp $dbsnp \\\n";
	print SH "1> $out_dir/$id/sh/$id.chr$chr.gvcf.gz.stdout ";
	print SH "2> $out_dir/$id/sh/$id.chr$chr.gvcf.gz.stderr ";
	print SH "\n";
	close SH;

  system ("sbatch $sh_file\n");
	}
}

exit;

sub check_dir	{
	my $temp_dir = shift;

	unless (-d $temp_dir)	{
		system ("mkdir $temp_dir");
	}
	return $temp_dir;
}

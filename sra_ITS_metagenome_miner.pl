#!/usr/bin/perl

# This program is used to search for specific ITS (or other) sequences in the SRA database using a mapping approach bwa
# Requires sratoolkit, bwa, samtools, ITSx, mafft

# Expected input is a fasta file of ITS sequences
# A list of SRA data in a text format with one experiments (acc. num.) per line.

# Make sure ITSx is in your $PATH or modify below

# In order not to fill up your home directory you need to redirect caching by sratoolkit or disable it 
# Code from https://standage.github.io/that-darn-cache-configuring-the-sra-toolkit.html
# mkdir -p ~/.ncbi
# echo '/repository/user/main/public/root = "/scratch/standage/sra-cache"' > ~/.ncbi/user-settings.mkfg
# Uncomment the next command if you want to disable network access altogether
# echo '/repository/user/cache-disabled = "true"' > ~/.ncbi/user-settings.mkfg
# I disabled caching completely

# open input and output files and 3rd argument which is y or n for whether it is a TLS (targeted locus study)
$inputfasta = $ARGV[0];
$inputSamples = $ARGV[1];

# Directory for the fastq data
$fastdirectory = '/scratch/tyjames_root/tyjames0/tyjames/Mala/03.07.23/temp3/data';

# Name of log file for fasterq-dump
$dumpfilelog = 'fasterq_dump.log';

if (@ARGV < 1) {
    die "Improper number of arguments. Usage: perl sra_ITS_miner.pl inputfastafile inputSRAlist\n";
}

##############################
######### RUN ITSX ########### 
##############################

# name for ITSx output file
$splitfasta = 'ITSx_out.concat.fasta'; 

my $cmd = 'ITSx -i ' . $inputfasta . ' -t Fungi';
print $cmd;
# ITSx -i targets_unsplit.fasta -o targets_ITSx -t Fungi

system($cmd);

# ITSx puts the sequences separately into two files ITSx_out.ITS1.fasta and ITSx_out.ITS2.fasta
$cmd = 'cat ITSx_out.ITS1.fasta ITSx_out.ITS2.fasta >ITSx_out.concat.fasta';
system($cmd);

# Next we will need to create an indexed file for bwa
$cmd = 'bwa index -p ITSx_bwa_reference ITSx_out.concat.fasta';
system($cmd); 

# Extract sequence names from fasta file
my @seqnames;
open (INPUTFASTA, $splitfasta) or die $!;

while (<INPUTFASTA>) {
	if ( /^>/ ) {
        # get the ID from the header
	my @header = split /\s/, $_;
        $header[0] =~ s/\>//;
	push (@seqnames, $header[0]);
    }
}

$numseqs = scalar(@seqnames);

print "\n$numseqs query sequences parsed. End ITSx.\n\n";
# close the file
close INPUTFASTA;


##############################
## PRINT OUT BIOSAMPLE DATA ##
##############################

open (OUTPUTTABLE, ">sra_ITS_miner_output.txt");

print OUTPUTTABLE "SRA_Expt\tData_size\tRead_length";
foreach (@seqnames) {
 	print OUTPUTTABLE "\t$_";
# also open file with the filename using the sequence name to dump the hits
# as well as a file with the metadata for each sequence search
	$outfilename = $_;
        $outfilename =~ s/\|/\_/g;
	$outfilename = $outfilename . '.fasta_list';
	open(FH, '>', $outfilename) or die $!;
}

print OUTPUTTABLE "\n";

# Open sample data set and read in line by line each study
open (SAMPLETABLE, $inputSamples) or die $!;

# Go through table of samples
while (<SAMPLETABLE>)
{
	$experiment = $_;
	chomp($experiment);
	$experiment =~ s/\s+$//;	
	print OUTPUTTABLE "$experiment\t";

##############################
###### DOWNLOAD FASTQ/A ######
##############################

# Use sratoolkit to download fastq file
# Reads are not treated as paired and separate, even if they are
# This may lead to double counting. May need to be dealt with later

# name for the prefetch file

	$srafile = $fastdirectory . '/' . $experiment . '/' . $experiment . '.sra'; 
	
# First download the data using prefetch

	until (-e $srafile)
	{
		$cmd = "prefetch -v -p -X 500g -f all -C yes " . $experiment . " -O " . $fastdirectory;
		print "$cmd\n";
		system($cmd);
	}

# name for the fastq/a file
	$fastfile = $experiment . '.fastq'; # currentlyused

# Next download the reads via fasterq-dump
# -e = num threads
# -f to overwrite existing file

	$cmd = "fasterq-dump --split-spot -e 8 -f -o " . $fastfile . ' ' . $srafile . " 2>" . $dumpfilelog;
	print "$cmd\n";
	system($cmd);			
	
# Now we will get the number of sequences and the length of each read

	open (DUMPLOGFILE, $dumpfilelog);
	while (<DUMPLOGFILE>)
	{
		print $_;
		if ($_ =~ /reads read/)
		{
			@data = split /\s+/, $_;
			$readnum = $data[3];
		} 
	}
	print "readnum=$readnum\n";

	$cmd = "head -n 2 " . $fastfile . " | tail -1 | wc -c";
	my $seqlength = `$cmd`;
	chomp($seqlength);
	print "seqlength=$seqlength\n";

# Check filesize which could prevent large files from being searched.  Uses a lot of memory	

	my $filesize = -s $fastfile;

# Make a hash of each seq name & set count for each to zero for each experiment
	foreach (@seqnames)
        {
                $count_hits{$_} = 0;
        }

# Now print the data size and read length in output table

    print OUTPUTTABLE "\t$readnum\t$seqlength";
    print "\nFinished downloading\n\n";
		
##############################
############ BWA #############
##############################

# Now use bwa to map reads to indexed file with the concatenated 5.8S removed sequences
# First checking if the file isn't massive, not used for this script so much
	if ($filesize < 1e12)
	{
# bwa -t = threads
# bwa -L = clipping penalty. High numbers mean that match needs to cover most of the read
# bwa -T = Threshold for reporting score.
# samtools -F 4 reports only matches
		$samfile = $experiment . '.sam';
		$cmd = 'bwa mem -t 8 -L 100 -T 32 -O 3 ITSx_bwa_reference ' . $fastfile . ' | samtools view -F 4 > ' . $samfile;
		print $cmd;
		system($cmd);
# Now open samfile output
		open (SAMF, $samfile);
# Read in each line and check if it matches a hit
		while (<SAMF>)
		{
			$line = $_;
			@targethit = split (/\t/, $_);
			foreach (@seqnames)
			{
                $tempoutfilename = $_;
                # Replace hashes with underscores
                $tempoutfilename =~ s/\|/\_/g;
                $outfastalistfilename = $tempoutfilename . '.fasta_list';
				open (FLFILE, '>>', $outfastalistfilename) or die $!;
				$name = $_;
# \Q allows match literal regardless of weird characters
				if ($line =~ /\Q$name/)
				{
					$count_hits{$_}++;
# Pulls fasta name from sam file, not used;
					$seqdesired = $targethit[1];
					print FLFILE "$experiment\t$line";
				} 				
			}
		}
	}
	else 
	{
		$analyzed = "no";
	}
	
	print "\nFinished bwa\n\n";

##############################
### PRINT COUNTS AND SEQS ####
##############################

# Print # observations of each sequence

	if ($analyze eq 'no')
	{
		foreach (@seqnames)
		{
			print OUTPUTTABLE "\tNA";
		}
	}
	else
	{
		foreach	(@seqnames)
		{
			print OUTPUTTABLE "\t$count_hits{$_}";
		}
# Close line in outputtable
        }
        print OUTPUTTABLE "\n";

# Can remove fastq and sra files to save space
	$cmd = "rm " . $fastfile;
	system ($cmd);
	$cmd = "rm " . $srafile;
	system ($cmd);
}


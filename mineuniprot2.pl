#! /usr/bin/perl -w

# mineuniprot.pl

# Copyright (C) 2023  John C. Newman
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA
#
# http://www.fsf.org/licenses/gpl.txt
#
######################################################

# This program mines and tabulates specific information from a uniprot database file (.dat)
# based on the uniprot accessions in an input file
# v1 mines interprot protein domain data
# v2 uses simpler input file format, runs a folder of input files, and integrates output

# Input directory is specified in the command-line argument
# Input file format:
# Column 0: unique uniprot accession
# Column 1: previously derived gene symbol (not used here)
# any other columns are ignored
$inputdir = $ARGV[0];

# Other necessary files specified here:
$uniprotdbfile = "../uniprot_sprot_rodents.dat";
#$uniprotdbfile = "../uniprot_sprot_dat_test";
#$uniprotdbfile = "../uniprot_sprot_dat_hmcs2.txt";

# Species ID for pulling out background list
# in sprot: "OX   NCBI_TaxID=10090;"
$speciesid = 10090;

# Read in the input
# get file list from the inputdir
@inputfiles = split(/[\t\n]/,`ls $inputdir`);

foreach $input (@inputfiles) {
	chomp $input;
	# we will use filename as an internal identifier, but convert any spaces to underscores first
	my $inputname = $input;
	$inputname =~ s/\s/\_/g;
	# Creates %input with key uniprotID
	open (INPUT,"$inputdir/$input") or die "Can't open input file $inputdir/$input\n";
	while (<INPUT>) {
		chomp $_;
		$line = $_;
		@elem = split(/\t/,$line);
		$input{$inputname}{$elem[0]} = 0;
		#print STDOUT ("Found $elem[0] $elem[1]\n"); # test code
	} # end while INPUT
	#print STDOUT ("Found file $input name $inputname\n"); #test code
	close INPUT;

} # end foreach @inputfiles


# Test code
#foreach $acc (keys %input) {
#	print STDOUT ("Found $acc\n"); # test code
#}

# Process Uniprot Database File
# to optimize speed, we'll first mine the entire dat file for the data of interest,
# then query our list against the internal data structrue
# I can reuse most of this code later for mining other data elements from the .dat

open (UNIPROTDB,$uniprotdbfile) or die "Can't open Uniprot dat file $uniprotdbfile\n";
my $readyfornewid = 1;
my @acclist = (); # temp array filled with each entry's set of accessions
my %updata = (); # keeps info with key accession -> subhash ("id" -> uniprot ID; "aa" -> aa length)
my $id = "";
my $aa = 0;
my %allinterpro = (); # a master list of all interpro ids in the database
while (<UNIPROTDB>) {
	chomp $_;
	$line = $_;
	# when we see the record separator line // we're ready for new uniprot accession(s)
	if ($line =~ /^\/\//) {
		$readyfornewid = 1;
		# test code
		#print STDOUT ("Found accessions:");
		#foreach $acc (@acclist) {
		#	print STDOUT ("\t$acc");
		#}
		#print STDOUT ("\n"); # end test code
	}
	if ($readyfornewid == 1) {
		# pull out the ID
		if ($line =~ /^ID\s+(\S+)/) {
			$id = $1;
			$readyfornewid = 0;
		}
		# separately pull out the aa length
		if ($line =~ /^ID.+\s(\d+)\s+AA/) {
			$aa = $1;
		}
	}
	
	# AC accessions
	# may have multiple accessions separated by semicolons and spaces like "AC   Q9CQV8; O70455; Q3TY33; Q3UAN6;"
	if ($line =~ /^AC\s+(.+)\;$/) {
		# note this tries to leave out the trailing semicolon to prevent a blank accession after split
		my $accentry = $1;
		# empty the accession list and refill
		@acclist = ();
		@acclist = split (/\;\s+/,$accentry);
		# backfill %updata by accession with id and aa from the prior ID line
		foreach $acc (@acclist) {
			$updata{$acc}{"id"} = $id;
			$updata{$acc}{"aa"} = $aa;
		}
	}
	
	# Species
	if ($line =~ /^OX.+\=(\d+)/) {
		#print STDOUT ("species $1\n"); # test code
		# if the species id matches...
		if ($speciesid == $1) {
			# ... add the first unniprot accession to the background subhash of input
			# note that from here the background is treated like other input lists
			$input{"_background"}{$acclist[0]} = 0;
		}
	}
	
	# DR interpro, add as subhash keys for each accession for this protein
	if ($line =~ /^DR\s+InterPro\;\s+(\S.+)\.$/) {
		my $ip = $1;
		foreach $acc (@acclist) {
			$updata{$acc}{"interpro"}{$ip} = 0;
		}
		# also add to the master list
		$allinterpro{$ip} = 0;
	}
	
}
close UNIPROTDB;

# test code
#foreach $acc (sort keys %updata) {
#	print STDOUT ("$acc\t" . $updata{$acc}{"id"} . "\t" . $updata{$acc}{"aa"} . "\n");
#	if (exists $updata{$acc}{"interpro"}) {
#		foreach $i (sort keys %{$updata{$acc}{"interpro"}}) {
#			print STDOUT ("$i\n");
#		}
#	}
#}

# test code
#foreach $ip (sort keys %allinterpro) {
#	print STDOUT ("$ip FOUND IN:");
#	foreach $acc (sort keys %updata) {
#		if (exists $updata{$acc}{"interpro"}{$ip}) {
#			print STDOUT (" $acc");
#		}
#	}
#	print STDOUT ("\n");
#}

# test code
#foreach $id (sort keys %mouse) {
#	print STDOUT ("Mouse protein $id\n");
#}

# Process input sequences
# Search domains
# %outputdata collects output data for all input lists (keys listnames)
# %outputcounts collects the total and found counts for each intput list (keys listnames, subkeys count and found)
# For now hard code the variables to search eg interpro (but could automate this if lots)
foreach $name (keys %input) {
	$outputcounts{$name} = ();
	$outputdata{$name} = ();
	foreach $id (sort keys %{$input{$name}}) {
		$outputcounts{$name}{"count"}++;
		# skip if not present
		unless (exists $updata{$id}) { 
			next; 
		}
		$outputcounts{$name}{"found"}++;
		foreach $dom (sort keys %{$updata{$id}{"interpro"}}) {
			$outputdata{$name}{"interpro"}{$dom}++;
		}
	}
}

# Integrated output
# Creates a single table with rows as interproIDs and columns as freq in bkgd and then input files

# Header rows
# Row 0: col names
# Row 1: first col is blank (interpro) and subsequent are alpha-sorted total counts
# Row 2: first col is blank (interpro) and subsequent are alpha-sorted found counts
print STDOUT ("Interpro");
foreach $name (sort keys %outputcounts) {
	print STDOUT ("\t$name");
}
print STDOUT ("\n");
print STDOUT ("Count");
foreach $name (sort keys %outputcounts) {
	if (exists $outputcounts{$name}{"count"}) {
		print STDOUT ("\t" . $outputcounts{$name}{"count"});
	} else {
		print STDOUT ("\t0");
	}
}
print STDOUT ("\n");
print STDOUT ("Found");
foreach $name (sort keys %outputcounts) {
	if (exists $outputcounts{$name}{"found"}) {
		print STDOUT ("\t" . $outputcounts{$name}{"found"});
	} else {
		print STDOUT ("\t0");
	}
}
print STDOUT ("\n");


# interpro rows
# this uses _background interpro domains as the reference for organizing the output
foreach $id (sort keys %{$outputdata{"_background"}{"interpro"}}) {
	# print the name of the interpro domain
	print STDOUT ("$id");
	# pull the relevant count from each input data list and print one per col in the output
	foreach $name (sort keys %outputdata) {
		if (exists $outputdata{$name}{"interpro"}{$id}) {
			print STDOUT ("\t" . $outputdata{$name}{"interpro"}{$id});
		} else {
			# catch if the domain didn't come up for that intput list
			print STDOUT ("\t0");
		}
	}
	# newline between interpro domains
	print STDOUT ("\n");
}
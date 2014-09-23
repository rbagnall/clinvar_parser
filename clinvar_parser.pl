#!/usr/bin/perl

# clinvar_parser.pl by Richard Bagnall
# use to format Clinvar pathogenicity predictions into a hash table
# puts mutliple predictions for the same nucleotide into one line
# usage clinvar_parser.pl variant_summary.txt
# ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

use strict; use warnings;

# make array of single nucleotide variants, sorted by position
my @snps = `grep \"GRCh37\" <$ARGV[0] \| grep \"single nucleotide variant\" | awk -F \"\t\" '{print\$14"\t"\$15"\t"\$9"\t"\$6}' | sort -k1,1V -k2,2n`;

# make array of indels, sorted by position
my @indels = `grep \"GRCh37\" <$ARGV[0] \| grep \"deletion\\|duplication\\|indel\\|insertion\" | awk -F \"\t\" '{print\$14"\t"\$15"\t"\$9"\t"\$6}' | sort -k1,1V -k2,2n`;

# combine multiple entries at same nucleotide into single entry

# make empty arrays
my @firstline = (); 
my @nextline = ();
my $chr;
my $current_start;

# make a snp_anno file
# merges events if they start at the same residue
open(SNP_ANNO, ">>clinvar_snp.anno") or die "error reading clinvar_snp.anno for reading";

# parse the snp array
foreach (@snps) {
	chomp;

# if @firstline is empty, populate @firstline with current line of the array and goto next
    if (!@firstline) {(@firstline = split("\t", $_));
        $chr=$firstline[0];
        $current_start = $firstline[1];
        next;
    }
    
# otherwise, populate nextline with <line>
    else {@nextline = split("\t", $_);
    }
    
# if $chr, $current_start and $current_end is the same as nextline, join them together
    if (($chr eq $nextline[0])&&($current_start == $nextline[1])) {
        $firstline[2]="$firstline[2];$nextline[2]";
        $firstline[3]="$firstline[3];$nextline[3]";
        next;
    }
    
#else populate snp_anno file and nextline becomes firstline
    
    else {print SNP_ANNO (join("\t", @firstline),"\n");
        @firstline=@nextline;
        $chr=$firstline[0];
        $current_start = $firstline[1];
        next;
    }
}

# make an indel_anno file
# merges different events if they start at the same nucleotide
open(INDEL_ANNO, ">>clinvar_indels.anno") or die "error reading clinvar_indel.anno for reading";

# parse the indel array
foreach (@indels) {
	chomp;
    
    # if @firstline is empty, populate @firstline with current line of the array and goto next
    if (!@firstline) {(@firstline = split("\t", $_));
        $chr=$firstline[0];
        $current_start = $firstline[1];
        next;
    }
    
    # otherwise, populate nextline with <line>
    else {@nextline = split("\t", $_);
    }
    
    # if $chr, $current_start and $current_end is the same as nextline, join them together
    if (($chr eq $nextline[0])&&($current_start == $nextline[1])) {
        $firstline[2]="$firstline[2];$nextline[2]";
        $firstline[3]="$firstline[3];$nextline[3]";
        next;
    }
    
    #else populate indel_anno file
    
    else {print INDEL_ANNO (join("\t", @firstline),"\n");
        @firstline=@nextline;
        $chr=$firstline[0];
        $current_start = $firstline[1];
        next;
    }
}

close SNP_ANNO;
close INDEL_ANNO;


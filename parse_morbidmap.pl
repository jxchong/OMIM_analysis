#!/usr/bin/env perl
#
# Description:
#
#
#
# Created by Jessica Chong on 2015-02-22.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($inputfile, $outputfile, $help);

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $inputfile) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
}


# For the file morbidmap, the fields are, in order:
# 1  - Disorder, <disorder MIM no.> (<phene mapping key>)
# 2  - Gene/locus symbols
# 3  - Gene/locus MIM no.
# 4  - cytogenetic location


open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
print $output_handle "phenoname\tphenoMIMnum\tphenoMappingKey\tLocusSymbols\tGeneMIMnum\tCytoLoc\tisComplex\n";
open (my $input_handle, "$inputfile") or die "Cannot read $inputfile: $!.\n";
while ( <$input_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
    if ($_ =~ /^#/) {
        next;
    }
	my ($phenoname, $locussymbol, $locusMIM, $cytoloc) = split(/\t/, $_);
	my ($phenoMIM, $phenomappingkey) = qw(NA -9);
	
	if ($phenoname =~ m/(\d{6})/) {
		$phenoMIM = $1;
	}
	if ($phenoname =~ m/\((\d)\)$/) {
		$phenomappingkey = $1;
	}	

	if ($phenomappingkey == 4) {
		# if a chromosomal del/dup syndrome, the gene MIM *is usually* the phenotype MIM
		if ($phenoMIM eq 'NA')  {
			$phenoMIM = $locusMIM;
		}
	} 
	
	# if no phenotype MIM number, the gene MIM *is usually* the phenotype MIM
	# if ($phenoMIM eq 'NA')  {
	# 	$phenoMIM = $locusMIM;
	# }
	
	my $isComplex = "no";
	# remove the following phenotypes:
	# QTL or quantitative trait locus
	# suscep* (susceptibility but susceptibility is spelled incorrectly in OMIM in a few entries)
	# risk
	# [] or {}
	# !!! do not remove somatic but should flag
	if ($phenoname =~ "risk" || $phenoname =~ m/QUANTITATIVE TRAIT LOCUS/i || $phenoname =~ m/QTL/i || $phenoname =~ /\[/ || $phenoname =~ /\{/ || $phenoname =~ m/suscep(\w+) to/) {
		$isComplex = "yes";
	}
	if ($phenoname =~ "somatic" && $phenoname =~ /carcinoma|cancer|tumor|leukemia|lymphoma|sarcoma|blastoma|adenoma|cytoma|myelodysplastic|Myelofibrosis|oma,/i ) {
		$isComplex = "cancer";
	} elsif ($phenoname =~ "somatic" && $isComplex ne 'yes') {
		$isComplex = "somatic";
	}
		
	print $output_handle "$phenoname\t$phenoMIM\t$phenomappingkey\t$locussymbol\t$locusMIM\t$cytoloc\t$isComplex\n";
}
close $input_handle;
close $output_handle;




################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


parse_morbidmap.pl - Parse morbidmap into a more useful form!


=head1 SYNOPSIS


perl B<xxxx.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--in> F<input file>	

	input file

=item B<--out> F<output file>

	name of output file

=item B<--help> I<help>

	print documentation

=back


=head1 FILES


xx


=head1 EXAMPLES


xxxxxxxx


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut

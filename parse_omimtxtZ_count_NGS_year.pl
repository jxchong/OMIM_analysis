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



open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
print $output_handle "MIMnum\tmimtitle\tMIMlink\tisComplex\tmentionsNGS\tmentionsNGSparagraph\tNGSyear\tyearDiscovered\n";

my ($mentionsNGS, $NGSparagraph, $NGSyear, $yeardiscovered) = ('no', 'NA', 'NA', 'NA');
my ($inrecord, $inMG) = (0, 0);
my $printentry = '';
my @MGtxt;
my $MGtxt_paragraph = '';
my ($mimnum, $mimtitle);

my $input_handle;
if ($inputfile =~ /\.Z$/) {
	open ($input_handle, "zcat $inputfile |") or die "Cannot read $inputfile: $!.\n";	
} else {
	open ($input_handle, "$inputfile") or die "Cannot read $inputfile: $!.\n";
}
while ( my $readline = <$input_handle> ) {
    # print STDERR "readline = $readline\t";
	if ($readline =~ /^\*RECORD\*/ && $inrecord == 1) {
		# if we are finishing a record/OMIM entry, then print stored info
		print $output_handle "$printentry\t";
		print $output_handle "$mentionsNGS\t$NGSparagraph\t$NGSyear\t$yeardiscovered\n";
		$printentry = '';
		$mentionsNGS = 'no';
		$NGSparagraph = 'NA';
		$NGSyear = 'NA';
		$yeardiscovered = 'NA';
		undef @MGtxt;
		$MGtxt_paragraph = "";
		$inMG = 0;
		$inrecord = 0;
		undef $mimnum;
		undef $mimtitle;
	}
	
	if ($readline =~ /^\*RECORD\*/) {			
		my $mimnum_expected = <$input_handle>;		# skip FIELD NO line
		my $mimnum = <$input_handle>;
		$mimnum =~ s/\s+$//;
		if ($mimnum_expected !~ /\*FIELD\* NO/) {
			die "line expected to contain *FIELD* NO, instead contains: $mimnum_expected\nnext line: $mimnum\n";
		}
		
		my $mimtitle_expected = <$input_handle>;		# skip FIELD TI line
		my $mimtitle = <$input_handle>;
		$mimtitle =~ s/\s+$//;
		if ($mimtitle_expected !~ /\*FIELD\* TI/) {
			die "line expected to contain *FIELD* TI, instead contains: $mimtitle_expected\nnext line: $mimtitle\n";
		}
		
		# An asterisk (*) before an entry number indicates a gene.
		#
		# A number symbol (#) before an entry number indicates that it is a descriptive entry, usually of a phenotype, and does not represent a unique locus. The reason for the use of the number symbol is given in the first paragraph of the entry. Discussion of any gene(s) related to the phenotype resides in another entry(ies) as described in the first paragraph.
		#
		# A plus sign (+) before an entry number indicates that the entry contains the description of a gene of known sequence and a phenotype.
		#
		# A percent sign (%) before an entry number indicates that the entry describes a confirmed mendelian phenotype or phenotypic locus for which the underlying molecular basis is not known.
		#
		# No symbol before an entry number generally indicates a description of a phenotype for which the mendelian basis, although suspected, has not been clearly established or that the separateness of this phenotype from that in another entry is unclear.
		#
		# A caret (^) before an entry number means the entry no longer exists because it was removed from the database or moved to another entry as indicated.
		# Brackets, "[ ]", indicate "nondiseases," mainly genetic variations that lead to apparently abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
		#
		# Braces, "{ }", indicate mutations that contribute to susceptibility to multifactorial disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).		
		if ($mimtitle =~ /^[\*\^]/) {					# only keep entries that begin with #,%,+,no symbol; exclude if beginning with *, +, or ^
			$inrecord = 0;
			next;
		} else {
			$inrecord = 1;
		}
		my $isComplex = "no";
		
		if ($mimtitle =~ m/QUANTITATIVE TRAIT LOCUS/i || $mimtitle =~ /risk/i || $mimtitle =~ m/QTL/i || $mimtitle =~ /\[/ || $mimtitle =~ /\{/ || $mimtitle =~ m/suscep(\w+) to/i) {
			$isComplex = "yes";
		}
		if ($mimtitle =~ /somatic/i && $isComplex ne 'yes') {
			$isComplex = "somatic";
		}
		$printentry =  "$mimnum\t$mimtitle\thttp://omim.org/entry/$mimnum\t$isComplex";
	}
	
	if ($readline =~ /^MOLECULAR GENETICS$/i && $inrecord == 1) {
		$inMG = 1;
		next;
	}
	
	if ($inMG == 1) {
		if ($readline =~ /Associations Pending Confirmation/i) { next; }
		if ($readline =~ /^- /i) { next; }
	
		if ($readline =~ /^\*FIELD\*/) {
			# print STDERR "detecting end of MG\n";
			if ($MGtxt_paragraph !~ /^\s*$/) {			# if paragraph/line break exists with no content or paragraph only contains a link to another entry, do not count as a paragraph
				push(@MGtxt, $MGtxt_paragraph);
			}
			
			# print "paragraphs=".join("####", @MGtxt)."\n";			
			if ($MGtxt[0] =~ /^\s*For [ \w]*[discussion|review] of [ ,\w]+, see [ \w(),]+\.$/) {		# if paragraph only contains a link to another entry, do not count as a paragraph
				my $test = shift(@MGtxt);
				# print "deleting $test\n";
			}
			
			for (my $i=0; $i<=$#MGtxt; $i++) {
				if ($mentionsNGS eq 'no') {
					# print STDERR "paragraph $i\n";
					# print STDERR "$MGtxt[$i]\n";
					($mentionsNGS, $NGSyear) = containsNGS($MGtxt[$i]);
					if ($mentionsNGS eq 'yes') {
						$NGSparagraph = "'".($i+1).'/'.($#MGtxt+1);	
					}
				}
			}
			$yeardiscovered = 0;
			if (scalar(@MGtxt) >= 1) {
				my @yearmatches = $MGtxt[0] =~ /\([a-z, \.;]*(\d{4})(?:[,;] [\w, \.;]*)*\)/ig;
				if (scalar(@yearmatches) == 1) {			# if length 1, then only one year
					$yeardiscovered = $yearmatches[0];
				} else {
					foreach my $match (@yearmatches) {
						if ($match > $yeardiscovered) {
							$yeardiscovered = $match;
						}
					}
				}
			} 
			$inMG = 0;
            # print STDERR "2: $printentry\tmentionsNGS = $mentionsNGS\tparagraph = $NGSparagraph\n";
		} elsif ($readline =~ /^\n/ && $MGtxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
            # print STDERR "saving new paragraph: $MGtxt_paragraph!!!!!!!!!!!!!! \n";
			push(@MGtxt, $MGtxt_paragraph);
			$MGtxt_paragraph = "";
		} else {
			$readline =~ s/\v//;
			$MGtxt_paragraph .= " $readline";
            # print STDERR "add line: $MGtxt_paragraph!!!!!!!!!!!!!!\n";
		}
	}
    # print STDERR "inrecord=$inrecord\tinMG=$inMG\n";
}
close $input_handle;
close $output_handle;






sub containsNGS {
	my $text = $_[0];
	my $NGS = 'no';
	my $NGSyear = 0;
	
	if ($text =~ /exome[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /genome[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /massively[\s\-]*parallel[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /next[\s\-]*generation[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /high[\s\-]*throughput[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /exome[\s\-]*capture/i) {
		$NGS = 'yes';
	}
	if ($text =~ /whole[\s\-]*exome/i) {
		$NGS = 'yes';
	}
	
	if ($NGS eq 'yes') {
		while ($text =~ /\((\d{4})\)/g) {
			if ($1 > $NGSyear) {
				$NGSyear = $1;
			}
		}
		if ($NGSyear < 2009) {
			# must be a false positive detection of NGS, like the mention of high-throughput sequencing in entry 243800 for Johanson-Blizzard Syndrome
			$NGS = 'no';
			$NGSyear = 'NA';
		}
	} else {
		$NGSyear = 'NA';
	}
	
	return ($NGS, $NGSyear);
}



################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


xxx.pl - 


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

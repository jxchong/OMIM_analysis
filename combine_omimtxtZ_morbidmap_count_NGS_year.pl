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


my ($morbidmap, $mim2gene, $omimtxt, $outputfile, $help);

GetOptions(
	'morbidmap=s' => \$morbidmap, 
	'mim2gene=s' => \$mim2gene, 
	'omimtxt=s' => \$omimtxt, 
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $morbidmap) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} 


# perl parse_morbidmap.pl --in raw_download_2014-02-20/morbidmap --out morbidmap.parsed.txt
# perl parse_omimtxtZ_count_NGS.pl --in raw_download_2014-02-20/omim.txt.Z --out omimtxtZ.parsed.NGS.txt
# perl combine_omimtxtZ_morbidmap_count_NGS.pl --morbidmap morbidmap.parsed.txt --mim2gene raw_download_2014-02-20/mim2gene.txt --omimtxt omimtxtZ.parsed.NGS.txt --out combinedOMIM.mentionsNGS.txt
 

# store phenotype MIM and gene MIM info
my $MIM_ref = {};

read_mim2gene($MIM_ref, $mim2gene);
read_omimtxt($MIM_ref, $omimtxt);
read_morbidmap($MIM_ref, $morbidmap);

open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
print $output_handle "MIMnum\torigMIMnum\tmorbidmapname\tOMIMtxtname\tOMIMlink\tphenotypecount\tphenomappingkey\tlocussymbol\taltlocussymbols\tgeneMIMnum\tMIMtype\tisComplex\tmentionsNGS\tNGSparagraph\tNGSyear\tyearDiscovered\n";
while (my($phenoMIMnum, $MIMdata_ref) = each %{$MIM_ref}) {
	# for unknown reasons, there are two entries that don't actually exist in OMIM despite having a number assigned
	# see: phenotypeMIM = 135853, geneMIM = 138160
	# see: phenotypeMIM = NA, geneMIM = 164350
	if (defined $MIM_ref->{$phenoMIMnum}{'Type'} && ($MIM_ref->{$phenoMIMnum}{'morbidmapname'} ne 'NA' || $MIM_ref->{$phenoMIMnum}{'OMIMtxtname'} ne 'NA')) {
		my @MIMnum = split("-", $phenoMIMnum);
		my $mimURL = "http://omim.org/entry/$MIMnum[0]";
		
		my @locussymbols = split(", ", $MIM_ref->{$phenoMIMnum}{'locussymbol'});
		print $output_handle "$phenoMIMnum\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'origMIMnum'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'morbidmapname'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'OMIMtxtname'}."\t";
		print $output_handle "$mimURL\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'phenocount'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'mappingkey'}."\t";
		print $output_handle "$locussymbols[0]\t";
		print $output_handle join(",", @locussymbols[1..$#locussymbols])."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'geneMIMnum'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'Type'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'isComplex'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'mentionsNGS'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'NGSparagraph'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'NGSyear'}."\t";
		print $output_handle $MIM_ref->{$phenoMIMnum}{'yearDiscovered'};
		print $output_handle "\n";
	}
}
close $output_handle;








sub read_mim2gene {
	open (my $input_handle, "$mim2gene") or die "Cannot read $mim2gene: $!.\n";
	<$input_handle>;
	while ( <$input_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($MIMnum, $Type, $GeneIDs, $ApprovedGeneSymbols) = split("\t", $_);
		if ($Type ne 'moved/removed' && $Type ne 'gene') {
			$MIM_ref->{$MIMnum}{'Type'} = $Type;
			$MIM_ref->{$MIMnum}{'morbidmapname'} = 'NA';
			$MIM_ref->{$MIMnum}{'OMIMtxtname'} = 'NA';
			$MIM_ref->{$MIMnum}{'phenocount'} = 0;
			$MIM_ref->{$MIMnum}{'mappingkey'} = 'NA';
			$MIM_ref->{$MIMnum}{'geneMIMnum'} = 'NA';
			$MIM_ref->{$MIMnum}{'isComplex'} = 'NA';
			$MIM_ref->{$MIMnum}{'mentionsNGS'} = 'NA';
			$MIM_ref->{$MIMnum}{'NGSparagraph'} = 'NA';
			$MIM_ref->{$MIMnum}{'NGSyear'} = 'NA';
			$MIM_ref->{$MIMnum}{'yearDiscovered'} = 'NA';
		}
	}
	close $input_handle;
}

sub read_omimtxt {
	open (my $input_handle, "$omimtxt") or die "Cannot read $omimtxt: $!.\n";
	<$input_handle>;
	while ( <$input_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($phenoMIMnum, $mimtitle, $mimlink, $isComplex, $mentionsNGS, $NGSparagraph, $NGSyear, $yeardiscovered) = split("\t", $_);
		
		if (!defined $MIM_ref->{$phenoMIMnum}) {
			print STDERR "$phenoMIMnum in omim.txt.Z but not in mim2gene\n";
		}
		$MIM_ref->{$phenoMIMnum}{'origMIMnum'} = $phenoMIMnum;		
		$MIM_ref->{$phenoMIMnum}{'morbidmapname'} = 'NA';
		$MIM_ref->{$phenoMIMnum}{'OMIMtxtname'} = $mimtitle;
		$MIM_ref->{$phenoMIMnum}{'phenocount'} = 0;
		$MIM_ref->{$phenoMIMnum}{'mappingkey'} = -9;
		$MIM_ref->{$phenoMIMnum}{'locussymbol'} = 'NA';
		$MIM_ref->{$phenoMIMnum}{'geneMIMnum'} = 'NA';
		$MIM_ref->{$phenoMIMnum}{'isComplex'} = $isComplex;
		$MIM_ref->{$phenoMIMnum}{'mentionsNGS'} = $mentionsNGS;
		$MIM_ref->{$phenoMIMnum}{'NGSparagraph'} = $NGSparagraph;
		$MIM_ref->{$phenoMIMnum}{'NGSyear'} = $NGSyear;
		$MIM_ref->{$phenoMIMnum}{'yearDiscovered'} = $yeardiscovered;
		
		## DEBUG CODE
		# if ($phenoMIMnum eq '615674') {
		# 	print STDERR "omimtxt: $_\n";
		# 	print STDERR "MIMnum=$phenoMIMnum\n";
		# 	print STDERR "morbidmapname=$MIM_ref->{$phenoMIMnum}{'morbidmapname'}\n";
		# 	print STDERR "OMIMtxtname=".$MIM_ref->{$phenoMIMnum}{'OMIMtxtname'}."\n";
		# 	print STDERR "phenomappingkey=$MIM_ref->{$phenoMIMnum}{'mappingkey'}\n";
		# 	print STDERR "geneMIMnum=$MIM_ref->{$phenoMIMnum}{'geneMIMnum'}\n";
		# 	print STDERR "MIMtype=".$MIM_ref->{$phenoMIMnum}{'Type'}."\n";
		# 	print STDERR "isComplex=$isComplex\n\n";
		# }	
		##
	}
	close $input_handle;
}

sub read_morbidmap {
	open (my $input_handle, "$morbidmap") or die "Cannot read $morbidmap: $!.\n";
	<$input_handle>;
	while ( <$input_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($phenoname, $phenoMIMnum, $phenoMappingKey, $LocusSymbols, $GeneMIMnum, $CytoLoc, $isComplex) = split("\t", $_);
		my $finalMIMid = $phenoMIMnum;

		if ($phenoMIMnum eq '135853') {
			# I don't know what the deal is with this entry; it's in morbidmap but not in genemap at all
			next;
		}
		
		## DEBUG CODE
		# print STDERR "original: $_\n";
		# print STDERR "MIMnum=$phenoMIMnum; final=$finalMIMid\n";
		# print STDERR "morbidmapname=$phenoname\n";
		# # print STDERR "OMIMtxtname=".$MIM_ref->{$finalMIMid}{'OMIMtxtname'}."\n";
		# print STDERR "phenomappingkey=$phenoMappingKey\n";
		# print STDERR "geneMIMnum=$GeneMIMnum\n";
		# # print STDERR "MIMtype=".$MIM_ref->{$finalMIMid}{'Type'}."\n";
		# print STDERR "isComplex=$isComplex\n\n";
		##	
		if ($phenoMIMnum eq 'NA' && defined $GeneMIMnum) {
			# handle special cases where phenotype doesn't have a phenotype number because it shares a MIM number with the gene
			$finalMIMid = $GeneMIMnum;
			$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;		
			$MIM_ref->{$finalMIMid}{'morbidmapname'} = $phenoname;
			$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
			$MIM_ref->{$finalMIMid}{'mappingkey'} = $phenoMappingKey;
			$MIM_ref->{$finalMIMid}{'locussymbol'} = $LocusSymbols;
			$MIM_ref->{$finalMIMid}{'geneMIMnum'} = $GeneMIMnum;
			$MIM_ref->{$finalMIMid}{'isComplex'} = $isComplex;
		} elsif (!defined $MIM_ref->{$finalMIMid}{'morbidmapname'} || $MIM_ref->{$finalMIMid}{'morbidmapname'} eq 'NA') {
			# if morbidmapname is still blank, store info
			# print STDERR "test2\n";
			$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;	
			$MIM_ref->{$finalMIMid}{'morbidmapname'} = $phenoname;
			$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
			$MIM_ref->{$finalMIMid}{'mappingkey'} = $phenoMappingKey;
			$MIM_ref->{$finalMIMid}{'locussymbol'} = $LocusSymbols;
			$MIM_ref->{$finalMIMid}{'geneMIMnum'} = $GeneMIMnum;
			$MIM_ref->{$finalMIMid}{'isComplex'} = $isComplex;
		} elsif ($MIM_ref->{$phenoMIMnum}{'geneMIMnum'} ne $GeneMIMnum) {
			# sometimes there are multiple phenotypes with same MIM number but linked to/caused by mutations in different genes, e.g. 602080 for Paget disease of bone
			# we need to split these into separate phenotypes
			# print STDERR "test3\n";
			$finalMIMid = "$phenoMIMnum-$GeneMIMnum";
			$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;	
			$MIM_ref->{$finalMIMid}{'morbidmapname'} = $phenoname;
			$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
			$MIM_ref->{$finalMIMid}{'OMIMtxtname'} = $MIM_ref->{$phenoMIMnum}{'OMIMtxtname'};
			$MIM_ref->{$finalMIMid}{'mentionsNGS'} = $MIM_ref->{$phenoMIMnum}{'mentionsNGS'};
			$MIM_ref->{$finalMIMid}{'NGSparagraph'} = $MIM_ref->{$phenoMIMnum}{'NGSparagraph'};
			$MIM_ref->{$finalMIMid}{'NGSyear'} = $MIM_ref->{$phenoMIMnum}{'NGSyear'};
			$MIM_ref->{$finalMIMid}{'yearDiscovered'} = $MIM_ref->{$phenoMIMnum}{'yearDiscovered'};
			$MIM_ref->{$finalMIMid}{'locussymbol'} = $LocusSymbols;
			$MIM_ref->{$finalMIMid}{'mappingkey'} = $phenoMappingKey;
			$MIM_ref->{$finalMIMid}{'geneMIMnum'} = $GeneMIMnum;
			$MIM_ref->{$finalMIMid}{'isComplex'} = $isComplex;
			$MIM_ref->{$finalMIMid}{'Type'} = "phenotype-added";
		} elsif ($MIM_ref->{$finalMIMid}{'geneMIMnum'} eq $GeneMIMnum && $MIM_ref->{$finalMIMid}{'mappingkey'} ne $phenoMappingKey) {
			# sometimes there are multiple phenotypes with same MIM number but linked to/caused by mutations at the same locus, e.g. 613443
			# where the phenotypes have different mapping keys (cleanest example is mappingkey=4 for chromosomal syndrome and =3 for mutations in one of the genes in the deldup)
			# we need to split these into separate phenotypes
			# print STDERR "test4\n";
			$finalMIMid = "$phenoMIMnum-$GeneMIMnum-$phenoMappingKey";
			$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;	
			$MIM_ref->{$finalMIMid}{'morbidmapname'} = $phenoname;
			$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
			$MIM_ref->{$finalMIMid}{'OMIMtxtname'} = $MIM_ref->{$phenoMIMnum}{'OMIMtxtname'};
			$MIM_ref->{$finalMIMid}{'mentionsNGS'} = $MIM_ref->{$phenoMIMnum}{'mentionsNGS'};
			$MIM_ref->{$finalMIMid}{'NGSparagraph'} = $MIM_ref->{$phenoMIMnum}{'NGSparagraph'};
			$MIM_ref->{$finalMIMid}{'NGSyear'} = $MIM_ref->{$phenoMIMnum}{'NGSyear'};
			$MIM_ref->{$finalMIMid}{'yearDiscovered'} = $MIM_ref->{$phenoMIMnum}{'yearDiscovered'};
			$MIM_ref->{$finalMIMid}{'mappingkey'} = $phenoMappingKey;
			$MIM_ref->{$finalMIMid}{'locussymbol'} = $LocusSymbols;
			$MIM_ref->{$finalMIMid}{'geneMIMnum'} = $GeneMIMnum;
			$MIM_ref->{$finalMIMid}{'isComplex'} = $isComplex;
			$MIM_ref->{$finalMIMid}{'Type'} = "phenotype-added-bykey";
		} elsif ($MIM_ref->{$finalMIMid}{'geneMIMnum'} eq $GeneMIMnum && $MIM_ref->{$finalMIMid}{'mappingkey'} eq $phenoMappingKey) {
			# sometimes there are multiple phenotypes using the same phenotype and gene MIM numbers, e.g. 113750 for skin/hair/eye pigmentation AND for OCA type VI
			if ($MIM_ref->{$finalMIMid}{'isComplex'} eq 'yes' && $isComplex ne 'yes') {
				# print STDERR "test5\n";
				# new phenotype is Mendelian, stored is complex
				# we care about OCA type VI as it is Mendelian, but not pigmentation, since it is complex
				$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;	
				$MIM_ref->{$finalMIMid}{'morbidmapname'} = "$phenoname && $MIM_ref->{$finalMIMid}{'morbidmapname'}";
				$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
				$MIM_ref->{$finalMIMid}{'OMIMtxtname'} = $MIM_ref->{$phenoMIMnum}{'OMIMtxtname'};
				$MIM_ref->{$finalMIMid}{'mentionsNGS'} = $MIM_ref->{$phenoMIMnum}{'mentionsNGS'};
				$MIM_ref->{$finalMIMid}{'NGSparagraph'} = $MIM_ref->{$phenoMIMnum}{'NGSparagraph'};
				$MIM_ref->{$finalMIMid}{'NGSyear'} = $MIM_ref->{$phenoMIMnum}{'NGSyear'};
				$MIM_ref->{$finalMIMid}{'yearDiscovered'} = $MIM_ref->{$phenoMIMnum}{'yearDiscovered'};
				$MIM_ref->{$finalMIMid}{'mappingkey'} = $phenoMappingKey;
				$MIM_ref->{$finalMIMid}{'locussymbol'} = $LocusSymbols;
				$MIM_ref->{$finalMIMid}{'geneMIMnum'} = $GeneMIMnum;
				$MIM_ref->{$finalMIMid}{'isComplex'} = $isComplex;
			} elsif ($MIM_ref->{$finalMIMid}{'isComplex'} ne 'yes' && $isComplex ne 'yes') {
				# both phenotypes are Mendelian
				# print STDERR "test6\n";
				$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;	
				$MIM_ref->{$finalMIMid}{'morbidmapname'} = "$MIM_ref->{$finalMIMid}{'morbidmapname'} && $phenoname";
				$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
			} elsif ($MIM_ref->{$finalMIMid}{'isComplex'} ne 'yes' && $isComplex eq 'yes') {
				# new phenotype is complex, stored is Mendelian
				# print STDERR "test7\n";
				$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;	
				$MIM_ref->{$finalMIMid}{'morbidmapname'} = "$MIM_ref->{$finalMIMid}{'morbidmapname'} && $phenoname";
				$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
			} else {
				# if the current and stored phenotypes are both complex
				# print STDERR "test8\n";
				$MIM_ref->{$finalMIMid}{'origMIMnum'} = $phenoMIMnum;	
				$MIM_ref->{$finalMIMid}{'morbidmapname'} = "$MIM_ref->{$finalMIMid}{'morbidmapname'} && $phenoname";
				$MIM_ref->{$finalMIMid}{'phenocount'} += 1;
				$MIM_ref->{$finalMIMid}{'OMIMtxtname'} = $MIM_ref->{$phenoMIMnum}{'OMIMtxtname'};
				$MIM_ref->{$finalMIMid}{'mentionsNGS'} = $MIM_ref->{$phenoMIMnum}{'mentionsNGS'};
				$MIM_ref->{$finalMIMid}{'NGSparagraph'} = $MIM_ref->{$phenoMIMnum}{'NGSparagraph'};
				$MIM_ref->{$finalMIMid}{'NGSyear'} = $MIM_ref->{$phenoMIMnum}{'NGSyear'};
				$MIM_ref->{$finalMIMid}{'yearDiscovered'} = $MIM_ref->{$phenoMIMnum}{'yearDiscovered'};
				$MIM_ref->{$finalMIMid}{'mappingkey'} = $phenoMappingKey;
				$MIM_ref->{$finalMIMid}{'locussymbol'} = $LocusSymbols;
				$MIM_ref->{$finalMIMid}{'geneMIMnum'} = $GeneMIMnum;
				$MIM_ref->{$finalMIMid}{'isComplex'} = $isComplex;
			}
		} else {
			print STDERR "yikes!\n";
			print STDERR "MIMnum=original=$phenoMIMnum;final=$finalMIMid\n";
			print STDERR "morbidmapname=".$MIM_ref->{$finalMIMid}{'morbidmapname'}."\n";
			print STDERR "OMIMtxtname=".$MIM_ref->{$finalMIMid}{'OMIMtxtname'}."\n";
			print STDERR "phenocount=".$MIM_ref->{$finalMIMid}{'phenocount'}."\n";
			print STDERR "phenomappingkey=".$MIM_ref->{$finalMIMid}{'mappingkey'}."\n";
			print STDERR "locussymbol=".$MIM_ref->{$finalMIMid}{'locussymbol'}."\n";
			print STDERR "geneMIMnum=".$MIM_ref->{$finalMIMid}{'geneMIMnum'}."\n";
			print STDERR "MIMtype=".$MIM_ref->{$finalMIMid}{'Type'}."\n";
			print STDERR "isComplex=".$MIM_ref->{$finalMIMid}{'isComplex'}."\n\n";
		}
		
		## DEBUG CODE
		# print STDERR "final: $_\n";
		# print STDERR "MIMnum=original=$phenoMIMnum; final=$finalMIMid\n";
		# print STDERR "morbidmapname=".$MIM_ref->{$finalMIMid}{'morbidmapname'}."\n";
		# print STDERR "OMIMtxtname=".$MIM_ref->{$finalMIMid}{'OMIMtxtname'}."\n";
		# print STDERR "phenocount=".$MIM_ref->{$finalMIMid}{'phenocount'}."\n";
		# print STDERR "phenomappingkey=".$MIM_ref->{$finalMIMid}{'mappingkey'}."\n";
		# print STDERR "geneMIMnum=".$MIM_ref->{$finalMIMid}{'geneMIMnum'}."\n";
		# print STDERR "MIMtype=".$MIM_ref->{$finalMIMid}{'Type'}."\n";
		# print STDERR "isComplex=".$MIM_ref->{$finalMIMid}{'isComplex'}."\n\n";
		##
	}
	close $input_handle;
}






################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


combine_omimtxtZ_morbidmap.pl - 


=head1 SYNOPSIS


perl B<combine_omimtxtZ_morbidmap.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--morbidmap> F<parsed morbidmap>	

	produced by parse_morbidmap.pl

=item B<--mim2gene> F<mim2gene>	

	from OMIM download

=item B<--morbidmap> F<parsed morbidmap>	

	produced by parse_morbidmap.pl

=item B<--omimtxt> F<parsed OMIM.txt.Z>	

	produced by parse_omimtxtZ_count_NGS_year.pl

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

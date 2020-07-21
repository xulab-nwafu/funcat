#!/usr/bin/perl

use strict;
use warnings;

my $input_fasta = $ARGV[0];
my $refseq = "data/ATprotein";
my $funcat_file = "data/athaliana_funcat2008";
my $funcat_scheme = "data/funcat-2.1_scheme";

funcat_classification();

sub funcat_classification {

	# perform blast search
	system("makeblastdb -dbtype prot -in $refseq");
	system("blastp -query $input_fasta -db $refseq -num_threads 8 -outfmt 6 -evalue 1e-5 -culling_limit 1 -out blastp.txt");

	# find best hit
	my %hit = ();
	open (IN, "blastp.txt") or die $!;
	while (<IN>) {
		chomp;
		next if (/^\s*$/);
		my @w = split /\t/;
		if (!exists $hit{$w[0]}) {
			$hit{$w[0]} = uc($w[1]);
		}
	}
	close IN;

	# assign catalogue
	my $funcat = load_athaliana_funcat();
	my %cat = ();
	foreach my $protein (sort keys %hit) {
		$cat{$protein}{"id"} = $hit{$protein};
		$cat{$protein}{"num"} = $funcat->{$hit{$protein}};
	}

	# count catalogue
	my %num = ();
	foreach my $protein (sort keys %hit) {
		next if (!exists $funcat->{$hit{$protein}});
		my @num = split /\,/, $cat{$protein}{"num"};
		foreach my $cat (@num) {
			$num{$cat}++;
		}
	}

	# output catalogue
	open (OUT, ">funcat.stat.txt") or die $!;
	my $scheme = load_funcat_scheme();
	foreach my $cat (sort keys $scheme) {
		if (!exists $num{$cat}) {
			$num{$cat} = 0;
		}
		print OUT $cat, "\t", $scheme->{$cat}, "\t", $num{$cat}, "\n";
	}
	close OUT;

	open (OUT, ">funcat.list.txt") or die $!;
	foreach my $protein (sort keys %hit) {
		next if (!exists $funcat->{$hit{$protein}});
		print OUT $protein, "\t", $cat{$protein}{"id"}, "\t", $cat{$protein}{"num"}, "\n";
	}
	close OUT;

}

sub load_athaliana_funcat {
	my %at_funcat = ();
	open (IN, $funcat_file) or die $!;
	while (<IN>) {
		chomp;
		next if (/^\s*$/);
		next if (/^concat/);
		my @w = split /\|/;
		next if ($w[1] =~ /\./);
		if (!exists $at_funcat{$w[0]}) {
			$at_funcat{uc($w[0])} = $w[1];
		} else {
			$at_funcat{uc($w[0])} .= "," . $w[1];
		}
	}
	close IN;
	return \%at_funcat;
}

sub load_funcat_scheme {
	my %scheme = ();
	open (IN, $funcat_scheme) or die $!;
	while (<IN>) {
		chomp;
		next if (/^\s*$/);
		next if (/^\#/);
		my @w = split /\s/, $_, 2;
		next if ($w[0] =~ /\./);
		$w[1] = ucfirst(lc($w[1]));
		$w[1] =~ s/dna/DNA/;
		$scheme{$w[0]} = $w[1];
	}
	close IN;
	return \%scheme;
}

#!/usr/bin/env perl
# function_score.pl
# Calculation of partial AUC (false positive rate from 0% to 1%) of coexpression matrix
# (c) 2022 ATTED-II
#
# Options:
# -d: Directory of coexpression data provided in the bulk download page in ATTED-II.
#     Example of data source
#     https://zenodo.org/record/4961962/files/Ath-m.v21-01.G20819-S12686.combat_pca_subagging.ls.d.zip
# -k: Annotation file (tab-separated: Pathway ID, Entrez Gene IDs)
#     Example of data source (KEGG FTP license is required):
#     ftp://ftp.bioinformatics.jp/kegg/genes/organisms/ath/ath_link.tar.gz
# -g: (option) Paralog genes (tab-separated: Paralog ID, Entrez Gene IDs)
#     Example of data source (KEGG FTP license is required):
#     ftp://ftp.kegg.net/kegg/genes/links/genes_ncbi-proteinid.list.gz
#     ftp://ftp.kegg.net/kegg/genes/links/genes_ko.list.gz
# -M: use a smaller-is-better coexpression index. (Default: a larger-is-better index)

use strict;
use File::Basename;
use Scalar::Util 'looks_like_number';
use Getopt::Std;

my %opt;
getopts('d:k:g:|M', \%opt);

my $coex_dir = $opt{d};
my $pathway_file = $opt{k};
my $paralog_file = $opt{g} // '';
my $coex_type = ($opt{'M'}) ? 'smaller' : 'larger';
my $FPR = 0.01;
my $max_genes_in_pathway = 50;

# (1) Load Pathway Data
my %is_test_gene;
my %is_same_pathway;
{
    open IN, $pathway_file or die $pathway_file, $!;
    while (my $ln = <IN>) {
	chomp $ln;
	my ($pathwayID, @genes) = split /\t/, $ln;
	next if scalar @genes > $max_genes_in_pathway;
	next if scalar @genes < 2;
	while (my $g0 = shift @genes){
	    $is_test_gene{$g0} = 1;
	    for my $g1 (@genes){
		$is_same_pathway{$g0}{$g1} = 1;
		$is_same_pathway{$g1}{$g0} = 1;
	    }
	}
    }
    close IN;
}
printf STDERR "[info] pathways: %s (%d genes)\n", $pathway_file, scalar keys %is_test_gene;

# (2) Load Paralog Data
my %gene2ko;
if ($paralog_file){
    open IN, $paralog_file or die;
    while (my $ln = <IN>){
	chomp $ln;
	my ($KO_ID, @genes) = split /\t/, $ln;
	my %v;
	for my $g (@genes){
	    $v{$g} = 1 if $is_test_gene{$g};
	}
	next if scalar keys %v < 2;
	for my $g (keys %v){
	    $gene2ko{$g} = $KO_ID;
	}
    }
    close IN;
}

printf STDERR "[info] paralogs: %s (%d genes)\n", $paralog_file, scalar keys %gene2ko;

# (3) Load Coex Data
my $total_false = 0;
my $total_true = 0;
my %already;   # to omit mirror gene relation: (gene A to B) and (gene B to A)
my @files = <$coex_dir/*>;

my %coex;
my %positive;
my %is_paralog;
my $i = 1;

foreach my $infile (@files){
    next unless $is_test_gene{basename $infile};

    open IN, $infile or die;
    while (my $ln = <IN>){
	chomp $ln;
	my $g0 = basename $infile;
	my ($g1, $coex_value) = split /\t/, $ln;
	die "[error] No gene ID: ", $! unless $g0 || $g1 || $coex_value;
	next unless $is_test_gene{$g0};
	next unless $is_test_gene{$g1};
	next if $g0 eq $g1;
	next if $already{$g0}{$g1} || $already{$g1}{$g0};
	die "[error] not number: $ln" unless Scalar::Util::looks_like_number $coex_value;

	$coex{$i} = $coex_value;
	$positive{$i} = ($is_same_pathway{$g0}{$g1}) ? 1 : 0;

	if ($gene2ko{$g0} && $gene2ko{$g0} eq $gene2ko{$g1}){
	    $is_paralog{$i} = 1;
	} else {
	    if ($is_same_pathway{$g0}{$g1}){
		$total_true ++;
	    }else{
		$total_false ++;
	    }
	}
	$i ++;
	$already{$g1}{$g0} = 1;
    }
    close IN;
}

# (4) ordered by coexpression strength (Stronger to Weaker)
my @coex_ordered_keys = sort {$coex{$b} <=> $coex{$a}} keys %coex;
@coex_ordered_keys = reverse @coex_ordered_keys if $coex_type eq 'smaller';

# (5) pAUC Calculation and Output
# AUC: running total of trapezoid areas
my $pAUC = 0;
my $coex_threshold = 0;
{
    my $current_x = 0; # running total of false
    my $current_y = 0; # running total of true
    my $previous_x = 0;
    my $previous_y = 0;
    my $trapezoid_area = 0;
    my $sum_area = 0;
    my $threshold_FP = $FPR * $total_false;

    for my $i (@coex_ordered_keys){
	next if $is_paralog{$i};
	if ($positive{$i}){
	    $current_y ++;
	}else{
	    $current_x ++;
	}

	if ($current_x < $threshold_FP){
	    $trapezoid_area = ($current_x - $previous_x) * ($current_y + $previous_y) / 2;
	    $sum_area += $trapezoid_area;
	    $previous_x = $current_x;
	    $previous_y = $current_y;
	}else{
	    my $slice_ratio = ($threshold_FP - $previous_x) / ($current_x - $previous_x);
	    my $sliced_x = $previous_x + ($current_x - $previous_x) * $slice_ratio;
	    my $sliced_y = $previous_y + ($current_y - $previous_y) * $slice_ratio;
	    $trapezoid_area = ($sliced_x - $previous_x) * ($sliced_y + $previous_y) / 2;
	    $sum_area += $trapezoid_area;
	    $coex_threshold = $coex{$i};
	    last;
	}
    }
    $pAUC = $sum_area / ($total_true * $total_false);
}

# (6) OUTPUT
{
    printf "%.3f\t", $pAUC / ($FPR ** 2);
    printf "%s\t", $coex_dir;
    printf "%s\t", $coex_type;
    printf "%s\t", $coex_threshold;
    printf "%d\n", scalar keys %is_test_gene;
}

exit;

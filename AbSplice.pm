=head1 LICENSE
TODO the licence

=head1 CONTACT
Christian Mertes <mertes ... in.tum.de>
=cut

=head1 NAME
 AbSplice

=head1 SYNOPSIS
 mv AbSplice.pm ~/.vep/Plugins
 ./vep -i variants.vcf --plugin AbSplice,scores=/path/to/absplice_scores.tsv.gz
=head1 DESCRIPTION
 A VEP plugin that retrieves pre-calculated annotations from AbSplice.
 
 AbSplice predicts the splicing impact of a variant on the given gene model. 
 It scores the variant per transcript in the rnage of 0-1, where the authors 
 provide the following cutoffs for the interpretation of the effect:
 TODO set correct cutoffs
   0.2 (high recall)
   0.5 (recommended)
   0.8 (high precision)
 This plugin is available for both GRCh37 and GRCh38.

 More information can be found at:
 https://github.com/gagneurlab/absplice/

 Please cite the AbSplice publication alongside VEP if you use this resource:
 TODO https://www.ncbi.nlm.nih.gov/pubmed/ TODO
 
 Running options:
 Output: 
 The following steps are necessary before running this plugin:
 The plugin can then be run:
 ./vep -i variants.vcf --plugin AbSplice,scores=/path/to/absplice_scores.tsv.gz,cutoff=0.4
 TODO set score correctly
 TODO set exclude tissue correctly
=cut

package AbSplice;

use strict;
use warnings;
use List::Util qw(max);
use List::MoreUtils qw(first_index indexes);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use Bio::EnsEMBL::Variation::VariationFeature;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $output_vcf;

sub new {
  my $class = shift;
  
  # init default Tabix Plugin
  my $self = $class->SUPER::new(@_);
  my $param_hash = $self->params_to_hash();
  $self->expand_left(0);
  $self->expand_right(0);
  $self->add_file($param_hash->{scores});
  
  # set cutoff if provided 
  if(defined($param_hash->{cutoff})) {
    my $cutoff = $param_hash->{cutoff};
    if($cutoff < 0 || $cutoff > 1) {
      die("ERROR: Cutoff score must be between 0 and 1!\n");
    }
    $self->{cutoff} = $cutoff;
  }

  # set tissues to annotate (all from file if not provided)
  my $inFile = IO::Uncompress::Gunzip->new($param_hash->{scores})
      or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  my @all_tissues = split("\t", $inFile->getline());
  @all_tissues = @all_tissues[6 .. $#all_tissues];
  $self->{all_tissues} = \@all_tissues;
  my @tissues = @all_tissues;
  if(defined($param_hash->{tissues})){
    @tissues = split "\\|", $param_hash->{tissues};
    foreach my $tissue (@tissues){
      if(!grep(/^$tissue$/, @all_tissues)){
        die("ERROR: Provided tissue '$tissue' not in the provided score file!");
      }
    }
  }
  $self->{tissues} = \@tissues;

  # define idx for the max score calculation
  my @max_score_idx;
  # TODO define default for exclusion testis
  my @max_score_exclusion = ("Testis");
  if(defined($param_hash->{maxScoreExclude})){
    @max_score_exclusion = split("\\|", $param_hash->{maxScoreExclude});
  }
  for(my $i = 0; $i < scalar @tissues; $i++){
    push(@max_score_idx, $i) if !grep( /^$tissues[$i]$/, @max_score_exclusion);
  }
  $self->{max_score_idx} = \@max_score_idx;

  # check output format
  if($self->{config}->{output_format} eq "vcf") {
    $output_vcf = 1;
  }
  my %self = %{$self};

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  my %header;

  if($output_vcf) {
    $header{'AbSplice_pred_ENSG'} = 'The gene symbol used for the  AbSplice score prediction';
    $header{'AbSplice_pred_MAX'}    = 'The maximum AbSplice score across tissues.';
    foreach my $tissue (@{$self->{tissues}}){
      $header{"AbSplice_pred_$tissue"} = "The AbSplice score for tissue $tissue"
    }
  } else {
    my $tsvHeader = join "|", map {"AbSplice_pred_$_"} ("ENSG", "MAX", @{$self->{tissues}});
    $header{'AbSplice_pred'} = "AbSplice predicted effect on splicing. The values are separated by | and the full annotation is: $tsvHeader";
  }

  if($self->{cutoff}) {
    $header{'AbSplice_cutoff'} = "Flag if prediction score pass the provided cutoff of '$self->{cutoff}' (PASS) or if it does not (FAIL).";
  }
  return \%header;
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;
  my $chr = $vf->{chr};

  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;

  # extracted data from parse_data function() aka a line from the AbSplice score file
  my @data = @{$self->get_data($chr, $start, $end)} if(defined $chr);

  return {} unless(@data);

  my $result_data = '';
  my $result_flag;

  # Store all AbSplice results
  # Iterate through the fetched results and keep it for the output
  my %hash_aux;
  foreach my $data_values (@data) {

    my $ref_allele;
    my $alt_allele;

    # get alt allele
    my $allele = $tva->variation_feature_seq;
    reverse_comp(\$allele) if $vf->{strand} < 0;
    my $new_allele_string = $vf->ref_allele_string.'/'.$allele;

    if($vf->ref_allele_string =~ /-/) {

      # convert to vcf format to compare the alt alleles
      my $vf_2 = Bio::EnsEMBL::Variation::VariationFeature->new(
          -start => $start, -end => $end, -strand => $vf->{strand},
          -allele_string => $new_allele_string);
      my $convert_to_vcf = $vf_2->to_VCF_record;
      $ref_allele = ${$convert_to_vcf}[3];
      $alt_allele = ${$convert_to_vcf}[4];
    }
    else {
      $ref_allele = $vf->ref_allele_string;
      $alt_allele = $allele;
    }

    my $matches = get_matched_variant_alleles(
        {
            ref => $ref_allele, alts => [$alt_allele],
            pos => $start, strand => $vf->strand},
        {
            ref => $data_values->{ref}, alts => [$data_values->{alt}],
            pos => $data_values->{start}});

    # if we found a matching entry report the results
    if (@$matches) {
      my %hash;
      my %data_values = %{$data_values};
      
      if($output_vcf || $self->{config}->{output_format} eq "json" || $self->{config}->{rest})  {
        my $prefix ="";
        $prefix = "AbSplice_pred_" if $output_vcf;
        $hash{$prefix. 'ENSG'} = $data_values->{gene};
        $hash{$prefix. 'MAX'}    = $data_values->{max};
        foreach my $tissue (@{$self->{tissues}}){
          my $tissue_index = first_index { $_ eq $tissue } @{$self->{all_tissues}};
          $hash{$prefix.$tissue} = @{$data_values->{scores}}[$tissue_index];
        }
      }
      else {
        $hash{'AbSplice_pred'} = join("|", ($data_values->{gene}, $data_values->{max}, @{$data_values->{scores}}));
      }
 
      # Add a flag if cutoff is used
      if($self->{cutoff}) {
        $hash{'AbSplice_cutoff'} = "FAIL";
        if($data_values->{max} ne "" && $data_values->{max} >= $self->{cutoff}) {
          $hash{'AbSplice_cutoff'} = "PASS";
        }
      }
      $hash_aux{$data_values->{gene}} = \%hash;
    }
  }
  return {} unless(%hash_aux);

  # find the AbSplice gene matching the variant gene symbol, if there is a match
  my $gene_symbol_id = $tva->transcript->{_gene_stable_id};
  if(($gene_symbol_id) && ($hash_aux{$gene_symbol_id})) {
    if($self->{config}->{output_format} eq "json" || $self->{config}->{rest}){
      return {AbSplice => $hash_aux{$gene_symbol_id}};
    }
    return $hash_aux{$gene_symbol_id};
  }
  return {};
}

# Parse data from AbSplice scoring file
sub parse_data {
  my ($self, $line) = @_;
  my ($chr, $start, $end, $ref, $alt, $gene, @scores) = split(/\t/, $line, -1);
  die("ERROR: Could not extract the scores. Please submit an issue.") if scalar @scores < 1;
  warn("ERROR: Could not extract the exact number of the scores. The offending line is: '$line'.") 
      if scalar @scores != scalar @{$self->{all_tissues}};
  
  my @tmp_scores = @scores[@{$self->{max_score_idx}}];
  my @non_empty_index = indexes { $_ ne "" } @tmp_scores; 
  my $max_score = max(@tmp_scores[@non_empty_index]);
  $max_score = "" if(!defined($max_score));

  return {
    chr    => $chr,
    start  => $start,
    ref    => $ref,
    alt    => $alt,
    gene   => $gene,
    max    => $max_score,
    scores => \@scores,
  };
}

sub myDebug {
  my ($self, $info, $data) = @_;
  print STDERR "\nDebug Info: $info,\t";
  if(ref($data) eq "HASH"){
    my %data_hash = %{$data};
    print STDERR "HASH:\n";
    print STDERR map { "\t$_ => $data_hash{$_}\n" } keys %data_hash;
  } elsif(ref($data) eq "ARRAY"){
    print STDERR "ARRAY:\n\t".join(", ", @{$data})."\n";
  } else {
    print STDERR "STRING:\n\t$data\n";
  }
  print STDERR "\n";
}

1;


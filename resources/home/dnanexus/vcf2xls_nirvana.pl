#!/usr/bin/perl 
# Kim Brugger (18 Jul 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use File::Spec;
use List::Util qw( max );
use POSIX;
use Spreadsheet::WriteExcel;


# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $path = $0;
  if ($path =~ /.*\//) {
    $path =~ s/(.*)\/.*/$1/;
    push @INC, "$path";
    $ENV{'PATH'} .= ":$path/";
  }
}

my $TABIX = "packages/htslib-1.7/tabix";

use Vcf;
use Getopt::Std;

my $opts = 'ti:v:a:R:g:e:o:u:p:MA:T:w:fHFNDC:I';
my %opts;
getopts($opts, \%opts);

my $samtools  = 'packages/samtools-1.7/samtools';
my $FLANK     = 250;

my $RARE_VARIANT_AF = $opts{A} || 0.02;
my $nb_usable_reads = $opts{"u"};
my $total_nb_reads = $opts{"T"};
my $workflow = $opts{"w"};
my $workflow_id = $opts{"i"};

my $manifest               = "BioinformaticManifest";
my $genes2transcripts_file = "nirvana_genes2transcripts";
my $genepanels_file        = "genepanels";
my $gemini_freq            = "gemini_freq.vcf.gz";
my $esp_vcf             = "esp_vcf.tab.gz";
my $kg_vcf              = "kg_vcf.tab.gz";
my $exac_vcf            = "exac_vcf.sites.vep.vcf.gz";
my $exon_file           = "exons_nirvana";

my %genes2transcripts;
my %transcript2gene;
my %VEP_transcripts;
my %genepanels;
my %panel_names;
my %panel_id;

readin_panels_n_manifest();

my $ID_QC_ONLY = $opts{ 'I' } || 0;
my $IGNORE_QC = $opts{ 'N' } || 0;

my @QC_sheets = (
              'Summary',
              'QC',  
              'Gene QC',
              );

my @variant_sheets = (
              'stop_gained',
              'frameshift_variant',
              'consensus splice',
              'missense_variant',
              'synonymous_variant',
              'other',
    );

my @sheets;
if ( $ID_QC_ONLY ) {
    @sheets = (@QC_sheets);
} else {
    @sheets = (@QC_sheets, @variant_sheets);
}

my %effect_levels = ('stop_gained'        => 6,
		     'frameshift_variant' => 5,
		     'consensus splice'   => 4,
		     'missense_variant'   => 3,
		     'synonymous_variant' => 2,
		     'other'              => 1);


my %field_index = ();
my %QCfield_index = ();
my %geneQCfield_index = ();

my $PIPELINE_VERSION = '1.4';
my $MIN_QUALITY_VAR  = 200;
my $LOW_COVERAGE_VAR = 10;
my $meta_only = $opts{'M'} || 0;
my $text_only = $opts{t} || 0;

my $homozygous_snps = $opts{ 'H' } || 0;
my $full_exome_snps = $opts{ 'F' } || 0;

my $coverage_file = $opts{ 'C' };
my $runfolder_coverage_file = $opts{ 'R'};

# Sometimes the sample is given before the parameters on the command
# line, so check and fail if this is the case
# foreach my $argv ( @ARGV) {
#   die "option parameter after the sample argument, please check you command" if ( $argv =~ /^-/);
# }

my $vcf_file = $opts{a} || shift || usage();
my $raw_vcf_file = $opts{v};
usage() if ( $opts{ 'h' });

check_vcf_integrity_after_annotation($vcf_file, $raw_vcf_file);

my $sample = find_sample_name( $vcf_file );

$sample =~ s/_.*//;

my %gene_list = readin_manifest( $manifest, $sample);
my %hotspots;

if ( $opts{ 'p' } ) {
  %gene_list = parameter_panels2genes($opts{ 'p' }, $sample);
}

die "No genes for $sample\n" if ( keys %gene_list == 0 );
my %transcript_list = gene_list_to_transcript_list( \%gene_list );
print "Gene list: ", join(",", sort keys %gene_list) , "\n";

my $excel_file = "/home/dnanexus/out/xls_reports/report.xls";

$excel_file = File::Spec->rel2abs( $excel_file );
print "Output excel file == $excel_file\n";
my ($workbook, $formatting) = setup_workbook( $excel_file );

my %added_worksheets;
my %worksheet_offset;
my %worksheet_last_gene;

my %csv_workbook = ();
my @csv_order;


my @panels_w_ids = map { $_ =~ s/^\ +//;
			 my $v = $panel_id{ uc($_) } || "NA"; 
			 $_ = "$_ ( $v )" }  split(",", $transcript_list{ 'PANEL'});

$gene_list{ 'PANEL_IDS'} = join(", ", @panels_w_ids );

my $sry;
$sry = check_sry($coverage_file);

setup_worksheets();

my %meta_stats = ();
$meta_stats{ 'PANEL'} = $gene_list{ 'PANEL'};
$meta_stats{ 'PANEL_IDS'} = $gene_list{ 'PANEL_IDS'};

if ( ! $ID_QC_ONLY ) {
  analyse_vcf_file( $vcf_file );
}

my ( $total_length, $total_plus20x) = (0,0);
print "Filling QC sheets\n";
fill_QC_sheets() if ( $IGNORE_QC == 0);

$workbook->close() if ( ! $text_only || !$meta_only);

print "SUCCESS\n";


sub check_sry {
  my ($coverage_file) = @_;
  my $sry_outcome;
  my $SRY_region = "Y:2655024-2655649";

  open(my $in, "$TABIX $coverage_file $SRY_region |") || die "Could not open '$coverage_file': $!\n";

  while ( my $line = <$in> ) {
    chomp $line;
    my @F = split("\t", $line);
    
    my ($region_chrom, $region_start, $region_end, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19);

    if ( $F[ 4 ] && $F[ 4 ] =~ m/^\d+\z/) {
      ($region_chrom, $region_start, $region_end, $min, $mean, $max, undef, $missing, $depth_1to5, $depth_6to9, $depth_10to19) = @F;
    } else {
      ($region_chrom, $region_start, $region_end, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19) = @F;
    }
      
    if ("Y" ne $region_chrom &&
        "2655024" ne $region_start &&
        "2655649" ne $region_end  ) {
          die "SRY region not found";
    } else {
      if ($mean >= 25) {
        $sry_outcome = 1
      } else {
        $sry_outcome = 0
      }
      return $sry_outcome;
    }
  }
  die "Issue in getting the SRY region";
}


sub check_vcf_integrity_after_annotation {
  my ( $vcf_file, $raw_vcf_file ) = @_;

  my $gvcf_file = $vcf_file;
  $gvcf_file =~ s/.annotated/.g/;

  if ( -e $gvcf_file ) {
    print "vcf check includes gvcf\n";
    system("./vcf_integrity_check.py -e -g $gvcf_file $vcf_file $raw_vcf_file");
  }
  else {
    print "./vcf_integrity_check.py -e $vcf_file $raw_vcf_file\n";
    system("./vcf_integrity_check.py -e $vcf_file $raw_vcf_file");
  }
      
  die "vcf file ( $vcf_file) have lost its integrity" 
      if ($?  != 0);

  print STDERR "vcf has kept its integrity after annotation\n";
}

# Kim Brugger (22 Oct 2015)
sub find_sample_name {
  my ( $vcf_file ) = @_;

  my $sample = "$vcf_file";
  $sample = "$vcf_file";
  $sample =~ s/.*\///;
  $sample =~ s/\..*//;

  return $sample ;
}

# Kim Brugger (20 May 2015)
sub fill_QC_sheets {
  my @gene_transcripts;

  foreach my $transcript ( sort keys %transcript_list ) {
    next if ( $transcript eq 'PANEL' );
    next if ( $transcript eq 'PANEL_IDS' );

    my $gene_name = $transcript_list{$transcript};

    push @gene_transcripts, "$transcript_list{$transcript} ( $transcript )";
    gene_performance( $gene_name, $transcript );
  }

  @gene_transcripts = sort @gene_transcripts;
  
  my $avg_plus_20x = "0.0 %";
  my $panel_coverage = "0.0 %";

  if ( $total_length ) {
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ 'Name' }, 'Total:', );
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ 'Region length' },  $total_length);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ '20+x'  }, $total_plus20x);
    $avg_plus_20x = sprintf("%.2f %%", $total_plus20x/$total_length*100);
    $panel_coverage = POSIX::floor($total_plus20x/$total_length*100);
    worksheet_write('QC', $worksheet_offset{ 'QC' },  $QCfield_index{ '20+x %'  }, $avg_plus_20x);

    worksheet_write('Summary', 0, 5, $avg_plus_20x);

    $meta_stats{ 'AVG_coverage'} = $avg_plus_20x;
  }

  my $report_blurb = "Next Generation Sequencing (NGS) of the coding region (+/-5 bp) of the following genes (reference sequences) using the Illumina TruSight One sequencing panel (NB. Whole exon deletions/duplications and other large rearrangements are not detected with this method) : \n\n";

  $report_blurb .= join("; ", @gene_transcripts ) . "\n\n";
  $report_blurb .= "$panel_coverage % of this panel was sequenced to a depth of 20X or greater (this includes homologous regions where reads do not map uniquely), with analytical sensitivity of 99.5-99.9% (95% confidence interval from benchmarking against GIAB HG001 reference material).  ";

  # Dementia needs different (additional) text:
  if (index($gene_list{ 'PANEL_IDS'}, "Dementia") != -1) {
    $report_blurb .= "Targeted analysis of the exon/intron boundary of exon 9 of the MAPT gene (NM_005910.5), a known hot-spot for pathogenic variants, has been performed.";
  }
  
  $report_blurb .= " \n\nThe presence of variants reported above, except for variants of unknown significance, has been confirmed by Sanger sequencing. Variants with a population frequency greater than 1 in 500 for dominant conditions, and 1 in 50 for recessive disorders have been deemed insignificant and are not reported. Variants are named using HGVS nomenclature, where nucleotide 1 is the A of the ATG-translation initiation codon. Identification of variants present in NGS data was performed using the Dias pipeline.";

  if ( $gene_list{ 'PANEL_IDS'} ) {
    worksheet_write('Summary', 1 ,  6 , "Panel(s) w/ id's", $$formatting{ 'bold' });
    worksheet_write('Summary', 1 ,  7 , $gene_list{ 'PANEL_IDS'});
  }
  worksheet_write('Summary', 6 ,  0 , $report_blurb);
}


# Kim Brugger (24 Jun 2015)
sub find_runfolder {
  my ( $in_file ) = @_;

  my $abs_path = File::Spec->rel2abs( $in_file );
  $abs_path =~ s/(.*\/).*/$1/;
  $abs_path =~ s/vcfs//;
  $abs_path =~ s/stats//;
  $abs_path =~ s/bams//;
  my @F = split("/", $abs_path);

  return $F[-1];
}


# Kim Brugger (13 May 2015)
sub print_meta_stats {
  $meta_stats{ 'stop_gained' }{ 'rare' } ||= 0;
  $meta_stats{ 'frameshift_variant' }{ 'rare' } ||= 0;
  $meta_stats{ 'missense_variant' }{ 'rare' } ||= 0;
  $meta_stats{ 'consensus splice' }{ 'rare' } ||= 0;
  $meta_stats{ 'synonymous' }{ 'rare' } ||= 0;

  print join(
    "\t", "Sample", 'Panel',   
    'AVG_coverage',
    'stop_gained', 
    'frameshift_variant',
    'consensus splice', 
    'missense_variant', 
    'synonymous' , 
    $excel_file || "N/A"
  )."\n";

  print join(
    "\t", $sample, $meta_stats{ 'PANEL'},   
    $meta_stats{ 'AVG_coverage'},
    $meta_stats{ 'stop_gained' }{ 'rare' }, 
    $meta_stats{  'frameshift_variant' }{ 'rare' },  
    $meta_stats{ 'consensus splice' }{ 'rare' }, 
    $meta_stats{ 'missense_variant' }{ 'rare' }, 
    $meta_stats{ 'synonymous' }{ 'rare' }, 
    $excel_file || "N/A"
  )."\n";
}


sub alts2csqallele {
  # Taken from Matts script and ported to perl.
  # """Convert vcf alt field values into equivalent VEP csq allele values
  
  # Args:
  #     vcf_record (_Record): a pyvcf _Record 
  
  # Returns:
  #     TYPE: An ordered list containing one allele for each alt
  # """

  my ($entry ) = @_;
  my @csq_alleles;
  
  for my $alt ( @{$$entry{ ALT }} ) {
    if (length($$entry{REF}) == 1) {
      if (length( $alt ) == 1) {
	      push @csq_alleles, $alt;
      }
      else {
	      # 1:many
	      push @csq_alleles, substr($alt, 1);
      }
    }
    else {
      if (length($alt) == 1) {
	      # many:1
	      push @csq_alleles, "-";
      }
      else {
	      # many:many
	      push @csq_alleles, substr($alt, 1);
      }
    }
  }
  return @csq_alleles;
}


# Kim Brugger (20 May 2015)
sub analyse_vcf_file {
  my ( $file ) = @_;

  my $vcf = Vcf->new(file=>$file);
  $vcf->parse_header();

  while (my $entry = $vcf->next_data_hash()) {
    my $CSQ_line = $$entry{INFO}{'CSQ'};
    next if ( ! $CSQ_line );
    my @CSQs;
    map {push @CSQs, [split(/\|/, $_)] } split(",", $CSQ_line);

    my $interesting_variant = 0;

    my $pos = "$$entry{CHROM}:$$entry{POS}";

    my @usable_CSQs;

    # Pull out the genotypes for the sample
    my ($gt1, $gt2) = split("/",$$entry{ gtypes }{ $sample }{ GT });

    # translate the csq-genotype to a base rather than a number
    if ( $gt1 == 0 ) {
      $gt1 = $$entry{REF};
    }
    else{
      $gt1 = $$entry{ALT}[$gt1 - 1];
    }

    if ( $gt2 == 0 ) {
      $gt2 = $$entry{REF};
    }
    else{
      $gt2 = $$entry{ALT}[$gt2 - 1];
    }
      
    if ( $homozygous_snps && $$entry{INFO}{AC} eq "2" ) {
      $interesting_variant++;
    }      

    $interesting_variant++ if ($full_exome_snps );

    foreach my $CSQ ( @CSQs) {
      my ($Allele,$ENS_gene, $HGNC,$RefSeq,$feature,$effects,$CDS_position,$Protein_position,$Amino_acid,$Existing_variation,$SIFT,$PolyPhen,$HGVSc,$Distance) = @$CSQ;

      next if (!$effects || $effects eq "upstream_gene_variant" || $effects eq "downstream_gene_variant");

      # Neither of the genotypes matches with our sample
      next if $Allele ne $gt1 && $Allele ne $gt2;

      $HGNC = uc ( $HGNC );

      next if ( ! $HGNC );

      if ( $transcript2gene{ $RefSeq } && $HGNC ne $transcript2gene{ $RefSeq } ) {
        $transcript2gene{ $RefSeq } = $HGNC;
      }

      if ( $homozygous_snps ) {
        if ( $$entry{INFO}{AC} eq "2" ){
          next if ( ! $genes2transcripts{ $HGNC } );
          $gene_list{ $HGNC } = $genes2transcripts{ $HGNC };
          $gene_list{ 'PANEL'} = 'Homozygous SNPs';
        }
        else {
          next;
        }
      }
      
      elsif ( $full_exome_snps ) {
        $gene_list{ $HGNC } = $genes2transcripts{ $HGNC };
        $gene_list{ 'PANEL'} = 'Full Exome';
      }
      elsif ( ! $transcript_list{ $RefSeq } ) {
        next;
      }

      if ($transcript_list{ uc ($RefSeq ) } && $transcript_list{ uc ($RefSeq ) } ne "" ) {
        $gene_list{ uc ($HGNC ) } = $RefSeq;
        $genes2transcripts{ $HGNC } = $RefSeq;

      }

      $effects = 'frameshift_variant' if ( $effects =~ /inframe_insertion/ ||
					   $effects =~ /inframe_deletion/  ||
					   $effects =~ /feature_elongation/ );

      # splice_region_variant = -8 to -3 and +1 to +3 therefore splice_region_variant not always intronic! 
      # Set as splice here, and intronic vars will be converted to other by write_variant if -6 to -8 (i.e. not +/-5)
      $effects = 'consensus splice' if ( $effects =~ /intron_variant/ && $effects =~ /splice_region_variant/);
      $effects = 'missense variant'   if ( $effects eq 'stop_lost');
      $effects = 'stop_gained'        if ( $effects eq 'start_lost');

      push @usable_CSQs, [$effects, $CSQ]
    }
     
    my ( $effects, $CSQ ) = (undef, undef);
    
    if ( @usable_CSQs == 0 ) {
      next;
    }

    my $comment = "";

    if ( @usable_CSQs > 1 ) {
      $comment = "Multi non-ref allelic site";
    }

    for (my $i = 0;  $i <@usable_CSQs; $i++) {
      ( $effects, $CSQ ) = @{$usable_CSQs[ $i ]};

      foreach my $effect ( split("&", $effects)) {
        if ( grep(/$effect/, @sheets )) {
          write_variant($effect, $entry, $CSQ, $comment);
            last;
        }
        # Keep any variants which don't match sheet names in other sheet
        else {
          write_variant('other', $entry, $CSQ, $comment);    
        }
      }
    }
  }
}


# Kim Brugger (20 May 2015)
sub setup_workbook {
  my ($outfile) = @_;

  my $book;
  my %formats;
  
  return if ( $text_only || $meta_only);
  
  my $workbook = Spreadsheet::WriteExcel->new( $excel_file );

  # default/ building colours
  #  8   =>   black
  #  9   =>   white
  # 10   =>   red
  # 11   =>   lime
  # 12   =>   blue
  # 13   =>   yellow
  # 14   =>   magenta
  # 15   =>   cyan
  # 16   =>   brown
  # 17   =>   green
  # 18   =>   navy
  # 20   =>   purple
  # 22   =>   silver
  # 23   =>   gray
  # 33   =>   pink
  # 53   =>   orange

  $formats{ 'bold'}           = $workbook->add_format(bold => 1);
  $formats{ 'red_cell'}       = $workbook->add_format(color => 'red'    );
  $formats{ 'purple_cell'}    = $workbook->add_format(color => 'purple' );
  $formats{ 'green_cell'}     = $workbook->add_format(color => 'green'  );
  $formats{ 'lightblue_cell'} = $workbook->add_format(color => 'lightblue'  );
  $formats{ 'orange_cell'}    = $workbook->add_format(color => 'orange'  );

  $formats{ 'dark_red_cell'}  = $workbook->set_custom_color(40, 216, 12, 12   ); # DarkRed
  $formats{ 'light_blue'}  = $workbook->set_custom_color(41, 78,  186,  255   ); 
  $formats{ 'magenta_cell'}  = $workbook->add_format(color => 'magenta'  );

  $formats{'table_head'}   = $workbook->add_format(border => 1, bold => 1, bg_color => $formats{ 'light_blue'});

  $formats{'table'}        = $workbook->add_format(top => 1, left=>1, right=>1, bottom => 1);

  $formats{'merged_table_head'}   = $workbook->add_format(bold => 1, bg_color => $formats{ 'light_blue'}, border => 1);
  $formats{'merged_table'}        = $workbook->add_format(top => 1, left=>1, right=>1, bottom => 1);

  return ($workbook, \%formats );
}


# Kim Brugger (18 Jan 2018)
sub gemini_af {
  my ( $chrom, $pos, $ref, $alt) = @_;

  if ( -e $gemini_freq ) {
    open(my $in, "$TABIX $gemini_freq $chrom:$pos-$pos |") || die "Could not open '$gemini_freq': $!\n";
    while (<$in>) {
      chomp;
      my ( $vcf_chrom, $vcf_pos, undef, $vcf_ref, $vcf_alt, undef, undef, $info ) = split("\t");

      if ( $vcf_pos == $pos && $vcf_ref eq $ref && $vcf_alt =~ /$alt/ ) {
        $info =~ /AF=(.*?);/;
        return $1;
      }
    }
    return 0;
  }
  else {
    die "Freq file '$gemini_freq' not found \n";

  }
}


# Kim Brugger (18 Jan 2018)
sub external_af {
  my ( $chrom, $pos, $ref, $alt) = @_;

  my %res;

  if ( -e $kg_vcf ) {
    open(my $in, "$TABIX $kg_vcf $chrom:$pos-$pos |") || die "Could not open '$kg_vcf': $!\n";

    while (<$in>) {
      chomp;
      my ( $vcf_chrom, $vcf_pos, $id, $vcf_ref, $vcf_alt, $AF_AFR, $AF_AMR, $AF_ASN, $AF_EUR, $AF_MAX ) = split("\t");

      if ( $vcf_pos == $pos && $vcf_ref eq $ref && $vcf_alt =~ /$alt/ ) {
        $res{ '1KG'} = $AF_MAX;
        last;
      }
    }
  }

  if ( -e $esp_vcf ) {
    open(my $in, "$TABIX $esp_vcf $chrom:$pos-$pos |") || die "Could not open '$esp_vcf': $!\n";
    
    while (<$in>) {
      chomp;
      my ( $vcf_chrom, $vcf_pos, $id, $vcf_ref, $vcf_alt, $AF_AA, $AF_AE, $AF_TA ) = split("\t");

      if ( $vcf_pos == $pos && $vcf_ref eq $ref && $vcf_alt =~ /$alt/ ) {
        $res{ 'ESP'} = max( $AF_AA, $AF_AE );
        last;
      }
    }
  }

  if ( -e $exac_vcf ) {
    open(my $in, "$TABIX $exac_vcf $chrom:$pos-$pos |") || die "Could not open '$exac_vcf': $!\n";

    while (<$in>) {
      chomp;
      my ( $vcf_chrom, $vcf_pos, undef, $vcf_ref, $vcf_alt, undef, undef, $info ) = split("\t");

      if ( $vcf_pos == $pos && $vcf_ref eq $ref && $vcf_alt =~ /$alt/ ) {
        my %fields;

        foreach my $field ( split(";",$info)) {
          $field =~ /(.*?)=(.*)/;
          $fields{ $1 } = $2;
        }
	
        if ( $fields{ AF } =~ s/,.*// ){
          $res{ ExAC } = $fields{ AF };
        }
        elsif ( $fields{ AF } ) {
          $res{ ExAC } = $fields{ AF };
        }
	
	      last;
      }
    }
  }
  return \%res;
}


# Kim Brugger (23 Aug 2013)
sub write_variant {
  my ($sheet_name, $entry, $CSQ, $comment) = @_;

  my ($Allele,$ENS_gene, $HGNC,$RefSeq,$feature,$Consequence,$CDS_position,$Protein_position,$Amino_acid,$Existing_variation,$SIFT,$PolyPhen,$HGVSc,$HGVSp) = @$CSQ;

  return if ($Consequence eq "upstream_gene_variant" || $Consequence eq "downstream_gene_variant");
  return if ( $HGVSc && $HGVSc eq "");

  if ( $Amino_acid ) {
    my ($aa1, $aa2) = split("/", $Amino_acid);
    $aa2 ||= $aa1;
    $sheet_name = 'synonymous_variant' if ( $aa1 eq $aa2 );
  }

  my $pure_cpos = pure_cpos($HGVSc);

  if ( $pure_cpos =~ /\+(\d+)/ || $pure_cpos =~ /\-(\d+)/) {
    $sheet_name = 'other' if ( $1 > 5 );
    $sheet_name = 'consensus splice' if ( $1 <= 5 );         
  }

  $sheet_name  = 'other' if ( $pure_cpos =~ /c\.\*/ || $pure_cpos =~ /c\.\-/ );

  my $change = $HGVSc || "";
  $change =~ s/.*://;
  $change =~ s/.*\d+(.*\>.*)/$1/;

  my @depths = split(",", $$entry{gtypes}{$sample}{AD});

  my ($gt1, $gt2) = split("/",$$entry{ gtypes }{ $sample }{ GT });
  my $genotype = "HET";
  $genotype = "HOMO" if ($gt1 == $gt2);

  my $depth = $$entry{ gtypes }{ $sample }{ DP };
  my $AAF = -1;

  #calculate the AAF for the variant, and transclate csq-gt into the right gt
  if ( $gt1 > 0 && $$entry{ALT}[ $gt1 - 1] eq $Allele ) {
    if ( $gt1 && $depth) {
      $AAF = sprintf("%.4f", $depths[ $gt1 ]/($depth));
    }
    else {
      $AAF = 1;
    }

    $Allele = $$entry{ ALT }[ $gt1 - 1];
  }
  
  elsif ( $gt2 > 0 && $$entry{ALT}[ $gt2 - 1] eq $Allele ) {
    if ($depth == 0) {
      $AAF = 0;
    }
    else {
      $AAF = sprintf("%.4f", $depths[ $gt2 ]/($depth));    
      $Allele = $$entry{ ALT }[ $gt2 - 1];
    }
  }


  my $AF_GEMINI = gemini_af( $$entry{CHROM}, $$entry{POS}, $$entry{REF}, $Allele);
  my $external_AFs = external_af($$entry{CHROM}, $$entry{POS}, $$entry{REF}, $Allele);

  foreach my $external_AF_source ( keys %$external_AFs  ) {
    if ( $external_AF_source =~ /ExAC/i ) {
      $$entry{ 'INFO'}{'AF_ExAC' } = $$external_AFs{ $external_AF_source };
    }
    elsif ( $external_AF_source =~ /ESP/i ) {
      $$entry{ 'INFO'}{'AF_ESP_MAX' } = $$external_AFs{ $external_AF_source };
    }
    elsif ( $external_AF_source =~ /1KG/i ) {
      $$entry{ 'INFO' }{ 'AF_1KG_MAX' } = $$external_AFs{ $external_AF_source };
    }
  }
  
  $HGNC   = $ENS_gene if ( !$HGNC   || $HGNC   eq "" );
  $RefSeq = $feature  if ( !$RefSeq || $RefSeq eq "" );

  $RefSeq = $genes2transcripts{ $HGNC } if ($genes2transcripts{ $HGNC });

  if ($sheet_name eq "hotspot" ) {
    $HGNC      = $$entry{ID};
    $RefSeq    = "";
    $pure_cpos = "";
    $SIFT      = "";
    $PolyPhen  = "";
    $change = "$$entry{REF}>". join(",", @{$$entry{ALT}});
  }

  if ($worksheet_last_gene{$sheet_name } && $worksheet_last_gene{$sheet_name } ne $HGNC ) {
    $worksheet_offset{ $sheet_name }++;
  }

  my $format = undef;
  my $sequence_neighbourhood = "";
  ###################################################################################################################
  if ( ($$entry{CHROM} eq "1" && $$entry{POS} eq "158597507" && $$entry{REF} eq "G" && ($$entry{ALT}[0] eq "C" || $$entry{ALT}[1] eq "C")) || ($$entry{CHROM} eq "1" && $$entry{POS} eq "158587858" && $$entry{REF} eq "G" && ( $$entry{ALT}[0] eq "A" || $$entry{ALT}[1] eq "A")) ) {
    $format = $$formatting{'red_cell'};
    $comment .= " Low expression allele LELY: confirmation needed";
  }

  if ( ($$entry{CHROM} eq "11" && $$entry{POS} eq "88911696" && $$entry{REF} eq "C" && ( $$entry{ALT}[0] eq "A" || $$entry{ALT}[1] eq "A"))) {
    $format = $$formatting{'red_cell'};
    $comment .= " Disease associated polymorphism causes mild disease allele when in cis with p.Arg402Gln";
  }

  if ( ($$entry{CHROM} eq "11" && $$entry{POS} eq "89017961" && $$entry{REF} eq "G" && ( $$entry{ALT}[0] eq "A" || $$entry{ALT}[1] eq "A"))) {
    $format = $$formatting{'red_cell'};
    $comment .= " Disease associated polymorphism causes mild disease allele when in cis with p.Ser192Tyr";
  }

  if ( ($$entry{CHROM} eq "17" && $$entry{POS} eq "78078341" && $$entry{REF} eq "T" && ( $$entry{ALT}[0] eq "G" || $$entry{ALT}[1] eq "G"))) {
    $format = $$formatting{'red_cell'};
    $comment .= " Pathogenic variant";
  }

  if ( ($$entry{CHROM} eq "16" && $$entry{POS} eq "89613145" && $$entry{REF} eq "C" && ( $$entry{ALT}[0] eq "T" || $$entry{ALT}[1] eq "T"))) {
    $format = $$formatting{'red_cell'};
    $comment .= " Pathogenic variant";
  }
  ###################################################################################################################
  if ( low_AF_variant($AF_GEMINI, $$entry{INFO}{ 'AF_1KG_MAX' }, $$entry{INFO}{ 'AF_ESP_MAX' }, $$entry{INFO}{ 'AF_ExAC' })) {
    $format = $$formatting{'red_cell'};

    if ( $$entry{'QUAL'} < $MIN_QUALITY_VAR || $$entry{INFO}{DP} < $LOW_COVERAGE_VAR ) {
      $format = $$formatting{'purple_cell'};
      $meta_stats{ $sheet_name }{ 'rare/LQ' }++;
    }
    else {
      $meta_stats{ $sheet_name }{ 'rare' }++;
    }
  }
  else {
    $meta_stats{ $sheet_name }{ 'common' }++;
  }

  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Change' }, "$change", $format );
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Score' }, $$entry{QUAL}, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Depth' }, $$entry{INFO}{DP}, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Genotype' }, $genotype, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Gene' }, $HGNC, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Transcript' }, $RefSeq, $format);

  my $position = "$$entry{'CHROM'}:$$entry{'POS'}";
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Position' }, "$position", $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Genomic Ref Allele' }, $$entry{'REF'}, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Genomic Alt Allele' }, $Allele, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'AAF' }, $AAF, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Nucleotide pos' }, pure_cpos($HGVSc), $format);

  if ( $Amino_acid && $sheet_name ne "hotspot" ) {
    my ($aa1, $aa2) = split("/", $Amino_acid);

    $aa2 ||= $aa1;

    $aa1 = one2three( $aa1 );
    $aa2 = one2three( $aa2 );

    worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'AA change' }, "p.$aa1$Protein_position $aa2", $format);
  }

  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'dbsnp' }, $$entry{ID}, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'PolyPhen' }, $PolyPhen, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'SIFT' }, $SIFT, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'AF_GEMINI' }, $AF_GEMINI, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'AF_1KG_MAX' }, $$entry{INFO}{ 'AF_1KG_MAX' }||= 0, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'AF_ESP_MAX' }, $$entry{INFO}{ 'AF_ESP_MAX' }||= 0, $format);
  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'AF_ExAC' }, $$entry{INFO}{ 'AF_ExAC'    }||= 0, $format);

  $comment ||= "";

  worksheet_write($sheet_name, $worksheet_offset{ $sheet_name }, $field_index{ 'Comment' }, "$comment", $format);

  $worksheet_offset{ $sheet_name }++;

  $worksheet_last_gene{$sheet_name } = $HGNC;
}


# Kim Brugger (20 May 2015)
sub low_AF_variant {
  my ( @freqs ) = @_;
  map { return 0 if ($_ && $_ > $RARE_VARIANT_AF) } @freqs;
  return 1;
}


# Kim Brugger (29 Jun 2012)
sub pure_cpos {
  my $pos = shift;

  return "" if (! $pos );

  $pos =~ s/ENST\d+.\d+://;
  $pos =~ s/NM_\d+.\d+://;

  $pos =~ s/[ACGT]+>[ACGT]//;
  return $pos;

  $pos =~ s/del.*\z//;
  $pos =~ s/ins.*\z//;
  $pos =~ s/\[\d+\][ACGT]//;

  $pos =~ s/\s//g;

  $pos =~ s/c\.//;
  $pos =~ s/^\*//;
  $pos =~ s/\+\d+\z//;
  $pos =~ s/\-\d+\z//;
  $pos =~ s/\s+//g;

  return 0 if ( $pos eq "");
  return $pos;
}


# Kim Brugger (17 Jan 2018)
sub region_coverage {
  my ($chrom, $start, $end) = @_;

  $end  ||= $start;

  my $region = "$chrom:$start-$end";

  if ($coverage_file) {
    open(my $in, "$TABIX $coverage_file $region |") || die "Could not open '$coverage_file': $!\n";

    while ( my $line = <$in> ) {
      chomp $line;
      my @F = split("\t", $line);
      
      my ($region_chrom, $region_start, $region_end, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19);
    
      if ( $F[ 4 ] && $F[ 4 ] =~ m/^\d+\z/) {
	      ($region_chrom, $region_start, $region_end, $min, $mean, $max, undef, $missing, $depth_1to5, $depth_6to9, $depth_10to19) = @F;
      }
      else {
	      ($region_chrom, $region_start, $region_end, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19) = @F;
      }
      
      if ( $chrom ne $region_chrom ||
           $start ne $region_start ||
           $end   ne $region_end  ) {
        next;
      }

      my %res;
      $res{ 'min' } = $min; 
      $res{ 'mean' } = $mean;  
      $res{ 'max' } = $max; 
      $res{ 'missing' } = $missing if $missing;
      $res{ '1to5' } = $depth_1to5 if $depth_1to5;
      $res{ '6to9' } = $depth_6to9 if $depth_6to9;
      $res{ '10to19' } = $depth_10to19 if $depth_10to19;
      
      return \%res;
    }

    print "coverage for $region not found in $coverage_file \n";

    return undef;
  }
  else {
    die "No coverage file - using coverage database..";
  }
}


# Kim Brugger (03 Apr 2018)
sub fetch_expected_coverage {
  my ($chrom, $start, $end) = @_;
  my $normalized_depth_factor = $nb_usable_reads/100000000;

  $end  ||= $start;

  my $region = "$chrom:$start-$end";

  open(my $in, "$TABIX $runfolder_coverage_file $region |") || die "Could not open '$runfolder_coverage_file': $!\n";
  my $line = <$in>;

  chomp $line;
  my @F = split("\t", $line);
  
  my ($fchrom, $fstart, $fend, $mean, $sd) = @F;

  my %res = ( mean => $mean*$normalized_depth_factor, 
	      sd => $sd);
  return \%res;
}


# Kim Brugger (09 Dec 2013)
sub gene_performance {
  my  ( $gene_name, $refseq ) = @_;

  my $refseq_basename = $refseq;
  $refseq_basename =~ s/\.\d+/./;

  my $python_stdout = qx(python get_transcript_regions.py -c $coverage_file -t $refseq -e $exon_file);
  my @exons = eval $python_stdout or warn "$@";
  
  my %exons;
  my %worst_exon;
  if ( ! @exons ) {
    die "Unknown or non-captures gene: $gene_name/$refseq for sample: $sample\n";
  }

  my $exons_qced = 0;

  foreach my $exon ( @exons ) {
    if ( $exon == -1 ) {
      print "Unknown gene: $gene_name/$refseq\n";
      next;
    }

    next if ( $refseq && $$exon{ refseq } !~ /$refseq_basename/);

    $exons_qced++;
    
    my $RefSeq = $$exon{ refseq } || "";
    $RefSeq    = $refseq || "";
    my $region = $$exon{ region };

    if (!$RefSeq ) {
      $RefSeq = $gene_name;
    }

    my $coverage = region_coverage($$region{ 'chrom' },$$region{ 'start' },$$region{ 'end' });

    if ( ! $coverage ) {
      $workbook->close() if ( ! $text_only && ! $meta_only);
      system "rm -f $excel_file";
      die "$gene_name exon:$$exon{ exon_nr } [$$region{ 'chrom' }:$$region{ 'start' }-$$region{ 'end' }] coverage information for $sample \"aid\" is not in the database!!! Bailing ...\n";
    }

    my $exon_name = $gene_name . "_Exon$$exon{ exon_nr }";

    if ($$exon{ exon_nr } < 0) {
      $exon_name = $gene_name . "_Intron". abs( $$exon{ exon_nr } );
    }

    if ( (! defined $worst_exon{ "$RefSeq" }{'min'}) || 
         (  $$coverage{min} == 0 ) ||
         ($worst_exon{ "$RefSeq" }{'min'} > $$coverage{min})) {

      $worst_exon{ "$RefSeq" }{'min'} = $$coverage{min};
      $worst_exon{ "$RefSeq" }{'position'} = "$$region{ 'chrom' }:$$region{ 'start' }-$$region{ 'end' }";
    }

    $worst_exon{ "$RefSeq" }{'missing'} += range2length($$coverage{'missing'});
    $worst_exon{ "$RefSeq" }{'1to5'}    += range2length($$coverage{'1to5'});
    $worst_exon{ "$RefSeq" }{'6to9'}    += range2length($$coverage{'6to9'});
    $worst_exon{ "$RefSeq" }{'10to19'}  += range2length($$coverage{'10to19'});
    $worst_exon{ "$RefSeq" }{'transcript_length'}  += $$region{ 'end' }-$$region{ 'start' } + 1; # 0 vs 1 based coordinates
    
    $exons{ "$RefSeq" }{ $exon_name }{ 'coverage' } = $coverage;
    $exons{ "$RefSeq" }{ $exon_name }{ "$RefSeq" }{ 'region' }   = $region;
    $exons{ "$RefSeq" }{ $exon_name }{ 'position' }   = "$$region{ 'chrom' }:$$region{ 'start' }-$$region{ 'end' }";
  }

  if ( !$exons_qced ) {
    print STDERR "Likely wrong transcript ($refseq) for $gene_name, sample: $sample\n";
    system "rm -f $excel_file";

    exit -1;
  }

  foreach my $transcript ( sort keys %worst_exon ) {
    my $plus20x = $worst_exon{$transcript}{'transcript_length'} -  $worst_exon{$transcript}{ 'missing'} - $worst_exon{$transcript}{ '1to5'} - $worst_exon{$transcript}{ '6to9'} - $worst_exon{$transcript}{ '10to19'};

    $total_length  += $worst_exon{$transcript}{'transcript_length'};
    $total_plus20x += $plus20x;

    my $format = undef;

    if ( $worst_exon{ $transcript }{'min'} && $worst_exon{ $transcript }{'min'} < 20  || $worst_exon{ $transcript }{'min'} == 0 ) {
      $format = $$formatting{ 'red_cell' };
    }

    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ 'Name' }, $gene_name, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ 'Transcript' }, $transcript, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ 'Min depth' },  $worst_exon{$transcript}{'min'}, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ 'Region length' },  $worst_exon{$transcript}{'transcript_length'}, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ 'Missing'  }, $worst_exon{$transcript}{ 'missing'}, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ '1-5x'  }, $worst_exon{$transcript}{ '1to5'}, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ '6-9x'  }, $worst_exon{$transcript}{ '6to9'}, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ '10-19x'  }, $worst_exon{$transcript}{ '10to19'}, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' }, $QCfield_index{ '20+x'  }, $plus20x, $format);
    worksheet_write('QC', $worksheet_offset{ 'QC' },  $QCfield_index{ '20+x %'  }, sprintf("%.2f %%", $plus20x/$worst_exon{$transcript}{'transcript_length'}*100), $format);

    $worksheet_offset{ 'QC' }++;
  }

  $worksheet_offset{ 'QC' }++;

  foreach my $transcript ( sort keys %exons ) {
    foreach my $exon ( sort {my $A = $a;
			     my $B = $b;
			     $A =~ s/.*exon//i;
			     $B =~ s/.*exon//i;
			     $A =~ s/.*intron//i;
			     $B =~ s/.*intron//i;
			     $A <=> $B } keys %{$exons{ $transcript }} ) {
      
      $exons{$transcript}{ $exon }{'coverage'}{'min'} ||= 0;
      my $format = undef;

      if ( $exons{$transcript}{ $exon }{'coverage'}{'min'} < 20 ) {
	      $format = $$formatting{ 'red_cell' };
      }
      
      my $region = $exons{$transcript}{$exon}{$transcript}{"region"};
      my $expected_coverage_hash = fetch_expected_coverage( $$region{ 'chrom' }, $$region{ 'start' }, $$region{ 'end' } );
	
      my $expected_formatted = 'NA';

      if ( $expected_coverage_hash ){
        my $expected_coverage =$$expected_coverage_hash{'mean'};
        my $expected_coverage_sd =$$expected_coverage_hash{'sd'};
	      $expected_formatted = sprintf("%.2f +/- %.2f", $expected_coverage, $expected_coverage_sd);
      }
      else {
	      die "no expected coverage for: $$region{ 'chrom' }, $$region{ 'start' }, $$region{ 'end' }\n";
      }

      # Disabled at scientists request
      #if ( $expected_coverage - 2 * $expected_coverage_sd > $exons{$transcript}{ $exon }{'coverage'}{'mean'} || 
        # $expected_coverage + 2 * $expected_coverage_sd < $exons{$transcript}{ $exon }{'coverage'}{'mean'} ) {
        # $format = $$formatting{ 'orange_cell' };
      #}

      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Name' }, $exon, $format );
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Transcript' }, $transcript, $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Min depth' },  $exons{$transcript}{ $exon }{'coverage'}{'min'}, $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Mean depth' },  $exons{$transcript}{ $exon }{'coverage'}{'mean'}, $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Max depth' }, $exons{$transcript}{ $exon }{'coverage'}{'max'}, $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Position'  }, $exons{$transcript}{ $exon }{'position'}, $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Exp mean depth'  }, $expected_formatted, $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ 'Missing'  }, 
		      coverage2cell_data($exons{$transcript}{ $exon }{'coverage'}{'missing'}), $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ '1-5x'  }, coverage2cell_data($exons{$transcript}{ $exon }{'coverage'}{'1to5'}), $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ '6-9x'  }, coverage2cell_data($exons{$transcript}{ $exon }{'coverage'}{'6to9'}), $format);
      worksheet_write('Gene QC', $worksheet_offset{ 'Gene QC' }, $geneQCfield_index{ '10-19x'  }, coverage2cell_data($exons{$transcript}{ $exon }{'coverage'}{'10to19'}), $format);

      $worksheet_offset{ 'Gene QC' }++;
    }
    $worksheet_offset{ 'Gene QC' }++;
  }
  $worksheet_offset{ 'Gene QC' }++;
}


# Takes coverage data and append a length of the range if applicable
# 
# Kim Brugger (20 May 2015)
sub coverage2cell_data {
  my ( $coverage ) = @_;
  my $length = range2length( $coverage );
  return "" if ( ! $length );
  return "$coverage/$length";
  
}


# Kim Brugger (24 Mar 2014)
sub range2length {
  my ($ranges ) = @_;

  return 0 if (! $ranges || $ranges eq "" || $ranges eq 'NULL');

  my $length = 0;

  foreach my $range (split(",", $ranges)) {
    $range =~ s/.*?://;
    my ( $start, $end) = split("-", $range);
    $length += $end - $start + 1;
  }

  return $length;
}


# Kim Brugger (23 Aug 2013)
sub setup_worksheets {
  foreach my $sheet ( @sheets ) {
    add_worksheet( $sheet );
  }
}


# Kim Brugger (23 Aug 2013)
sub usage {
  $0 =~ s/.*\///;
  print STDERR "USAGE :: $0 is filtering an annotated vcf file into an xls file for interpretation\n";
  print STDERR "USAGE :: Example cmd line:\n";
  print STDERR "USAGE :: $0 -a sample.annotated.vcf -v sample.vcf -R sample.refseq_5bp.gz -C runfolder.refseq_5bp.gz ";
  print STDERR "-w workflow_name -i workflow_id -u 10000000 -T 10500000 \n";
  print STDERR "USAGE :: $0 -p \"Demantia\" -a sample.annotated.vcf -v sample.vcf -R sample.refseq_5bp.gz -C runfolder.refseq_5bp.gz ";
  print STDERR "-w workflow_name -i workflow_id -u 10000000 -T 10500000 \n";
  print STDERR "USAGE :: Additional (historical/advanced) options:\n";
  print STDERR "USAGE :: -t[ext output], default is 0\n";
  print STDERR "USAGE :: -M[eta sheet only info], default is 0\n";
  print STDERR "USAGE :: -F[ull exome report], defaut is 0\n";
  print STDERR "USAGE :: -N[o QC data in report]\n";
  print STDERR "USAGE :: -I[D QC only report]\n";
  print STDERR "USAGE :: -H (homozygous_snps)\n";
  print STDERR "USAGE :: -A (RARE_VARIANT_AF), default is 0.02\n";
  print STDERR "USAGE :: -u (nb_usable_reads)\n";
  print STDERR "USAGE :: -T (total_nb_reads)\n";
  print STDERR "USAGE :: -w (workflow)\n";
  print STDERR "USAGE :: -i (workflow_id)\n";
  exit 1;
}

 
# Kim Brugger (09 Jul 2013)
sub add_worksheet {
  my ( $sheet_name ) = @_;

  return if ( $added_worksheets{ $sheet_name } );

  my @QCfields = ('Name', 'Transcript', 'Region length', 'Min depth', 'Missing', '1-5x', '6-9x', '10-19x', '20+x', '20+x %');
  my @geneQCfields = ('Name','Transcript', 'Position', 'Min depth', 'Max depth', 'Mean depth', 'Exp mean depth', 
		      'Missing', '1-5x', '6-9x', '10-19x');
  my @fields = ('Gene', 'Transcript', 
		'Position', 'Genomic Ref Allele', 'Genomic Alt Allele', 
		'Nucleotide pos','Change', 'AA change', 'Score', 'Depth', 'AAF', 'Genotype',       
		'dbsnp', 'PolyPhen', 'SIFT',
		'AF_GEMINI',
		'AF_1KG_MAX',
		'AF_ESP_MAX',
		'AF_ExAC',
		'Comment',
  );

  $added_worksheets{ $sheet_name } = $workbook->add_worksheet( $sheet_name ) if ( !$text_only && !$meta_only);

  my $i = 0;

  # add some std cells to all work sheets.

  if ( $sheet_name =~ /Gene QC/) {
    foreach my $QCfield ( @geneQCfields ) {
      $geneQCfield_index{ $QCfield } = $i;
      worksheet_write($sheet_name, 0, $i++, $QCfield, $$formatting{ 'bold' });
    }  
  }
  elsif ( $sheet_name =~ /QC/) {
    foreach my $QCfield ( @QCfields ) {
      $QCfield_index{ $QCfield } = $i;
      worksheet_write($sheet_name, 0, $i++, $QCfield, $$formatting{ 'bold' });
    }   
    $worksheet_offset{ 'QC' } = 2;
  }
  elsif ( $sheet_name =~ /Summary/) {
    worksheet_write($sheet_name,  0, 0, "Gemini ID:", $$formatting{ 'bold' });
    worksheet_write($sheet_name,  0, 1, $sample);
    worksheet_write($sheet_name, 0, 2, "Inferred gender", $$formatting{ 'bold' });
    worksheet_write($sheet_name, 1, 2, "SRY present", $$formatting{ 'bold' });

    if ($sry) {
      worksheet_write($sheet_name, 1, 3, 'Yes');
    } else {
      worksheet_write($sheet_name, 1, 3, 'No');
    }

    worksheet_write($sheet_name,  0, 4, "Panel coverage", $$formatting{ 'bold' });
    worksheet_write($sheet_name,  0, 6, "Panel(s):", $$formatting{ 'bold' });
    worksheet_write($sheet_name, 0, 7, $gene_list{ 'PANEL'});
    worksheet_write($sheet_name,  1, 0, "GM number:", $$formatting{ 'bold' });
    worksheet_write($sheet_name,  2, 0, "Name:", $$formatting{ 'bold' });
    worksheet_write($sheet_name,  5, 0, "Report text:", $$formatting{ 'bold' });

    my $offset = 8;

    worksheet_cell_width( $sheet_name, 0, 0, 12);
    worksheet_cell_width( $sheet_name, 1, 7, 20);

    worksheet_write($sheet_name,  $offset, 1, "Phenotype:", $$formatting{'table_head'} );
    worksheet_write($sheet_name,  $offset, 2, "", $$formatting{'table_head'} );
    worksheet_write($sheet_name,  $offset, 3, "", $$formatting{'table_head'} );
    worksheet_write($sheet_name,  $offset, 4, "", $$formatting{'table_head'} );

    for( my $i = 0; $i < 4; $i++ ) {
      $offset += 1;

      for( my $j = 0; $j < 4; $j++ ) {
	      worksheet_write($sheet_name,  $offset, 1 + $j, "", $$formatting{'table'} );
      }
    }

    $offset += 3;

    worksheet_write($sheet_name,  $offset, 1, "Panels", $$formatting{ 'table_head' } );
    worksheet_write($sheet_name,  $offset, 2, "Excel file", $$formatting{'table_head'} );
    worksheet_merge_cells_and_write( $sheet_name, $offset, 3, $offset, 4, 'Comments', $$formatting{'merged_table_head'});
    worksheet_write($sheet_name,  $offset, 5, "Analysis by", $$formatting{'table_head'} );
    worksheet_write($sheet_name,  $offset, 6, "Date", $$formatting{'table_head'} );
    worksheet_write($sheet_name,  $offset, 7, "Checked by", $$formatting{'table_head'} );
    worksheet_write($sheet_name,  $offset, 8, "Date", $$formatting{'table_head'} );

    for( my $i = 0; $i < 2; $i++ ) {
      $offset += 1;

      worksheet_write($sheet_name,  $offset, 1, "", $$formatting{'table'} );
      worksheet_write($sheet_name,  $offset, 2, "", $$formatting{'table'} );
      worksheet_merge_cells_and_write( $sheet_name, $offset, 3, $offset, 4, '', $$formatting{'merged_table'});
      worksheet_write($sheet_name,  $offset, 5, "", $$formatting{'table'});
      worksheet_write($sheet_name,  $offset, 6, "", $$formatting{'table'});
      worksheet_write($sheet_name,  $offset, 7, "", $$formatting{'table'});
      worksheet_write($sheet_name,  $offset, 8, "", $$formatting{'table'});
    }

    $offset += 3;

    worksheet_merge_cells_and_write( $sheet_name, $offset, 1, $offset, 7, 'Sanger sequencing confirmation', $$formatting{'merged_table_head'});

    $offset += 1;

    worksheet_write($sheet_name,  $offset, 1, "gene", $$formatting{ 'table_head' } );
    worksheet_write($sheet_name,  $offset, 2, "NM_#", $$formatting{ 'table_head' } );
    worksheet_write($sheet_name,  $offset, 3, "coordinate", $$formatting{ 'table_head' } );
    worksheet_write($sheet_name,  $offset, 4, "cDNA", $$formatting{ 'table_head' } );
    worksheet_write($sheet_name,  $offset, 5, "protein change", $$formatting{ 'table_head' } );
    worksheet_write($sheet_name,  $offset, 6, "WS#", $$formatting{ 'table_head' } );
    worksheet_write($sheet_name,  $offset, 7, "confirmed (y/n)", $$formatting{ 'table_head' } );

    for( my $i = 0; $i < 3; $i++ ) {
      $offset += 1;
      
      for( my $j = 0; $j < 7; $j++ ) {
	      worksheet_write($sheet_name,  $offset, 1 + $j, "", $$formatting{'table'} );
      }
    }
    $offset += 3;

    worksheet_merge_cells_and_write( $sheet_name, $offset, 1, $offset, 2, 'GEM comments summary', $$formatting{'merged_table_head'});
    worksheet_merge_cells_and_write( $sheet_name, $offset, 3, $offset, 5, 'date', $$formatting{'merged_table_head'});

    for( my $i = 0; $i < 4; $i++ ) {
      $offset += 1;

      for( my $j = 0; $j < 4; $j++ ) {
        worksheet_merge_cells_and_write( $sheet_name, $offset, 1, $offset, 2, '', $$formatting{'merged_table'});
        worksheet_merge_cells_and_write( $sheet_name, $offset, 3, $offset, 5, '', $$formatting{'merged_table'});
      }
    }
    $offset += 2;
    $offset += 5;

    worksheet_write($sheet_name,  $offset, 0, "Bics/Seq QC", $$formatting{ 'bold' });
    $offset += 1;

    worksheet_write($sheet_name,  $offset, 0, "Reads:", $$formatting{ 'bold' });
    worksheet_write($sheet_name, $offset, 1, $total_nb_reads, undef);
    $offset += 1;

    worksheet_write($sheet_name, $offset, 0, "Usable Reads", $$formatting{ 'bold' });
    worksheet_write($sheet_name, $offset, 1, $nb_usable_reads, undef);
    $offset += 2;

    worksheet_write($sheet_name, $offset, 0, "Workflow", $$formatting{ 'bold' });
    worksheet_write($sheet_name, $offset, 1, $workflow, undef);
    $offset += 1;

    worksheet_write($sheet_name, $offset, 0, "Workflow id", $$formatting{ 'bold' });
    worksheet_write($sheet_name, $offset, 1, $workflow_id, undef);
    $offset += 1;
  }
  else {
    foreach my $field ( @fields ) {
      $field_index{ $field } = $i;
      worksheet_write($sheet_name, 0, $i++, $field, $$formatting{ 'bold' });    
    }    
  }

  $worksheet_offset{ $sheet_name } = 1;
}

  
# Kim Brugger (27 Aug 2013)
sub readin_gene_list {
  my ($gene_list) = @_;

  my %gene_list = ();

  open( my $in, $gene_list ) || die "Could not open '$gene_list': $!\n";

  while(<$in>) {
    chomp;
    s/\ //g;
    s/\r//g;
    my ( $gene, $transcript) = split("\t");
    next if ( /^\z/ );
    next if ( $gene eq "ALL:FULLCLINICALEXOME");

    $transcript ||= $genes2transcripts{uc($gene)};

    if ( ! $transcript || $transcript eq "") {
      print STDERR  "No transcript for $gene \n";
      system "rm -f $excel_file";

      die;
      next;
    }

    $gene_list{ uc($gene) } = uc($transcript);
  }

  return %gene_list;
}


# Kim Brugger (08 Jan 2018)
sub gene_list_to_transcript_list {
  my ($gene_list) = @_;
  my %res;
  map{ $res{ $gene_list{ $_ } } = $_ if $_ ne "PANEL"} keys %$gene_list;

  $res{ 'PANEL' } = $$gene_list{ 'PANEL'};
  
  return %res;
}


# Kim Brugger (27 Aug 2013)
sub readin_manifest {
  my ($manifest, $sample) = @_;
  my %gene_list = ();
  my %panels;

  open( my $in, $manifest ) || die "Could not open '$manifest': $!\n";

  my $sample_in_manifest = 0;

  while(<$in>) {
    chomp;
    s/\r//g;
    next if (/^\z/);
    my ( $gemini, $panel, $panel_id, $gene, $transcript ) = split("\t", $_);
    $gene =~ s/ //g;

    next if ($panel eq "BLANK");
    next if ($gene  eq "BLANK");
    next if ($gene  eq "");
    next if ( $gene eq "ALL:FULLCLINICALEXOME");
    next if ( $gemini ne $sample );

    $sample_in_manifest = 1;

    $panel_id{ uc ( $panel ) } = $panel_id;

    $panels{ $panel }++;

    $transcript ||= $genes2transcripts{uc($gene)};

    if ( ! $transcript || $transcript eq "") {
      print STDERR  "No transcript for $gene  ( $sample )\n";
      die;
      next;
    }
    
    $gene_list{ uc($gene) } = uc($transcript);
  }

  # Require either sample is in manifest or panels are specified with -p
  if (! $opts{ 'p' } ) {
    if (! $sample_in_manifest) {
      die "Sample name not found in the manifest";
    } 
  }

  $gene_list{ 'PANEL' } = join(", ", sort keys %panels );
  
  return %gene_list;
}


# Kim Brugger (03 Feb 2016)
sub worksheet_cell_width {
  my ( $id, $x, $y, $value ) = @_;
  $added_worksheets{ $id }->set_column($x, $y,  $value);
}


# Kim Brugger (01 Jul 2013)
sub worksheet_write {
  my ( $id, $x, $y, $value, $formatting ) = @_;

  return if ( not defined $value );

  $added_worksheets{ $id }->write($x, $y, $value, $formatting)  if ( !$text_only && !$meta_only);

  if ( ! $csv_workbook{ $id } ) {
    push @csv_order, $id;
  }

  $csv_workbook{ $id }[$x][$y]=$value;
}


# Kim Brugger (01 Jul 2013)
sub worksheet_write_url {
  my ( $id, $x, $y, $link, $value, $formatting ) = @_;

  return if ( not defined $value );

  $added_worksheets{ $id }->write_url($x, $y, $link, $value, $formatting)  if ( !$text_only && !$meta_only);

  if ( ! $csv_workbook{ $id } ) {
    push @csv_order, $id;
  }

  $csv_workbook{ $id }[$x][$y]=$value;
}


# Kim Brugger (03 Feb 2016)
sub worksheet_merge_cells_and_write {
  my ( $id, $x_start, $y_start, $x_end, $y_end, $value, $formatting ) = @_;

  my $range = sprintf("%s%d:%s%d", chr( ord('A') + $y_start), $x_start + 1, chr( ord('A') + $y_end), $x_end + 1);

  $added_worksheets{ $id }->merge_range($range, $value, $formatting);
}


# Kim Brugger (02 Jun 2010)
#
# Now handles more than a single AA
#
# Kim Brugger (14 Apr 2016)
sub one2three {
  my ( $aminoacid) = @_;

  return "" if ( ! defined $aminoacid);

  my %trans = ('A' => 'Ala',
	       'R' => 'Arg',
	       'N' => 'Asn',
	       'D' => 'Asp',
	       'C' => 'Cys',
	       'E' => 'Glu',
	       'Q' => 'Gln',
	       'G' => 'Gly',
	       'H' => 'His',
	       'I' => 'Ile',
	       'L' => 'Leu',
	       'K' => 'Lys',
	       'M' => 'Met',
	       'F' => 'Phe',
	       'P' => 'Pro',
	       'S' => 'Ser',
	       'T' => 'Thr',
	       'W' => 'Trp',
	       'Y' => 'Tyr',
	       'V' => 'Val',
	       '*' => 'Ter');

  my $three = "";

  foreach my $aa ( split("", $aminoacid)) {
    if ( ! $trans{ $aa } ) {
      $three .= "Xxx";
    }
    else {
      $three .= $trans{ $aa }
    }
  }

  return ( $three );
  return $trans{ $aminoacid } if ($trans{ $aminoacid });
  return $aminoacid;
}


# Kim Brugger (01 Dec 2010)
sub grantham_score {
  my ($aa1, $aa2) = @_;

  $aa1 = uc($aa1);
  $aa2 = uc($aa2);

  return 0 if ( $aa1 eq $aa2);
  
  ($aa1, $aa2) = ($aa2, $aa1) if ($aa1 gt $aa2);
  
  my %grantham;
  $grantham{ALA}{ARG}=112; 
  $grantham{ALA}{ASN}=111; 
  $grantham{ALA}{ASP}=126; 
  $grantham{ALA}{CYS}=195; 
  $grantham{ALA}{GLN}=91; 
  $grantham{ALA}{GLU}=107; 
  $grantham{ALA}{GLY}=60; 
  $grantham{ALA}{HIS}=86; 
  $grantham{ALA}{ILE}=94; 
  $grantham{ALA}{LEU}=96; 
  $grantham{ALA}{LYS}=106; 
  $grantham{ALA}{MET}=84; 
  $grantham{ALA}{PHE}=113; 
  $grantham{ALA}{PRO}=27; 
  $grantham{ALA}{SER}=99; 
  $grantham{ALA}{THR}=58; 
  $grantham{ALA}{TRP}=148; 
  $grantham{ALA}{TYR}=112; 
  $grantham{ALA}{VAL}=64;

  $grantham{ARG}{ASN}=86; 
  $grantham{ARG}{ASP}=96; 
  $grantham{ARG}{CYS}=180; 
  $grantham{ARG}{GLN}=43; 
  $grantham{ARG}{GLU}=54; 
  $grantham{ARG}{GLY}=125; 
  $grantham{ARG}{HIS}=29; 
  $grantham{ARG}{ILE}=97; 
  $grantham{ARG}{LEU}=102; 
  $grantham{ARG}{LYS}=26; 
  $grantham{ARG}{MET}=91; 
  $grantham{ARG}{PHE}=97; 
  $grantham{ARG}{PRO}=103; 
  $grantham{ARG}{SER}=110; 
  $grantham{ARG}{THR}=71; 
  $grantham{ARG}{TRP}=101; 
  $grantham{ARG}{TYR}=77; 
  $grantham{ARG}{VAL}=96; 

  $grantham{ASN}{ASP}=23; 
  $grantham{ASN}{CYS}=139; 
  $grantham{ASN}{GLN}=46; 
  $grantham{ASN}{GLU}=42; 
  $grantham{ASN}{GLY}=80; 
  $grantham{ASN}{HIS}=68; 
  $grantham{ASN}{ILE}=149; 
  $grantham{ASN}{LEU}=153; 
  $grantham{ASN}{LYS}=94; 
  $grantham{ASN}{MET}=142; 
  $grantham{ASN}{PHE}=158; 
  $grantham{ASN}{PRO}=91; 
  $grantham{ASN}{SER}=46; 
  $grantham{ASN}{THR}=65;
  $grantham{ASN}{TRP}=174;
  $grantham{ASN}{TYR}=143;
  $grantham{ASN}{VAL}=133;

  $grantham{ASP}{CYS}=154; 
  $grantham{ASP}{GLN}=61; 
  $grantham{ASP}{GLU}=45; 
  $grantham{ASP}{GLY}=94; 
  $grantham{ASP}{HIS}=81; 
  $grantham{ASP}{ILE}=168; 
  $grantham{ASP}{LEU}=172; 
  $grantham{ASP}{LYS}=101; 
  $grantham{ASP}{MET}=160; 
  $grantham{ASP}{PHE}=177; 
  $grantham{ASP}{PRO}=108; 
  $grantham{ASP}{SER}=65; 
  $grantham{ASP}{THR}=85; 
  $grantham{ASP}{TRP}=181; 
  $grantham{ASP}{TYR}=160; 
  $grantham{ASP}{VAL}=152;

  $grantham{CYS}{GLN}=154; 
  $grantham{CYS}{GLU}=170; 
  $grantham{CYS}{GLY}=159; 
  $grantham{CYS}{HIS}=174; 
  $grantham{CYS}{ILE}=198; 
  $grantham{CYS}{LEU}=198; 
  $grantham{CYS}{LYS}=202; 
  $grantham{CYS}{MET}=196; 
  $grantham{CYS}{PHE}=205; 
  $grantham{CYS}{PRO}=169; 
  $grantham{CYS}{SER}=112; 
  $grantham{CYS}{THR}=149; 
  $grantham{CYS}{TRP}=215; 
  $grantham{CYS}{TYR}=194; 
  $grantham{CYS}{VAL}=192;

  $grantham{GLN}{GLU}=29; 
  $grantham{GLN}{GLY}=87; 
  $grantham{GLN}{HIS}=24; 
  $grantham{GLN}{ILE}=109; 
  $grantham{GLN}{LEU}=113; 
  $grantham{GLN}{LYS}=53; 
  $grantham{GLN}{MET}=101; 
  $grantham{GLN}{PHE}=116; 
  $grantham{GLN}{PRO}=76;
  $grantham{GLN}{SER}=68; 
  $grantham{GLN}{THR}=42; 
  $grantham{GLN}{TRP}=130; 
  $grantham{GLN}{TYR}=99; 
  $grantham{GLN}{VAL}=96;

  $grantham{GLU}{GLY}=98; 
  $grantham{GLU}{HIS}=40; 
  $grantham{GLU}{ILE}=134; 
  $grantham{GLU}{LEU}=138; 
  $grantham{GLU}{LYS}=56; 
  $grantham{GLU}{MET}=126; 
  $grantham{GLU}{PHE}=140; 
  $grantham{GLU}{PRO}=93; 
  $grantham{GLU}{SER}=80; 
  $grantham{GLU}{THR}=65; 
  $grantham{GLU}{TRP}=152; 
  $grantham{GLU}{TYR}=122; 
  $grantham{GLU}{VAL}=121;

  $grantham{GLY}{HIS}=89; 
  $grantham{GLY}{ILE}=135; 
  $grantham{GLY}{LEU}=138; 
  $grantham{GLY}{LYS}=127; 
  $grantham{GLY}{MET}=127; 
  $grantham{GLY}{PHE}=153; 
  $grantham{GLY}{PRO}=42; 
  $grantham{GLY}{SER}=56; 
  $grantham{GLY}{THR}=59; 
  $grantham{GLY}{TRP}=184; 
  $grantham{GLY}{TYR}=147; 
  $grantham{GLY}{VAL}=109;
  
  $grantham{HIS}{ILE}=94; 
  $grantham{HIS}{LEU}=99; 
  $grantham{HIS}{LYS}=32; 
  $grantham{HIS}{MET}=87; 
  $grantham{HIS}{PHE}=100; 
  $grantham{HIS}{PRO}=77; 
  $grantham{HIS}{SER}=89; 
  $grantham{HIS}{THR}=47; 
  $grantham{HIS}{TRP}=115; 
  $grantham{HIS}{TYR}=83; 
  $grantham{HIS}{VAL}=84; 
  
  $grantham{ILE}{LEU}=5; 
  $grantham{ILE}{LYS}=102; 
  $grantham{ILE}{MET}=10; 
  $grantham{ILE}{PHE}=21; 
  $grantham{ILE}{PRO}=95; 
  $grantham{ILE}{SER}=142; 
  $grantham{ILE}{THR}=89; 
  $grantham{ILE}{TRP}=61; 
  $grantham{ILE}{TYR}=33; 
  $grantham{ILE}{VAL}=29;
  
  $grantham{LEU}{LYS}=107; 
  $grantham{LEU}{MET}=15; 
  $grantham{LEU}{PHE}=22; 
  $grantham{LEU}{PRO}=98; 
  $grantham{LEU}{SER}=145; 
  $grantham{LEU}{THR}=92; 
  $grantham{LEU}{TRP}=61; 
  $grantham{LEU}{TYR}=36; 
  $grantham{LEU}{VAL}=32;
  
  $grantham{LYS}{MET}=95; 
  $grantham{LYS}{PHE}=102; 
  $grantham{LYS}{PRO}=103; 
  $grantham{LYS}{SER}=121; 
  $grantham{LYS}{THR}=78; 
  $grantham{LYS}{TRP}=110; 
  $grantham{LYS}{TYR}=85; 
  $grantham{LYS}{VAL}=97;
  
  $grantham{MET}{PHE}=28; 
  $grantham{MET}{PRO}=87; 
  $grantham{MET}{SER}=135; 
  $grantham{MET}{THR}=81; 
  $grantham{MET}{TRP}=67; 
  $grantham{MET}{TYR}=36; 
  $grantham{MET}{VAL}=21;
  
  $grantham{PHE}{PRO}=114; 
  $grantham{PHE}{SER}=155; 
  $grantham{PHE}{THR}=103; 
  $grantham{PHE}{TRP}=40; 
  $grantham{PHE}{TYR}=22; 
  $grantham{PHE}{VAL}=50;
  
  $grantham{PRO}{SER}=74; 
  $grantham{PRO}{THR}=38; 
  $grantham{PRO}{TRP}=147; 
  $grantham{PRO}{TYR}=110; 
  $grantham{PRO}{VAL}=68;
  
  $grantham{SER}{THR}=58; 
  $grantham{SER}{TRP}=177; 
  $grantham{SER}{TYR}=144; 
  $grantham{SER}{VAL}=124;
  
  $grantham{THR}{TRP}=128; 
  $grantham{THR}{TYR}=92; 
  $grantham{THR}{VAL}=69;
  
  $grantham{TRP}{TYR}=37; 
  $grantham{TRP}{VAL}=88;
  
  $grantham{TYR}{VAL}=55;

  return $grantham{$aa1}{$aa2} if ($grantham{$aa1}{$aa2});
  return "NA";
}


# Kim Brugger (23 Apr 2015)
sub readin_panels_n_manifest {
  undef %genepanels;
  
  open(my $in, $genes2transcripts_file);

  map{ 
    my @F = split(/\s+/, $_); 
    $genes2transcripts{ uc($F[0])} = uc($F[1]) if ( $F[1] );
    $transcript2gene{uc($F[1])} = $F[0] if ( $F[1] );
  } <$in>;

  close($in);

  open( $in, $genepanels_file);

  map{ 
    chomp; 
    my ( $panel_name, $panel_id, $gene) = split(/\t/, $_); 

    if ( $gene && ($gene ne "BLANK" or $gene ne "")) {
      push @{$genepanels{ uc( $panel_name )}}, uc( $gene );
      $panel_names{ uc($panel_name) } = $panel_name;
      $panel_id{ uc ( $panel_name ) } = $panel_id;
    }
  } <$in>;

  close( $in );
}


# Kim Brugger (23 Apr 2015)
sub parameter_panels2genes {
  my ( $param ) = @_;
  my %gene_list;
  my %panels;

  $param =~ s/^\s+//;
  $param =~ s/^\s\z//;

  foreach my $panel (split(",", $param)) {
    $panel =~ s/^\s+//;
    $panel =~ s/\s+\z//g;

    $panel = $panel_names{ uc( $panel )} if ( $panel_names{ uc( $panel )} );

    $panels{ $panel }++;
    # Panels with an _ at the start are single genes. So make an ad hoc single gene panel when needed.
    if ( $panel =~ /^_/ ) {
      my $gene = $panel;
      $gene =~ s/^_//;
      push @{$genepanels{ uc( "$panel" )}}, uc( $gene );
    }

    if ( ! $genepanels{ uc( $panel )} ) {
      print STDERR  "Unknown panel requested: $panel\n";
      find_closest_panel( $panel );
      exit -1;
    }

    foreach my $gene ( @{$genepanels{ uc( $panel )}} ) {
      next if ( $gene eq 'BLANK');
      my $transcript ||= $genes2transcripts{uc($gene)};
      
      if ( ! $transcript || $transcript eq "") {
        print STDERR  "No transcript for $gene ( $sample )\n";
        die;
        next;
      }
      $gene_list{ uc($gene) } = uc($transcript);
    }
  }

  $gene_list{ 'PANEL' } = join(", ", sort keys %panels);

  return %gene_list;
}


# If the panel name is not found, we do a binary search through the panel array to find the best candidates.
# 
# Kim Brugger (24 Apr 2015)
sub find_closest_panel {
  my ( $panel_name ) = @_;

  $panel_name =~ s/^\s+//;
  $panel_name =~ s/\s+\z//;

  use POSIX qw(ceil floor);

  my @panels;

  map { push @panels, $panel_names{ $_ }}  sort { $a cmp $b } keys %panel_names;

  $| = 1;
  # set the start and end of the array and find the 
  # the middle of the array
  my ( $left, $right ) = (0, int(@panels) - 1);
  my $middle = floor(($right - $left)/2);
  
  # Flush the buffer constantly
  $| = 1;
  
  my $loop_counter = 0;

  while (1) {
    # The new block is to the left of the middle.
    if ( uc($panel_name) lt uc($panels[ $middle ]) ) {
      $right = $middle;
      $middle = $left + floor(($right - $left)/2);
      last if ( $right <= $left || $middle == $left || $middle == $right);
    }
    # The new block is to the right of the middle.
    elsif ( uc($panel_name) gt uc($panels[ $middle ]) ) {
      $left = $middle;
      $middle = $left + floor(($right - $left)/2);
      last if ( $right <= $left || $middle == $left || $middle == $right);
    }
    
    # Now things gets interesting, we here start to calculate
    # overlapping and contained panels.
    #
    # this is a contained snp, exactly what we want!!!!
    elsif ( uc($panel_name) ge uc($panels[ $middle ])  &&
            uc($panel_name) le uc($panels[ $middle ]) ) {
      return 1;
      last;
    }
    else {
      last;
    }
  }

  my @syndromes;

  my $start = $middle - 2;
  $start = 0 if ( $start < 0 );

  for(my $i = $start; $i < $start + 5; $i++) {
    push @syndromes, $panels[ $i ];
  }

  print "Closest syndromes matching :  \n";
  print "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n";

  print join("\n", @syndromes ) . "\n";

  return 0;
}


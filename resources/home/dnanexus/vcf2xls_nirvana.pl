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

my $opts = 'p:a:v:u:T:w:i:c:h';
my %opts;
getopts($opts, \%opts);

my $samtools  = 'packages/samtools-1.7/samtools';

my $RARE_VARIANT_AF = 0.02;
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

my %genes2transcripts;
my %transcript2gene;
my %VEP_transcripts;
my %genepanels;
my %clinical_ind_panels;
my %panel_names;
my %panel_id;

readin_panels_n_manifest();

my @QC_sheets = ('Summary');

my @variant_sheets = (
              'stop_gained',
              'frameshift_variant',
              'consensus splice',
              'missense_variant',
              'synonymous_variant',
              'other',
    );

my @sheets;
@sheets = (@QC_sheets, @variant_sheets);

my %effect_levels = (
    'stop_gained'        => 6,
    'frameshift_variant' => 5,
    'consensus splice'   => 4,
    'missense_variant'   => 3,
    'synonymous_variant' => 2,
    'other'              => 1
);

my %field_index = ();

my $MIN_QUALITY_VAR  = 200;
my $LOW_COVERAGE_VAR = 10;

my $homozygous_snps = 0;
my $full_exome_snps = 0;

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

my %gene_list;
my %hotspots;

if ( $opts{ 'p' } ) {
  %gene_list = parameter_panels2genes($opts{ 'p' }, $sample);
} else {
  %gene_list = readin_manifest( $manifest, $sample);
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
my @clinical_inds = map {
  $_ =~ s/^\ +//;
  my $v = $panel_id{ uc($_) } || "NA"; 
  $_ = "$_"
}  split(",", $transcript_list{ 'PANEL'});

my @panel_w_version = map {
  $_ =~ s/^\ +//;
  my $v = $panel_id{ uc($_) } || "NA"; 
  $_ = "$v"
}  split(",", $transcript_list{ 'PANEL'});

$gene_list{ 'Clinical indication(s)'} = join(", ", @clinical_inds );
$gene_list{ 'Panel id(s)'} = join(", ", @panel_w_version );

my $coverage_file = $opts{c};

my $sry;
$sry = check_sry($coverage_file);

setup_worksheets();

my %meta_stats = ();
$meta_stats{ 'PANEL'} = $gene_list{ 'PANEL'};
$meta_stats{ 'PANEL_IDS'} = $gene_list{ 'PANEL_IDS'};

analyse_vcf_file( $vcf_file );

print "Filling summary sheet\n";
fill_summary_sheet();

$workbook->close();

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
  $sample =~ m/^X[0-9]+/;

  return $& ;
}

# Kim Brugger (20 May 2015)
sub fill_summary_sheet {
  if ( $gene_list{ 'Panel id(s)'} ) {
    worksheet_write('Summary', 1 ,  4 , "Panel(s)", $$formatting{ 'bold' });
    worksheet_write('Summary', 1 ,  5 , $gene_list{ 'Panel id(s)'});
  }
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

    my $gt1;
    my $gt2;

    if ( $$entry{ gtypes }{ $sample }) {
      # Pull out the genotypes for the sample
      ($gt1, $gt2) = split("/",$$entry{ gtypes }{ $sample }{ GT });
    }
    else {
      my $sample_in_vcf;

      for my $sub_hash ( $$entry{gtypes} ) {
        my %sample2alleles = %$sub_hash;
        $sample_in_vcf = join "\n", keys %sample2alleles;
      }

      die "ERROR: \"$sample\" extracted from the file name is different than \"$sample_in_vcf\" found in the vcf file";
    }

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

  $added_worksheets{ $sheet_name } = $workbook->add_worksheet( $sheet_name );

  my $i = 0;

  # add some std cells to all work sheets.

  if ( $sheet_name =~ /Summary/) {
    worksheet_write($sheet_name,  0, 0, "Gemini ID:", $$formatting{ 'bold' });
    worksheet_write($sheet_name,  0, 1, $sample);
    worksheet_write($sheet_name, 1, 2, "SRY present", $$formatting{ 'bold' });

    if ($sry) {
      worksheet_write($sheet_name, 1, 3, 'Yes');
    } else {
      worksheet_write($sheet_name, 1, 3, 'No');
    }

    worksheet_write($sheet_name,  0, 4, "Clinical indication(s):", $$formatting{ 'bold' });
    worksheet_write($sheet_name, 0, 5, $gene_list{ 'Clinical indication(s)'});
    worksheet_write($sheet_name,  1, 0, "GM number:", $$formatting{ 'bold' });
    worksheet_write($sheet_name,  2, 0, "Name:", $$formatting{ 'bold' });

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
  my %clinical_inds;

  open( my $in, $manifest ) || die "Could not open '$manifest': $!\n";

  my $sample_in_manifest = 0;

  while(<$in>) {
    chomp;
    s/\r//g;
    next if (/^\z/);
    my ( $gemini, $clinical_ind, $panel_id, $gene, $transcript) = split("\t", $_);
    $gene =~ s/ //g;

    next if ($clinical_ind eq "BLANK");
    next if ($gene  eq "BLANK");
    next if ($gene  eq "");
    next if ( $gene eq "ALL:FULLCLINICALEXOME");
    next if ( $gemini ne $sample );

    $sample_in_manifest = 1;

    $panel_id{ uc ( $clinical_ind ) } = $panel_id;

    $clinical_inds{ $clinical_ind }++;

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

  $gene_list{ 'PANEL' } = join(", ", sort keys %clinical_inds );
  
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

  $added_worksheets{ $id }->write($x, $y, $value, $formatting);

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


# Kim Brugger (23 Apr 2015)
sub readin_panels_n_manifest {
  undef %genepanels;
  undef %clinical_ind_panels;
  
  open(my $in, $genes2transcripts_file);

  map{ 
    my @F = split(/\s+/, $_); 
    $genes2transcripts{ uc($F[0])} = uc($F[1]) if ( $F[1] && $F[2] eq "clinical_transcript" );
    $transcript2gene{uc($F[1])} = $F[0] if ( $F[1] && $F[2] eq "clinical_transcript" );
  } <$in>;

  close($in);

  open( $in, $genepanels_file);

  map{ 
    chomp; 
    my ( $clinical_ind, $panel_id, $gene) = split(/\t/, $_); 

    if ( $gene && ($gene ne "BLANK" or $gene ne "")) {
      push @{$genepanels{ uc( $clinical_ind )}}, uc( $gene );
      push @{$clinical_ind_panels{ uc( $clinical_ind )}}, uc( $panel_id );
      $panel_names{ uc($clinical_ind) } = $clinical_ind;
      $panel_id{ uc ( $clinical_ind ) } = $panel_id;
    }
  } <$in>;

  close( $in );
}


# Kim Brugger (23 Apr 2015)
sub parameter_panels2genes {
  my ( $param ) = @_;
  my %gene_list;
  my %clinical_inds;

  $param =~ s/^\s+//;
  $param =~ s/^\s\z//;

  foreach my $panel (split(";", $param)) {
    $panel =~ s/^\s+//;
    $panel =~ s/\s+\z//g;

    # Panels with an _ at the start are single genes. So make an ad hoc single gene panel when needed.
    if ( $panel =~ /^_/ ) {
      my $gene = $panel;
      $gene =~ s/^_//;
      push @{$genepanels{ uc( "$panel" )}}, uc( $gene );
    }

    if ( $panel_names{ uc( $panel )} ) {
      my $clinical_ind = $panel_names{ uc( $panel )};
      my $panels = $panel_id{ uc($clinical_ind) };
    }

    $clinical_inds{ $panel }++;

    if ( ! $genepanels{ uc( $panel )} ) {
      print STDERR  "Unknown panel requested: $panel\n";
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

  $gene_list{ 'PANEL' } = join(", ", sort keys %clinical_inds);

  return %gene_list;
}

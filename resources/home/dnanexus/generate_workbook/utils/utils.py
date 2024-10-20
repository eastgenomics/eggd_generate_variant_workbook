from csv import Sniffer

import pandas as pd


def is_numeric(value:str) -> bool:
    """
    Returns true if given value is in some form numeric

    Parameters
    ----------
    value : str
        string to check if can be cast to numeric type

    Returns
    -------
    bool
        True if value can be numeric
    """
    try:
        if str(value)[0] == '0' and not str(value)[1] == '.':
            # catch things such as 01 which aren't real numbers
            return False
    except Exception:
        # could raise lots of errors (i.e IndexErrors from slicing)
        # lazily catch all and continue
        pass

    try:
        float(value)
        return True
    except ValueError:
        return False


def determine_delimiter(data, suffixes) -> None:
    """
    Attempt to determine delimiter from a given string and list of
    file suffixes.

    Will check for tsv or csv in suffixes and return accordingly, if not
    will infer from data using csv.Sniffer, defaults to tabs if it can't
    be determined.

    Parameters
    ----------
    data : str
        data to check for delimiter
    suffixes : list
        list of file suffixes

    Returns
    -------
    delimiter : str
        delimiter inferred from given data
    """
    if '.tsv' in suffixes:
        return '\t'

    if '.csv' in suffixes:
        return ','

    try:
        delimiter = Sniffer().sniff(str(data)).delimiter
    except Exception as error:
        print(
            "Error in determing delimiter from given data. Will default "
            f"to using tabs.\n\nError: {error}\n\n"
        )
        delimiter = '\t'

    return delimiter


def map_chr_to_nc(chrom, build) -> str:
    """
    Maps given chromosome to NC value for generating MasterMind URL, taken
    from here: https://www.genomenon.com/mastermind-faq/#How%20can%20I%20search%20for%20variants%20in%20Mastermind

    Parameters
    ----------
    chrom : str
        chromosome to return NC value of
    build : int
        reference build inferred from reference stored in vcf header,
        will be either 37 or 38

    Returns
    -------
    nc_id : str
        mapped NC value from chromosome
    """
    mapping = {
        "1": {37: "NC_000001.10", 38: "NC_000001.11"},
        "2": {37: "NC_000002.11", 38: "NC_000002.12"},
        "3": {37: "NC_000003.11", 38: "NC_000003.12"},
        "4": {37: "NC_000004.11", 38: "NC_000004.12"},
        "5": {37: "NC_000005.9", 38: "NC_000005.10"},
        "6": {37: "NC_000006.11", 38: "NC_000006.12"},
        "7": {37: "NC_000007.13", 38: "NC_000007.14"},
        "8": {37: "NC_000008.10", 38: "NC_000008."},
        "9": {37: "NC_000009.11", 38: "NC_000009.12"},
        "10": {37: "NC_000010.10", 38: "NC_000010.11"},
        "11": {37: "NC_000011.9", 38: "NC_000011.10"},
        "12": {37: "NC_000012.11", 38: "NC_000012.12"},
        "13": {37: "NC_000013.10", 38: "NC_000013.11"},
        "14": {37: "NC_000014.8", 38: "NC_000014.9"},
        "15": {37: "NC_000015.9", 38: "NC_000015.10"},
        "16": {37: "NC_000016.9", 38: "NC_000016.10"},
        "17": {37: "NC_000017.10", 38: "NC_000017.11"},
        "18": {37: "NC_000018.9", 38: "NC_000018.10"},
        "19": {37: "NC_000019.9", 38: "NC_000019.10"},
        "20": {37: "NC_000020.10", 38: "NC_000020.11"},
        "21": {37: "NC_000021.8", 38: "NC_000021.9"},
        "22": {37: "NC_000022.10", 38: "NC_000022.11"},
        "X": {37: "NC_000023.10", 38: "NC_000023.11"},
        "Y": {37: "NC_000024.9", 38: "NC_000024.10"}
    }

    nc_id = mapping.get(chrom, None)

    if nc_id:
        nc_id = nc_id.get(build, None)

    return nc_id


class buildHyperlink():
    """
    Functions for generating annotation resource specific hyperlinks to
    display in the workbook. Each link is formatted slightly different, so
    they get their own hand crafted artisanal functions to generate them in
    the correct format from the below base URLs.
    """
    def __init__(self) -> None:
        self.urls = {
            "existing_variation": "https://www.ncbi.nlm.nih.gov/snp/",
            "clinvar": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
            "cosmic": "https://cancer.sanger.ac.uk/cosmic/search?genome=BUILD&q=",
            "hgmd": "https://my.qiagendigitalinsights.com/bbp/view/hgmd/pro/mut.php?acc=",
            "mastermind": "https://mastermind.genomenon.com/detail?mutation=",
            "gnomad_base_url": "https://gnomad.broadinstitute.org/variant/CHROM-POS-REF-ALT",
            "decipher": "https://www.deciphergenomics.org/sequence-variant/CHROM-POS-REF-ALT",
            "oncokb": "https://www.oncokb.org/gene/",
            "cbioportal": "https://www.cbioportal.org/results/mutations?case_set_id=all&gene_list=SYMBOL&cancer_study_list=5c8a7d55e4b046111fee2296",
            "pecan": "https://pecan.stjude.cloud/variants/proteinpaint/?gene=SYMBOL"
        }

    def build(self, column, value, build) -> str:
        """
        Build annotation resource specific hyperlinks from the above given URLs

        Parameters
        ----------
        column : str
            current column for which to add URL
        value : pd.Series
            series of current row data, used to get column value from as well
            as any other required fields (i.e. chrom, pos, ref, alt etc.)
        build : int
            either 37 or 38, inferred perviously from reference parsed from
            vcf header, controls URL to use for those that are build specific

        Returns
        -------
        str
            Excel formatted hyperlink if a URL is present for that resource,
            else just the original value is returned
        """
        if (
            not value[column] or
            value[column] == '.' or
            value[column] == 'nan' or
            pd.isna(value[column])
        ):
            return value[column]

        url = None

        # partially match against column names and add appropriate hyperlink
        # for clinvar and hgmd, only match where the column will be named
        # with CSQ_ (i.e. CSQ_HGMD) prefix and not capture other columns such
        # as CSQ_HGMD_CLASS which do not need linking and would not
        # build a valid URL
        if 'gnomad' in column.lower():
            url = self.gnomad(value, build)
        elif 'cosmic' in column.lower():
            url = self.cosmic(value, build, column)
        elif 'existing_variation' in column.lower():
            url = self.existing_variation(value[column])
        elif 'mastermind' in column.lower() or 'mmid3' in column.lower():
            url = self.mastermind(value, build, column)
            value[column] = url.split('=')[1]
        elif column.lower().endswith('clinvar'):
            url = self.clinvar(value[column])
        elif column.lower().endswith('hgmd'):
            url = self.hgmd(value[column])

        # below will be exact matches as columns are from --additional_columns
        # for each we are setting the value to just be the gene symbol since
        # they aren't 'real' annotation from the vcf and therefore are empty
        if column.lower() == 'decipher':
            url = self.decipher(value, build, column)
            value[column] = url.split('/')[-1]
        elif column.lower() == 'oncokb':
            url = self.oncokb(value)
            value[column] = value.CSQ_SYMBOL
        elif column.lower() == 'cbioportal':
            url = self.cbioportal(value)
            value[column] = value.CSQ_SYMBOL
        elif column.lower() == 'pecan':
            url = self.pecan(value)
            value[column] = value.CSQ_SYMBOL

        if not url:
            # URL not generated, likely just a regular vcf column
            return value[column]

        if len(url) > 255:
            # Excel has a string limit of 255 characters inside a formula,
            # if URL is too long just display the value
            return value[column]

        if is_numeric(value[column]):
            # return numeric values not wrapped in quotes
            return f'=HYPERLINK("{url}", {value[column]})'
        else:
            # values for everything else which is hyperlinked
            # needs to be cast to string
            return f'=HYPERLINK("{url}", "{value[column]}")'

    def gnomad(self, value, build) -> str:
        # gnomad URL has build specific suffix
        if build == 37:
            url = f"{self.urls['gnomad_base_url']}?dataset=gnomad_r2_1"
        elif build == 38:
            url = f"{self.urls['gnomad_base_url']}?dataset=gnomad_r3"
        else:
            return value[column]

        url = url.replace('CHROM', str(value.CHROM).replace('chr', ''))
        url = url.replace('POS', str(value.POS))
        url = url.replace('REF', str(value.REF))
        url = url.replace('ALT', str(value.ALT))

        return url

    def cosmic(self, value, build, column) -> str:
        url = self.urls.get('cosmic').replace('BUILD', str(build))
        url = f'{url}{value[column]}'

        return url

    def existing_variation(self, rsid) -> str:
        if not rsid.startswith('rs'):
            # non-rsID in Existing_variation column => return without URL
            return None

        return f"{self.urls.get('existing_variation')}{rsid}"

    def mastermind(self, value, build, column) -> str:
        if not build:
            # no reference build, can't generate URL
            return value[column]

        # get NC value for chromosome
        nc_id = map_chr_to_nc(str(value.CHROM).replace('chr', ''), build)

        # build URL and set value to display equal to what is in the URL
        url = self.urls.get('mastermind')
        url = f'{url}{nc_id}:g.{value.POS}{value.REF}%3E{value.ALT}'

        return url

    def clinvar(self, clinvar_id) -> str:
        return f"{self.urls.get('clinvar')}{clinvar_id}"

    def hgmd(self, hgmd_id) -> str:
        return f"{self.urls.get('hgmd')}{hgmd_id}"

    def decipher(self, value, build, column) -> str:
        url = self.urls.get('decipher')
        url = url.replace('CHROM', str(value.CHROM).replace('chr', ''))
        url = url.replace('POS', str(value.POS))
        url = url.replace('REF', str(value.REF))
        url = url.replace('ALT', str(value.ALT))

        return url

    def oncokb(self, value) -> str:
        return f"{self.urls.get('oncokb')}{value.CSQ_SYMBOL}"

    def cbioportal(self, value) -> str:
        url = self.urls.get('cbioportal')
        return f"{url.replace('SYMBOL', value.CSQ_SYMBOL)}"

    def pecan(self, value) -> str:
        url = self.urls.get('pecan')
        return f"{url.replace('SYMBOL', value.CSQ_SYMBOL)}"


def parse_cvo(cvo_df) -> pd.DataFrame:
    """
    Parse MSI, TMB and Gene Amplifications section from Illumina TSO500
    local app CombinedVariantOutput.tsv file. Rest of file contains file
    metrics and SNVs, but as we already have the VCF of variants we will
    exclude this to just keep the additional metrics.

    This file is structured with the following sections:

        [Analysis Details]
        ...

        [TMB]
        ...

        [MSI]
        ...

        [Gene Amplifications]
        ...

        [Splice Variants]
        ...


    Parameters
    ----------
    cvo_df : pd.DataFrame
        Dataframe of CombinedVariantOutput file to parse from

    Returns
    -------
    pd.DataFrame
        MSI, TMB and Gene Amplifications sections of dataframe
    """
    # get row indexes of TMB and Splice Variants to slice from dataframe
    tmb_idx = cvo_df[0].eq('[TMB]').idxmax()
    splice_idx = cvo_df[0].eq('[Splice Variants]').idxmax()

    if tmb_idx==0 or splice_idx==0:
        # didn't correctly parse out fields, return original dataframe to write
        print(
            'Could not parse TMB and Splice Variants fields from '
            'CombinedVariantOutput file, full file will be written to the '
            'additional sheet'
        )
        return cvo_df

    return cvo_df.iloc[tmb_idx:splice_idx]


def parse_metrics_output(metrics_df, sample_vcf) -> pd.DataFrame:
    """
    Parse out sample metrics from run level MetricsOutput.tsv TSO500 local
    app output file if provided to --additional_files.

    Uses first VCF file prefix to assume as sample name for parsing out
    of metrics file.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        DataFrame of all metrics from MetricsOutput.tsv
    sample_vcf : str
        name of sample vcf file

    Returns
    -------
    pd.DataFrame
        parsed MetricsOutput dataframe for given sample
    """
    # get index of where per sample metrics start in file, -2 to drop footer
    metrics_idx = metrics_df[0].eq('[DNA Library QC Metrics]').idxmax()

    if not metrics_idx:
        # file doesn't have expected field for MetricsOutput
        print(
            'WARNING: Could not parse "[DNA Library QC Metrics]" from '
            'MetricsOutput.tsv. Writing whole file to sheet.'
        )
        return metrics_df

    # select metrics from file without top section of run header and footer
    metrics_df = metrics_df.iloc[metrics_idx+1:-2].reset_index(drop=True)

    # get column of given samples metrics
    vcf_prefix = sample_vcf.split('_')[0]

    # get column index of our sample
    idx = [
        idx for idx, column in enumerate(metrics_df.iloc[0])
        if column.startswith(vcf_prefix)
    ]

    if not idx:
        # can't find sample, return whole file
        print(
            f'WARNING: Could not parse sample from MetricsOutput.tsv with '
            f'prefix: {vcf_prefix}. Writing whole file to sheet.'
        )
        return metrics_df

    if len(idx) > 1:
        # found more than one match for given prefix, return whole file
        print(
            'WARNING: Found more than one sample match in MetricsOutput.tsv '
            f'for prefix: {vcf_prefix}. Writing whole file to sheet.'
        )
        return metrics_df

    # select first 3 cols with labels and our sample column
    metrics_df = metrics_df.iloc[:, [0, 1, 2, idx[0]]]

    return metrics_df

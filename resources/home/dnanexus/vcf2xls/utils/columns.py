import sys
from typing import Union
import pandas as pd


class splitColumns():
    """
    Functions for spliting and checking of variant annotation and
    attribute columns (FORMAT, SAMPLE, INFO, CSQ), called during reading
    of VCFs from file in vcf.read()
    """
    @staticmethod
    def format_fields(vcf_df) -> pd.DataFrame:
        """
        Get format fields from FORMAT column to split out sample values to
        individual columns, this transforms the data as such:

        -----------------------------------------------------------------------
        FORMAT                 | SAMPLE
        -----------------------------------------------------------------------
        GT:AD:DP:GQ:PL         | 0/1:83,88:171:99:2136,0,1880
        GT:AD:DP:GQ:PL         | 0/1:128,128:256:99:3163,0,3015
        GT:AD:DP:GQ:PGT:PID:PL | 0/1:120,93:213:99:0|1:46896303_C_G:3367,0,6381
        -----------------------------------------------------------------------
                                      |
                                      |
                                      ▼
        -----------------------------------------------------------------------
          GT  | AD      |  DP   |  GQ  |  PGT  |   PID         |  PL
        -----------------------------------------------------------------------
          0/1 | 83,88   |  171  |  99  |  na   |  na           |  2136,0,1880
          0/1 | 128,128 |  256  |  99  |  na   |  na           |  3163,0,3015
          0/1 | 120,93  |  213  |  99  |  0|1  |  46896303_C_G |  3367,0,6381
        -----------------------------------------------------------------------

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf with split out FORMAT fields
        """
        # get unique list of FORMAT fields from all rows
        fields = list(set(':'.join(vcf_df.FORMAT.tolist()).split(':')))

        # join respective pairs from FORMAT and SAMPLE columns to list of
        # '=' joined pairs
        # before -> GT:AD:DP:GQ:PL  1/1:0,41:41:99:1442,123,0
        # after  -> [GT=1/1, AD=0,41, DP=41, GQ=99, PL=1442,123,0]
        tmp_df = pd.DataFrame()
        tmp_df['tmp'] = vcf_df.apply(lambda x: [
            "=".join(x) for x in zip(x.FORMAT.split(':'), x.SAMPLE.split(':'))
        ], axis=1)

        # join everything with a '|', this allows us to split and cast the
        # pairs to dict to turn into a list of values => pd.DataFrame
        # can't use ',' since some SAMPLE values may contain them
        tmp_df['tmp'] = tmp_df['tmp'].apply(lambda x: '|'.join(x))

        # split every row to a dict of key value pairs and build df
        format_cols = tmp_df.join(pd.DataFrame([
            dict(val) for val in tmp_df.pop('tmp').str.findall((
                r'(\w+)=([^|\}]+)'
            ))
        ]))

        # add split out columns to main df
        vcf_df = vcf_df.join(format_cols, rsuffix=" (FMT)")
        vcf_df.drop(['FORMAT', 'SAMPLE'], axis=1, inplace=True)

        return vcf_df


    @staticmethod
    def info(vcf_df) -> pd.DataFrame:
        """
        Splits out the INFO column of vcf to all separate values, excluding
        the CSQ values. This transforms as:

        -----------------------------------------------------------------------
        INFO
        -----------------------------------------------------------------------
        AC=1;AF=0.5;AN=2;BaseQRankSum=-0.244;DB;DP=215;......
        AC=1;AF=0.5;AN=2;BaseQRankSum=-0.291;;DB;DP=207;.....
        AC=2;AF=1;AN=2;DB;DP=310;ExcessHet=3.0103;FS=0;......
        -----------------------------------------------------------------------
                                      |
                                      |
                                      ▼
        -----------------------------------------------------------------------
         AC  |  AF  |  AN  |  DB  |  DP  |  ExcessHet  |  BaseQRankSum  |  FS
        -----------------------------------------------------------------------
         1   |  0.5 |  2   | True |  215 |             |  -0.244        |
        -----------------------------------------------------------------------
         1   |  0.5 |  2   | True |  207 |             |  -0.291        |
        -----------------------------------------------------------------------
         2   |  1   |  2   | True |  310 |  3.0103     |                |  0
        -----------------------------------------------------------------------


        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf with separated INFO fields
        """
        # get the key value pairs of INFO data
        info_pairs = [
            x.split('CSQ=')[0].split(';') for x in vcf_df['INFO'].tolist()
        ]
        info_pairs = [[x for x in lst if x] for lst in info_pairs]

        # get unique list of keys from key=value pairs in INFO from all rows
        info_keys = sorted(list(set([
            x.split('=')[0] if '=' in x else x for pair in info_pairs for x in pair
        ])))
        info_keys = [x for x in info_keys if x]  # can end up with empty string

        info_values = []

        # info_pairs -> list of list of pairs, one list per variant
        for variant_pairs in info_pairs:
            # for every variants values, split them out to dict to add to df
            pair_values = {}

            for pair in variant_pairs:
                if '=' in pair:
                    # key value pair
                    key, value = pair.split('=')
                else:
                    # Flag value present (e.g STR)
                    key, value = pair, True

                pair_values[key] = value

            info_values.append(pair_values)

        # build df of values to add to main df
        info_df = pd.DataFrame(
            info_values, columns=info_keys
        )


        for col in info_keys:
            # add all info values to main vcf df
            vcf_df[col] = info_df[col]

        # drop INFO and CSQ as we fully split them out
        vcf_df.drop(['INFO'], axis=1, inplace=True)

        return vcf_df


    @staticmethod
    def csq(vcf_df, csq_fields) -> Union[pd.DataFrame, list]:
        """
        Split out CSQ string from other values in the INFO column to separate
        fields to get annotation.  Column headers taken from format stored by
        VEP in the header, read from vcf.parse_csq_fields().

        Variants with multiple transcript annotation will have duplicate CSQ
        data that is comma sepparated => expand this to multiple rows, if no
        ',' present rows will remain unaffacted (i.e. one transcript)

        Transforms as:

        -----------------------------------------------------------------------
        INFO
        -----------------------------------------------------------------------
        ...;CSQ=SIK1|SNV|missense_variant|13/14|NM_173354.5:c.1844C>T|...
        ...;CSQ=COL18A1|SNV|synonymous_variant|6/42|NM_0013750.1:c.846G>T...
        ...;CSQ=NSD1|SNV|missense_variant|6/24|NM_001384.1:c.1369T>C|...
        -----------------------------------------------------------------------
                                      |
                                      |
                                      ▼
        -----------------------------------------------------------------------
        SYMBOL | VAR_CLASS | Consequence        | EXON  | HGVSc
        -----------------------------------------------------------------------
        SIK1   | SNV       | missense_variant   | 13/14 | NM_173354.5:c.1844C>T
        -----------------------------------------------------------------------
        OL18A1 | SNV       | synonymous_variant | 6/42  | NM_0013750.1:c.846G>T
        -----------------------------------------------------------------------
        NSD1   | SNV       | missense_variant   | 6/24  | NM_001384.1:c.1369T>C
        -----------------------------------------------------------------------


        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf
        csq_fields : list
            list of CSQ field names parsed from vcf header, used to assign as
            column headings when splitting out CSQ data for each variant

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf with separated CSQ fields
        expanded_vcf_rows : int
            total number of extra rows from expanding multiple transcript
            annotations to individual rows
        """
        df_rows = len(vcf_df.index)
        columns = list(vcf_df.columns)

        # split out CSQ string from INFO col from each row
        vcf_df['CSQ'] = vcf_df['INFO'].apply(lambda x: x.split('CSQ=')[-1])

        # set index to be everything except CSQ, expand this to multiple rows
        # then set index back
        vcf_df = vcf_df.set_index(columns).apply(
            lambda x: x.str.split(',').explode()
        ).reset_index()

        # split each CSQ value to own columns
        vcf_df[csq_fields] = vcf_df.CSQ.str.split('|', expand=True)

        # drop INFO and CSQ as we fully split them out
        vcf_df.drop(['CSQ'], axis=1, inplace=True)

        if df_rows != len(vcf_df.index):
            # total rows has changed => we must have multiple transcripts
            print(f"Total rows of VCF changed on splitting CSQ values")
            print(f"Total rows before: {df_rows}")
            print(f"Total rows after: {len(vcf_df.index)}")
            expanded_vcf_rows = len(vcf_df.index) - df_rows
        else:
            expanded_vcf_rows = 0

        return vcf_df, expanded_vcf_rows

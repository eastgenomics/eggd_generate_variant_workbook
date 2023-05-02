import sys
from typing import Union
import pandas as pd


class splitColumns():
    """
    Functions for spliting and checking of variant annotation and
    attribute columns (FORMAT, SAMPLE, INFO, CSQ), called during reading
    of VCFs from file in vcf.read()
    """
    def split(self, vcf_df) -> Union[pd.DataFrame, int]:
        """
        Call methods to split out INFO and FORMAT columns of given dataframe

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of variants to split

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of variants
        expanded_vcf_rows : int
            total number of extra rows from expanding multiple transcript
            annotations to individual rows
        """
        vcf_df = self.info(vcf_df)
        vcf_df = self.format_fields(vcf_df)
        vcf_df = self.unique_vep(vcf_df)

        return vcf_df


    def unique_vep(self, vcf_df) -> pd.DataFrame:
        """
        Handle known bug in VEP annotation where it duplicates COSMIC IDs
        This creates a

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of variants

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of variants
        """

        # Find all columns that start with 'csq'
        csq_columns = [col for col in vcf_df.columns if col.lower().startswith('csq')]

        # Join the 'csq' columns using '&' and remove duplicates
        for col in csq_columns:
            vcf_df[col] = vcf_df[col].apply(
                lambda x: '&'.join(sorted(set(x.split('&')))) if isinstance(x, str) else x
            )

        return vcf_df


    def format_fields(self, vcf_df) -> pd.DataFrame:
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

        Raises
        ------
        AssertionError
            Raised when SAMPLE column contains '##' that is used to split out
            data to separte columns on

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf with split out FORMAT fields
        """
        # join respective pairs from FORMAT and SAMPLE columns to '##'
        # separated '=' joined pairs
        # before -> GT:AD:DP:GQ:PL  1/1:0,41:41:99:1442,123,0
        # after  -> GT=1/1##AD=0,41##DP=41##GQ=99##PL=1442,123,0
        # n.b. using '##' to join since VCF spec for what is allowed in
        # sample strings is not strict, but ## is very unlikely to be used
        tmp_df = pd.DataFrame()

        # sense check what we use to split on does not already exist
        assert vcf_df[vcf_df['SAMPLE'].str.contains('##')].empty, (
            'Error attempting to split SAMPLE column, one or more rows '
            '  contain "##" which is used to separate the columns.'
        )

        tmp_df['tmp'] = vcf_df.apply(lambda x: '##'.join([
            "=".join(x) for x in zip(x.FORMAT.split(':'), x.SAMPLE.split(':'))
        ]), axis=1)

        # split every row to a dict of key value pairs and build df
        format_cols = tmp_df.join(pd.DataFrame([
            dict(val) for val in tmp_df.pop('tmp').str.findall((
                r'(\w+)=([^##\}]+)'
            ))
        ]))

        # add split out columns to main df
        vcf_df = vcf_df.join(format_cols, rsuffix="_FMT")
        vcf_df.drop(['FORMAT', 'SAMPLE'], axis=1, inplace=True)

        return vcf_df


    def info(self, vcf_df) -> pd.DataFrame:
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
        info_df = pd.DataFrame(info_values, columns=info_keys)
        vcf_df = pd.concat([vcf_df, info_df], axis=1)

        # drop INFO and CSQ as we fully split them out
        vcf_df.drop(['INFO'], axis=1, inplace=True)


        return vcf_df

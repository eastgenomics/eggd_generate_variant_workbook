import re

import numpy as np
import pandas as pd


class filter():
    """"
    Functions related to filtering dataframes of variants from args.filter

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    vcfs : list of pd.DataFrame
        list of dataframes read in from self.args.vcfs
    filtered_rows : pd.DataFrame
        dataframe of all rows dropped from all vcfs
    """
    def __init__(self, args, vcfs):
        self.args = args
        self.vcfs = vcfs
        self.filtered_rows = pd.DataFrame()


    def validate(self) -> None:
        """
        Validate filters passed for filtering variants down.

        These come from args.filter and wil be in the format
        <column><operator><value> (i.e. "Consequence!=synonymous"), currently
        supports the operators >, <, >=, <=, == and !=

        Method will check that a valid operator has been passed, and that the
        given column is present in all columns of all vcfs
        """
        for filter in self.args.filter:
            # check a valid operand passed
            assert len(re.findall('>|<|>=|<=|==|!=', filter)) == 1, (
                f"invalid operand passed in filter: {filter}"
            )

        for vcf in self.vcfs:
            # for each filter, check the specified column is in all the vcfs
            for filter in self.args.filter:
                assert re.split(r'>|<|>=|<=|==|!=', filter)[0] in vcf.columns, (
                    f"Column specified in filter '{filter}' not in vcf "
                    f"columns: {vcf.columns}"
                )


    def build(self) -> None:
        """
        Formats cmd line passed filters as list of lists of each filter
        expression for filtering df of variants.

        This will separate out the passed filters from the format
        "Consequence!=synonymous" -> ["Consequence", "!=", "synonymous"]
        """
        field_value = [
            re.split(r'>|<|>=|<=|==|!=', x) for x in self.args.filter
        ]
        operator = [
            re.findall('>|<|>=|<=|==|!=', x) for x in self.args.filter
        ]

        self.filters = [
            [x[0], y[0], x[1]] for x, y in zip(field_value, operator)
        ]


    def filter(self) -> None:
        """
        Apply filters passed to ech dataframe of variants.

        Filters first checked in self.validate_filters then formatted in
        self.build filters ready to be passed to np.where()

        Currently wrapped in an eval() call which is not ideal but works for
        interpreting the operator passed as a string from the cmd line
        """
        # build list of indices of variants to filter out against specified
        # filters, then apply filter to df, retain filtered rows if --keep set
        for idx, vcf in enumerate(self.vcfs):
            all_filter_idxs = []
            for filter in self.filters:
                column, operator, value = filter[0], filter[1], filter[2]
                if pd.api.types.is_numeric_dtype(vcf[column]):
                    # check column we're filtering is numeric and set types
                    value = float(value)
                else:
                    # string values have to be wrapped in quotes from np.where
                    value = f"'{value}'"

                # get row indices to filter
                filter_idxs = eval((
                    f"np.where(vcf['{column}'].apply("
                    f"lambda x: x {operator} {value}))[0]"
                ))
                all_filter_idxs.extend(filter_idxs)

            # get unique list of indexes matching filters
            all_filter_idxs = sorted(list(set(all_filter_idxs)))

            if not self.args.always_keep.empty:
                # a list of variants / positions specified to never filter,
                # check for these in our filter list and remove if any present
                all_filter_idxs = self.retain_variants(vcf, all_filter_idxs)

            # apply the filter, assign back to the filtered df
            self.filtered_rows = self.filtered_rows.append(
                    vcf.loc[all_filter_idxs], ignore_index=True
                )

            # drop from current vcf dataframe
            self.vcfs[idx] = vcf.drop(all_filter_idxs)

            filter_string = ', '.join([' '.join(x) for x in self.filters])

            print((
                f"\nApplied the following filters: {filter_string} to vcf(s)\n"
                f"Filtered out {len(all_filter_idxs)} rows\n"
                f"Total rows remaining: {len(self.vcfs[idx].index)}\n\n"
            ))

        self.filtered_rows = self.filtered_rows.reset_index()

        if self.args.keep:
            # keeping filtered variants to write to file
            self.vcfs.append(self.filtered_rows)
            self.args.sheets.append("filtered")


    def retain_variants(self, vcf_df, filter_idxs) -> list:
        """
        Given a vcf and list of indices, check if any to filter are in the
        given list of variants to never filter and remove if so.

        Regions to not apply filters to are taken from tsv file passed with
        self.args.always_keep, and coordinates in file are inclusive.

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of variants from VCF
        filter_idxs : list
            list of df indices selected to filter out

        Returns:
        filter_idxs : list
            list of df indices selected to filter out
        """
        retain_idxs = []

        for _, row in self.args.always_keep.iterrows():
            # get indices of df of any variants to always retain, then
            # drop these from the list of filter indices
            retain_idxs.extend(
                np.where((
                        vcf_df['CHROM'] == row['chrom']
                    ) & (
                        vcf_df['POS'] >= row['start']
                    ) & (
                        vcf_df['POS'] <= row['end']
                ))[0]
            )

        filter_idxs = list(set(filter_idxs) - set(retain_idxs))

        return filter_idxs


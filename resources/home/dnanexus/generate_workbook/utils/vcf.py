from cmath import nan
import gzip
import os
from pathlib import Path, PurePath
import re
import shutil
import subprocess
import sys
from typing import Union
import urllib.parse

import pandas as pd
import numpy as np

from .columns import splitColumns
from .filters import filter
from .utils import buildHyperlink, is_numeric, determine_delimiter, \
    parse_cvo, parse_metrics_output


class vcf():
    """
    Functions to handle reading and manipulating vcf data

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    vcfs : list
        list of dataframe(s) of vcf data read in and formatted
    additional_files : dict
        dict of dataframes from (optionaly) passed additional files of
        tabulated data to write to additional sheets in workbook
    refs : list
        list of genome reference files used for given VCFs
    filtered_rows : list
        list of dataframes of rows filtered out from vcfs
    urls : dict
        mapping dictionary of column name to URLs, used for adding hyperlinks
            to column values before writing to file
    """

    def __init__(self, args) -> None:
        self.args = args
        self.vcfs = []
        self.additional_files = {}
        self.refs = []
        self.filtered_vcfs = []


    def process(self) -> None:
        """
        Function to call all methods in vcf() for processing given VCFs and
        formatting ready to write to output file

        Calls methods in following order:

            - self.filter() (optional with --filter)
            - self.read()
            - splitColumns.split()
            - self.merge()
            - self.drop_columns()
            - self.reorder()
            - self.format_strings()
            - self.add_hyperlinks()
            - self.rename_columns()
        """
        filters = filter(self.args)

        if self.args.additional_files:
            # additional non VCF files given, try read these in to dataframe(s)
            self.read_additional_files()

        # read in the each vcf, optionally filter, and then apply formatting
        for vcf in self.args.vcfs:
            # names for intermediary vcfs
            vcf_stem = Path(vcf).stem.replace('.vcf', '')
            split_vcf = f"{vcf_stem}.split.vcf"
            split_vcf_gz = f"{vcf_stem}.split.vcf.gz"
            filter_vcf = f"{vcf_stem}.filter.vcf"
            filter_vcf_gz = f"{vcf_stem}.filter.vcf.gz"

            # if VCF annotated with VEP and has not been split with bcftools,
            # first split multiple transcript annotation to separate VCF
            # records, and separate CSQ fields to separate INFO fields
            if self.check_vep_vcf(vcf, split_vcf):
                self.bcftools_pre_process(vcf, split_vcf)
                self.bgzip(split_vcf)

            if self.args.filter:
                # filter vcf against specified filters using bcftools
                _, columns = self.parse_header(vcf)
                filters.filter(split_vcf, filter_vcf, columns)
                self.bgzip(filter_vcf)

                # filters.filter() writes temp vcf with modified FILTER column
                variant_df = self.read(filter_vcf, Path(vcf).stem)

                # get filtered rows and read back to new dfs
                keep_df, filtered_df = filters.split_include_exclude(variant_df)

                # check we haven't dropped any variants
                filters.verify_total_variants(split_vcf_gz, keep_df, filtered_df)

                # split out INFO and FORMAT column values to individual
                # columns in dataframe, need to handle cases where one or both
                # may be empty to not raise errors downstream
                if not keep_df.empty and not filtered_df.empty:
                    # both have variants
                    keep_df = splitColumns().split(keep_df)
                    filtered_df = splitColumns().split(filtered_df)
                elif keep_df.empty and not filtered_df.empty:
                    # everything excluded, make empty keep df with same columns
                    # as those in excluded/filtered df
                    filtered_df = splitColumns().split(filtered_df)
                    keep_df = filtered_df.copy().drop(filtered_df.index)
                elif not keep_df.empty and filtered_df.empty:
                    # nothing filtered out, make empty filtered df with same
                    # columns as included variants
                    keep_df = splitColumns().split(keep_df)
                    filtered_df = keep_df.copy().drop(keep_df.index)
                else:
                    # both empty, we can't magic up columns names so just
                    # continue and workbook will have standard VCF columns
                    pass

                self.vcfs.append(keep_df)
                self.filtered_vcfs.append(filtered_df)

                if not self.args.keep_tmp:
                    os.remove(filter_vcf_gz)

                # clean up some memory since big dataframes can use a lot
                del keep_df, filtered_df
            else:
                # not filtering vcf, read in full vcf and split out INFO and
                # FORMAT/SAMPLE column values to individual columns in df
                vcf_df = self.read(split_vcf, Path(vcf).stem)
                if not vcf_df.empty:
                    vcf_df = splitColumns().split(vcf_df)
                self.vcfs.append(vcf_df)

                del vcf_df

            # delete tmp vcf from splitting CSQ str in bcftools_pre_process()
            os.remove(split_vcf)

            if self.args.filter:
                os.remove(filter_vcf)

            if not self.args.keep_tmp:
                os.remove(split_vcf_gz)

        if self.args.merge:
            self.vcfs = self.merge(self.vcfs)

        if self.args.filter and self.args.keep:
            # merge all filtered dataframes to one and add to list of vcfs for
            # doing column operations and writing to Excel file
            if not all(x.empty for x in self.filtered_vcfs):
                self.filtered_vcfs = self.merge(self.filtered_vcfs)
            self.vcfs.append(self.filtered_vcfs[0])
            self.args.sheets.append('excluded')

        if self.args.summary == 'dias':
            # if it is dias pipeline, add the empty col
            # named Interpreted in the first variant sheet
            self.vcfs[0]['Interpreted'] = ''

        if self.args.split_hgvs:
            self.split_hgvs(self.vcfs)

        if self.args.add_raw_change:
            self.add_raw_change()

        if self.args.print_columns:
            self.print_columns()

        if self.args.additional_columns:
            self.add_additional_columns()

        if self.args.af_format == "percent":
            self.percent_af(self.vcfs)

        if self.args.join_columns:
            self.joining_columns(self.vcfs)

        if self.args.report_text:
            self.make_report_text(self.vcfs)

        self.vcfs = self.format_strings(self.vcfs)
        self.vcfs = self.add_hyperlinks(self.vcfs)

        if self.args.exclude or self.args.include:
            self.drop_columns()

        if self.args.reorder:
            self.order_columns(self.vcfs)

        self.vcfs = self.rename_columns(self.vcfs)

        print("\nSUCCESS: Finished munging variants from vcf(s)\n")


    def bcftools_pre_process(self, vcf, output_vcf) -> None:
        """
        Decompose multiple transcript annotation to individual records, and
        split VEP CSQ string to individual INFO keys. Adds a 'CSQ_' prefix
        to each field extracted from the CSQ string to stop potential conflicts
        with existing INFO fields, which is then stripped before writing
        to the Excel file

        Parameters
        ------
        vcf : str
            path to vcf file to use
        output_vcf : str
            name for output vcf

        Outputs
        -------
        {vcf}.split.vcf : file
            vcf file output from bcftools

        Raises
        ------
        AssertionError
            Raised when non-zero exit code returned by bcftools
        """
        print(f"Splitting {vcf} with bcftools +split-vep")

        # check total rows before splitting
        pre_split_total = subprocess.run(
            f"zgrep -v '^#' {vcf} | wc -l", shell=True,
            capture_output=True, check=True
        )

        cmd = (
            f"bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_' -d {vcf} | "
            f"bcftools annotate -x INFO/CSQ -o {output_vcf}"
        )

        output = subprocess.run(cmd, shell=True, capture_output=True, check=False)

        assert output.returncode == 0, (
            f"\n\tError in splitting VCF with bcftools +split-vep. VCF: {vcf}"
            f"\n\tExitcode:{output.returncode}"
            f"\n\t{output.stderr.decode()}"
        )

        # check total rows after splitting
        post_split_total = subprocess.run(
            f"zgrep -v '^#' {output_vcf} | wc -l", shell=True,
            capture_output=True, check=True
        )

        print(
            f"Total lines before splitting: {pre_split_total.stdout.decode()}"
            f"Total lines after splitting: {post_split_total.stdout.decode()}"
        )


    def bgzip(self, file) -> None:
        """
        Call bgzip on given file

        Parameters
        ----------
        file : file to compress

        Outputs
        -------
        input file, but compressed

        Raises
        ------
        AssertionError
            Raised when non-zero exit code returned by bgzip
        """
        output = subprocess.run(
            f"bgzip --force {file} -c > {file}.gz",
            shell=True, capture_output=True, check=False
        )

        assert output.returncode == 0, (
            f"\n\tError in compressing file with bgzip. File: {file}"
            f"\n\tExitcode:{output.returncode}"
            f"\n\t{output.stderr.decode()}"
        )


    def check_vep_vcf(self, vcf, split_vcf) -> bool:
        """
        Checks if VCF is annotated with VEP and has not been previously
        split with bcftools +split-vep.

        This is required to control if to call self.bcftools_pre_process(),
        if it is not annotated or has been previously split, then we
        create a tmp vcf with the same name as when split to use downstream.

        Parameters
        ----------
        vcf : str
            filename of vcf to check
        split_vcf : str
            filename of tmp vcf to create if not annotated / already split

        Returns
        -------
        bool
            True if annotated and not already split, False if not

        Outputs
        -------
        file
            vcf file with name from split_vcf input
        file
            bgzipped equivalent of above
        """
        header, _ = self.parse_header(vcf)

        if not any([
            x.startswith('##VEP') for x in header
        ]) or any([
            x.startswith('##bcftools_split-vep') for x in header
        ]):
            # VCF not annotated with VEP or annotated with VEP but
            # already split => create tmp vcf to continue processing
            print('VCF not annotated with VEP or already split, continuing...')
            if vcf.endswith('.gz'):
                input_file = gzip.open(vcf, 'rb')
            else:
                input_file = open(vcf, 'rb')

            with open(split_vcf, 'wb') as output_file:
                shutil.copyfileobj(input_file, output_file)
                input_file.close()

            self.bgzip(split_vcf)

            return False
        else:
            return True


    def read(self, vcf, sample=None) -> pd.DataFrame:
        """
        Reads given vcf into pd.DataFrame and parses header

        Parameters
        ------
        vcf : str
            path to vcf file to use
        sample : str
            name of vcf, used for adding name to df if --add_name passed

        Returns
        -------
        vcf_df : pandas.DataFrame
            dataframe of all variants
        """
        print(f"\n\nReading in vcf {vcf} for sample {sample}\n")

        if sample:
            sample = sample.replace('.vcf', '').replace('.gz', '')
            if '_' in sample:
                sample = sample.split('_')[0]

        header, columns = self.parse_header(vcf)
        self.parse_reference(header)

        if self.args.print_header:
            self.print_header(header)

        # read vcf into pandas df
        vcf_df = pd.read_csv(
            vcf, sep='\t', comment='#', names=columns, compression='infer'
        )

        if self.args.add_name:
            # add sample name from filename as 1st column
            vcf_df.insert(loc=0, column='sampleName', value=sample)

        if self.args.add_comment_column:
            # add empty 'Comment' column to end of df
            vcf_df['Comment'] = ''

        if self.args.add_classification_column:
            # add empty 'Classification' column to end of df
            vcf_df['Classification'] = ''

        return vcf_df


    def read_additional_files(self):
        """
        Attempt to read in additional tabulated files to dataframes for
        writing as additional sheets to the output workbook

        Updates self.additional_files dictionary with file_prefix : dataframe
        """
        for idx, file in enumerate(self.args.additional_files):
            # get prefix from filename for naming sheet if not specified
            if self.args.additional_sheets:
                prefix = self.args.additional_sheets[idx]
            else:
                prefix = PurePath(file).name.replace(
                    ''.join(PurePath(file).suffixes), ''
                )

            # Excel has a limit of 31 characters for sheet name -> trim
            if len(prefix) > 31:
                prefix = prefix[:31]
                print(
                    f"Prefix of additional file {file} is >31 character "
                    "limit for an Excel worksheet. Name will be trimmed to "
                    f"maximum length: {prefix}"
                )

            # read file contents in to list
            if file.endswith('.gz'):
                with gzip.open(file) as fh:
                    file_contents = [
                        x.decode() for x in fh.read().splitlines()
                    ]
            else:
                with open(file) as fh:
                    file_contents = fh.read().splitlines()

            if file.endswith('vcf') or file.endswith('vcf.gz'):
                # vcf passed => process and format nicer for displaying
                split_additional_vcf=file.replace('.vcf', '_split.vcf')

                if self.check_vep_vcf(file, split_additional_vcf):
                    self.bcftools_pre_process(
                        vcf=file,
                        output_vcf=split_additional_vcf
                    )
                    file_df = self.read(split_additional_vcf, Path(split_additional_vcf).stem)
                else:
                    file_df = self.read(file)

                # call some of the formatting methods for regular vcfs
                # to get things like split INFO columns and hyperlinks
                if not file_df.empty:
                    file_df = splitColumns().split(file_df)

                if self.args.split_hgvs:
                    file_df = self.split_hgvs([file_df])[0]

                if self.args.af_format == "percent":
                    file_df = self.percent_af([file_df])[0]
                if self.args.report_text:
                    file_df = self.make_report_text([file_df])[0]

                file_df = self.format_strings([file_df])[0]
                file_df = self.add_hyperlinks([file_df])[0]
                if self.args.join_columns:
                    file_df = self.joining_columns([file_df])[0]

                if self.args.exclude or self.args.include:
                    self.drop_columns([file_df])

                if self.args.reorder:
                    file_df = self.order_columns([file_df])[0]

                file_df = self.rename_columns([file_df])[0]
                # force header to also be first line of df so it is written
                # to the Excel sheet
                file_df = pd.DataFrame(
                    [file_df.columns], columns=file_df.columns).append(file_df)
            else:
                # check what delimiter the data uses
                # check end of file to avoid potential headers causing issues
                delimiter = determine_delimiter(
                    '\n'.join(file_contents[-5:]), PurePath(file).suffixes
                )

                file_df = pd.DataFrame(
                    [line.split(delimiter) for line in file_contents]
                )

            if file.endswith('_CombinedVariantOutput.tsv'):
                # file passed is a CombinedVariantOutput file from Illumina
                # TSO500 app, just parse out TMB, MSI and Amplifications
                print(
                    'CombinedVariantOutput file passed to --additional_files, '
                    'parsing out TMB, MSI and Amplifications'
                )
                file_df = parse_cvo(cvo_df=file_df)

            if file.endswith('MetricsOutput.tsv'):
                # file passed is run level MetricsOutput.tsv from Illumina
                # TSO500 app, attempt to parse out just sample metrics to display
                print(
                    'TSO500 MetricsOutput passed to --additional_files, '
                    'attempting to parse sample metrics from file'
                )
                file_df = parse_metrics_output(
                    metrics_df=file_df,
                    sample_vcf=Path(self.args.vcfs[0]).name
                )

            self.additional_files[prefix] = file_df


    def parse_header(self, vcf) -> Union[list, list]:
        """
        Read in header lines of given vcf to list, returning the list and the
        vcf column names

        Parameters
        ----------
        vcf : str
            vcf filename to read header from

        Returns
        -------
        header : list
            list of header lines read from vcf
        columns : list
            column names from vcf

        Raises
        ------
        AssertionError
            Raised when header looks to be malformed and column names incorrect
        """
        if vcf.endswith('.gz'):
            fh = gzip.open(vcf)
        else:
            fh = open(vcf)

        # read in header of vcf
        header = []
        for line in fh.readlines():
            if vcf.endswith('.gz'):
                line = line.decode()

            if line.startswith('#'):
                header.append(line.rstrip('\n'))
            else:
                break

        fh.close()

        columns = [x.strip('#') for x in header[-1].split()]
        columns[-1] = 'SAMPLE'

        assert columns[0] == 'CHROM', (
            "Parsed header appears to be malformed, column names parsed as: "
            f"{columns}"
        )

        return header, columns


    def parse_reference(self, header) -> None:
        """
        Parse reference file used from VCF header

        Parameters
        ----------
        header : list
            lines of vcf header
        """
        ref = ''

        # first check if we can get reference build from VEP command string
        vep = [x for x in header if x.startswith('##VEP')]
        if vep:
            assembly = re.search(r'assembly="[\w\d\.]+"', vep[0])
            if assembly:
                ref = assembly.group(0).replace('assembly=', '').strip('"\'')

        if not ref:
            ref = next(
                iter([x for x in header if x.startswith('##reference')]), None
            )

        if ref:
            if ref not in self.refs:
                # add reference file if found and same not already in list
                self.refs.append(Path(ref).name)

                # check we don't have a mix of 37 and 38
                assert not ('37' in str(self.refs) and '38' in str(self.refs)), (
                    'References from vcfs appear to be a mix of reference '
                    f'builds.\n References parsed: {self.refs}'
                )


    def add_hyperlinks(self, vcfs) -> list:
        """
        Format column value as an Excel hyperlink if URL for column specified

        Parameters
        ----------
        vcfs : list
            list of pd.DataFrames of vcfs to add hyperlinks to
        Returns
        -------
        list
            list of dataframes with added hyperlinks
        """
        # some URLs are build specific, infer which to use from build in header
        build = None
        reference = ''

        if self.refs:
            reference = self.refs[0].lower()
            if '37' in reference or 'hg19' in reference:
                build = 37
            elif '38' in reference:
                build = 38

        for idx, vcf in enumerate(vcfs):
            if vcf.empty:
                # empty dataframe => nothing to add links to
                continue
            for column in vcf.columns:
                vcfs[idx][column] = vcfs[idx].apply(
                    lambda x: buildHyperlink().build(
                        column=column,
                        value=x,
                        build=build
                    ), axis=1
                )

        return vcfs


    def format_strings(self, vcfs) -> list:
        """
        Fix formatting of string values with different encoding and nans

        Parameters
        ----------
        vcfs : list
            list of pd.DataFrames of vcfs to fix strings for
        Returns
        -------
        list
            list of dataframes with fixed strings
        """
        for idx, vcf in enumerate(vcfs):
            # pass through urllib unqoute and UTF-8 to fix any weird symbols
            vcf = vcf.applymap(
                lambda x: urllib.parse.unquote(x).encode('UTF-8').decode()
                if type(x) == str else x
            )

            # remove any nans that are strings
            vcf = vcf.applymap(
                lambda x: x.replace('nan', '')
                if x == 'nan' and type(x) == str else x
            )

            vcfs[idx] = vcf

        return vcfs


    def print_columns(self) -> None:
        """
        Simple method to just print the columns from each vcf and exit.

        Useful for identify what columns are present in INFO and CSQ fields
        for using --include, --exclude and --reorder arguments
        """
        for name, vcf in zip(self.args.vcfs, self.vcfs):
            print(f"Columns for {Path(name).name}: ")
            print(f"\n\t{list(vcf.columns)}\n\n")

        sys.exit(0)


    def print_header(self, header) -> None:
        """
        Simple method to print vcf header(s) after splitting CSQ string with
        bcftools to show all available fields and their types

        Parameters
        ----------
        header : list
            vcf header as list read from file
        """
        [print(x) for x in header]
        sys.exit(0)


    def drop_columns(self, vcfs=None) -> None:
        """
        If `--exclude` or `--include` passed, drop given columns
        (or inverse of) from vcf data if they exist.

        If `--include` passed will take the given list of columns and drop the
        remaining columns not specified from all dataframes

        If `--exclude` passed will take the given list of columns and drop
        from all dataframes

        Parameters
        ----------
        vcfs : pandas.DataFrame
            pandas dataframe containing all rows and column from VCF

        Raises
        ------
        AssertionError
            Raised when columns specified with --include / --exclude are not
            present in one or more of the dataframes
        """

        if vcfs is None:
            vcfs = self.vcfs

        for idx, vcf in enumerate(vcfs):
            if self.args.include:
                # include passed => select all columns not specified to drop
                columns = self.args.include
                to_drop = list(
                    set(vcf.columns.tolist()) - set(columns)
                )
            elif self.args.exclude:
                columns = self.args.exclude
                to_drop = self.args.exclude
            else:
                continue

            # get any columns passed that aren't present in vcf
            invalid = list(set(columns) - set(vcf.columns))

            if invalid:
                print(
                    f"WARNING: Columns passed with `--include / --exclude not "
                    f"present in vcf ({invalid}), skipping these columns..."
                )
                if self.args.exclude:
                    # only need to remove in the case of excluding since
                    # include has already selected valid columns
                    for col in invalid:
                        to_drop.remove(col)

            if self.args.report_text:
                to_drop.remove("Report_text")

            vcfs[idx].drop(to_drop, axis=1, inplace=True, errors='ignore')


    def order_columns(self, vcfs) -> list:
        """
        Reorder columns by specified order from `--reorder` argument, any not
        specified will retain original order after reorder columns

        Parameters
        ----------
        vcfs : list
            list of pd.DataFrames of vcfs to order columns from
        Returns
        -------
        list
            list of dataframes with ordered columns

        Raises
        ------
        AssertionError
            Raised when columns specified with --reorder are not
            present in one or more of the dataframes
        """
        for idx, vcf in enumerate(vcfs):
            vcf_columns = list(vcf.columns)

            # check columns given are present in vcf
            invalid = list(
                set(self.args.reorder) - set(vcf_columns) -
                set(self.args.additional_columns)
            )
            if invalid:
                print(
                    f"WARNING: columns passed to --reorder not present in vcf:"
                    f" {invalid}. Skipping these columns and continuing..."
                )
                for col in invalid:
                    self.args.reorder.remove(col)

            [vcf_columns.remove(x) for x in self.args.reorder]
            column_order = self.args.reorder + vcf_columns

            vcfs[idx] = vcf[column_order]

        return vcfs


    def add_additional_columns(self) -> None:
        """
        Append empty columns specified from --additional_columns for adding
        additional hyperlinks to external resources (e.g. decipher, oncoKB etc.)
        """
        for column in self.args.additional_columns:
            if column.lower() in ['decipher']:
                # column is only for b38, check if vcf also is
                if not self.refs:
                    print(
                        f'WARNING: {column} specified to --additional_columns '
                        'but no reference could be parsed from vcf header. '
                        f'Continuing without adding {column} column.'
                    )
                    continue

                if 'hg19' in self.refs[0] or '37' in self.refs[0]:
                    print(
                        f'WARNING: {column} specified to `--additional_columns '
                        'but VCF appears to be for b37. Continuing without '
                        f'adding {column} column.'
                    )
                    continue

            for idx, vcf in enumerate(self.vcfs):
                vcf[column] = column
                self.vcfs[idx] = vcf


    def rename_columns(self, vcfs) -> list:
        """
        Rename columns from key value pairs passed from --rename argument,
        also remove underscores from all names for nicer reading

        Parameters
        ----------
        vcfs : list
            list of pd.DataFrames of vcfs to rename columns from
        Returns
        -------
        list
            list of dataframes with renamed columns

        Raises
        ------
        AssertionError
            Raised when columns specified with --rename do not exist in more
            or more of the vcfs columns

        AssertionError
            Raised when new column names specified are already present in the
            vcf
        """
        for idx, vcf in enumerate(vcfs):
            if self.args.rename:
                # check the given new name(s) not already a column name
                assert all(
                    x not in vcf.columns for x in self.args.rename.values()
                ), (
                    f"Column(s) specified with --rename already present in "
                    f"one or more of the given vcfs. \n\ Column names: "
                    f"\n\n\t{vcf.columns}. \n\nNew column names passed to "
                    f"--rename: \n\n\t{list(self.args.rename.values())}"
                )

                # check specified columns are present in vcf, if not print
                # warning, remove and continue
                new_names_dict = self.args.rename.copy()

                invalid = list(
                    set(new_names_dict.keys()) - set(vcf.columns.tolist())
                )

                if invalid:
                    print(
                        f"WARNING: columns passed to --rename not present in vcf:"
                        f" {invalid}. Skipping these columns and continuing..."
                    )
                    for key in invalid:
                        new_names_dict.pop(key)

                vcfs[idx].rename(
                    columns=dict(new_names_dict.items()), inplace=True
                )

            # strip prefix from column name if present and not already a column
            vcfs[idx].columns = self.strip_csq_prefix(vcfs[idx])

            # remove underscores from all column names
            vcfs[idx].columns = [
                x.replace('_', ' ') for x in vcfs[idx].columns
            ]

        return vcfs


    def strip_csq_prefix(self, vcf) -> list:
        """
        Strip CSQ prefix added by bcftools -split-vep from column names

        Any conflicts in names with already present columns will retain prefix

        Parameters
        ----------
        vcf : pd.DataFrame
            dataframe to modify column names of

        Returns
        -------
        list
            list of column names with CSQ_ prefixes removed
        """
        return [
            x.replace('CSQ_', '', 1) if (
                x.startswith('CSQ_') and x.replace('CSQ_', '') not in vcf.columns
            ) else x for x in vcf.columns
        ]


    def merge(self, vcfs) -> None:
        """
        Merge all variants into one big dataframe, should be used with
        --add_name argument if provenance of variants in merged dataframe
        is important
        """
        # don't attempt to merge empty vcfs as likely to have diff. columns
        vcfs = [x for x in vcfs if not x.empty]

        return [pd.concat(vcfs).reset_index(drop=True)]


    def split_hgvs(self, vcfs) -> pd.DataFrame:
        """
        If --split_hgvs specified, attempt to split HGVSc and HGVSp columns
        into 2 separate ones: c. change (DNA) and p. change (Protein).
        """
        for idx, vcf in enumerate(vcfs):
            # check required columns are in the dataframe
            if not all(col in vcf.columns for col in ['CSQ_HGVSc', 'CSQ_HGVSp']):
                print(
                    'WARNING: --split_hgvs specified but CSQ_HGVSc and/or '
                    'CSQHGVSp not present in VCF fields. Continuing without '
                    'splitting HGVS.'
                )
                continue

            # ensure columns we're going to create don't already exist
            if any(col in vcf.columns for col in ['DNA', 'Protein']):
                print(
                    'WARNING: --split_hgvs specified but DNA and/or Protein '
                    'already exist in the vcf. Continuing without splitting HGVS.'
                )
                continue

            vcf['DNA'] = vcf['CSQ_HGVSc'].str.split(':').str[1]
            vcf['Protein'] = vcf['CSQ_HGVSp'].str.split(':').str[1]

            vcfs[idx] = vcf

        return vcfs


    def add_raw_change(self) -> None:
        """
        Adds a column named 'rawChange' of the 'raw' genomic change, formatted
        as {CHROM}:g.{POS}{REF}>{ALT}. This is analogous to HGVSg output by
        VEP, but will not include text such as INV and DEL.
        """
        for idx, vcf in enumerate(self.vcfs):
            if vcf.empty:
                self.vcfs[idx]['rawChange'] = ''
                continue

            if not all(
                col in vcf.columns for col in ['CHROM', 'POS', 'REF', 'ALT']
            ):
                # one or more required columns missing => skip
                print(
                    "WARNING: one or more required columns missing for "
                    "add_raw_change, continuing without adding"
                )
                continue

            self.vcfs[idx]['rawChange'] = vcf.agg(
                '{0[CHROM]}:g.{0[POS]}{0[REF]}>{0[ALT]}'.format, axis=1)


    def percent_af(self, vcfs) -> list:
        """
        Finds the column with "AF" and will convert the number format
        to percent
        Parameters
        ----------
        vcfs : list
            list of pd.DataFrames of vcfs to change the AF columns type
        Returns
        -------
        list
            list of dataframes with AF changed to percent type
        """
        # find the sheets and apply to all sheets
        for idx, vcf in enumerate(vcfs):
            if 'AF' not in vcf.columns:
                continue
            vcf['AF'] = vcf['AF'].astype(np.float16)
            vcf['AF'] = vcf['AF'].map(lambda n: '{:,.1%}'.format(n))

            vcfs[idx] = vcf

        return vcfs


    def make_report_text(self, vcfs):
        """
        Makes a report text that follows the has the details per row
        gene_symbol consequence, hgvsc, hgvsp, cosmic, dbsnp and
        allele frequency

        Parameters
        ----------
        vcfs : list
            list of vcf dataframes to which to add report text column to

        Returns
        -------
        list
            list of vcf dataframes with added column
        """
        for idx, vcf in enumerate(vcfs):
            vcf['Report_text'] = vcf.apply(self.format_report_text, axis=1)

            vcfs[idx] = vcf

        return vcfs


    @staticmethod
    def format_report_text(row) -> str:
        """
        Format the report text for a given row

        Parameters
        ----------
        row : pd.Series
            Series of variant row

        Returns
        -------
        str
            Report text formatted as a single string
        """
        def add_none(value):
            """
            Replaces absent VCF values ('.') as 'None'
            """
            return value if value != '.' else 'None'

        # force the field names to be lower case to handle differences
        # in case and remove CSQ prefix
        row.index = [
            re.sub(r'^csq_', '', x.lower()) for x in row.index.tolist()
        ]
        text = ""
        if row.get('symbol') and row.get('consequence'):
            text += f"{row.get('symbol', '')} {row.get('consequence')} "

        if row.get('exon', '').replace('.', ''):
            text += f"in exon {str(row.get('exon', '')).split('/')[0]}\n"

        if row.get('intron', '').replace('.', ''):
            text += f"in intron {str(row.get('intron', '')).split('/')[0]}\n"

        if row.get('dna'):
            text += f"HGVSc: {add_none(row.get('dna', ''))}\n"
        elif row.get('hgvsc'):
            text += f"HGVSc: {add_none(row.get('hgvsc', ''))}\n"

        if row.get('protein'):
            text += f"HGVSp: {add_none(row.get('protein', ''))}\n"
        elif row.get('hgvsp'):
            text += f"HGVSp: {add_none(row.get('hgvsp', ''))}\n"

        if row.get('cosmiccmuts', '').replace('.', ''):
            text += f"COSMIC coding ID: {row.get('cosmiccmuts')}\n"

        if row.get('cosmicncmuts', '').replace('.', ''):
            text += f"COSMIC non-coding ID: {row.get('cosmicncmuts')}\n"

        if row.get('cosmic', '').replace('.', ''):
            text += f"COSMIC ID: {row.get('cosmic')}\n"

        if row.get('existing_variation', '').replace('.', ''):
            text += f"dbSNP: {row.get('existing_variation', '')}\n"

        text += f"Allele Frequency (VAF): {add_none(str(row.get('af', '')))}"

        return text


    def joining_columns(self, vcfs) -> list:
        """
        Joins two exisiting columns into a new column. The user will provide the
        columns with a choice of seperator

        Parameters
        ----------
        vcfs : list
            list of vcf dataframes to which to add report text column to

        Returns
        -------
        list
            list of vcf dataframes with new column containing info of two columns
        """
        list_of_joins=self.args.join_columns

        # if the list has multiple joins, seperate the string by space
        #### check that the lists are seperated by columns
        list_of_joins = list_of_joins[0].split(" ")

        list_join_dicts = []
        for join in list_of_joins:
            joins_dict={}
            # assumption is user provides header with "=" as seperator, so there
            # should be one equal sign
            if join.count('=') not in (1, 2):
                raise TypeError (
                    'WARNING: --join_columns requires the header to be stated '
                    'with an equal symbol(e.g --join_columns="Prev_Count=CSQ_Prev_Count_AC,/,CSQ_Prev_Count_NS"'
                )

            joins_dict["header"], column_to_join = join.split("=", 1)

            # assumption is user provides "," as seperator, so there
            # should be two or three commas (three due to column seperator being chosen as comma)
            if join.count(',') not in (2, 3):
                raise TypeError (
                     'WARNING: --join_columns requires the seperator of columns '
                    'to be comma (e.g --join_columns="Prev_Count=CSQ_Prev_Count_AC,/,CSQ_Prev_Count_NS"). '
                    'Expects either two or three minimum commas in string.'
                )
            # split differs depending on column seperator, so check this
            # first and last element should be strings and no
            if len(column_to_join.split(',')) == 3 :
                # user has not provided comma as the column seperator
                # e.g Prev_Count=CSQ_Prev_Count_AC,/,CSQ_Prev_Count_NS
                # so the split becomes
                # ['CSQ_Prev_Count_AC', '/', 'CSQ_Prev_Count_NS']

                # check that a blank space was not provided
                if column_to_join.split(',')[1] == '':
                    raise TypeError (
                        'WARNING: --join_columns requires the seperator of columns '
                        'to be provided.'
                    )
                joins_dict['join_1'], joins_dict['seperator'], joins_dict['join_2'] = column_to_join.split(',')
            else:
                # if user provides a comma as the seperator for the
                # column the split wont work cleanly
                # e.g Prev_Count=CSQ_Prev_Count_AC,,,CSQ_Prev_Count_NS
                # so the split becomes
                # ['CSQ_Prev_Count_AC', '', '', 'CSQ_Prev_Count_NS']
                joins_dict['join_1'], joins_dict['join_2'] = column_to_join.split(',')[0::3]
                joins_dict['seperator'] = ','
            list_join_dicts.append(joins_dict)


        for idx, vcf in enumerate(vcfs):
            for join_dicts in list_join_dicts:
                columns_requested = [joins_dict['join_1'], joins_dict['join_2']]
                # check that the headers exist in VCF, if not raise warning
                if not all(x in vcf.columns for x in columns_requested):
                    print (
                    f"WARNING: Column {columns_requested} requested "
                    f"to join do not exist in VCF."
                )
                else:
                    # sometimes we may want to combine two int columns (e.g chrom,pos)
                    # consider that it will need to be sring type for ease
                    vcf[join_dicts['header']] = vcf[join_dicts['join_1']].astype(str) + join_dicts['seperator']  + vcf[join_dicts['join_2']].astype(str)
                    vcfs[idx] = vcf


        return vcfs


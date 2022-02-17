from utils.arguments import arguments
from utils.excel import excel
from utils.vcf import vcf


def main():
    parser = arguments()

    # read in and process vcf(s)
    vcf_handler = vcf(parser.args)
    vcf_handler.process()

    # generate output Excel file
    excel_handler = excel(parser.args, vcf_handler.vcfs, vcf_handler.refs)
    excel_handler.generate()


if __name__ == "__main__":
    main()

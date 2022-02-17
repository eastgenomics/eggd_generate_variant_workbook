from arguments import arguments
from excel import excel
from vcf import vcf


def main():
    parser = arguments()
    vcf_handler = vcf(parser.args)
    excel(parser.args, vcf_handler.vcfs, vcf_handler.refs)


if __name__ == "__main__":
    main()

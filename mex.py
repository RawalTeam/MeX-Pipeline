import argparse

from mex_pipeline import Mex


def cli_options():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    required_parser = parser.add_argument_group("required arguments")
    required_parser.add_argument("-1", "--fq1", required=True, help="FASTQ Read 1")
    required_parser.add_argument("-g", "--genome", required=True, help="Genome FASTA")
    required_parser.add_argument("-te", "--te", required=True, help="TE FASTA")
    required_parser.add_argument("-O", "--outdir", required=True, help="Output Directory")

    optional_parser = parser.add_argument_group("optional arguments")
    optional_parser.add_argument("-h", "--help", action="help", help="show this help message and exit")
    optional_parser.add_argument("-2", "--fq2", required=False, help="FASTQ Read 2", default=None)
    optional_parser.add_argument("-p", "--processes", required=False, help="Number of processes for multiprocessing",
                                 default=2, type=int)
    optional_parser.add_argument("--force", required=False, help="Rerun entire MeX pipeline",
                                 action='store_true')

    ngstemapper_parser = parser.add_argument_group("ngs_te_mapper2 arguments",
                                                   description="https://github.com/bergmanlab/ngs_te_mapper2"
                                                               "#command-line-help-page")
    ngstemapper_parser.add_argument("--annotation", required=False, help="reference TE annotation in GFF3 "
                                                                         "format (must have 'Target' attribute"
                                                                         " in the 9th column)",
                                    default=None)
    ngstemapper_parser.add_argument("--window", required=False, help="merge window for identifying TE clusters",
                                    type=int, default=10)
    ngstemapper_parser.add_argument("--min_mapq", required=False, help="minimum mapping quality of alignment",
                                    type=int, default=20)
    ngstemapper_parser.add_argument("--min_af", required=False, help="minimum allele frequency",
                                    type=float, default=0.1)
    ngstemapper_parser.add_argument("--tsd_max", required=False, help="maximum TSD size", type=int, default=25)
    ngstemapper_parser.add_argument("--gap_max", required=False, help="maximum gap size", type=int, default=5)
    ngstemapper_parser.add_argument("--keep_files", required=False, help="If provided then all ngs_te_mapper2 "
                                                                         "intermediate files will be kept",
                                    action='store_true')

    vep_parser = parser.add_argument_group("Ensembl Variant Effect Predictor (VEP) arguments",
                                           description="https://asia.ensembl.org/info/docs/tools/vep/script/"
                                                       "vep_options.html#basic")
    vep_parser.add_argument("--assembly", required=False, help="Genome assembly ex., GRCh38, GRCh37, and other. "
                                                               "See VEP docs (https://www.ensembl.org/info/docs/"
                                                               "tools/vep/script/vep_other.html#assembly)",
                            default="GRCh38", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = cli_options()
    Mex(args).execute_pipeline()

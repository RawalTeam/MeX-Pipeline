import argparse
import os

from mex_pipeline import Mex


def cli_options():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-1", "--fq1", required=True, help="FASTQ Read 1 *")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA *")
    parser.add_argument("-te", "--te", required=True, help="TE FASTA *")
    parser.add_argument("-O", "--outdir", required=True, help="Output Directory *")
    parser.add_argument("-2", "--fq2", required=False, help="FASTQ Read 2", default=None)
    parser.add_argument("-p", "--processes", required=False, help="Number of processes for multiprocessing",
                        default=2, type=int)
    parser.add_argument("-a", "--assembly", required=False, help="Genome assembly ex., GRCh38, GRCh37, and other. "
                                                                 "See VEP docs (https://www.ensembl.org/info/docs/"
                                                                 "tools/vep/script/vep_other.html#assembly)",
                        default="GRCh38", type=str)
    parser.add_argument("--force", required=False, help="Rerun entire MeX pipeline",
                        action='store_true')
    return parser.parse_args()


if __name__ == "__main__":
    args = cli_options()
    Mex(args).execute_pipeline()

import os
import sys

# noinspection PyUnresolvedReferences
sys.path.append(snakemake.params.mex_path)
from mex_pipeline import MexUtils


def execute(vcf, assembly, output, log_file, threads):
    if MexUtils.check_dir_exists(os.path.dirname(output)):
        MexUtils.rm_dirs(os.path.dirname(output))

    if not MexUtils.make_dirs(os.path.dirname(output)):
        sys.exit(1)

    cmd = [f"vep -i {vcf} -o {output} -e --warning_file {output + '_warn.txt'} --fork {threads} --cache --refseq"]
    cmd = " ".join(cmd)
    MexUtils.run_subprocess(cmd, log_file)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    execute(vcf=snakemake.input.vcf,
            output=snakemake.output[0],
            threads=snakemake.threads,
            log_file=snakemake.params.log_file,
            assembly=snakemake.params.assembly)

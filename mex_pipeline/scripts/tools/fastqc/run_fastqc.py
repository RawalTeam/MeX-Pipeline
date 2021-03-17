import os
import sys

# noinspection PyUnresolvedReferences
sys.path.append(snakemake.params.mex_path)
from mex_pipeline import MexUtils


def execute(fq1, fq2, paired, output_dir, thread, log_file):

    if MexUtils.check_dir_exists(output_dir):
        MexUtils.rm_dirs(output_dir)

    if not MexUtils.make_dirs(output_dir):
        sys.exit(1)

    cmd = ["fastqc " + fq1]
    if paired:
        cmd.append(fq2)
    cmd.append(f"--outdir={output_dir}")
    cmd.append(f"--threads={thread}")

    cmd = " ".join(cmd)
    MexUtils.run_subprocess(cmd, log_file)

    with open(os.path.join(output_dir, "fastqc_done.txt"), "w") as f:
        f.write("done")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    execute(fq1=snakemake.input.fq1,
            fq2=snakemake.params.fq2,
            paired=snakemake.params.paired,
            output_dir=snakemake.params.output_dir,
            thread=snakemake.threads,
            log_file=snakemake.params.log_file)

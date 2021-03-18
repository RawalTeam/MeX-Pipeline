import sys
import os

import time

# noinspection PyUnresolvedReferences
sys.path.append(snakemake.params.mex_path)
from mex_pipeline import MexUtils


def execute(input1_file, input2_file, output1_file, output2_file, log_file):
    output_dir = os.path.dirname(output1_file)

    if MexUtils.check_dir_exists(output_dir):
        MexUtils.rm_dirs(output_dir)

    if not MexUtils.make_dirs(output_dir):
        sys.exit(1)

    cmd = ["fastp", f"-i {input1_file}", f"-o {output1_file}"]

    if input2_file is not None and output2_file is not None:
        cmd.append(f"-I {input2_file} -O {output2_file}")

    cmd.append(f"-j {os.path.dirname(output1_file)}/fastp.json -h {os.path.dirname(output1_file)}/fastp.html")

    cmd = " ".join(cmd)
    MexUtils.run_subprocess(cmd, log_file)


if __name__ == "__main__":
    start = time.time()
    # noinspection PyUnresolvedReferences
    execute(input1_file=snakemake.input.fq1,
            input2_file=snakemake.params.fq2,
            output1_file=snakemake.params.output1_file,
            output2_file=snakemake.params.output2_file,
            log_file=snakemake.params.log_file)
    exec_time = (time.time() - start) / 60
    MexUtils.console_print("info", f"FASTP completed in {exec_time} minutes")

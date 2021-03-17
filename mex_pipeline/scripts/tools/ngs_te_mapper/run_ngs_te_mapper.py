import os
import sys
import json
from glob import glob

# noinspection PyUnresolvedReferences
sys.path.append(snakemake.params.mex_path)
from mex_pipeline import MexUtils
from mex_pipeline.scripts.utils.output import NgsTeMapper, make_vcf_files


def execute(fq1, fq2, paired, genome, te_file, output_dir, tool_config, log_file, threads, tsd=20):
    if len(glob(os.path.join(output_dir, "bed_tsd", "*.bed"))) > 0:
        bed_file = glob(os.path.join(output_dir, "bed_tsd", "*.bed"))[0]
        inserts = NgsTeMapper.get_inserts(bed_file)
        make_vcf_files(inserts, output_dir)
        return

    if MexUtils.check_dir_exists(output_dir):
        MexUtils.rm_dirs(output_dir)

    if not MexUtils.make_dirs(output_dir):
        sys.exit(1)

    with open(tool_config) as f:
        tool_config = json.load(f)

    samples = [fq1]
    if paired:
        samples.append(fq2)
    samples = "\;".join(samples)

    cmd = ['Rscript',
           "--vanilla",
           tool_config["ngs_te_mapper"]['install_path'] + "/sourceCode/ngs_te_mapper.R",
           "genome=" + genome,
           "teFile=" + te_file,
           "tsd=" + str(tsd),
           "thread=" + str(threads),
           "output=" + output_dir,
           "sourceCodeFolder=" + tool_config["ngs_te_mapper"]['install_path'] + "/sourceCode/",
           "sample=" + samples]
    cmd = " ".join(cmd)
    MexUtils.run_subprocess(cmd, log_file)

    if len(glob(os.path.join(output_dir, "bed_tsd", "*.bed"))) > 0:
        bed_file = glob(os.path.join(output_dir, "bed_tsd", "*.bed"))[0]
        inserts = NgsTeMapper.get_inserts(bed_file)
        make_vcf_files(inserts, output_dir)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    execute(fq1=snakemake.input.fq1,
            genome=snakemake.input.genome,
            te_file=snakemake.input.te_fasta,
            fq2=snakemake.params.fq2,
            paired=snakemake.params.paired,
            tool_config=snakemake.params.tool_config,
            output_dir=snakemake.params.output_dir,
            threads=snakemake.threads,
            log_file=snakemake.params.log_file)

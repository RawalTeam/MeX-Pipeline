import os
import sys
import json
from glob import glob
import time

# noinspection PyUnresolvedReferences
sys.path.append(snakemake.params.mex_path)
from mex_pipeline import MexUtils
from mex_pipeline.scripts.utils.output import NgsTeMapper2, make_vcf_file


def execute(fq1, fq2, paired, genome, te_file, output_dir, tool_config, log_file, threads,
            annotation, window, min_mapq, min_af, tsd_max, gap_max, keep_files):
    if len(glob(os.path.join(output_dir, "*.bed"))) > 0:
        bed_nonref_file = glob(os.path.join(output_dir, "*.nonref.bed"))[0]
        bed_ref_file = glob(os.path.join(output_dir, "*.ref.bed"))[0]
        inserts_nonref = NgsTeMapper2.get_inserts(bed_nonref_file, "non-reference")
        inserts_ref = NgsTeMapper2.get_inserts(bed_ref_file, "reference")
        make_vcf_file(inserts_nonref, output_dir, "non_reference")
        make_vcf_file(inserts_ref, output_dir, "reference")
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
    samples = ",".join(samples)

    cmd = ['python',
           tool_config["ngs_te_mapper2"]['install_path'] + "/sourceCode/ngs_te_mapper2.py",
           "-f " + str(samples),
           "-r " + str(genome),
           "-l " + str(te_file),
           "-t " + str(threads),
           "--window " + str(window),
           "--min_mapq " + str(min_mapq),
           "--min_af " + str(min_af),
           "--tsd_max " + str(tsd_max),
           "--gap_max " + str(gap_max),
           "--out " + str(output_dir)]
    if annotation is not None:
        cmd.append("--annotation " + str(annotation))
    if keep_files:
        cmd.append("--keep_files")
    cmd = " ".join(cmd)
    MexUtils.console_print("cmd", cmd)
    MexUtils.run_subprocess(cmd, log_file)

    if len(glob(os.path.join(output_dir, "*.bed"))) > 0:
        bed_nonref_file = glob(os.path.join(output_dir, "*.nonref.bed"))[0]
        bed_ref_file = glob(os.path.join(output_dir, "*.ref.bed"))[0]
        inserts_nonref = NgsTeMapper2.get_inserts(bed_nonref_file, "non-reference")
        inserts_ref = NgsTeMapper2.get_inserts(bed_ref_file, "reference")
        make_vcf_file(inserts_nonref, output_dir, "non_reference")
        make_vcf_file(inserts_ref, output_dir, "reference")


if __name__ == "__main__":
    start = time.time()
    # noinspection PyUnresolvedReferences
    execute(fq1=snakemake.input.fq1,
            genome=snakemake.input.genome,
            te_file=snakemake.input.te_fasta,
            fq2=snakemake.params.fq2,
            paired=snakemake.params.paired,
            tool_config=snakemake.params.tool_config,
            output_dir=snakemake.params.output_dir,
            threads=snakemake.threads,
            log_file=snakemake.params.log_file,
            annotation=snakemake.params.annotation,
            window=snakemake.params.window,
            min_mapq=snakemake.params.min_mapq,
            min_af=snakemake.params.min_af,
            tsd_max=snakemake.params.tsd_max,
            gap_max=snakemake.params.gap_max,
            keep_files=snakemake.params.keep_files)
    exec_time = (time.time() - start) / 60
    MexUtils.console_print("info", f"NGS_TE_MAPPER2 completed in {exec_time} minutes")

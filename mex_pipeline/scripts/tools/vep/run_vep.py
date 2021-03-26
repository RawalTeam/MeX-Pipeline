import os
import sys
import time

# noinspection PyUnresolvedReferences
sys.path.append(snakemake.params.mex_path)
from mex_pipeline import MexUtils


def execute(ngs_te_mapper2_vcfs, assembly, outputs, log_file, threads, cache_dir):
    ngs_te_mapper2_non_ref_vcf = ngs_te_mapper2_vcfs[0]
    ngs_te_mapper2_ref_vcf = ngs_te_mapper2_vcfs[1]
    ngs_te_mapper2_non_ref_ann_vcf = outputs[0]
    ngs_te_mapper2_ref_ann_vcf = outputs[1]

    if MexUtils.check_dir_exists(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf)):
        MexUtils.rm_dirs(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf))

    if not MexUtils.make_dirs(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf)):
        sys.exit(1)

    if MexUtils.check_file_exists(ngs_te_mapper2_non_ref_vcf):
        if not MexUtils.is_file_empty(ngs_te_mapper2_non_ref_vcf):
            cmd = f"vep -i {ngs_te_mapper2_non_ref_vcf} -o {ngs_te_mapper2_non_ref_ann_vcf} -e --warning_file " \
                  f"{ngs_te_mapper2_non_ref_ann_vcf + '_warn.txt'} --fork {threads} --cache --refseq -a {assembly} " \
                  f"--dir {cache_dir}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
        else:
            cmd = f"touch {ngs_te_mapper2_non_ref_ann_vcf}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
    else:
        cmd = f"touch {ngs_te_mapper2_non_ref_ann_vcf}"
        MexUtils.run_subprocess(cmd, log_file, log_mode="a+")

    if MexUtils.check_file_exists(ngs_te_mapper2_ref_vcf):
        if not MexUtils.is_file_empty(ngs_te_mapper2_ref_vcf):
            cmd = f"vep -i {ngs_te_mapper2_ref_vcf} -o {ngs_te_mapper2_ref_ann_vcf} -e --warning_file " \
                  f"{ngs_te_mapper2_ref_ann_vcf + '_warn.txt'} --fork {threads} --cache --refseq -a {assembly} " \
                  f"--dir {cache_dir}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
        else:
            cmd = f"touch {ngs_te_mapper2_ref_ann_vcf}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
    else:
        cmd = f"touch {ngs_te_mapper2_ref_ann_vcf}"
        MexUtils.run_subprocess(cmd, log_file, log_mode="a+")


if __name__ == "__main__":
    start = time.time()
    # noinspection PyUnresolvedReferences
    execute(ngs_te_mapper2_vcfs=snakemake.input.ngs_te_mapper2_vcfs,
            outputs=snakemake.output,
            threads=snakemake.threads,
            log_file=snakemake.params.log_file,
            assembly=snakemake.params.assembly,
            cache_dir=snakemake.params.cache_dir)
    exec_time = (time.time() - start) / 60
    MexUtils.console_print("info", f"VEP completed in {exec_time} minutes")

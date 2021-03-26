import os
import sys
import time

# noinspection PyUnresolvedReferences
sys.path.append(snakemake.params.mex_path)
from mex_pipeline import MexUtils


def execute(ngs_te_mapper2_vcfs, assembly, outputs, log_file, threads, cache_dir):
    ngs_te_mapper2_non_ref_vcf = ngs_te_mapper2_vcfs[0]
    ngs_te_mapper2_ref_vcf = ngs_te_mapper2_vcfs[1]
    ngs_te_mapper2_non_ref_ann_vcf = outputs[1]
    ngs_te_mapper2_ref_ann_vcf = outputs[0]

    ngs_te_mapper2_non_ref_ann_vcf_name = os.path.basename(outputs[1]).split(".")[0]
    ngs_te_mapper2_ref_ann_vcf_name = os.path.basename(outputs[0]).split(".")[0]

    out_dir = os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf)

    if MexUtils.check_dir_exists(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf)):
        MexUtils.rm_dirs(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf))

    if not MexUtils.make_dirs(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf)):
        sys.exit(1)

    warn_dir = os.path.join(out_dir, "warn")
    html_dir = os.path.join(out_dir, "html")

    if not MexUtils.make_dirs(warn_dir):
        sys.exit(1)

    if not MexUtils.make_dirs(html_dir):
        sys.exit(1)

    if MexUtils.check_file_exists(ngs_te_mapper2_non_ref_vcf):
        if not MexUtils.is_file_empty(ngs_te_mapper2_non_ref_vcf):
            cmd = [
                "vep",
                f"-i {ngs_te_mapper2_non_ref_vcf}",
                f"-o {ngs_te_mapper2_non_ref_ann_vcf}",
                f"--warning_file {os.path.join(warn_dir, ngs_te_mapper2_non_ref_ann_vcf_name + '_warn.txt')}",
                "-e"
                "--cache",
                "--refseq",
                f"-a {assembly}",
                f"--dir {cache_dir}",
                "--tab",
                "--fields",
                "\"Uploaded_variation,Location,SYMBOL,ZYG,CANONICAL,Allele,Gene,Feature,Feature_type,Consequence,"
                "Amino_acids,Codons,Existing_variation,IMPACT,VARIANT_CLASS,HGVSc,HGVSp,HGVSg,INTRON,EXON,BIOTYPE,"
                "AF,PUBMED,PHENO,SIFT,PolyPhen\"",
                f"--stats_file {os.path.join(html_dir, ngs_te_mapper2_non_ref_ann_vcf_name + '.html')}"
            ]
            cmd = " ".join(cmd)
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
        else:
            cmd = f"touch {ngs_te_mapper2_non_ref_ann_vcf}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
    else:
        cmd = f"touch {ngs_te_mapper2_non_ref_ann_vcf}"
        MexUtils.run_subprocess(cmd, log_file, log_mode="a+")

    if MexUtils.check_file_exists(ngs_te_mapper2_ref_vcf):
        if not MexUtils.is_file_empty(ngs_te_mapper2_ref_vcf):
            cmd = [
                "vep",
                f"-i {ngs_te_mapper2_ref_vcf}",
                f"-o {ngs_te_mapper2_ref_ann_vcf}",
                f"--warning_file {os.path.join(warn_dir, ngs_te_mapper2_ref_ann_vcf_name + '_warn.txt')}",
                "-e"
                "--cache",
                "--refseq",
                f"-a {assembly}",
                f"--dir {cache_dir}",
                "--tab",
                "--fields",
                "\"Uploaded_variation,Location,SYMBOL,ZYG,CANONICAL,Allele,Gene,Feature,Feature_type,Consequence,"
                "Amino_acids,Codons,Existing_variation,IMPACT,VARIANT_CLASS,HGVSc,HGVSp,HGVSg,INTRON,EXON,BIOTYPE,"
                "AF,PUBMED,PHENO,SIFT,PolyPhen\"",
                f"--stats_file {os.path.join(html_dir, ngs_te_mapper2_ref_ann_vcf_name + '.html')}"
            ]
            cmd = " ".join(cmd)
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
        else:
            cmd = f"touch {ngs_te_mapper2_ref_ann_vcf}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
    else:
        cmd = f"touch {ngs_te_mapper2_ref_ann_vcf}"
        MexUtils.run_subprocess(cmd, log_file, log_mode="a+")

    if not MexUtils.is_file_empty(ngs_te_mapper2_ref_ann_vcf):
        cmd = [
            "filter_vep",
            f"-i {ngs_te_mapper2_ref_ann_vcf}",
            "-filter \"CANONICAL\"",
            f"-o {os.path.join(out_dir, 'canonical_' + ngs_te_mapper2_ref_ann_vcf_name)}"
        ]
        cmd = " ".join(cmd)
        MexUtils.run_subprocess(cmd, log_file, log_mode="a+")

    if not MexUtils.is_file_empty(ngs_te_mapper2_non_ref_ann_vcf):
        cmd = [
            "filter_vep",
            f"-i {ngs_te_mapper2_non_ref_ann_vcf}",
            "-filter \"CANONICAL\"",
            f"-o {os.path.join(out_dir, 'canonical_' + ngs_te_mapper2_non_ref_ann_vcf_name)}"
        ]
        cmd = " ".join(cmd)
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

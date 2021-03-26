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

    if MexUtils.check_dir_exists(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf)):
        MexUtils.rm_dirs(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf))

    if not MexUtils.make_dirs(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf)):
        sys.exit(1)

    warn_dir = os.path.join(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf), "warn")
    html_dir = os.path.join(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf), "html")
    canonical_dir = os.path.join(os.path.dirname(ngs_te_mapper2_non_ref_ann_vcf), "canonical")

    if not MexUtils.make_dirs(os.path.dirname(warn_dir)):
        sys.exit(1)

    if not MexUtils.make_dirs(os.path.dirname(html_dir)):
        sys.exit(1)

    if not MexUtils.make_dirs(os.path.dirname(canonical_dir)):
        sys.exit(1)

    if MexUtils.check_file_exists(ngs_te_mapper2_non_ref_vcf):
        if not MexUtils.is_file_empty(ngs_te_mapper2_non_ref_vcf):
            cmd = f"vep -i {ngs_te_mapper2_non_ref_vcf} -o {ngs_te_mapper2_non_ref_ann_vcf} -e --warning_file " \
                  f"{os.path.join(warn_dir, ngs_te_mapper2_non_ref_ann_vcf + '_warn.txt')} --fork {threads} " \
                  f"--cache --refseq -a {assembly}  --dir {cache_dir} --vcf --fields " \
                  f"\"Uploaded variation,Location,SYMBOL,ZYG,CANONICAL,Allele,Gene,Feature,Feature type,Consequence," \
                  f"Amino acid change,Codon change,Co-located variation,IMPACT,VARIANT_CLASS,HGVSc,HGVSp,HGVSg," \
                  f"INTRON,EXON,AF,PUBMED,PHENO,SIFT,PolyPhen\" --stats_file " \
                  f"{os.path.join(html_dir, ngs_te_mapper2_non_ref_ann_vcf + '.html')}"
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
                  f"{os.path.join(warn_dir, ngs_te_mapper2_ref_vcf + '_warn.txt')} --fork {threads} " \
                  f"--cache --refseq -a {assembly}  --dir {cache_dir} --vcf --fields " \
                  f"\"Uploaded variation,Location,SYMBOL,ZYG,CANONICAL,Allele,Gene,Feature,Feature type,Consequence," \
                  f"Amino acid change,Codon change,Co-located variation,IMPACT,VARIANT_CLASS,HGVSc,HGVSp,HGVSg," \
                  f"INTRON,EXON,AF,PUBMED,PHENO,SIFT,PolyPhen\" --stats_file " \
                  f"{os.path.join(html_dir, ngs_te_mapper2_ref_ann_vcf + '.html')}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
        else:
            cmd = f"touch {ngs_te_mapper2_ref_ann_vcf}"
            MexUtils.run_subprocess(cmd, log_file, log_mode="a+")
    else:
        cmd = f"touch {ngs_te_mapper2_ref_ann_vcf}"
        MexUtils.run_subprocess(cmd, log_file, log_mode="a+")

    if not MexUtils.is_file_empty(ngs_te_mapper2_ref_ann_vcf):
        cmd = f"filter_vep -i {ngs_te_mapper2_ref_ann_vcf} -filter \"CANONICAL\" " \
              f"-o {os.path.join(canonical_dir, ngs_te_mapper2_ref_ann_vcf)}"
        MexUtils.run_subprocess(cmd, log_file, log_mode="a+")

    if not MexUtils.is_file_empty(ngs_te_mapper2_non_ref_ann_vcf):
        cmd = f"filter_vep -i {ngs_te_mapper2_non_ref_ann_vcf} -filter \"CANONICAL\" " \
              f"-o {os.path.join(canonical_dir, ngs_te_mapper2_non_ref_ann_vcf)}"
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

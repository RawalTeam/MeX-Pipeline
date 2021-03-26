import subprocess


def execute(log, target, cache_dir, assembly):
    with open(log, "a") as f:
        subprocess.check_call([f"vep_install -a cf -s homo_sapiens_refseq -y {assembly} --CONVERT -c {cache_dir}"],
                              shell=True, stderr=f, stdout=f)

    with open(target, "w") as f:
        subprocess.check_call([f"vep --help"], shell=True, stderr=f, stdout=f)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    execute(log=snakemake.params.log, target=snakemake.output[0], assembly=snakemake.params.assembly,
            cache_dir=snakemake.params.cache_dir)

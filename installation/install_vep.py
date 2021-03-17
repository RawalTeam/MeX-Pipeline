import subprocess


def execute(log, target):
    with open(log, "w") as f:
        subprocess.check_call([f"vep_install -a cf -s homo_sapiens_refseq -y GRCh38 --CONVERT"],
                              shell=True, stderr=f, stdout=f)

    with open(target, "w") as f:
        subprocess.check_call([f"vep --help"],
                              shell=True, stderr=f, stdout=f)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    execute(log=snakemake.params.log, target=snakemake.output[0])

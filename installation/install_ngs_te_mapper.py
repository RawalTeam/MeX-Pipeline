import os
import shutil
import subprocess


def execute(url, commit, log, install_path):
    if os.path.exists(install_path):
        shutil.rmtree(install_path)

    try:
        os.makedirs(install_path)
    except OSError:
        pass

    with open(log, "w") as f:
        subprocess.check_call([f"git clone {url} {install_path}"], shell=True, stderr=f, stdout=f)

    with open(log, "a+") as f:
        subprocess.check_call(['git checkout ' + commit], stdout=f, stderr=f, shell=True, cwd=install_path)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    execute(url=snakemake.params.url,
            commit=snakemake.params.commit,
            log=snakemake.params.log,
            install_path=snakemake.params.install_path)

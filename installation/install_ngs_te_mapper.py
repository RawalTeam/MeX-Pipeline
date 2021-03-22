import os
import shutil
import subprocess


def execute(url, commit, log, install_path):
    if os.path.exists(install_path):
        shutil.rmtree(install_path) #If install path exists then remove the directory entirely 

    try:
        os.makedirs(install_path) # Make directory for installation
    except OSError:
        pass #If there's an OS error , output nothing. 

    with open(log, "w") as f:
        subprocess.check_call([f"git clone {url} {install_path}"], shell=True, stderr=f, stdout=f) #Clone the repository

    with open(log, "a+") as f:
        subprocess.check_call(['git checkout ' + commit], stdout=f, stderr=f, shell=True, cwd=install_path) # git checkout the repository

if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    #Execute Snakemake for install
    execute(url=snakemake.params.url,
            commit=snakemake.params.commit,
            log=snakemake.params.log,
            install_path=snakemake.params.install_path)

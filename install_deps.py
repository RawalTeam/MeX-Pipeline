import argparse
import os
import json
import subprocess

INSTALL_BASE = os.path.dirname(os.path.abspath(__file__)) + "/installation"

INSTALL_DIR = os.path.dirname(os.path.abspath(__file__)) + "/installation/tools"

ENV_DIR = os.path.dirname(os.path.abspath(__file__)) + "/envs"

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

DEP = {
    "ngs_te_mapper": {
        "url": "https://github.com/bergmanlab/ngs_te_mapper.git",
        "commit": "f9f48996ac346ac86d57edbd00534aa1227b753e",
        "install_path": os.path.join(INSTALL_DIR, "ngs_te_mapper"),
        "type": "github",
        "env": os.path.join(ENV_DIR, "ngstemapper.yaml"),
        "log": INSTALL_DIR + "/ngs_te_mapper.log",
        "target": os.path.join(INSTALL_DIR, "ngs_te_mapper/sourceCode/ngs_te_mapper.R"),
        "script": os.path.join(INSTALL_BASE, "install_ngs_te_mapper.py")
    },
    "ngs_te_mapper2": {
        "url": "https://github.com/bergmanlab/ngs_te_mapper2.git",
        "commit": "59d20a4bf1b2491193383293240c1c13dcfaa0d3",
        "install_path": os.path.join(INSTALL_DIR, "ngs_te_mapper2"),
        "type": "github",
        "env": os.path.join(INSTALL_DIR, "ngs_te_mapper2/envs/ngs_te_mapper2.yml"),
        "env_download": os.path.join(ENV_DIR, "ngstemapper2.yaml"),
        "log": INSTALL_DIR + "/ngs_te_mapper2.log",
        "target": os.path.join(INSTALL_DIR, "ngs_te_mapper2/README.md"),
        "script": os.path.join(INSTALL_BASE, "install_ngs_te_mapper2.py")
    },
    "vep": {
        "type": "conda",
        "env": os.path.join(ENV_DIR, "vep.yaml"),
        "log": INSTALL_DIR + "/vep.log",
        "target": os.path.join(INSTALL_DIR, "vep/version.txt"),
        "script": os.path.join(INSTALL_BASE, "install_vep.py")
    },
    "preprocess": {
        "type": "conda",
        "env": os.path.join(ENV_DIR, "mexpreprocess.yaml")
    }
}


def cli_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", required=False, help="Number of processes used", default=2, type=int)
    return parser.parse_args()


def execute(cores):
    try:
        os.makedirs(INSTALL_DIR)
    except OSError:
        pass

    params = []
    for tool, info in DEP.items():
        params.append({"tool": tool, "info": info})

    with open(os.path.join(INSTALL_BASE, "config.json"), "w") as f:
        json.dump(DEP, f, indent=4)

    cmd = [f"snakemake --configfile {os.path.join(INSTALL_BASE, 'config.json')}",
           "--use-conda",
           "--conda-prefix " + os.path.join(ENV_DIR, "conda"),
           f"--cores {cores}",
           DEP['vep']['target'],
           DEP['ngs_te_mapper']['target'],
           DEP['ngs_te_mapper2']['target']]
    subprocess.check_call([" ".join(cmd)], cwd=os.path.join(PROJECT_DIR, "installation"), shell=True)


if __name__ == "__main__":
    options = cli_options()
    execute(options.processes)

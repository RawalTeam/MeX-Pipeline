import argparse
import os
import json
import subprocess
import pathlib
import sys

INSTALL_BASE = os.path.dirname(os.path.abspath(__file__)) + "/installation"

INSTALL_DIR = os.path.dirname(os.path.abspath(__file__)) + "/installation/tools"

ENV_DIR = os.path.dirname(os.path.abspath(__file__)) + "/envs"

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))


def create_dep_config(opts):
    dep = {
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
            "cache_dir": os.path.abspath(opts.cachedir),
            "script": os.path.join(INSTALL_BASE, "install_vep.py"),
            "assemblies": [opts.assembly],
            "only_assembly": opts.assembly
        },
        "preprocess": {
            "type": "conda",
            "env": os.path.join(ENV_DIR, "mexpreprocess.yaml")
        }
    }
    return dep


def cli_options():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-p", "--processes", required=False, help="Number of processes used", default=2, type=int)
    parser.add_argument("-a", "--assembly", required=False, help="Genome assembly ex., GRCh38, GRCh37, and other. "
                                                                 "See VEP docs (https://www.ensembl.org/info/docs/"
                                                                 "tools/vep/script/vep_other.html#assembly)",
                        default="GRCh38", type=str)
    parser.add_argument("-d", "--cachedir", required=False, help="VEP Data directory",
                        default=os.path.join(str(pathlib.Path.home()), ".vep"), type=str)
    parser.add_argument("-oa", "--only-assembly", required=False,
                        help="Download Genome assembly ex., GRCh38, GRCh37, and other. "
                             "See VEP docs (https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#assembly"
                             ") in existing VEP cache directory. Requires config.json in installation directory",
                        type=str)
    return parser.parse_args()


def execute(opts):
    try:
        os.makedirs(INSTALL_DIR)
    except OSError:
        pass

    dep = create_dep_config(opts)

    params = []
    for tool, info in dep.items():
        params.append({"tool": tool, "info": info})

    with open(os.path.join(INSTALL_BASE, "config.json"), "w") as f:
        json.dump(dep, f, indent=4)

    cmd = [f"snakemake --configfile {os.path.join(INSTALL_BASE, 'config.json')}",
           "--use-conda",
           "--conda-prefix " + os.path.join(ENV_DIR, "conda"),
           f"--cores {opts.processes}",
           dep['vep']['target'],
           dep['ngs_te_mapper2']['target']]
    subprocess.check_call([" ".join(cmd)], cwd=os.path.join(PROJECT_DIR, "installation"), shell=True)


def add_assembly(opts):
    assembly = opts.only_assembly
    if not os.path.exists(os.path.join(INSTALL_BASE, "config.json")):
        sys.stderr.write("config.json not found in installation directory. Run install_deps.py without -oa or "
                         "--only-assembly first")
        sys.exit(1)
    with open(os.path.join(INSTALL_BASE, "config.json")) as f:
        dep = json.load(f)

    assemblies_installed = dep['vep']['assemblies']
    if assembly in assemblies_installed:
        sys.stdout.write(f"{assembly} is already added")
        sys.exit(0)

    dep['vep']['only_assembly'] = assembly
    with open(os.path.join(INSTALL_BASE, "config.json"), "w") as f:
        json.dump(dep, f, indent=4)

    cmd = [f"snakemake --configfile {os.path.join(INSTALL_BASE, 'config.json')}",
           "--use-conda",
           "--conda-prefix " + os.path.join(ENV_DIR, "conda"),
           f"--cores {opts.processes}",
           dep['vep']['target'] + f".{assembly}.txt"]
    try:
        subprocess.check_call([" ".join(cmd)], cwd=os.path.join(PROJECT_DIR, "installation"), shell=True)
        dep['vep']['assemblies'].append(assembly)
        with open(os.path.join(INSTALL_BASE, "config.json"), "w") as f:
            json.dump(dep, f, indent=4)
    except subprocess.CalledProcessError:
        raise


if __name__ == "__main__":
    options = cli_options()
    if options.only_assembly is not None:
        add_assembly(options)
    else:
        execute(options)

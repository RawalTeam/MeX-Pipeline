import shutil
import subprocess
import os
import sys
import json
import traceback
import time
from datetime import datetime
from pprint import pprint


class MexUtils:

    @staticmethod
    def check_fastq(fastq_file):
        if MexUtils.check_file_exists(fastq_file):
            if not MexUtils.is_file_empty(fastq_file):
                return True
        return False

    @staticmethod
    def check_fasta(fasta_file):
        if MexUtils.check_file_exists(fasta_file):
            if not MexUtils.is_file_empty(fasta_file):
                return True
        return False

    @staticmethod
    def check_file_exists(f):
        if os.path.exists(f):
            if os.path.isfile(f):
                return True
        return False

    @staticmethod
    def check_dir_exists(f):
        if os.path.exists(f):
            if os.path.isdir(f):
                return True
        return False

    @staticmethod
    def is_file_empty(f):
        if os.stat(f).st_size == 0:
            return True
        return False

    @staticmethod
    def make_dirs(f):
        if not MexUtils.check_dir_exists(f):
            try:
                os.makedirs(f)
                return True
            except OSError:
                return False
        return True

    @staticmethod
    def rm_dirs(f):
        if MexUtils.check_dir_exists(f):
            try:
                shutil.rmtree(f)
                return True
            except OSError:
                return False
        return True

    @staticmethod
    def run_subprocess(cmd, log=None, fatal=True, cwd=None, log_mode="w"):
        if cwd is None:
            if log is not None:
                with open(log, log_mode) as f:
                    try:
                        f.write(cmd)
                        subprocess.check_call(cmd, stdout=f, stderr=f, shell=True)
                    except subprocess.CalledProcessError:
                        track = traceback.format_exc()
                        error = f"\nERROR **** {cmd} **** ERROR\n{track}"
                        MexUtils.write_into_log(f, error)
                        sys.stderr.write(track)
                        if fatal:
                            sys.exit(1)
            else:
                try:
                    subprocess.check_call(cmd, shell=True)
                except subprocess.CalledProcessError:
                    track = traceback.format_exc()
                    sys.stderr.write(track)
                    if fatal:
                        sys.exit(1)
        else:
            if log is not None:
                with open(log, log_mode) as f:
                    try:
                        f.write(cmd)
                        subprocess.check_call(cmd, stdout=f, stderr=f, shell=True, cwd=cwd)
                    except subprocess.CalledProcessError:
                        track = traceback.format_exc()
                        error = f"\nERROR **** {cmd} **** ERROR\n{track}"
                        MexUtils.write_into_log(f, error)
                        sys.stderr.write(track)
                        if fatal:
                            sys.exit(1)
            else:
                try:
                    subprocess.check_call(cmd, shell=True, cwd=cwd)
                except subprocess.CalledProcessError:
                    track = traceback.format_exc()
                    sys.stderr.write(track)
                    if fatal:
                        sys.exit(1)

    @staticmethod
    def write_into_log(log, message):
        log.write(message)

    @staticmethod
    def is_paired_reads_provided(fq2):
        paired = False
        if fq2 is not None:
            paired = True
        return paired

    @staticmethod
    def console_print(typ, msg):
        print(f"[{datetime.now().strftime('%d-%b-%Y %H:%M:%S')}] [{typ.upper()}] {msg}")

    @staticmethod
    def calculate_max_threads_per_rule(max_threads, multi_methods_used):
        if multi_methods_used > 1:
            is_even = False
            if max_threads % 2 == 0:
                is_even = True

            max_threads = max_threads // 2
            if is_even:
                max_threads = max_threads - 1
        return max_threads


class Mex:

    def __init__(self, cli_args):
        self.PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

        self.fastq1 = os.path.abspath(cli_args.fq1)
        if cli_args.fq2 is not None:
            self.fastq2 = os.path.abspath(cli_args.fq2)
        else:
            self.fastq2 = None
        self.genome = os.path.abspath(cli_args.genome)
        self.tefasta = os.path.abspath(cli_args.te)
        self.out_dir = os.path.abspath(cli_args.outdir)
        self.threads = cli_args.processes
        self.force = cli_args.force
        self.assembly = cli_args.assembly

        self.genome_name = os.path.basename(self.genome).split(".")[0]

        MexUtils.make_dirs(self.out_dir)
        MexUtils.make_dirs(os.path.join(self.out_dir, "logs"))
        self.check_input_files()
        self.place_genome_fasta()
        self.place_te_fasta()

        self.tools_config = os.path.join(self.PROJECT_DIR, "installation/config.json")
        with open(self.tools_config) as f:
            self.tools_config_obj = json.load(f)
        self.config = self.create_config()

        with open(os.path.join(self.out_dir, "config.json"), "w") as config_json:
            json.dump(self.config, config_json, indent=4)

        self.snakefile_path = os.path.join(self.PROJECT_DIR, "Snakefile")

    def place_genome_fasta(self):
        genome_dir = os.path.join(self.out_dir, self.genome_name, "genome")
        MexUtils.make_dirs(genome_dir)
        try:
            os.symlink(self.genome, os.path.join(genome_dir, self.genome_name + ".fasta"))
        except FileExistsError:
            pass
        self.genome = os.path.join(genome_dir, self.genome_name + ".fasta")

    def place_te_fasta(self):
        te_dir = os.path.join(self.out_dir, self.genome_name, "te")
        MexUtils.make_dirs(te_dir)
        try:
            os.symlink(self.tefasta, os.path.join(te_dir, os.path.basename(self.tefasta)))
        except FileExistsError:
            pass
        self.tefasta = os.path.join(te_dir, os.path.basename(self.tefasta))

    def create_config(self):
        scripts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")
        is_paired = MexUtils.is_paired_reads_provided(self.fastq2)

        fqp1_name = f"{os.path.basename(self.fastq1).split('.')[0]}_fq1_preprocessed"
        fqp2_name = ""

        sample_name = fqp1_name
        if is_paired:
            fqp2_name = f"{os.path.basename(self.fastq2).split('.')[0]}_fq2_preprocessed"
            sample_name += fqp2_name

        config_obj = {
            "args": {
                "fastq1": self.fastq1,
                "fastq2": self.fastq2,
                "genome": self.genome,
                "te": self.tefasta,
                "outdir": self.out_dir,
                "force": self.force,
                "processes": self.threads,
                "threads": max(1, MexUtils.calculate_max_threads_per_rule(int(self.threads), multi_methods_used=5)),
                "assembly": self.assembly
            },
            "params": {
                "paired": is_paired,
                "tools_config": self.tools_config,
                "mex_path": os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                "vep_cache_dir": self.tools_config_obj['vep']['cache_dir']
            },
            "outputs": {
                "fastp": {
                    "fqp1": os.path.join(self.out_dir, "outputs/fastp", fqp1_name + ".fastq"),
                    "fqp2": os.path.join(self.out_dir, "outputs/fastp", fqp2_name + ".fastq") if is_paired else None
                },
                "fastqc": os.path.join(self.out_dir, "outputs/fastqc"),
                "ngs_te_mapper2": os.path.join(self.out_dir, "outputs/ngs_te_mapper2"),
                "vep": os.path.join(self.out_dir, "outputs/vep")
            },
            "targets": {
                "fastp": os.path.join(self.out_dir, "outputs/fastp", "fastp.html"),
                "fastqc": os.path.join(self.out_dir, "outputs/fastqc", "fastqc_done.txt"),
                "ngs_te_mapper2": [
                    os.path.join(self.out_dir, "outputs/ngs_te_mapper2", "non_reference.vcf"),
                    os.path.join(self.out_dir, "outputs/ngs_te_mapper2", "reference.vcf")
                ],
                "vep": [
                    os.path.join(self.out_dir, "outputs/vep/ngs_te_mapper2", sample_name + "_reference_ann.vcf"),
                    os.path.join(self.out_dir, "outputs/vep/ngs_te_mapper2", sample_name + "_non_reference_ann.vcf"),
                ]
            },
            "envs": {
                "ngs_te_mapper2": self.tools_config_obj['ngs_te_mapper2']['env'],
                "preprocessing": self.tools_config_obj['preprocess']['env'],
                "vep": self.tools_config_obj['vep']['env']
            },
            "logs": {
                "fastp": os.path.join(self.out_dir, "logs/fastp.log"),
                "fastqc": os.path.join(self.out_dir, "logs/fastqc.log"),
                "ngs_te_mapper2": os.path.join(self.out_dir, "logs/ngs_te_mapper2.log"),
                "vep": os.path.join(self.out_dir, "logs/vep.log")
            },
            "scripts": {
                "fastp": os.path.join(scripts_dir, "tools/fastp/run_fastp.py"),
                "fastqc": os.path.join(scripts_dir, "tools/fastqc/run_fastqc.py"),
                "ngs_te_mapper2": os.path.join(scripts_dir, "tools/ngs_te_mapper2/run_ngs_te_mapper2.py"),
                "vep": os.path.join(scripts_dir, "tools/vep/run_vep.py")
            }
        }
        return config_obj

    def check_input_files(self):
        if not MexUtils.check_fasta(self.genome):
            sys.exit(1)
        if not MexUtils.check_fasta(self.tefasta):
            sys.exit(1)
        if not MexUtils.check_fastq(self.fastq1):
            sys.exit(1)
        if self.fastq2 is not None and not MexUtils.check_fastq(self.fastq2):
            sys.exit(1)

    def execute_pipeline(self):
        start = time.time()
        print(f"{'*'*25} Executing MeX (v1) {'*'*25}")
        print("Parameters info:")
        print(f"FASTQ 1 = {self.fastq1}")
        print(f"FASTQ 2 = {self.fastq2}")
        print(f"GENOME = {self.genome}")
        print(f"TE FASTA = {self.tefasta}")
        print(f"ASSEMBLY = {self.assembly}")
        print(f"OUTPUT DIRECTORY = {self.out_dir}/outputs")
        print(f"LOG DIRECTORY = {self.out_dir}/logs")
        print(f"NUMBER OF PROCESS = {self.threads}")
        print(f"CONFIG FILE = {os.path.join(self.out_dir, 'config.json')}")
        print(f"FORCE = {self.force}")

        shutil.copy(self.snakefile_path, os.path.join(self.out_dir, "Snakefile"))

        cmd = [f"snakemake --configfile {os.path.join(self.out_dir, 'config.json')}"]
        if self.force:
            cmd.append("--force")
        cmd.append("--use-conda")
        cmd.append("--conda-prefix " + os.path.join(self.PROJECT_DIR, "envs/conda"))
        cmd.append(f"--cores {self.threads}")
        cmd.append(self.config['targets']['fastqc'])
        cmd.append(self.config['targets']['vep'][0])
        MexUtils.run_subprocess(" ".join(cmd), fatal=False, cwd=self.out_dir)
        exec_time = (time.time() - start) / 60
        MexUtils.console_print("info", f"MeX completed in {exec_time} minutes")
        MexUtils.console_print("info", "Generating Workflow Report")
        cmd.append(f"--report {self.out_dir + '/workflow.html'}")
        MexUtils.run_subprocess(" ".join(cmd), fatal=False, cwd=self.out_dir)
        end = (time.time() - start) / 60
        MexUtils.console_print("info", f"Done in {end} minutes")

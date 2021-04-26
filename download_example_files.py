import os
import shutil
import subprocess
import sys
import urllib.request
from multiprocessing import Pool

EXAMPLE_FOLDER = os.path.dirname(os.path.abspath(__file__)) + "/example"
MAX_ATTEMPTS = 5


# noinspection PyBroadException
def download(url, out_file, attempts=1):
    tmp_file = out_file + ".tmp"
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    if attempts > MAX_ATTEMPTS:
        print("Maximum attempts reached. Try again later")
        sys.exit(1)

    try:
        print("Downloading", url, "to", out_file)
        urllib.request.urlretrieve(url, tmp_file)
    except KeyboardInterrupt:
        sys.exit(1)
    except:
        print(sys.exc_info())
        print("Download Failed... Retrying...")
        download(url, out_file, attempts + 1)

    shutil.move(tmp_file, out_file)


def uncompress(out):
    subprocess.check_call([f"gzip -k -d {out}"])


def download_files():
    fastq1_url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/sequence_read/" \
                 "SRR622461_1.filt.fastq.gz"
    fastq2_url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/sequence_read/" \
                 "SRR622461_2.filt.fastq.gz"
    genome_url = "https://drive.google.com/file/d/15Z1hzyYzfsVfIqzqAWelZBP1Fxw91Xl7/view?usp=sharing"
    te_url = "https://drive.google.com/file/d/1NYGSSuafwQx7l_cm6h0nVS_tck3uvOLX/view?usp=sharing"

    try:
        os.makedirs(EXAMPLE_FOLDER)
    except OSError:
        pass

    download(fastq1_url, os.path.join(EXAMPLE_FOLDER, "SRR622461_1.filt.fastq.gz"))
    download(fastq2_url, os.path.join(EXAMPLE_FOLDER, "SRR622461_2.filt.fastq.gz"))
    download(genome_url, os.path.join(EXAMPLE_FOLDER, "hg38_chr123.fa.gz"))
    download(te_url, os.path.join(EXAMPLE_FOLDER, "RMRBSeqs_Original_Alu.fasta.gz"))

    with Pool() as p:
        p.map(uncompress, [
            os.path.join(EXAMPLE_FOLDER, "SRR622461_1.filt.fastq.gz"),
            os.path.join(EXAMPLE_FOLDER, "SRR622461_2.filt.fastq.gz"),
            os.path.join(EXAMPLE_FOLDER, "hg38_chr123.fa.gz"),
            os.path.join(EXAMPLE_FOLDER, "RMRBSeqs_Original_Alu.fasta.gz")
        ])


if __name__ == "__main__":
    download_files()

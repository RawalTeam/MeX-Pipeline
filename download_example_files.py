import os
import subprocess
import sys
import urllib.request
import hashlib


EXAMPLE_FOLDER = os.path.dirname(os.path.abspath(__file__)) + "/example"
MAX_ATTEMPTS = 2


# noinspection PyBroadException
def download(url, out_file, hash_md5, attempts=1):
    if attempts > MAX_ATTEMPTS:
        print("Maximum attempts reached. Try again later")
        sys.exit(1)

    try:
        print("Downloading ", url, " to ", out_file)
        urllib.request.urlretrieve(url, out_file)
    except KeyboardInterrupt:
        sys.exit(1)
    except:
        print(sys.exc_info())
        print("Download Failed... Retrying...")
        download(url, out_file, hash_md5, attempts + 1)

    print("Checking MD5 hash of ", out_file)
    with open(out_file, "rb") as f:
        current_md5 = hashlib.md5(f.read()).hexdigest()

        if current_md5 != hash_md5:
            print("MD5 hash of", out_file, " : ", current_md5, "does not match current MD5 hash: ", hash_md5)
            print("Retrying")
            download(url, out_file, hash_md5, attempts + 1)
        else:
            print("MD5 hash of ", out_file, "matches")
            print("Download complete")

    if out_file.endswith(".fa.gz") or out_file.endswith(".fasta.gz"):
        subprocess.check_call([f"gzip {out_file}"], shell=True)


def download_files():
    fastq1_url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/sequence_read/" \
                 "SRR622461_1.filt.fastq.gz"
    fastq2_url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/sequence_read/" \
                 "SRR622461_2.filt.fastq.gz"
    genome_url = "https://drive.google.com/file/d/15Z1hzyYzfsVfIqzqAWelZBP1Fxw91Xl7/view?usp=sharing"
    te_url = "https://drive.google.com/file/d/1NYGSSuafwQx7l_cm6h0nVS_tck3uvOLX/view?usp=sharing"

    fastq1_hash = "d9d0d334bd142a08e1992c0d568213ae"
    fastq2_hash = "f2fad79d5e73ec745fa5880f299dc810"
    genome_hash = "61bb21477145ebacece13c53093054d3"
    te_hash = "bfedd1a50690d9f86fe062e6e70bb2dd"

    try:
        os.makedirs(EXAMPLE_FOLDER)
    except OSError:
        pass

    download(fastq1_url, os.path.join(EXAMPLE_FOLDER, "SRR622461_1.filt.fastq.gz"), fastq1_hash)
    download(fastq2_url, os.path.join(EXAMPLE_FOLDER, "SRR622461_2.filt.fastq.gz"), fastq2_hash)
    download(genome_url, os.path.join(EXAMPLE_FOLDER, "hg38_chr123.fa.gz"), genome_hash)
    download(te_url, os.path.join(EXAMPLE_FOLDER, "RMRBSeqs_Original_Alu.fasta.gz"), te_hash)


if __name__ == "__main__":
    download_files()

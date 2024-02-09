#!/usr/bin/env python3
import argparse
import subprocess
from Bio import Entrez
from Bio import SeqIO
import logging
from tqdm import tqdm

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Fetch sequences from NCBI databases.')
    parser.add_argument('--db', type=str, choices=['prot', 'nucl'], required=True, help='Database: prot (protein) or nucl (nucleotide).')
    parser.add_argument('--email', type=str, required=True, help='Email for using NCBI API.')
    parser.add_argument('--infile', type=str, required=True, help='Input file with accession numbers.')
    parser.add_argument('--outfile', type=str, required=True, help='Output FASTA file.')
    return parser.parse_args()

def check_conda_environment(env_name='fetchthat'):
    # Simplified check; consider improving for real-world usage
    try:
        subprocess.run(['conda', 'env', 'list'], check=True)
        logging.info(f"Conda environment '{env_name}' is assumed to be set up.")
    except subprocess.CalledProcessError:
        logging.error(f"Conda environment '{env_name}' is not setup correctly. Please check.")

def fetch_sequences(email, db, infile, outfile):
    Entrez.email = email
    accession_numbers = []

    with open(infile, 'r') as f:
        accession_numbers = [line.strip() for line in f.readlines()]

    sequences = []
    db_name = 'protein' if db == 'prot' else 'nucleotide'

    for accession in tqdm(accession_numbers, desc="Fetching sequences"):
        try:
            handle = Entrez.efetch(db=db_name, id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            sequences.append(record)
            handle.close()
        except Exception as e:
            logging.error(f"Failed to fetch {accession}: {e}")

    SeqIO.write(sequences, outfile, "fasta")
    logging.info(f"Sequences saved to {outfile}")


if __name__ == "__main__":
    args = parse_arguments()
    check_conda_environment()
    fetch_sequences(args.email, args.db, args.infile, args.outfile)

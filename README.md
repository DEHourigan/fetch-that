## fetch that
A tool to fetch protein and nucleotide sequences using Entrez from Biopython

![Example of](data/smash1.jpeg "Example Plot1")



## Requirements

- Python 3.12.1
- tqdm

## Installation

1. Clone the repository to your local machine:
    ```bash
    git clone https://github.com/DEHourigan/fetch-that.git
    ```
2. Navigate to the cloned directory:
    ```bash
    cd fetch-that
    ```
3. Install the required dependencies:
    ```bash
    conda env create -f fetchthat.yml
	conda activate fetchthat
    ```

## Usage
python fetchthat.py --infile test_data/protein_acc.txt --db prot --email [email here] --outfile test_data/out.faa

1. Prepare your data:
    - input is a text file with 1 accession per line


## Output

fasta file containing desired sequences

## Contributing

DEHourigan
 
## Contact

For any queries, please reach out via GitHub issues or directly to `114402828@umail.icc.ie`.

---


---


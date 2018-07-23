# searchPepepe
Align amino acid sequences to genome with Pepepeptide

## Requirement
- Python3
- BioPython

## Usage
usage: searchPepepe.py [-h] -i INPUT [-l LENGTH] [-m MISMATCH]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        a fasta file of amino acids
  -l LENGTH, --length LENGTH
                        default length of query
  -m MISMATCH, --mismatch MISMATCH
                        number of mismatches

## Example
python3 searchPepepe.py -i test.fa -l 13 -m 1

Output is a bed file (test_pepepe.bed)

![This program was developed with ข้าวซอย (Thai curry noodle)](curry_noodle.JPG)

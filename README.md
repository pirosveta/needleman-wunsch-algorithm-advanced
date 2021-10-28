# Needleman-Wunsch algorithm advanced affine gap penalty
## Options:
+ `-g, --gap`       Penalty for the gap   
+ `-a`              Alphabet
+ \* `-i`           Paths to input sequences
+ `-o`              Path to output file


## Examples:
+ `-i ./seq1.fasta ./seq2.fasta -a BLOSUM62 --gap -10`

+ `-i ./seq1.fasta ./seq2.fasta -a DNAFull -g -10 -o ./out.txt`

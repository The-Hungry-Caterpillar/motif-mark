# motif-mark

The purpose of `motif-mark-oop.py` is to generate an image displaying motifs of interest around exons of a given contig.

To use, first download `env.yml`, `bioinfo.py`, and `motif-mark-oop.py` to your working directory.

Then run the command:
```conda env create -f env.yml && conda activate cairo```

Now you're ready to use motif mark.

## Motif-mark Input
* `-f`: input fasta file with lowercase representing introns and uppercase represting exons.
* `-m`: input sequences of motifs, separated by new lines. Can be upper or lower case. Can be [ambiguously labeled](https://en.wikipedia.org/wiki/Nucleic_acid_notation) (note that if your motif contains "U" then this code will look for "U" in the input fasta file, as noted in the linked wiki table.)
* `-d`: pass this option to have the output suited to a dark background

## Output
A really nice looking .png image of all input contigs and features. The motifs will be randomly colored and stacked, so if you're having trouble with the colors simply run the script again.


## Example Usage
```python motif-mark-oop.py -f data/Figure_1.fasta -m data/Fig_1_motifs.txt```

Out:
`Figure_1.fasta.png`:
![Figure_1.fasta.png](figures/Figure_1.fasta.png)

or,

`Figure_1.fasta.dark.png`:
![Figure_1.fasta.dark.png](figures/Figure_1.fasta.dark.png)

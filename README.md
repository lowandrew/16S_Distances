# 16S Distances

Want to know how close some 16S regions tend to be between two different genera/families/phyla? This is the repo for
you.

### Installation

- Clone this repository
- Make sure you have the BBMap suite installed and accessible on your $PATH
- Make sure you have Biopython installed

### Usage

- Download the SILVA 16S database in FASTA format. Make sure you choose the ungapped (aka unaligned) version.
- Extract your region of interest using your favourite primers using the `extract_primer_regions.py` script.
(This will probably take a while)
- With the output FASTA file you get from the previous step, run `16s_distances.py` on your genera of interest.

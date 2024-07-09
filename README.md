# HaploCrackling

**Computational design of allele-specific CRISPR guide RNAs for phased genomes**

Bo Zhou, Steve S. Ho, Louis C. Leung, Thomas R. Ward, Marcus Ho, Melanie J. Plastini, Scott C. Vermilyea, Marina E. Emborg, Thaddeus G. Golos, Megan A. Albertelli, Philippe Mourrain, Dimitri Perrin, Karen J. Parker, Alexander E. Urban.  Haplotype-phased common marmoset embryonic stem cells for genome editing using CRISPR/Cas9. bioRxiv, 21 October 2020, 373886; doi: https://doi.org/10.1101/373886 

## Preamble

This tool has been forked from our open-source tool called Crackling, which is available at https://github.com/bmds-lab/Crackling under the Berkeley Software Distribution (BSD) 3-Clause license.

## Dependencies

- ISSL-based search off-target sites (included)

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)

- sgRNAScorer 2.0 model (included)

- Python v3.8+

## Input files

- FASTA sequence file

- VCF file

- Bowtie2 index of FASTA files (see below for instructions)

- List of off-target sites (see below for instructions)

## Installation & Usage

1. Clone or [download](https://github.com/bmds-lab/Crackling-phased/archive/master.zip) the source.

    ```bash
    git clone https://github.com/bmds-lab/Crackling-phased.git ~/Crackling-phased/
    cd ~/Crackling-phased
    ```

2. Be sure that your current working directory is `Crackling-phased`, containing the `setup.py` file. Install using pip:

    ```bash
    python3.9 -m pip install .
    ```

2. Configure the pipeline. See `config.ini`.

4. Ensure Bowtie2 and RNAfold are reachable system-wide, by adding them to your environments *PATH* variable.

    Check these are reachable by typing (the version numbers and directories may differ slightly):

    ```bash
    $ bowtie2 --version
    /home/<user>/bowtie2-2.3.4.1/bowtie2-align-s version 2.3.4.1
    64-bit
    Built on UbuntuDesktopMachine
    Monday 25 June  09:17:27 AEST 2018
    Compiler: gcc version 5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.9)
    Options: -O3 -m64 -msse2 -funroll-loops -g3 -std=c++98 -DPOPCNT_CAPABILITY
    Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
    ```
    
    ```bash
    $ RNAfold --version
    RNAfold 2.4.14
    ```

5. Compile the off-target indexing and scoring functions. An index of off-targets is required: to prepare this, read in the *Utilities* section (*Off-target Indexing*).

    ```bash
    make
    ```

5. Create a Bowtie2 index

    The Bowtie2 manual can be found [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).
    
    Our recommended usage:
    
    ```bash
    bowtie2-build --threads 128 input-file output-file
    ```
    
    For example:
    
    ```bash
    bowtie2-build --threads 128 ~/genomes/mouse.fa ~/genomes/mouse.fa.bowtie2
    ```
    
    Bowtie2 produces multiple files for its index. When referring to the index, use the base-name (i.e. `output-file`) that you provided `bowtie2-build`.
    
5. Configure the Crackling pipeline by editing `config.ini`.

    You will most likely need to update:

       - `general` -> `name`

       - `input` -> `exon-sequences`

       - `input` -> `vcf-file`

       - `input` -> `offtarget-sites`

       - `input` -> `bowtie2-index`

       - `output` -> `dir`

       - `sgrnascorer2` -> `model`

5. Run the pipeline: 

    ```bash
    haplocrackling -c config.ini
    ```

# Utilities

The package provides a number of utilities:

- Off-target indexing (including extracting target sites and generating the ISSL index)
- Counting targeted transcripts per guide RNA
- Retraining the provided sgRNAScorer 2.0 model (if needed)

## Off-target Indexing

This utility does not use `config.ini` but rather, command-line arguments are used to specify inputs and outputs.

1. Extract off-target sites:

   ```bash
    extractOfftargets [-h] --vcf VCF [--maxOpenFiles MAXOPENFILES] [--threads THREADS] output inputs [inputs ...]
   ```

   For example:

   ```bash
   extractOfftargets --vcf ~/genomes/mouse.vcf ~/genomes/mouse_offtargets.txt ~/genomes/mouse.fa
   ```

   The input provided can be:

   - A single, or a space sperated list, of multi-FASTA formatted files

   - A directory, for which we scan every file by parsing, using [glob](https://docs.python.org/3/library/glob.html): `<input-dir>/*`

2. Generate the index:

   ```bash
   $ bin/isslCreateIndex
   Usage: bin/isslCreateIndex [offtargetSites.txt] [sequence length] [slice width (bits)] [sissltable]
   ```

   For example:

   *For a 20bp sgRNA where up to four mismatches are allowed, use a slice width of eight (4 mismatches \* 2 bits per mismatch)*

   ```bash
   bin/isslCreateIndex ~/genomes/mouse_offtargets.txt 20 8 ~/genomes/mouse_offtargets.txt.issl
   ```

   A progress indicator is printed to *stderr*, formatted as `<current line of input file> / <number of lines in input file> : <running total of distinct sites>`.

## Counting targeted transcripts per guide RNA

This utility is useful after you have ran HaploCrackling.

Using the CLI command `countHitTranscripts`:

```bash
usage: countHitTranscripts [-h] [-a ANNOTATION] [-c CRACKLING] [-o OUTPUT]
                           [-s]

optional arguments:
  -h, --help            show this help message and exit
  -s, --sample          Run sample

group:
  -a ANNOTATION, --annotation ANNOTATION
                        The GFF3 annotation file
  -c CRACKLING, --crackling CRACKLING
                        The HaploCrackling output file
  -o OUTPUT, --output OUTPUT
                        The output file
```

For example, two guides, *A* and *B*, have been selected by Crackling as safe and efficient. How many transcripts of a gene do each guide target?

Exons are presented by `|||||`.

    Chromosome 1:
    
     (Target A)    (Target B)  (Target C)    (Target D)
    		 *             *           *             *                           
       ----||*|||-------|||*|||------||*|||----------*--- (Gene 1 - Transcript 1)
       ----||*|||----------*---------||*|||----------*--- (Gene 1 - Transcript 2)
       ------*----------|||*|||------||*|||----------*--- (Gene 1 - Transcript 3)
       ------*-----------------------||*|||----------*--- (Gene 1 - Transcript 4)
    		 *             *           *             *      

Use `--sample` to run the utility for the example above:

```bash
$ countHitTranscripts --sample
Writing test data to file.
The expected results from the test are:
AAAA 2/4
AAAT 2/4
AATA 4/4
ATAA 0/0

Pickled to: /tmp/tmp68qd5n6y.p
['seq', 'bowtieChr', 'bowtieStart', 'bowtieEnd', 'hits']
['AAAA', 'Chr1', '60', '83', '2/4']
['AAAT', 'Chr1', '200', '223', '2/4']
['AATA', 'Chr1', '320', '343', '4/4']
['ATAA', 'Chr1', '460', '483', '0/0']
```

## Training the sgRNAScorer 2.0 model (if needed)

We have provided a pre-trained model, however, dependent on your environment (Python and package versions), you may need to retrain it, using the CLI command `trainModel`. The model is saved to file and may break when newer versions of sklearn try to load it. All arguments to this command are optional, as the utility will compute the default values for you.

You can run the command without any arguments:

```bash
trainModel
```

or, by specifying arguments:

```bash
Using user specified arguments
usage: trainModel [-h] -g GOOD -b BAD -s SPACERLENGTH -p PAMORIENTATION -l
                  PAMLENGTH -o SVMOUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -g GOOD, --good GOOD
  -b BAD, --bad BAD
  -s SPACERLENGTH, --spacerLength SPACERLENGTH
  -p PAMORIENTATION, --pamOrientation PAMORIENTATION
  -l PAMLENGTH, --pamLength PAMLENGTH
  -o SVMOUTPUT, --svmOutput SVMOUTPUT
```

## References

Ben Langmead and Steven L Salzberg. Fast gapped-read alignment with Bowtie2. Nature Methods, 9(4):357, 2012.

Bradford, J., & Perrin, D. (2019). A benchmark of computational CRISPR-Cas9 guide design methods. PLoS Computational Biology, 15(8), e1007274.

Bradford, J., & Perrin, D. (2019). Improving CRISPR guide design with consensus approaches. BMC Genomics, 20(9), 931.

Bradford, J., Chappell, T., & Perrin, D. (2022). Rapid whole-genome identification of high quality CRISPR guide RNAs with the Crackling method. The CRISPR Journal, 5(3), 410-421.

Chari, R., Yeo, N. C., Chavez, A., & Church, G. M. (2017). sgRNA Scorer 2.0: a species-independent model to predict CRISPR/Cas9 activity. ACS synthetic biology, 6(5), 902-904.

Lorenz, R., Bernhart, S. H., Zu  Siederdissen, C. H., Tafer, H., Flamm, C., Stadler, P. F., &  Hofacker, I. L. (2011). ViennaRNA Package 2.0. *Algorithms for molecular biology*, *6*(1), 1-14.

Montague, T. G., Cruz, J. M., Gagnon, J. A., Church, G. M., & Valen, E. (2014). CHOPCHOP: a CRISPR/Cas9 and TALEN web tool for genome editing. Nucleic acids research, 42(W1), W401-W407.

Sunagawa, G. A., Sumiyama, K., Ukai-Tadenuma, M., Perrin, D., Fujishima, H., Ukai, H., ... & Shimizu, Y. (2016). Mammalian reverse genetics without crossing reveals Nr3a as a short-sleeper gene. Cell reports, 14(3), 662-677.

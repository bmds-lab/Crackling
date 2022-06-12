# Crackling

**Rapid Whole-Genome Identification of High Quality CRISPR Guide RNAs with the Crackling Method**

Jacob Bradford, Timothy Chappell, and Dimitri Perrin. The CRISPR Journal. Jun 2022.410-421. http://doi.org/10.1089/crispr.2021.0102

## Preamble

> The design of CRISPR-Cas9 guide RNAs is not trivial and is a computationally demanding task. Design tools need to identify target sequences that will maximize the likelihood of obtaining the desired cut, while minimizing off-target risk. There is a need for a tool that can meet both objectives while remaining practical to use on large genomes. 
> 
> In this study, we present Crackling, a new method that is more suitable for meeting these objectives. We test its performance on 12 genomes and on data from validation studies. Crackling maximizes guide efficiency by combining multiple scoring approaches. On experimental data, the guides it selects are better than those selected by others. It also incorporates Inverted Signature Slice Lists (ISSL) for faster off-target scoring. ISSL provides a gain of an order of magnitude in speed compared with other popular tools, such as Cas-OFFinder, Crisflash, and FlashFry, while preserving the same level of accuracy. Overall, this makes Crackling a faster and better method to design guide RNAs at scale. 
>
> Crackling is available at https://github.com/bmds-lab/Crackling under the Berkeley Software Distribution (BSD) 3-Clause license.

## Dependencies

- ISSL-based search off-target sites (included)

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)

- sgRNAScorer 2.0 model (included)

- Python v3.6+

## Installation & Usage

1. Clone or [download](https://github.com/bmds-lab/Crackling/archive/master.zip) the source.

    ```bash
    git clone https://github.com/bmds-lab/Crackling.git ~/Crackling/
    cd ~/Crackling
    ```

2. Install using pip

    ```bash
    python3.6 -m pip install -e .
    ```

    Important: the dot `.` indicates that *pip* will run `setup.py` from the current working directory.

    The `-e` flag is for *editable*,

    > -e	Install a project in editable mode (i.e. setuptools "develop mode") from a local project path or a VCS url.

2. Configure the pipeline. See `config.ini`.

4. Ensure Bowtie2 and RNAfold are reachable system-wide, by adding them to your environments *PATH* variable.

    Check these are reachable by typing (the version numbers and directories may differ slightly):

    ```
    $ bowtie2 --version
    /home/<user>/bowtie2-2.3.4.1/bowtie2-align-s version 2.3.4.1
    64-bit
    Built on UbuntuDesktopMachine
    Monday 25 June  09:17:27 AEST 2018
    Compiler: gcc version 5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.9)
    Options: -O3 -m64 -msse2 -funroll-loops -g3 -std=c++98 -DPOPCNT_CAPABILITY
    Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
    
    
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
    
    ```
    bowtie2-build --threads 128 input-file output-file
    ```
    
    For example:
    
    ```bash
    bowtie2-build --threads 128 ~/genomes/mouse.fa ~/genomes/mouse.fa.bowtie2
    ```
    
    Bowtie2 produces multiple files for its index. When referring to the index, use the base-name (i.e. `output-file`) that you provided `bowtie2-build`.
    
5. Configure the Crackling pipeline by editing `config.ini`.

5. Run the pipeline: 

    ```bash
    Crackling -c config.ini
    ```

# Utilities

The Crackling package provides a number of utilities:

- Off-target indexing (including extracting target sites and generating the ISSL index)
- Counting targeted transcripts per guide RNA
- Retraining the provided sgRNAScorer 2.0 model (if needed)

## Off-target Indexing

1. Extract off-target sites:

   ```bash
   extractOfftargets <output-file>  {<input-files>... | input-dir>}
   ```

   For example:

   ```
   extractOfftargets ~/genomes/mouse_offtargets.txt ~/genomes/mouse.fa
   ```

   The input provided can be:

   - A single, or a space sperated list, of multi-FASTA formatted files

   - A directory, for which we scan every file by parsing, using [glob](https://docs.python.org/3/library/glob.html): `<input-dir>/*`

   Note: Unlike previous versions, sorting the extracted off-targets is no longer required as extractOfftargets.py completes this automatically now.

2. Generate the index:

   ```
   usage: createIsslIndex [-h] -t OFFTARGETS -l GUIDELENGTH -w SLIDEWIDTH -o
                          OUTPUT [-b BINARY]
   
   optional arguments:
     -h, --help            show this help message and exit
     -t OFFTARGETS, --offtargets OFFTARGETS
                           A text file containing off-target sites
     -l GUIDELENGTH, --guidelength GUIDELENGTH
                           The length of an off-target site
     -w SLIDEWIDTH, --slidewidth SLIDEWIDTH
                           The ISSL slice width in bits
     -o OUTPUT, --output OUTPUT
                           A filepath to save the ISSL index
     -b BINARY, --binary BINARY
                           A filepath to the createIsslIndex binary (optional)
   ```

   For example:

   *For a 20bp sgRNA where up to four mismatches are allowed, use a slice width of eight (4 mismatches \* 2 bits per mismatch)*

   ```
   createIsslIndex -t ~/genomes/mouse_offtargets.txt -l 20 -w 8 - o ~/genomes/mouse_offtargets-sorted.txt.issl
   ```

   A progress indicator is printed to *stderr*, like so:

   > 8576/8583 : 6548
   >
   > 8577/8583 : 6549
   >
   > 8578/8583 : 6549
   >
   > 8579/8583 : 6549
   >
   > 8580/8583 : 6549
   >
   > 8581/8583 : 6549
   >
   > 8582/8583 : 6549
   >
   > 8583/8583 : 6550

   formatted as `<current line of input file> / <number of lines in input file> : <running total of distinct sites>`.

   This is indicating that the 6549'th distinct site has been seen on lines 8577 through 8582.

   The indicator is provided for every 10,000 input lines that are processed, and for every of the last 100 input lines.


## Counting targeted transcripts per guide RNA

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
                        The Crackling output file
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

We provided a pre-trained model, however, dependent on your environment (Python and package versions), you may need to retrain it, using the CLI command `trainModel`. All arguments to this command are optional, as the utility will compute the default values for you.

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

Bradford, J., & Perrin, D. (2019). A benchmark of computational CRISPR-Cas9 guide design methods. PLoS computational biology, 15(8), e1007274.

Bradford, J., & Perrin, D. (2019). Improving CRISPR guide design with consensus approaches. BMC genomics, 20(9), 931.

Chari, R., Yeo, N. C., Chavez, A., & Church, G. M. (2017). sgRNA Scorer 2.0: a species-independent model to predict CRISPR/Cas9 activity. ACS synthetic biology, 6(5), 902-904.

Lorenz, R., Bernhart, S. H., Zu  Siederdissen, C. H., Tafer, H., Flamm, C., Stadler, P. F., &  Hofacker, I. L. (2011). ViennaRNA Package 2.0. *Algorithms for molecular biology*, *6*(1), 1-14.

Montague, T. G., Cruz, J. M., Gagnon, J. A., Church, G. M., & Valen, E. (2014). CHOPCHOP: a CRISPR/Cas9 and TALEN web tool for genome editing. Nucleic acids research, 42(W1), W401-W407.

Sunagawa, G. A., Sumiyama, K., Ukai-Tadenuma, M., Perrin, D., Fujishima, H., Ukai, H., ... & Shimizu, Y. (2016). Mammalian reverse genetics without crossing reveals Nr3a as a short-sleeper gene. Cell reports, 14(3), 662-677.

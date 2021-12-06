# Crackling

**Faster and better CRISPR guide RNA design with the Crackling method**

Jacob Bradford, Timothy Chappell, Dimitri Perrin

bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261

## Preamble

> The design of CRISPR-Cas9 guide RNAs is not trivial, and is a computationally demanding task. Design tools need to identify target sequences that will maximise the likelihood of obtaining the desired cut, whilst minimising off-target risk. There is a need for a tool that can meet both objectives while remaining practical to use on large genomes.
>
> Here, we present Crackling, a new method that is more suitable for meeting these objectives. We test its performance on 12 genomes and on data from validation studies. Crackling maximises guide efficiency by combining multiple scoring approaches. On experimental data, the guides it selects are better than those selected by others. It also incorporates Inverted Signature Slice Lists (ISSL) for faster off-target scoring. ISSL provides a gain of an order of magnitude in speed compared to other popular tools, such as Cas-OFFinder, Crisflash and FlashFry, while preserving the same level of accuracy. Overall, this makes Crackling a faster and better method to design guide RNAs at scale.
>
> Crackling is available at https://github.com/bmds-lab/Crackling under the Berkeley Software Distribution (BSD) 3-Clause license. 

## Dependencies

- ISSL-based search off-target sites (included)

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)

- sgRNAScorer 2.0 model (included)

- Python v3.6+

## Installation

1. Clone or [download](https://github.com/bmds-lab/Crackling/archive/master.zip) the source.

    ```
    git clone https://github.com/bmds-lab/Crackling.git ~/Crackling/
    cd ~/Crackling
    ```

2. Install using pip

    ```
    python3.6 -m pip install -e ~/Crackling/
    ```

    The `-e` flag is for *editable*,

    > -e	Install a project in editable mode (i.e. setuptools "develop mode") from a local project path or a VCS url.

2. Configure the pipeline. See `config.ini`.

3. Ensure Bowtie2 and RNAFold are reachable system-wide, by adding them to your environments *PATH* variable.

5. Compile the off-target indexing and scoring functions. An index of off-targets is required: to prepare this, read the next section (*Off-target Indexing*).

    ```
    make
    ```

    or

    ```
    g++ -o ./bin/isslScoreOfftargets ./src/isslScoreOfftargets.cpp -O3 -std=c++11 -fopenmp -mpopcnt -Iparallel_hashmap
    ```

5. Run the pipeline: 

    ```
    python Crackling.py -c config.ini
    ```

## Off-target Indexing

1. Extract off-target sites:

    ```
    python bin/extractOfftargets.py <output-file>  {<input-files>... | input-dir>}
    ```

    For example:

    ```
    python bin/extractOfftargets.py ~/genomes/mouse_offtargets.txt ~/genomes/mouse.fa
    ```

   The input provided can be:

   - A single, or a space sperated list, of multi-FASTA formatted files

   - A directory, for which we scan every file by parsing, using [glob](https://docs.python.org/3/library/glob.html): `<input-dir>/*`

   Note: Unlike previous versions, sorting the extracted off-targets is no longer required as extractOfftargets.py completes this automatically now.



2. Compile the ISSL indexer

    Compile the indexer first: 
    
    ```
	make
	```
	
	or 
	
	```
    g++ -o ./bin/isslCreateIndex ./src/isslCreateIndex.cpp -O3 -std=c++11 -fopenmp -mpopcnt
    ```
    
3. Generate the index:
   
    ```
    ./bin/isslCreateIndex <offtargets-sorted> <guide-length> <slice-width-bits> <index-name>
    ```
    
    For example:
    
    *For a 20bp sgRNA where up to four mismatches are allowed, use a slice width of eight (4 mismatches \* 2 bits per mismatch)*
    
    ```
    ./bin/isslCreateIndex ~/genomes/mouse_offtargets-sorted.txt 20 8 ~/genomes/mouse_offtargets-sorted.txt.issl
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

## Bowtie2 index

The Bowtie2 manual can be found [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

Crackling requires a Bowtie2 index to be provided.

Our recommended usage:

```
bowtie2-build --threads 128 input-file output-file
```

For example:

```bash
bowtie2-build --threads 128 ~/genomes/mouse.fa ~/genomes/mouse.fa.bowtie2
```

Bowtie2 produces multiple files for its index. When referring to the index, use the base-name (i.e. `output-file`) that you provided `bowtie2-build`.

## References

Ben Langmead and Steven L Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4):357, 2012.

Bradford, J., & Perrin, D. (2019). A benchmark of computational CRISPR-Cas9 guide design methods. PLoS computational biology, 15(8), e1007274.

Bradford, J., & Perrin, D. (2019). Improving CRISPR guide design with consensus approaches. BMC genomics, 20(9), 931.

Chari, R., Yeo, N. C., Chavez, A., & Church, G. M. (2017). sgRNA Scorer 2.0: a species-independent model to predict CRISPR/Cas9 activity. ACS synthetic biology, 6(5), 902-904.

Montague, T. G., Cruz, J. M., Gagnon, J. A., Church, G. M., & Valen, E. (2014). CHOPCHOP: a CRISPR/Cas9 and TALEN web tool for genome editing. Nucleic acids research, 42(W1), W401-W407.

Sunagawa, G. A., Sumiyama, K., Ukai-Tadenuma, M., Perrin, D., Fujishima, H., Ukai, H., ... & Shimizu, Y. (2016). Mammalian reverse genetics without crossing reveals Nr3a as a short-sleeper gene. Cell reports, 14(3), 662-677.

# Crackling

**Faster and better CRISPR guide RNA design with the Crackling method**

Jacob Bradford, Timothy Chappell, Dimitri Perrin

bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261

## Preamble

We present Crackling, a new method for whole-genome identification of suitable CRISPR targets. The method maximises the efficiency of the guides by combining the results of multiple scoring approaches. On experimental data, the set of guides it selects are better than those produced by existing tools. The method also incorporates a new approach for faster off-target scoring, based on Inverted Signature Slice Lists (ISSL). This approach provides a gain of an order of magnitude in speed, while preserving the same level of accuracy.

## Dependencies

- ISSL-based search off-target sites (included)

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)

- sgRNAScorer 2.0 model (included)

## Installation

1. Clone or [download](https://github.com/bmds-lab/Crackling/archive/master.zip) the source.

2. Configure the pipeline. See `config.py`.

3. Ensure Bowtie2 and RNAFold are reachable from the installation directory.

4. Compile the off-target scoring function. An index of off-targets is required: to prepare this, read the next section (*Off-target Indexing*).

    ```
    g++ -o search_ots_score search_ots_score.cpp -O3 -std=c++11 -fopenmp -mpopcnt
    ```

5. Run the pipeline: 

    ```
    python3 process.py -c config
    ```

## Off-target Indexing

1. Prepare a file for your genome which contains one line per chromosome.

2. Extract off-target sites:

    ```
    python prepareListOfftargetSites.py <input-file> <output-file>
    ```
    
    For example:
    
    ```
    python prepareListOfftargetSites.py ~/genomes/mouse.txt ~/genomes/mouse_offtargets.txt
    ```

3. Sort the off-target sites. 

    On Linux:
    
    ```
    sort --parallel=64 ~/genomes/mouse_offtargets.txt > ~/genomes/mouse_offtargets-sorted.txt
    ```

4. Build the ISSL index

    Compile the indexer first: 
    
    ```
    g++ -o index_offtargetSites index_offtargetSites.cpp -O3 -std=c++11 -fopenmp -mpopcnt
    ```
    
    Generate the index:
    
    *For a 20bp sgRNA where up to four mismatches are allowed, use a slice width of eight*
    
    ```
    ./index_offtargetSites <offtargets-sorted> <guide-length> <slice-width-bits> <index-name>
    ```
    
    For example:
    
    ```
    ./index_offtargetSites ~/genomes/mouse_offtargets-sorted.txt 20 8 ~/genomes/mouse_offtargets-sorted.txt.issl
    ```


## References

Ben Langmead and Steven L Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4):357, 2012.

Bradford, J., & Perrin, D. (2019). A benchmark of computational CRISPR-Cas9 guide design methods. PLoS computational biology, 15(8), e1007274.

Bradford, J., & Perrin, D. (2019). Improving CRISPR guide design with consensus approaches. BMC genomics, 20(9), 931.

Chari, R., Yeo, N. C., Chavez, A., & Church, G. M. (2017). sgRNA Scorer 2.0: a species-independent model to predict CRISPR/Cas9 activity. ACS synthetic biology, 6(5), 902-904.

Montague, T. G., Cruz, J. M., Gagnon, J. A., Church, G. M., & Valen, E. (2014). CHOPCHOP: a CRISPR/Cas9 and TALEN web tool for genome editing. Nucleic acids research, 42(W1), W401-W407.

Sunagawa, G. A., Sumiyama, K., Ukai-Tadenuma, M., Perrin, D., Fujishima, H., Ukai, H., ... & Shimizu, Y. (2016). Mammalian reverse genetics without crossing reveals Nr3a as a short-sleeper gene. Cell reports, 14(3), 662-677.
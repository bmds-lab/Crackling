# Crackling

*CRISPR, faster, better – The Crackling method for whole-genome target detection*

Jacob Bradford, Dimitri Perrin

## Preamble

We present Crackling, a new method for whole-genome identification of suitable CRISPR targets. The method maximises the efficiency of the guides by combining the results of multiple scoring approaches. On experimental data, the set of guides it selects are better than those produced by existing tools. The method also incorporates a new approach for faster off-target scoring, based on Inverted Signature Slice Lists (ISSL). This approach provides a gain of an order of magnitude in speed, while preserving the same level of accuracy.

## Dependencies

- ISSL-based search off-target sites (included)

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)

- [sgRNAScorer 2.0 model (included)]()

## Installation

1. Clone or [download](https://github.com/bmds-lab/Crackling/archive/master.zip) the source.

2. Configure the pipeline. See `config.py`.

3. Ensure Bowtie2 and RNAFold are reachable from the installation directory.

4. Compile the off-target scoring function

```
g++ -o search_ots_score search_ots_score.cpp -std=c++11 -fopenmp -mpopcnt
```

5. Run the pipeline: `python3 process.py -c config`

## References

Ben Langmead and Steven L Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4):357, 2012.

Bradford, J., & Perrin, D. (2019). A benchmark of computational CRISPR-Cas9 guide design methods. PLoS computational biology, 15(8), e1007274.

Bradford, J., & Perrin, D. (2019). Improving CRISPR guide design with consensus approaches. BMC genomics, 20(9), 931.

Chari, R., Yeo, N. C., Chavez, A., & Church, G. M. (2017). sgRNA Scorer 2.0: a species-independent model to predict CRISPR/Cas9 activity. ACS synthetic biology, 6(5), 902-904.
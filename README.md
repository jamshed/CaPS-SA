# themis

themis is a simple parallel suffix array construction algorithm

# installation

to compare correctness against Shun et. al.'s 2014 parallel algorithm, install via:

```
tar xzvf suffixTree.tar
cd suffixTree.tar
make suffixArray
```

then compare the programs via:

```
cd utils
bash test_correctness.sh <N>
```

This will generate a random DNA string of length of length N and run both programs on it, comparing their output and running time. N=1000 by default.

# Themis

Themis is a simple and parallel Suffix Array and LCP array construction algorithm.

## Installation

- From source:

  ```bash
  git clone https://github.com/jamshed/themis.git
  cd themis/
  mkdir build && cd build/
  cmake -DCMAKE_INSTALL_PREFIX=../ ..
  make install
  cd ..
  ```

  This installs Themis in a sub-directory named `bin`, inside the project root directory.

## Usage

```bash
themis <input_file> <output_file> <(optional)-subproblem-count>
```

## Miscellaneous

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

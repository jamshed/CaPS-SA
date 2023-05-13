# CaPS-SA

CaPS-SA is a simple, parallel, and cache-friendly Suffix Array and LCP array construction algorithm.

## Installation

From source:

```bash
git clone https://github.com/jamshed/CaPS-SA.git
cd CaPS-SA/
mkdir build && cd build/
cmake -DCMAKE_INSTALL_PREFIX=../ ..
make install
cd ..
```

This installs `caps_sa` in a sub-directory named `bin`, inside the project root directory.
If the suffix and the LCP arrays are to be constructed for texts with sizes $\geq 2^{32}$, then please replace the `cmake` command in the installation instructions with the following one: `cmake -DCMAKE_INSTALL_PREFIX=../ -DLARGE_IDX=ON` .

## Usage

  ```bash
  export PARLAY_NUM_THREADS=<thread-count>
  caps_sa <input_file> <output_file> <(optional)-subproblem-count>
  ```

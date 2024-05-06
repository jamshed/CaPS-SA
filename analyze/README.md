## launching runs

### adversarial

to launch a set of adversarial runs:
```
python launch.py --vary <var>
```

where `<var>` is {`wc` or `wl`}, this'll vary the word count / word length.

by default, we'll vary the word counts from 2^1 - 2^8 and then the word length from 2^1 - 2^10

you can change these values by adding the `--powers` flag, e.g.

```
python launch.py --vary wc --powers 8 10
```

will vary the word counts from 2^1 - 2^8 for a total of 8 runs and have word length set to 2^10 for each run. default gen length is 2^27, if you want to change that, then add the `--gen_len <v>` flag and type in the actual value, e.g. `--gen_len 1e9` or `--gen_len 2**20` or `--gen_len 123456789`,

so an adversarial launch script using all params would look like:

```
python launch.py --vary wc --powers 8 10 --gen_len 2**27
```
(these are the values by default)

---

### uniform

to launch a set of uniform runs, just do

```
python launch.py --gen unif
```

and add the `--powers <lower> <upper>` flag to specify the lower and upper powers of 2 to test for the *generated text length*, e.g.

```
python launch.py --gen unif --powers 25 30
```
will test gen length from 2**25 to 2**30 (inclusive), and by default it tests from 2^range(18,27) (inclusive) for a total of 10 runs

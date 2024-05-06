import os, subprocess, argparse

wc_pow = 8
wl_pow = 10

num_procs = 128

parser = argparse.ArgumentParser()
parser.add_argument('--vary', type=str, default='wc')
parser.add_argument('--gen', type=str, default='adv') # adv(ersarial) or unif(orm)
parser.add_argument('--powers', nargs='+', type=int, default=None)
parser.add_argument('--gen_len', type=int, default=2**27)
args = parser.parse_args()

if args.vary not in ['wc', 'wl']:
    print(f"--vary has to be one of (wc, wl), not {args.vary}")
    exit()
if args.powers and len(args.powers) != 2:
    print('need to have 2 ints for the --powers range')
    exit()
if args.powers:
    wc_pow, wl_pow = args.powers[0], args.powers[1]

zaratan_template = """#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -A bhatele-lab-cmsc
#SBATCH -N 1
#SBATCH --ntasks-per-node={0}
#SBATCH -o {3}

export PARLAY_NUM_THREADS={0}

../bin/caps_sa {1} {2}
"""

# if uniform input
if args.gen.lower() == 'unif':
    for power in range(args.powers[0], args.powers[1]+1):
        gen_length = 2**power
        input_file = f"unif-{power}.out"
        SA_file = "SA-" + input_file
        times_file = "times-" + input_file

        # generate adversarial input
        subprocess.run(['../gen/unif_gen', input_file, str(gen_length)])

        # run the job, capture output
        bash_text = zaratan_template.format(num_procs, input_file, SA_file, times_file)

        with open(f'run{power}.sh', 'w') as f:
            f.write(bash_text)
        subprocess.run(['sbatch', f'run{power}.sh'])
    exit('done launching uniform runs')

# else adversarial input
# fix word count or word length, and vary the other.
rango, wc, wl = [0]*3
if args.vary == 'wc':
    rango = wc_pow
else:
    rango = wl_pow

for i in range(rango):
    if args.vary == 'wc':
        wc = 2**(i+1)
        wl = 2**wl_pow
    else:
        wc = 2**wc_pow
        wl = 2**(i+1)

    input_file = f"adv-{args.gen_len}-{wc}-{wl}.out"
    SA_file = "SA-" + input_file
    times_file = "times-" + input_file

    # generate adversarial input
    subprocess.run(['../gen/adv_gen', input_file, str(args.gen_len), str(wc), str(wl)])
        
    # run the job, capture output
    bash_text = zaratan_template.format(num_procs, input_file, SA_file, times_file)
        
    with open(f'run{i}.sh', 'w') as f:
        f.write(bash_text)
    subprocess.run(['sbatch', f'run{i}.sh'])

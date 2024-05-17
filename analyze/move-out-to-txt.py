import glob,os,subprocess, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-dir', '-d', type=str, default='')
args = parser.parse_args()

for f in glob.glob(os.path.join(args.dir, 'times*.out')):
    a,b = f.split('.')
    subprocess.run(['mv', f, f"{a}.txt"])

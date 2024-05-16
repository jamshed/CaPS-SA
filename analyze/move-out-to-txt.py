import glob,os,subprocess

for f in glob.glob('times*.out'):
    a,b = f.split('.')
    subprocess.run(['mv', f, f"{a}.txt"])

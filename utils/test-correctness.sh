#/bin/bash

N=1000  # specify the default size of the random sequence
# check if user provided a value for N
if [ ! -z "$1" ]; then
    N=$1
fi

#TRUE_PROGRAM="../simple_baseline"
TRUE_PROGRAM="../suffixTree/suffixArray"
TEST_PROGRAM="../simple_baseline"

# Step 1: generate a random string of length N
seed=$RANDOM  # generate a random seed
string=$(python3 gen_rand_seq.py $seed $N)

# Step 2: write the random string to a temporary file
tempfile=$(mktemp)
echo $string > $tempfile

# Step 3: run both programs on the temporary file

tempout1=$(mktemp)
tempout2=$(mktemp)
program1_time=$({ time $TRUE_PROGRAM $tempfile $tempout1 &> /dev/null; } 2>&1 1>/dev/null)
program2_time=$({ time $TEST_PROGRAM $tempfile $tempout2 &> /dev/null; } 2>&1 1>/dev/null)

program1_output=$(cat $tempout1)
program2_output=$(cat $tempout2)
#IFS has to be empty for newlines to stick around...
IFS=
#echo $program1_output
#echo $program2_output

# Step 4: verify that both programs have the same output
if [ "$program1_output" = "$program2_output" ]; then
  result="\033[32mcorrect\033[0m"
else
  result="\033[31mincorrect\033[0m"
fi

# Step 5: delete the temporary file
rm $tempfile

# Step 6: write to stdout the random seed used and whether the output was correct or not
echo "Random seed: $seed"
echo "True program runtime: $program1_time seconds"
echo "Test program runtime: $program2_time seconds"
echo -e "Output correctness: $result"


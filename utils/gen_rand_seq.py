import sys
import random

# Get the random seed and sequence length from command-line arguments
seed = int(sys.argv[1])
N = int(sys.argv[2])

# Set the random seed
random.seed(seed)

# Generate a random DNA sequence of length N
bases = ['A', 'C', 'G', 'T']
seq = ''.join(random.choice(bases) for _ in range(N))

# Print the sequence
print(seq)


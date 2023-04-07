import sys

def suffix_array(s):
    """Compute the suffix array of a string."""
    n = len(s)
    sa = list(range(n))
    sa.sort(key=lambda i: s[i:])
    return sa

def lcp_array(s, sa):
    """Compute the longest common prefix array of a string and its suffix array."""
    n = len(s)
    rank = [0] * n
    for i in range(n):
        rank[sa[i]] = i
    lcp = [0] * (n-1)
    h = 0
    for i in range(n):
        if rank[i] == 0:
            continue
        j = sa[rank[i]-1]
        while i+h < n and j+h < n and s[i+h] == s[j+h]:
            h += 1
        lcp[rank[i]-1] = h
        if h > 0:
            h -= 1
    return lcp


def main(argv):
    inp = argv[1]
    out = argv[2]
    with open(inp, 'r') as f:
        s = f.read().strip()

    sa = suffix_array(s)
    lcp = lcp_array(s, sa)

    with open(out,'w') as f:
        f.write(" ".join([str(x) for x in sa])+"\n")
        f.write("- " + " ".join([str(x) for x in lcp]))
if __name__ == "__main__":
    main(sys.argv)

import sys

# From ChatGPT Mar 23 Version.

def suffix_array(s):
    """Generate the suffix array of a string"""
    n = len(s)
    suffixes = [(s[i:], i) for i in range(n)]
    suffixes.sort()
    return [suf[1] for suf in suffixes]

def lcp_array(s, sa):
    """Generate the LCP array of a string given its suffix array"""
    n = len(s)
    rank = [0] * n
    lcp = [0] * n
    for i in range(n):
        rank[sa[i]] = i
    h = 0
    for i in range(n):
        if rank[i] > 0:
            j = sa[rank[i]-1]
            while i+h < n and j+h < n and s[i+h] == s[j+h]:
                h += 1
            lcp[rank[i]] = h
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
        f.write(" ".join([str(x) for x in lcp])+"\n")
if __name__ == "__main__":
    main(sys.argv)

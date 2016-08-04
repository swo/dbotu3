def kmer_dict(nt, k):
    dat = {}
    for kmer in [nt[i: i + k] for i in range(len(nt) - k)]:
        if kmer in dat:
            dat[kmer] += 1
        else:
            dat[kmer] = 1

    return dat

def distance(kmers1, kmers2):
    assert isinstance(kmers1, dict)
    assert isinstance(kmers2, dict)
    keys1 = set(kmers1.keys())
    keys2 = set(kmers2.keys())
    common_keys = keys1 & keys2
    unique_keys1 = keys1 - keys2
    unique_keys2 = keys2 - keys1
    out = sum([(kmers1[x] - kmers2[x]) ** 2 for x in common_keys])
    out += sum([kmers1[x] ** 2 for x in unique_keys1])
    out += sum([kmers2[x] ** 2 for x in unique_keys2])
    return out

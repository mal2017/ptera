import pandas as pd
import sys


"""

This script ripped from https://github.com/COMBINE-lab/terminus/blob/master/scripts/python/extract_txp_group.py

terminus commit e16fc4f

"""


def main():
    salmon_quant_file = sys.argv[1]
    terminus_cluster_file = sys.argv[2]
    txp_group_tsv = sys.argv[3]

    qdf_transcripts = pd.read_csv(salmon_quant_file, sep = '\t')['Name'].values
    remap = {}
    for t in qdf_transcripts:
        remap[t] = t
    with open(terminus_cluster_file) as fp:
        for line in fp:
            tokens = line.strip().split(',')
            for t in tokens[1:]:
                remap[t] = tokens[0]
    assert(len(remap) == len(qdf_transcripts))
    pd.DataFrame.from_dict(remap, orient='index').reindex(qdf_transcripts).reset_index().to_csv(
        txp_group_tsv,
        index=False,
        header = None
    )


if __name__ == '__main__':
    main()

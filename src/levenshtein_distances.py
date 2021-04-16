import logging
import multiprocessing
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from diff_match_patch import diff_match_patch

logging.basicConfig(format='%(asctime)s %(message)s',
                    level=logging.INFO,
                    datefmt="%Y-%m-%d %H:%M:%S")

df = pd.read_csv("../data/ncbi_sgene_good_unique.csv")

N = df.shape[0]     # distance matrix size
N = 100
M = np.arange(N)

accessions = df["accession"].tolist()
accessions = accessions[:N]


def levenshtein_metric(x, y):
    dmp = diff_match_patch()
    diffs = dmp.diff_main(df.iloc[int(x)]['sgene_nucleotide'],
                          df.iloc[int(y)]['sgene_nucleotide'])
    return dmp.diff_levenshtein(diffs)


def worker(index):
    if index % 10 == 0:
        # Tarpine informacija
        logging.info(f"Current index {index}")

    d = list(map(levenshtein_metric, repeat(index), list(range(index))))
    d = np.pad(d, (0, N - len(d)))
    return d


def main():

    nproc = multiprocessing.cpu_count()
    # nproc = 4
    logging.info(f"Calculating distances. Using {nproc} CPUs")

    with ProcessPoolExecutor(nproc) as pool:
        results = list(pool.map(worker, M))
        distances = np.stack(results)
    distances += distances.T

    df_dist = pd.DataFrame(distances, columns=accessions, index=accessions)
    df_dist.to_csv("distances_lev.csv", index=False)

    logging.info("Done")


if __name__ == "__main__":
    main()

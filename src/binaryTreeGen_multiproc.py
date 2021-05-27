import logging
import multiprocessing

from pathlib import Path
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from Bio import AlignIO
from Bio import Phylo
from anytree import Node, RenderTree, AsciiStyle
from anytree.exporter import DotExporter

logging.basicConfig(format='%(asctime)s %(message)s',
                    level=logging.INFO,
                    datefmt="%Y-%m-%d %H:%M:%S")

tree_filename = "../data/RAxML_bestTree.sgene_good_unique.tre"
tree = Phylo.read(tree_filename, "newick")
terminals = tree.get_terminals()

reference_name = "MN908947_3" # 2246_6593532f926e48cc68421ef20a33018c
names = [terminal.name for terminal in terminals]
# names = names[:500]

if reference_name not in names:
    names.append(reference_name)

def _calc_distance(leaf):
    """ Helper function to calculate distances between 
        given leaf and other leaves."""

    distances = np.zeros(len(names))
    leaf_index = names.index(leaf)

    if leaf_index % 10 == 0:
        # Tarpine informacija
        logging.info(f"Current leaf index {leaf_index}")

    if leaf_index == 0:
        return distances

    for index in range(leaf_index):
        distances[index] = tree.distance(leaf, names[index])

    return distances


def calc_distances(reuse_distances=True):
    """ Calculate distance matrix and save dataframe to file.
    """

    distances_path = Path("distances.csv")

    if distances_path.is_file() and reuse_distances:
        logging.info("Loading distances from file")
        return pd.read_csv(distances_path, index_col=0)

    nproc = multiprocessing.cpu_count()
    # nproc = 4
    logging.info(f"Calculating distances. Using {nproc} CPUs")

    with ProcessPoolExecutor(nproc) as pool:
        results = list(pool.map(_calc_distance, names))
        distances = np.stack(results)

    distances += distances.T

    df = pd.DataFrame(distances, columns=names, index=names)
    df.to_csv(distances_path)

    return df


def build_tree(reuse_distances=True):
    MinPair = namedtuple("MinPair", 'distance target leaf')

    distances = calc_distances(reuse_distances)
    distances = distances[distances.gt(0)]  # pavercia nulius i NaN

    # distances = distances[names]
    # print(distances.shape)
    # distances = distances.loc[names]
    # print(distances.shape)

    # Drop ref sgene row to prevent child referencing back to ref sgene
    distances.drop(reference_name, inplace=True, axis=0)

    rooted_tree = Node(reference_name)
    leaves = [rooted_tree]

    logging.info("Building tree")
    while distances.shape[0] > 0:
        min_pair = None

        for leaf in leaves: 

            # Suranda maziausia atstuma didesni uz nuli
            distance = distances[leaf.name].min()
            # ir tos sekos ID. Indexas yra seku ID
            target = distances[leaf.name].idxmin()

            if min_pair is None or min_pair.distance > distance:
                min_pair = MinPair(distance, target, leaf)

        leaves.append(Node(min_pair.target, parent=min_pair.leaf))

        # distances.drop(min_pair.target, inplace=True, axis=1)
        distances.drop(min_pair.target, inplace=True, axis=0)

        # checking whether node has 2 leaves
        if len(min_pair.leaf.children) == 2:
            leaves.remove(min_pair.leaf)

        if distances.shape[0] % 10 == 0:
            # Tarpine informacija
            logging.info(f"Remaining leaves {distances.shape[0]}")

    # print(RenderTree(rooted_tree, style=AsciiStyle()).by_attr())
    DotExporter(rooted_tree).to_dotfile("tree.dot")
    logging.info("Done")

if __name__ == "__main__":
    build_tree(reuse_distances=True)

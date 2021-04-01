from Bio import AlignIO
from Bio import Phylo
from anytree import Node, RenderTree, AsciiStyle
from operator import itemgetter

file = "data/speed_test.raxml.bestTree"
# file = "data/RAxML_bestTree.sgene_good_unique.tre"
tree = Phylo.read(file, "newick")
terminals = tree.get_terminals()

reference_name = terminals[5].name
names = [terminal.name for terminal in terminals if terminal.name != reference_name]
rooted_tree = Node(reference_name)

# Phylo.draw_ascii(tree)

def getDistances(parent, children):
    return list(map(lambda x: tree.distance(x, parent), children))

def getMin(parent, children):
    return min(
        enumerate(
            getDistances(parent, children)
        ),
        key=itemgetter(1)
    )

leafs = [rooted_tree]
while len(names):
    min_distance = None
    current_leaf = None
    current_target = None
    for leaf in leafs:
        target = getMin(leaf.name, names)
        if min_distance == None or min_distance > target[1]:
            min_distance = target[1]
            current_target = target
            current_leaf = leaf
    # adding target to leaf
    leafs.append(Node(names[current_target[0]], parent=current_leaf))
    del names[current_target[0]]
    # checking whether node has 2 leafs
    if len(current_leaf.children) == 2:
        leafs.remove(current_leaf)
print(RenderTree(rooted_tree, style=AsciiStyle()).by_attr())

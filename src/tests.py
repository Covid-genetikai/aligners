from Bio import AlignIO
from Bio import Phylo
####################################################################
# CLUSTALW2                                                        #
# Ref: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec93 #
####################################################################

# from Bio.Align.Applications import ClustalwCommandline
# cline = ClustalwCommandline("clustalw2", infile="data/genomic.fna")
# cline()

####################################################################
# MUSCLE                                                           #
# Ref: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec94 #
####################################################################

# from Bio.Align.Applications import MuscleCommandline
# muscle_cline = MuscleCommandline(input="data/genomic.fna", clwout="data/genomic.aln")
# muscle_cline()

####################################################################

def printAln(file):
    align = AlignIO.read(file, "clustal")
    for record in align:
        print("%s - %s" % (record.seq, record.id))


# tree = Phylo.read("data/seqdump.dnd", "newick")
# Phylo.draw_ascii(tree)
#!/usr/bin/env python2
'''Converting a phylogenetic tree in Newick format into input for
Hybrid-Lambda (Zhu et al. 2013, arXiv:1303.0673)'''
__author__ = "Michael Gruenstaeudl, PhD"
__copyright__ = "Copyright (C) 2014 Michael Gruenstaeudl"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.07.31.1400"
__status__ = "Testing"

#####################
# IMPORT OPERATIONS #
#####################

from collections import OrderedDict
from termcolor import colored
import argparse
import dendropy
import GeneralStringOperations as GSO
import os
import sys

###############
# DEFINITIONS #
###############


def addNodeNumb(intree):
    ''' Adding node numbers to tree.
    Args:
        intree:     a left-ladderized tree string without node labels
    Returns:
        outtree:    a left-ladderized tree string with node labels
    '''
    outtree = []
    nNodes = intree.count(")")*2
    for c in intree:
        if c == ")":
            outtree.append(c + "s" + str(nNodes))
            nNodes -= 1
        else:
            outtree.append(c)
    return ''.join(outtree)


def addHybrDict(intree, hybrDict):
    ''' Adds hybrid information to tree string in accordance with the
        requirements of hybrid-Lambda.
    Args:
        intree:     a tree string without hybrid info
        hybrDict:   a dictionary with hybrid parent info of
                    format {"B":"0.6", "C":"0.4"}

    Returns:
        outtree:    a tree string with hybrid info
    '''
    # Looping over the specified parental taxa and incorporating the hybrid
    # info for each taxon in list
    for key, val in hybrDict.iteritems():

        # 1. Replace parent's brlen with adjusted brlen
        if not intree.find(key):
            print colored("Error: Specified hybrid parent not present in \
                           input tree.", 'red')
            break
        else:
            # Parse out the branch length immediately following the key,
            # save as "brlen"; then replace said branch length with
            # half of value
            try:
                brlen = float(GSO.exstr(intree, key+":", ")"))
                intree = GSO.replstr(intree, key+":", ")", str(brlen*0.5))
            except ValueError:
                brlen = float(GSO.exstr(intree, key+":", ","))
                intree = GSO.replstr(intree, key+":", ",", str(brlen*0.5))

        # 2. Split intree into three sections by keywords
        split1 = GSO.csplit(intree, key, rightflag=True)
        aList = [split1[0]] + GSO.csplit(split1[1], key+":"+str(brlen*0.5),
                                         rightflag=False)

        # 3. Compile hybrid info in a string
        hybStr = "(h#"+val+":"+str(brlen*0.5)+","+aList[1]+"):"+str(brlen*0.5)

        # 4. Replace second intree element with hybrid string
        aList[1] = hybStr
        # Note: intree in TFL needs to be inside the loop!
        intree = ''.join(aList)

    return intree


def main(treeName, parentInfo):
    # Reading tree as string
    # treeStr = open(treeName, "r").read()
    treelines = open(treeName, "r").readlines()
    treeStrs = []
    for line in treelines:
        l = line.strip()
        if len(l) > 0:
            if l[0] != "#":
                treeStrs.append(line)

    for treeStr in treeStrs:
        # Reading tree by DendroPy
        tree = dendropy.Tree.get_from_string(treeStr, "newick")
        # DEBUGLINE: print(tree.as_ascii_plot())
        # Placeholder to potentially modify tree further

        # Left-ladderize tree
        tree.ladderize(ascending=True)
        treeStr = tree.as_string('newick')

        # Parsing parentInfo into dictionary
        aDict = {}
        for i in parentInfo.split(","):
            aDict[i.split(":")[0]] = i.split(":")[1]
        # Hypothetical content of aDict at this point:
        #    {"B":"0.6", "C":"0.4"}

        # Adding hybrid information to intree
        outtree = addHybrDict(treeStr, aDict)

        # Adding nodes to tree
        outtree = addNodeNumb(outtree)

        # Saving output to file
        outf = open(GSO.rmext(treeName)+".hybrLmbda.INPUT", "a")
        outf.write(outtree)
        outf.close()

###########
# EXECUTE #
###########

print ""
print colored("  Script name: "+sys.argv[0], 'cyan')
print colored("  Author: "+__author__, 'cyan')
print colored("  Version: "+__version__, 'cyan')
print ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converting a phylogenetic \
    tree in Newick format into input for Hybrid-Lambda (Zhu et al. 2013, \
    arXiv:1303.0673); '+__copyright__)
    parser.add_argument('-t', '--tree', help='name of input tree',
                        default="infile.tre", required=True)
    parser.add_argument('-p', '--parentinfo', help='info on parental taxa of \
    hybrids; format: <parent1>:<likelih.parent1>,<parent2>:<likelih.parent2>',
                        default="A:0.6,B:0.4", required=True)
    args = parser.parse_args()

main(args.tree, args.parentinfo)

print colored("  Done.", 'cyan')
print ""

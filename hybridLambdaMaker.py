#!/usr/bin/env python2
'''Converting a phylogenetic tree in Newick format into input for
Hybrid-Lambda (Zhu et al. 2013, arXiv:1303.0673)'''
__author__ = "Michael Gruenstaeudl, PhD"
__copyright__ = "Copyright (C) 2014 Michael Gruenstaeudl"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.09.07.2200"
__status__ = "Testing"

#####################
# IMPORT OPERATIONS #
#####################

import argparse
from collections import OrderedDict
import dendropy
import GeneralStringOperations as GSO
import os
import re
import sys
from termcolor import colored

###############
# DEFINITIONS #
###############


def addNodeNumb(inTree):
    ''' Adding node numbers to tree.
    Args:
        inTree:     a left-ladderized tree string without node labels
    Returns:
        outTree:    a left-ladderized tree string with node labels
    '''
    outTree = []
    nNodes = inTree.count(")")*2
    for c in inTree:
        if c == ")":
            outTree.append(c + "s" + str(nNodes))
            nNodes -= 1
        else:
            outTree.append(c)
    return ''.join(outTree)


def addHybrDict(inTree, hybrDict):
    ''' Adds hybrid information to tree string in accordance with the
        requirements of hybrid-Lambda.
    Args:
        inTree:     a tree string without hybrid info
        hybrDict:   a dictionary with hybrid parent info of
                    format {"B":"0.6", "C":"0.4"}

    Returns:
        outTree:    a tree string with hybrid info
    '''
    # Looping over the specified parental taxa and incorporating the hybrid
    # info for each taxon in list
    for key, val in hybrDict.iteritems():

        # 1. Replace parent's brlen with adjusted brlen
        if key not in inTree:
            print colored(" ERROR: Specified hybrid parent not present in \
input tree.\n  Quitting ...\n", 'red')
            sys.exit()
        else:
            # Parse out the branch length immediately following the key,
            # save as "brlen"; then replace said branch length with
            # half of value
            try:
                brlen = float(GSO.exstr(inTree, key+":", ")"))
                inTree = GSO.replstr(inTree, key+":", ")", str(brlen*0.5))
            except ValueError:
                brlen = float(GSO.exstr(inTree, key+":", ","))
                inTree = GSO.replstr(inTree, key+":", ",", str(brlen*0.5))

        # 2. Split inTree into three sections by keywords
        split1 = GSO.csplit(inTree, key, rightflag=True)
        aList = [split1[0]] + GSO.csplit(split1[1], key+":"+str(brlen*0.5),
                                         rightflag=False)

        # 3. Compile hybrid info in a string
        hybStr = "(h#"+val+":"+str(brlen*0.5)+","+aList[1]+"):"+str(brlen*0.5)

        # 4. Replace second inTree element with hybrid string
        aList[1] = hybStr
        # Note: inTree in TFL needs to be inside the loop!
        inTree = ''.join(aList)

    return inTree


########
# MAIN #
########

def main(treeName, parentInfo):

    # If hyrbid mode: confirm that parent likelihoods unequal
    if "," in parentInfo:
        search = re.search('\w+:(\d*\.\d+),\w+:(\d*\.\d+)', parentInfo,
                           re.IGNORECASE)
        if search.group(1) == search.group(2):
            print colored("  ERROR: Parent likelihoods must not be the \
same.\n  Quitting ...\n", "red")
            sys.exit()

    # Reading tree as string
    inData = open(treeName, "r").readlines()
    trees = []
    for line in inData:
        l = line.strip()
        if len(l) > 0:
            if l[0] != "#":
                trees.append(line)

    outList = []
    for treeStr in trees:
        # Reading tree by DendroPy
        tree = dendropy.Tree.get_from_string(treeStr, "newick")
        # DEBUGLINE: print(tree.as_ascii_plot())
        # Placeholder to potentially modify tree further

        # Left-ladderize tree
        tree.ladderize(ascending=False)
        treeStr = tree.as_string('newick')

        # remove [&U]
        treeStr = re.sub(r'(\[.*?\])\s*?', '', treeStr)
        treeStr = treeStr.lstrip().rstrip()

        # Parsing parentInfo into dictionary
        aDict = {}
        for i in parentInfo.split(","):
            aDict[i.split(":")[0]] = i.split(":")[1]
        # Hypothetical content of aDict at this point: {"B":"0.6", "C":"0.4"}

        # Adding hybrid information to inTree
        outTree = addHybrDict(treeStr, aDict)

        # Adding nodes to tree
        outTree = addNodeNumb(outTree)

        # Adding tree to outList
        outList.append(outTree)

    # Saving output to file
    if "," in parentInfo:
        outf = open(GSO.rmext(treeName)+".wHybrid.tre", "w")
    else:
        outf = open(GSO.rmext(treeName)+".wSister.tre", "w")
    outf.write("\n".join(outList))
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

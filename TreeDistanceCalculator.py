#! /usr/bin/env python3

desc = """This script check several trees (in the working folder) against a simulation tree.
For equal topologies, it also outputs weighted Robinson-Foulds and Felsenstein branch length distances. Output to stdout or a csv file.

João Brazão version 1.1 (September 2021)
"""

import os
import argparse
from p4 import *
import numpy as np
import textwrap

def tree_check(simulation_tree,file_format,dirpath,output_file):

    if not os.path.exists(simulation_tree):
        print("Error: %s does not exist" % (simulation_tree))
        sys.exit()

    if not os.path.exists(dirpath):
        print("Error: %s path does not exist"%(dirpath))
        sys.exit()
    else:
        startdir = os.getcwd()
        os.chdir(dirpath)

    read(simulation_tree)
    m = var.trees[0]
    m.taxNames = list([i.name for i in m.iterLeavesNoRoot()])

    tree_names = [] #trees file name
    wrf_list = []   #wrf_distances
    bld_list = []   #bld_distances
    tlen_list = []   #tree length
    sym_list = []   #unweighted Robinson Foulds distance, symmetric difference
    xtrees = [] #trees (name files) with different topology
    table = []

    n = 1
    for tree in sorted(os.listdir(dirpath)):
        if tree.endswith(file_format):	# treefile termination is for IQ-Tree outputs
            try:
                read(tree)
                t = var.trees[n]
                t.taxNames = m.taxNames
                tree_names.append(t.fName)
                wrf = m.topologyDistance(t,'wrf')
                wrf_list.append(float(wrf))
                bld = m.topologyDistance(t,'bld')
                bld_list.append(float(bld))
                tlen = t.getLen()
                tlen_list.append(float(tlen))
                sym = m.topologyDistance(t) #RF
                sym_list.append(float(sym))
                if sym != 0:
                    xtrees.append(t.fName)
                n+=1
            except:
                print ("\n... P4 error. file: %s ...\n" %(tree))

    wavg = np.mean(wrf_list)
    bavg = np.mean(bld_list)
    tavg = np.mean(tlen_list)
    savg = np.mean(sym_list)
    wstdev = np.std(wrf_list)
    bstdev = np.std(bld_list)
    tstdev = np.std(tlen_list)
    sstdev = np.std(sym_list)

    wrf_list.append(wavg)
    bld_list.append(bavg)
    tlen_list.append(tavg)
    sym_list.append(savg)
    wrf_list.append(wstdev)
    bld_list.append(bstdev)
    tlen_list.append(tstdev)
    sym_list.append(sstdev)
    table += [tree_names, wrf_list, bld_list, tlen_list, sym_list, xtrees]

    print (f"\n--> {n-1} tree(s) checked and {len(xtrees)} different topologies. The Robinson-foulds mean is {round(savg,2)}.\n")
    print ("--> The Weighted Robinson-Foulds mean is %f +- %f and the Felsensteins branch-length mean is %f +- %f \n" %(wavg,wstdev,bavg,bstdev))
    print ("--> The tree length mean is %f and difference is %f +-% f" %(tavg,abs(m.getLen()-tavg),tstdev))

    if output_file:
        table[0].append("Average")
        table[0].append("Standard_Deviation")
        with open(output_file+".csv", "w") as f:
            f.write("Tree_file, Weighted_Robinson-Foulds, Felsensteins_branch-length, Tree_length, Unweighted_RF")
            f.write("\n")
            for k in range(len(table[1])):
                for i in table[:-1]:
                    f.write(str(i[k]))
                    f.write(",")
                f.write("\n")
            f.write("\nTrees with different topologies\n")
            f.write(",".join(table[-1]))
        print (f"\n--> Output saved to {dirpath}{output_file}.csv")

###################################################################################

def main(simulation_tree,file_format,dirpath,output_file):

    tree_check(simulation_tree,file_format,dirpath,output_file)

###################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(desc),
        )
    parser.add_argument("-s","--simtree",
                    help="The simulation tree.",
                    )
    parser.add_argument("-f","--file_format",
                    help="The extension of the tree files format (eg.'nwk', 'treefile', etc)",
                    default='treefile'
                    )
    parser.add_argument("-d", "--dirpath",
                    dest="dirpath",
                    help="Directory containing analysis " +\
                         "Default: current directory",
                    default="./"
                    )
    parser.add_argument("-o","--outputfile",
                    help="The output file name. Default: None",
                    default=None
                    )

    args = parser.parse_args()
    main(args.simtree,args.file_format,args.dirpath,args.outputfile)

#! /usr/bin/env python3

desc="""Alignment data simulator. It generates several data sets in phylip format using a given model and tree.
This script is derived from the 'Scripts for common tasks' of the P4 phyloinformatic toolkit (vers. 1.3; Foster, 2004; https://p4.nhm.ac.uk/scripts.html)

"""

__author__='Joao Brazao'

import os
import re
from p4 import *
import argparse
import textwrap


def main(tree,length,number,model,rates,compositions,nexus,dirpath):

    if not os.path.exists(dirpath):
        print("Error: %s path does not exist"%(dirpath))
        sys.exit()
    else:
        startdir=os.getcwd()
        os.chdir(dirpath)

    read(tree)
    t=var.trees[0]
    taxNames=list([i.name for i in t.iterLeavesNoRoot()])

    #create empty data
    a=func.newEmptyAlignment(dataType='protein', taxNames=taxNames, length=length)
    d=Data([a])
    t.taxNames=taxNames

    t.data=d
    optimise=0
    if model: #a model present in P4
        t.newComp(free=0, spec=model)
        t.newRMatrix(free=0, spec=model)
    else: #Rates and compositions from the user.
        compositions_list=[float(n) for n in re.split(r'[,\s\t]+', compositions)]
        compositions_list[0]=compositions_list[0] + (1 - sum(compositions_list)) #Some times, due to rounded numbers, the sum != 1
        rates_list=[float(n) for n in re.split(r'[,\s\t]+', rates)]
        t.newComp(free=0, spec='specified', val=compositions_list)
        t.newRMatrix(free=0, spec='specified', val=rates_list)

    t.setNGammaCat(nGammaCat=4)
    t.newGdasrv(free=0, val=0.75)
    t.setPInvar(free=0, val=0.0)

    #Generate the data sets
    func.reseedCRandomizer(os.getpid())
    for n in range(1,number+1):
        t.simulate()
        if nexus:
            d.writeNexus(f"dataset{length}_{n}.nex", writeDataBlock=True)
        d.alignments[0].writePhylip(f"dataset{length}_{n}.phy", interleave=False, flat=True)
    print ("Simulated alignments created.")

####################################################################################################

if __name__ == "__main__":

    parser=argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument("tree",
            help="The simulation tree."
            )
    parser.add_argument("length",type=int,
            help="Simulated alignment length."
            )
    parser.add_argument("number",type=int,
            help="Number of replicates."
            )
    parser.add_argument("-m","--model",
            help="The amino-acid substitution model to use (it has to be available in P4). Or specify the rates (--rates) and composition (--compositions)."
            )
    parser.add_argument("-r","--rates",
            help="Substitution rates matrix. All 189 parameters separated by space or comma. Eg. '0.05 0.89 0.33 ...'"
            )
    parser.add_argument("-c","--compositions",
            help="Composition frequencies."
            )
    parser.add_argument("-n", "--nex",
                        help="Save also in nexus format.",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-d", "--dirpath",
                    dest="dirpath",
                    help="Favourite directory." +\
                         "Default: current directory",
                    default="./"
                    )
    args=parser.parse_args()
    main(args.tree, args.length, args.number, args.model, args.rates, args.compositions,args.nex, args.dirpath)

#! /usr/bin/env python3

desc = """Convert MCMC parameters (inferred using the P4 package (Foster,2004)), paml, and RAxML model formats to paml, RAxML or Phylobayes model formats.
(Part of this script is derived from Cymon J. Cox and P4 tutorial scripts)

\ne.g.\n$ aminoacids_model_format_convertor.py paml phylobayes -m gcpREV.dat -o gcpREV.pb\n

João Brazão version 1.1
"""

import os
import sys
import numpy
import argparse
import textwrap
from p4 import *
from itertools import zip_longest

def import_data(input_data_format, burnin, model, dirpath):

    if input_data_format not in ["p4","paml", "raxml"]:
        print ("Error: expecting either 'p4', 'paml' or 'raxml'")
        sys.exit()

    if not os.path.exists(dirpath):
        print("Error: %s path does not exist"%(dirpath))
        sys.exit()
    else:
        startdir = os.getcwd()
        os.chdir(dirpath)

    if input_data_format == "p4":
        if burnin:
            try:
                burnin=int(burnin)/100
            except ValueError:
                print("Error: expecting interger for burnin value")
                sys.exit()
        else:
            burnin=0

        if not os.path.exists("mcmc_pramsProfile_0.py"):
            print("Error: mcmc_pramsProfile_0.py does not exist")
            sys.exit()

        if not os.path.exists("mcmc_prams_0"):
            print("Error: mcmc_prams_0 does not exist")
            sys.exit()

        p = func.summarizeMcmcPrams(skip=burnin, makeDict=True)
        """{'part0': [['comp[0]', '0.100268', '0.000007 ', '101.1'],
                      ['comp[1]', '0.069685', '0.000005 ', '394.1'],
                      ['comp[2]', '0.028398', '0.000002 ', '383.6'],
        etc
        """
        for part in p.keys():
            if part != "part0":
                print("Error: More than one partition detected! This script only works for one partition")
                sys.exit()

        pramValues=[float(v[1])*10e4 for v in p["part0"][20:-1]]    #substitution rates
        compositions = [float(v[1]) for v in p["part0"][0:20]]

    if input_data_format == "paml":
        with open (model,'r') as f:
            matrix_values = [n for n in f.read().split()]
        if len(matrix_values) >= 210:
            matrixValues = [float(n) for n in matrix_values[0:190]]
            compositions = [float(n) for n in matrix_values[190:210]]
        else:
            sys.exit(f'{len(matrix_values)}: Exit. Your matrix has not the correct number of parameters')

        #Creation of a second list with the proper order of the rate values
        start_pst=[0, 2, 5, 9, 14, 20, 27, 35, 44, 54, 65, 77, 90, 104, 119, 135, 152, 170, 189]
        pramValues=[]
        c=0
        for k in start_pst: #  k is the top position of which pamL matrix column
            a=0+c
            e=k
            for i in range(19-c):
                pramValues.append(matrixValues[e])
                a+=1
                e+=a
            c+=1

    if input_data_format == "raxml":
        pramValues = []
        with open(model,"r") as f:
            fullmatrix=[float(value) for value in f.readlines()]# if value.strip() != "0.0"]
            assert len(fullmatrix) == 420, "Error: RAxML model does not have 420 input values, but %s" % len(fullmatrix)
        m=19
        for i in range(1,400,21):
            pramValues += fullmatrix[i:i+m]
            m -= 1
        compositions = fullmatrix[-20:]

    if sum(compositions) != 1.0:
        #print("Composition frequences added to > 1.0 - actual value: %s - adjusting" % sum(compositions))
        newvalue = 1.0 - (sum(compositions)-compositions.pop())
        compositions.append(newvalue)

    assert len(pramValues) == 190, "Incorrect number of Substitution rates %s" % len(pramValues)
    assert len(compositions) == 20, "Incorrect number of composition frequencies %s" % len(compositions)
    return pramValues, compositions

def output_model(pramValues, compositions, output_model_format, output_name):

    if output_model_format not in ["paml", "raxml", "phylobayes"]:
        print ("Error: expecting either 'paml', 'raxml' or 'phylobayes'")
        sys.exit()

    #Write 19-1 top triangle
    rows = []
    i_count = 0
    for n in range(19,0,-1):
        i = i_count
        j = i + n
        #print "[%i:%i]" % (i, j)
        values = []
        for value in pramValues[i:j]:
            values.append(value)
        rows.append(values)
        i_count += n

    if output_model_format == "phylobayes":
        doPHYLOBAYES(rows, output_name)
    else:
        nr = []
        for r in rows:
            r.reverse()
            nr.append(r)

        #zip cols into rows remove Nones from short cols
        nc = []
        for col in zip_longest(*nr):
            values = [v for v in col if v != None]
            nc.append(values)

        #Turn upside down
        nc.reverse()
        if output_model_format == "paml":
            doPAML(nc,compositions,output_name)
        elif output_model_format == "raxml":
            doRAxML(nc,compositions,output_name)

def doPHYLOBAYES(rows, output_name):
    if output_name:
        with open(output_name,"w") as f:
            f.write("A R N D C Q E G H I L K M F P S T W Y V\n")
            for row in rows:
                f.write (" ".join([str(v) for v in row]))
                f.write('\n')
        print("Your data was converted and saved in %s"%(output_name))
    else:
        print("A R N D C Q E G H I L K M F P S T W Y V")
        for row in rows:
            print(" ".join([str(v) for v in row]))

def doPAML(nc, compositions, output_name):
    if output_name:
        with open(output_name,"w") as f:
            for rates in nc:
                f.write("   ".join([str(round(v,3)) for v in rates]))
                f.write("\n")
            f.write("\n")
            f.write("  ".join([str(round(v,6)) for v in compositions]))
        print("Your data was converted and saved in %s"%(output_name))
    else:
        for rates in nc:
            print("  ".join([str(v) for v in rates]))
        print("\n")
        print(" ".join([str(v) for v in compositions]))

def doRAxML(nc,compositions,output_name):

    m = numpy.zeros((20, 20), numpy.float)
    lineNum = 0
    for row in nc:
        for i in range(len(row)):
            m[lineNum +1][i] = row[i]
            m[i][lineNum +1] = row[i]
        lineNum += 1

    if output_name:
        with open(output_name,"w") as f:
            for row in m:
                for value in row:
                    f.write (str(value) + '\n')
            for value in compositions:
                f.write (str(value))
                f.write ('\n')
        print("Your data was converted and saved in %s"%(output_name))
    else:
        for row in m:
            for value in row:
                print(value)
        for value in compositions:
            print(value)
##################################################################################
def main(input_data_format,burnin, model, dirpath, output_model_format, output_name):

    rates, compositions = import_data(input_data_format, burnin, model, dirpath)

    output_model(rates, compositions, output_model_format, output_name)

##################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument("input_data_format",
            help="Input format options: p4, paml or raxml.",
            )
    parser.add_argument("-b","--burnin",
            help="Number of samples to burnin of a P4 mcmc analysis.",
            )
    parser.add_argument("-m","--model",
            help="Amino-acid substitution model in paml or RAxML format.",
            )
    parser.add_argument("output_model_format",
            help="Output format options: paml, raxml, or phylobayes.",
            )
    parser.add_argument("-o","--outputfile", help="Output model name.")
    parser.add_argument("-d", "--dirpath",
                        dest="dirpath",
                        help="Directory containing analysis " +\
                             "Default: current directory",
                        default="./"
                        )

    args = parser.parse_args()
    main(args.input_data_format, args.burnin, args.model, args.dirpath, args.output_model_format, args.outputfile)

#aminoacids_model_format_convertor.py

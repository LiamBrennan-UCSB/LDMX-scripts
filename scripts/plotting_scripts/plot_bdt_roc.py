import argparse
import os
import sys

##################
# Set up
##################

# Set up the argument parser
parser = argparse.ArgumentParser(description = 'A python script to handle/pass options for ROC plots to C plotting macros')

parser.add_argument('--bdtName1',
                    dest = 'bdtName1',
                    default = 'bdt',
                    help = 'Name of the BDT')

parser.add_argument('--bdtName2',
                    dest = 'bdtName2',
                    default = '',
                    help = 'Name of the BDT')

parser.add_argument('--flatDir1',
                    dest = 'flatDir1',
                    default = '/nfs/slac/g/ldmx/users/aechavez/ldmx-sw-v3.0.0/LDMX-scripts/pyEcalVeto/training/eval_trees',
                    help = 'Directory where the flat evaluation trees are saved')

parser.add_argument('--flatDir2',
                    dest = 'flatDir2',
                    default = '',
                    help = 'Directory where the flat evaluation trees are saved')

parser.add_argument('--outDir',
                    dest = 'outDir',
                    default = '/nfs/slac/g/ldmx/users/aechavez/ldmx-sw-v3.0.0/plots',
                    help = 'Directory where the plots will be saved')

parser.add_argument('--withSelect',
                    dest = 'withSelect',
                    default = 'base',
                    help = 'Selection to make')

parser.add_argument('--bkgEff',
                    dest = 'bkgEff',
                    default = '0.0001',
                    help = 'Background efficiency')

parser.add_argument('--isLog',
                    dest = 'isLog',
                    default = False,
                    help = 'Whether or not to make the plot in log scale',
                    action = 'store_true')

parser.add_argument('--noZoom',
                    dest = 'noZoom',
                    default = False,
                    help = 'Whether or not to zoom in on the plot',
                    action = 'store_true')

args = parser.parse_args()

# Initialize the arguments from the parser
bdtName1 = args.bdtName1
bdtName2 = args.bdtName2
flatDir1 = args.flatDir1
flatDir2 = args.flatDir2
outDir = args.outDir
withSelect = args.withSelect
bkgEff = args.bkgEff
isLog = args.isLog
noZoom = args.noZoom

###########################
# Main subroutine
###########################

def main():

    # Convert boolean values into integers
    log = 1 if isLog else 0
    zoom = 0 if noZoom else 1

    # Pass the options to the plotting macro
    if bdtName2 == '' or flatDir2 == '':
        os.system("root -l -q -b 'plot_bdt_roc_0.C+("
            + '"{}", "{}", "{}", "{}", {}, {}, {}'.format(bdtName1, flatDir1, outDir, withSelect, bkgEff, log, zoom)
            + ")'"
        )
    else:
        os.system("root -l -q -b 'plot_bdt_roc_1.C+("
            + '"{}", "{}", "{}", "{}", "{}", "{}", {}, {}, {}'.format(bdtName1, bdtName2, flatDir1, flatDir2, outDir, withSelect, bkgEff, log, zoom)
            + ")'"
        )

###############
# RUN
###############

if __name__ == '__main__':
    main()


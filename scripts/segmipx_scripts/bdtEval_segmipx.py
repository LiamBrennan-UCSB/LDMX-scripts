import os
import sys
import numpy as np
import pickle as pkl
import xgboost as xgb
import mods.ROOTmanager as manager

pkl_file = '/nfs/slac/g/ldmx/users/aechavez/ldmx-sw-v3.0.0/segmipx_bdt/training/segmipx_train_out_0/segmipx_train_out_0_weights.pkl'
model = pkl.load(open(pkl_file,'rb'))

def main():

    # Inputs and their trees and stuff
    pdict = manager.parse()
    inlist = pdict['inlist']
    outlist = pdict['outlist']
    group_labels = pdict['groupls']
    maxEvent = pdict['maxEvents']

    # TreeModel to build here
    branches_info = {
            # Base Fernand variables
            'nReadoutHits':                 {'rtype': int,   'default': 0  },
            'summedDet':                    {'rtype': float, 'default': 0. },
            'summedTightIso':               {'rtype': float, 'default': 0. },
            'maxCellDep':                   {'rtype': float, 'default': 0. },
            'showerRMS':                    {'rtype': float, 'default': 0. },
            'xStd':                         {'rtype': float, 'default': 0. },
            'yStd':                         {'rtype': float, 'default': 0. },
            'avgLayerHit':                  {'rtype': float, 'default': 0. },
            'stdLayerHit':                  {'rtype': float, 'default': 0. },
            'deepestLayerHit':              {'rtype': int,   'default': 0  },
            'ecalBackEnergy':               {'rtype': float, 'default': 0. },
            # MIP tracking variables
            'straight4':                    {'rtype': int,   'default': 0  },
            'firstNearPhLayer':             {'rtype': int,   'default': 33 },
            'nNearPhHits':                  {'rtype': int,   'default': 0  },
            'photonTerritoryHits':          {'rtype': int,   'default': 0  },
            'epSep':                        {'rtype': float, 'default': 0. },
            'epDot':                        {'rtype': float, 'default': 0. },
            # Longitudinal segment variables
            'energy_s1':                    {'rtype': float, 'default': 0. },
            'xMean_s1':                     {'rtype': float, 'default': 0. },
            'yMean_s1':                     {'rtype': float, 'default': 0. },
            'layerMean_s1':                 {'rtype': int,   'default': 0  },
            'energy_s2':                    {'rtype': float, 'default': 0. },
            'yMean_s3':                     {'rtype': float, 'default': 0. },
            # Electron RoC variables
            'eContEnergy_x1_s1':            {'rtype': float, 'default': 0. },
            'eContEnergy_x2_s1':            {'rtype': float, 'default': 0. },
            'eContYMean_x1_s1':             {'rtype': float, 'default': 0. },
            'eContEnergy_x1_s2':            {'rtype': float, 'default': 0. },
            'eContEnergy_x2_s2':            {'rtype': float, 'default': 0. },
            'eContYMean_x1_s2':             {'rtype': float, 'default': 0. },
            # Photon RoC variables
            'gContNHits_x1_s1':             {'rtype': int,   'default': 0  },
            'gContYMean_x1_s1':             {'rtype': float, 'default': 0. },
            'gContNHits_x1_s2':             {'rtype': int,   'default': 0  },
            # Outside RoC variables
            'oContEnergy_x1_s1':            {'rtype': float, 'default': 0. },
            'oContEnergy_x2_s1':            {'rtype': float, 'default': 0. },
            'oContEnergy_x3_s1':            {'rtype': float, 'default': 0. },
            'oContNHits_x1_s1':             {'rtype': int,   'default': 0  },
            'oContXMean_x1_s1':             {'rtype': float, 'default': 0. },
            'oContYMean_x1_s1':             {'rtype': float, 'default': 0. },
            'oContYMean_x2_s1':             {'rtype': float, 'default': 0. },
            'oContYStd_x1_s1':              {'rtype': float, 'default': 0. },
            'oContEnergy_x1_s2':            {'rtype': float, 'default': 0. },
            'oContEnergy_x2_s2':            {'rtype': float, 'default': 0. },
            'oContEnergy_x3_s2':            {'rtype': float, 'default': 0. },
            'oContLayerMean_x1_s2':         {'rtype': int,   'default': 0  },
            'oContLayerStd_x1_s2':          {'rtype': float, 'default': 0. },
            'oContEnergy_x1_s3':            {'rtype': float, 'default': 0. },
            'oContLayerMean_x1_s3':         {'rtype': int,   'default': 0  },
            # Quantities needed for BDT analysis
            'discValue_segmipx':            {'rtype': float, 'default': 0.5},
            'isAtTSP':                      {'rtype': int,   'default': 0  },
            'isAtESP':                      {'rtype': int,   'default': 0  },
            'recoilPT':                     {'rtype': float, 'default': 0. }
            }

    # Construct tree processes
    procs = []
    for gl, group in zip(group_labels, inlist):
        procs.append( manager.TreeProcess(event_process, group, ID=gl, tree_name='EcalVeto',
            pfreq=100) )

    # Process jobs
    for proc in procs:

        print('\nRunning %s'%(proc.ID))
        
        # Move into appropriate scratch dir
        os.chdir(proc.tmp_dir)

        # Make an output file and new tree (copied from input + discValue)
        proc.tfMaker = manager.TreeMaker(group_labels[procs.index(proc)]+'.root',\
                                         "EcalVeto",\
                                         branches_info,\
                                         outlist[procs.index(proc)]
                                         )

        # RUN
        proc.extrafs = [ proc.tfMaker.wq ] # Gets executed at the end of run()
        proc.run(maxEvents=maxEvent)

    # Remove scratch directory if there is one
    manager.rmScratch()

    print('\nDone!\n')


def event_process(self):

    # Feature list from input tree
    # Exp: feats = [ feat_value for feat_value in self.tree~ ]
    feats = [
            # Base variables
            self.tree.nReadoutHits             , # 0
            self.tree.summedDet                , # 1
            self.tree.summedTightIso           , # 2
            self.tree.maxCellDep               , # 3
            self.tree.showerRMS                , # 4
            self.tree.xStd                     , # 5
            self.tree.yStd                     , # 6
            self.tree.avgLayerHit              , # 7
            self.tree.stdLayerHit              , # 8
            self.tree.deepestLayerHit          , # 9
            self.tree.ecalBackEnergy           , # 10
            # MIP Tracking variables
            self.tree.straight4                , # 11
            self.tree.firstNearPhLayer         , # 12
            self.tree.nNearPhHits              , # 13
            self.tree.photonTerritoryHits      , # 14
            self.tree.epSep                    , # 15
            self.tree.epDot                    , # 16
            # Longitudinal segment variables
            self.tree.energy_s1                , # 17
            self.tree.xMean_s1                 , # 18
            self.tree.yMean_s1                 , # 19
            self.tree.layerMean_s1             , # 20
            self.tree.energy_s2                , # 21
            self.tree.yMean_s3                 , # 22
            # Electron RoC variables
            self.tree.eContEnergy_x1_s1        , # 23
            self.tree.eContEnergy_x2_s1        , # 24
            self.tree.eContYMean_x1_s1         , # 25
            self.tree.eContEnergy_x1_s2        , # 26
            self.tree.eContEnergy_x2_s2        , # 27
            self.tree.eContYMean_x1_s2         , # 28
            # Photon RoC variables
            self.tree.gContNHits_x1_s1         , # 29
            self.tree.gContYMean_x1_s1         , # 30
            self.tree.gContNHits_x1_s2         , # 31
            # Outside RoC variables
            self.tree.oContEnergy_x1_s1        , # 32
            self.tree.oContEnergy_x2_s1        , # 33
            self.tree.oContEnergy_x3_s1        , # 34
            self.tree.oContNHits_x1_s1         , # 35
            self.tree.oContXMean_x1_s1         , # 36
            self.tree.oContYMean_x1_s1         , # 37
            self.tree.oContYMean_x2_s1         , # 38
            self.tree.oContYStd_x1_s1          , # 39
            self.tree.oContEnergy_x1_s2        , # 40
            self.tree.oContEnergy_x2_s2        , # 41
            self.tree.oContEnergy_x3_s2        , # 42
            self.tree.oContLayerMean_x1_s2     , # 43
            self.tree.oContLayerStd_x1_s2      , # 44
            self.tree.oContEnergy_x1_s3        , # 45
            self.tree.oContLayerMean_x1_s3       # 46
            ]

    # Copy input tree feats to new tree
    for feat_name, feat_value in zip(self.tfMaker.branches_info, feats):
        self.tfMaker.branches[feat_name][0] = feat_value

    # Add the quantities needed for BDT analysis to the tree
    evtarray = np.array([feats])
    pred = float(model.predict(xgb.DMatrix(evtarray))[0])
    self.tfMaker.branches['discValue_segmipx'][0] = pred
    self.tfMaker.branches['isAtTSP'][0] = self.tree.isAtTSP
    self.tfMaker.branches['isAtESP'][0] = self.tree.isAtESP
    self.tfMaker.branches['recoilPT'][0] = self.tree.recoilPT

    # Fill new tree with current event values
    self.tfMaker.tree.Fill()

if __name__ == "__main__":
    main()

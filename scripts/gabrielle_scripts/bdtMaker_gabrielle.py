import os
import sys
import logging
import argparse
import ROOT as r
import numpy as np
import pickle as pkl
import xgboost as xgb
import matplotlib as plt
from array    import array
from optparse import OptionParser

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
plt.use('Agg')

class sampleContainer:

    def __init__(self,filename,maxEvts,trainFrac,isSig):

        print("Initializing container!")
        self.tree = r.TChain("EcalVeto")
        self.tree.Add(filename)
        self.maxEvts = maxEvts
        self.trainFrac = trainFrac
        self.isSig   = isSig

    def root2PyEvents(self):
        self.events =  []
        for event in self.tree:

            if len(self.events) >= self.maxEvts:
                continue

            # Pre-selection to filter out anomalous events
            if event.isAtTSP == 0 and event.isAtESP == 1:
                continue

            evt = [
                    # Base Fernand variables
                    event.nReadoutHits                , # 0
                    event.summedDet                   , # 1
                    event.summedTightIso              , # 2
                    event.maxCellDep                  , # 3
                    event.showerRMS                   , # 4
                    event.xStd                        , # 5
                    event.yStd                        , # 6
                    event.avgLayerHit                 , # 7
                    event.stdLayerHit                 , # 8
                    event.deepestLayerHit             , # 9
                    event.ecalBackEnergy              , # 10
                    # Electron RoC variables
                    event.electronContainmentEnergy_x1, # 11
                    event.electronContainmentEnergy_x2, # 12
                    event.electronContainmentEnergy_x3, # 13
                    event.electronContainmentEnergy_x4, # 14
                    event.electronContainmentEnergy_x5, # 15
                    # Photon RoC variables
                    event.photonContainmentEnergy_x1  , # 16
                    event.photonContainmentEnergy_x2  , # 17
                    event.photonContainmentEnergy_x3  , # 18
                    event.photonContainmentEnergy_x4  , # 19
                    event.photonContainmentEnergy_x5  , # 20
                    # Outside RoC variables
                    event.outsideContainmentEnergy_x1 , # 21
                    event.outsideContainmentEnergy_x2 , # 22
                    event.outsideContainmentEnergy_x3 , # 23
                    event.outsideContainmentEnergy_x4 , # 24
                    event.outsideContainmentEnergy_x5 , # 25
                    event.outsideContainmentNHits_x1  , # 26
                    event.outsideContainmentNHits_x2  , # 27
                    event.outsideContainmentNHits_x3  , # 28
                    event.outsideContainmentNHits_x4  , # 29
                    event.outsideContainmentNHits_x5  , # 30
                    event.outsideContainmentXStd_x1   , # 31
                    event.outsideContainmentXStd_x2   , # 32
                    event.outsideContainmentXStd_x3   , # 33
                    event.outsideContainmentXStd_x4   , # 34
                    event.outsideContainmentXStd_x5   , # 35
                    event.outsideContainmentYStd_x1   , # 36
                    event.outsideContainmentYStd_x2   , # 37
                    event.outsideContainmentYStd_x3   , # 38
                    event.outsideContainmentYStd_x4   , # 39
                    event.outsideContainmentYStd_x5     # 40
            ]

            self.events.append(evt)

        new_idx=np.random.permutation(np.arange(np.shape(self.events)[0]))
        self.events = np.array(self.events)
        np.take(self.events, new_idx, axis=0, out=self.events)
        print("Final Event Shape" + str(np.shape(self.events)))

    def constructTrainAndTest(self):
        self.train_x = self.events[0:int(len(self.events)*self.trainFrac)]
        self.test_x = self.events[int(len(self.events)*self.trainFrac):]

        self.train_y = np.zeros(len(self.train_x)) + (self.isSig == True)
        self.test_y = np.zeros(len(self.test_x)) + (self.isSig == True)

class mergedContainer:
    def __init__(self, sigContainer,bkgContainer):
        self.train_x = np.vstack((sigContainer.train_x,bkgContainer.train_x))
        self.train_y = np.append(sigContainer.train_y,bkgContainer.train_y)
        
        self.train_x[np.isnan(self.train_x)] = 0.000
        self.train_y[np.isnan(self.train_y)] = 0.000

        self.test_x  = np.vstack((sigContainer.test_x,bkgContainer.test_x))
        self.test_y  = np.append(sigContainer.test_y,bkgContainer.test_y)
        
        self.dtrain = xgb.DMatrix(self.train_x,self.train_y)
        self.dtest  = xgb.DMatrix(self.test_x,self.test_y)

if __name__ == "__main__":
    
    # Parse
    parser = OptionParser()
    parser.add_option('--seed', dest='seed',type="int",  default=2, help='Numpy random seed.')
    parser.add_option('--max_evt', dest='max_evt',type="int",  default=1500000, help='Max Events to load')
    parser.add_option('--train_frac', dest='train_frac',  default=.8, help='Fraction of events to use for training')
    parser.add_option('--eta', dest='eta',type="float",  default=0.023, help='Learning Rate')
    parser.add_option('--tree_number', dest='tree_number',type="int",  default=1000, help='Tree Number')
    parser.add_option('--depth', dest = 'depth',type = 'int',  default = 10, help = 'Max tree depth')
    parser.add_option('-b', dest = 'bkg_file', default = './bdt_0/bkg_train.root', help = 'Name of background file')
    parser.add_option('-s', dest = 'sig_file', default = './bdt_0/sig_train.root', help = 'Name of signal file')
    parser.add_option('-o', dest = 'out_name',  default = 'bdt_test', help = 'Output pickle name')
    (options, args) = parser.parse_args()

    # Seed numpy's randomness
    np.random.seed(options.seed)
   
    # Get BDT num
    bdt_num = 0
    Check = True
    while Check:
        if not os.path.exists(options.out_name + '_' + str(bdt_num)):
            try:
                os.makedirs(options.out_name + '_' + str(bdt_num))
                Check = False
            except:
               Check = True
        else:
            bdt_num += 1

    # Print run info
    print( 'Random seed is = {}'.format(options.seed)             )
    print( 'You set max_evt = {}'.format(options.max_evt)         )
    print( 'You set tree number = {}'.format(options.tree_number) )
    print( 'You set max tree depth = {}'.format(options.depth)    )
    print( 'You set eta = {}'.format(options.eta)                 )

    # Make Signal Container
    print( 'Loading sig_file = {}'.format(options.sig_file) )
    sigContainer = sampleContainer(options.sig_file, options.max_evt, options.train_frac, True)
    sigContainer.root2PyEvents()
    sigContainer.constructTrainAndTest()

    # Make Background Container
    print( 'Loading bkg_file = {}'.format(options.bkg_file) )
    bkgContainer = sampleContainer(options.bkg_file, options.max_evt, options.train_frac, False)
    bkgContainer.root2PyEvents()
    bkgContainer.constructTrainAndTest()

    # Merge
    eventContainer = mergedContainer(sigContainer,bkgContainer)

    params = {
               'objective': 'binary:logistic',
               'eta': options.eta,
               'max_depth': options.depth,
               'min_child_weight': 20,
               # 'silent': 1,
               'subsample': 0.9,
               'colsample_bytree': 0.85,
               # 'eval_metric': 'auc',
               'eval_metric': 'error',
               'seed': 1,
               'nthread': 1,
               'verbosity': 1
               # 'early_stopping_rounds' : 10
    }

    # Train the BDT model
    evallist = [(eventContainer.dtrain, 'train'), (eventContainer.dtest, 'eval')]
    gbm = xgb.train(params, eventContainer.dtrain, num_boost_round = options.tree_number, evals = evallist, early_stopping_rounds = 10)

    # Store BDT
    output = open(options.out_name + '_' + str(bdt_num) + '/' + \
            options.out_name + '_' + str(bdt_num) + '_weights.pkl', 'wb')
    pkl.dump(gbm, output)

    # Plot feature importances
    xgb.plot_importance(gbm)
    plt.pyplot.savefig(options.out_name + '_' + str(bdt_num) + '/' + \
            options.out_name + '_' + str(bdt_num) + '_fimportance.png', # png file name
            dpi = 500, bbox_inches = 'tight', pad_inches = 0.5) # png parameters
    
    # Closing statment
    print('Files saved in: ', options.out_name + '_' + str(bdt_num))


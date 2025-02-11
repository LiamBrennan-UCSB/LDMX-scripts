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

        print("Initializing Container!")
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
                    event.nReadoutHits             , # 0
                    event.summedDet                , # 1
                    event.summedTightIso           , # 2
                    event.maxCellDep               , # 3
                    event.showerRMS                , # 4
                    event.xStd                     , # 5
                    event.yStd                     , # 6
                    event.avgLayerHit              , # 7
                    event.stdLayerHit              , # 8
                    event.deepestLayerHit          , # 9
                    event.ecalBackEnergy           , # 10
                    # MIP Tracking variables
                    event.straight4                , # 11
                    event.firstNearPhLayer         , # 12
                    event.nNearPhHits              , # 13
                    event.photonTerritoryHits      , # 14
                    event.epSep                    , # 15
                    event.epDot                    , # 16
                    # Longitudinal segment variables
                    event.energy_s1                , # 17
                    event.xMean_s1                 , # 18
                    event.yMean_s1                 , # 19
                    event.layerMean_s1             , # 20
                    event.energy_s2                , # 21
                    event.yMean_s3                 , # 22
                    # Electron RoC variables
                    event.eContEnergy_x1_s1        , # 23
                    event.eContEnergy_x2_s1        , # 24
                    event.eContYMean_x1_s1         , # 25
                    event.eContEnergy_x1_s2        , # 26
                    event.eContEnergy_x2_s2        , # 27
                    event.eContYMean_x1_s2         , # 28
                    # Photon RoC variables
                    event.gContNHits_x1_s1         , # 29
                    event.gContYMean_x1_s1         , # 30
                    event.gContNHits_x1_s2         , # 31
                    # Outside RoC variables
                    event.oContEnergy_x1_s1        , # 32
                    event.oContEnergy_x2_s1        , # 33
                    event.oContEnergy_x3_s1        , # 34
                    event.oContNHits_x1_s1         , # 35
                    event.oContXMean_x1_s1         , # 36
                    event.oContYMean_x1_s1         , # 37
                    event.oContYMean_x2_s1         , # 38
                    event.oContYStd_x1_s1          , # 39
                    event.oContEnergy_x1_s2        , # 40
                    event.oContEnergy_x2_s2        , # 41
                    event.oContEnergy_x3_s2        , # 42
                    event.oContLayerMean_x1_s2     , # 43
                    event.oContLayerStd_x1_s2      , # 44
                    event.oContEnergy_x1_s3        , # 45
                    event.oContLayerMean_x1_s3       # 46
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
    parser.add_option('--depth', dest='depth',type="int",  default=10, help='Max Tree Depth')
    parser.add_option('-b', dest='bkg_file', default='./bdt_0/bkg_train.root', help='name of background file')
    parser.add_option('-s', dest='sig_file', default='./bdt_0/sig_train.root', help='name of signal file')
    parser.add_option('-o', dest='out_name',  default='bdt_test', help='Output Pickle Name')
    (options, args) = parser.parse_args()

    # Seed numpy's randomness
    np.random.seed(options.seed)
   
    # Get BDT num
    bdt_num=0
    Check=True
    while Check:
        if not os.path.exists(options.out_name+'_'+str(bdt_num)):
            try:
                os.makedirs(options.out_name+'_'+str(bdt_num))
                Check=False
            except:
               Check=True
        else:
            bdt_num+=1

    # Print run info
    print( 'Random seed is = {}'.format(options.seed)             )
    print( 'You set max_evt = {}'.format(options.max_evt)         )
    print( 'You set tree number = {}'.format(options.tree_number) )
    print( 'You set max tree depth = {}'.format(options.depth)    )
    print( 'You set eta = {}'.format(options.eta)                 )

    # Make Signal Container
    print( 'Loading sig_file = {}'.format(options.sig_file) )
    sigContainer = sampleContainer(options.sig_file,options.max_evt,options.train_frac,True)
    sigContainer.root2PyEvents()
    sigContainer.constructTrainAndTest()

    # Make Background Container
    print( 'Loading bkg_file = {}'.format(options.bkg_file) )
    bkgContainer = sampleContainer(options.bkg_file,options.max_evt,options.train_frac,False)
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
               'subsample':.9,
               'colsample_bytree': .85,
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
    output = open(options.out_name+'_'+str(bdt_num)+'/' + \
            options.out_name+'_'+str(bdt_num)+'_weights.pkl', 'wb')
    pkl.dump(gbm, output)

    # Plot feature importances
    xgb.plot_importance(gbm)
    plt.pyplot.savefig(options.out_name+'_'+str(bdt_num)+"/" + \
            options.out_name+'_'+str(bdt_num)+'_fimportance.png', # png file name
            dpi=500, bbox_inches='tight', pad_inches=0.5) # png parameters
    
    # Closing statment
    print("Files saved in: ", options.out_name+'_'+str(bdt_num))

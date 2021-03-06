import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da

#eosDir='/eos/uscms/store/user/mba2012/IMGs/HighLumi_ROOTv2'
#eosDir='/eos/uscms/store/user/mba2012/IMGs/h24gamma_eta14'
eosDir='/eos/cms/store/user/mandrews/OPENDATA/IMGs/MGG90_Eta14v3'
#decays = ['h22gammaSM_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
#decays = ['SM2gamma_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
#decays = ['SM2gamma_1j_1M_noPU', 'h22gammaSM_1j_1M_noPU']
#decays = ['SM2gamma_1j_1M_noPU', 'h22gammaSM_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
decays = ['DiPhotonBorn_MGG90_Eta14v3', 'GluGluHToGG_MGG90_Eta14v3', 'GJet_MGG90_Eta14v3']

#chunk_size_ = 250
chunk_size_ = 200
#chunk_size_ = 273
scale = 1.

@delayed
def load_X(tree, start_, stop_, branches_, readouts, scale):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    # Convert the object array X to a multidim array:
    # 1: for each event x in X, concatenate the object columns (branches) into a flat array of shape (readouts*branches)
    # 2: reshape the flat array into a stacked array: (branches, readouts)
    # 3: embed each stacked array as a single row entry in a list via list comprehension
    # 4: convert this list into an array with shape (events, branches, readouts) 
    X = np.array([np.concatenate(x).reshape(len(branches_),readouts[0]*readouts[1]) for x in X])
    #print "X.shape:",X.shape
    X = X.reshape((-1,len(branches_),readouts[0],readouts[1]))
    X = np.transpose(X, [0,2,3,1])
    # Rescale
    X /= scale 
    return X

@delayed
def load_single(tree, start_, stop_, branches_):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    X = np.array([x[0] for x in X])

    return X

for j,decay in enumerate(decays):

    #if j == 2 or j == 1:
    if j != 2:
        pass
        #continue

    tfile_str = '%s/%s_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_nXXX_IMG.root'%(eosDir,decay)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    #tree = ROOT.TChain("fevt/RHTree")
    #tree.Add('%s/DiPhotonBorn_MGG90All_Eta14_IMG.root'%eosDir)
    #tree.Add('%s/DiPhotonBox_MGG90All_Eta14_IMG.root'%eosDir)
    nevts = tree.GetEntries()
    #neff = (nevts//1000)*1000
    #neff = 200
    #neff = 63000
    neff = 57400
    #neff = 84600
    #neff = 110000 
    #neff = 81900 
    chunk_size = chunk_size_
    if neff > nevts:
        neff = int(nevts)
        chunk_size = int(nevts)
    #neff = 1000 
    #neff = 233000
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    # EB
    readouts = [170,360]
    #branches = ["TracksPt_EB","EB_energy"]
    branches = ["EB_energyT", "EB_energyZ"]
    X = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X.shape

    ## eventId
    #branches = ["eventId"]
    #eventId = da.concatenate([\
    #            da.from_delayed(\
    #                load_single(tree,i,i+chunk_size, branches),\
    #                shape=(chunk_size,),\
    #                dtype=np.int32)\
    #            for i in range(0,neff,chunk_size)])
    #print " >> Expected shape:", eventId.shape

    # runId
    branches = ["runId"]
    runId = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.int32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", runId.shape

    # m0
    branches = ["m0"]
    m0 = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", m0.shape

    # diPhoE
    branches = ["diPhoE"]
    diPhoE = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", diPhoE.shape

    # diPhoPt
    branches = ["diPhoPt"]
    diPhoPt = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", diPhoPt.shape

    # FC inputs 
    branches = ["FC_inputs"]
    X_FC = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,5),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_FC.shape

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(X.shape[0], label, dtype=np.float32),\
            chunks=(chunk_size,))

    file_out_str = "%s/%s_IMG_RH%d+FC_n%d_label%d.hdf5"%(eosDir,decay,int(scale),neff,label)
    #file_out_str = "%s/DiPhotonAll_MGG90All_Eta14_IMG_RH%d_n%d_label%d.hdf5"%(eosDir,int(scale),neff,label)
    #file_out_str = "%s/%s_IMG_RH%d_n%dk_label%d.hdf5"%(eosDir,decay,int(scale),neff//1000.,label)
    #file_out_str = "test.hdf5"
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y, 'runId': runId, 'm0': m0, 'diPhoE': diPhoE, 'diPhoPt': diPhoPt}, compression='lzf')
    da.to_hdf5(file_out_str, {
                              '/X': X,
                              '/X_FC': X_FC,
                              '/y': y, 
                              'runId': runId, 
                              'm0': m0, 
                              'diPhoE': diPhoE, 
                              'diPhoPt': diPhoPt}, compression='lzf')

    print " >> Done.\n"

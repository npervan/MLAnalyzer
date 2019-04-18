import glob, os
import ROOT
import numpy as np
import dask.array as da
from convert_Tree2Dask_utils import *
import itertools

import argparse
parser = argparse.ArgumentParser(add_help=True, description='Process some integers.')
parser.add_argument('-l', '--label', required=True, type=int, help='Decay label.')
parser.add_argument('-of', '--nOutFiles', type=int, default=0, help='Number of files the inputs will be merged into. If this equals 0, # of outfiles = # of in files')
parser.add_argument('-c', '--chooseFile', action='store_true', help='Choose whether or not you wish to only run on one file, if you select this argument, must also use -n')
parser.add_argument('-if', '--fileNum', type=int, help='Choose which file you are processing. If you choose this option, you must also use -c')
parser.add_argument('-j', '--nJets', type=int, default=1, help='How many are in your root file. For ttbar, select 1 or 2. For QCD, this selection does not matter (it will be overwritten as 1)')
parser.add_argument('-s', '--file_idx_start', default=1, type=int, help='File index start.')
parser.add_argument('-e', '--file_idx_end', default=9999999, type=int, help='File index end.')
parser.add_argument('-g', '--granularity', default=1, type=int, help='Increased Pixel Granularity')
args = parser.parse_args()
if not ( (args.chooseFile) != (args.fileNum is None) ):
    parser.error('Cannot use option -if or -c without also using the other')

#outDir='test_output'
outDir='/home/bjornb/Granularity_tests/hdf5'
xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
#decays = ['QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU', 'QCDToQQ_Pt_80_120_13TeV_TuneCUETP8M1_noPU']
#decays = ['QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8','TT_Mtt1500toInf_TuneCUETP8M1_14TeV-powheg-pythia8']
decays = ['QCD-Granularity-Tests', 'TTbar-Granularity-Tests']
#decays=['test']
#test_file=xrootd+'//store/user/bburkle/E2E/ML_TopTagging/TT_Mtt1500toInf_TuneCUETP8M1_14TeV-powheg-pythia8/TTbar_Image_NTuplizer_TopTagged/181110_043333/0000/TT_Mtt1500toInf_TuneCUETP8M1_14TeV-powheg-pythia8_976.root'

scale = [1., 1.]
chunk_size = 50
jet_shape = 125 * args.granularity
njets = args.nJets
nOutFiles = args.nOutFiles
tot_evts = 0

#qcd_300to600 = open('filelists/qcd_300to600_noshift').read().splitlines()
#for i in range(len(qcd_300to600)):
#    qcd_300to600[i] = '/mnt/hadoop/store/user/bburkle/E2E/opendata_ntuples/CRAB_UserFiles/QCD_300to600_njets1_opendata_shifted-pixel/190328_214655/0000/'+qcd_300to600[i]
#for i in range(len(qcd_600to3000)):
#    qcd_600to3000[i] = '/mnt/hadoop/store/user/bburkle/E2E/opendata_ntuples/CRAB_UserFiles/QCD_600to3000_njets1_opendata_shifted-pixel/190328_214806/0000/'+qcd_600to3000[i]

ttbar_files = glob.glob('/home/bjornb/Granularity_tests/ntuples/ttbar_fixSig/*')

qcd_300to600 = glob.glob('/home/bjornb/Granularity_tests/ntuples/qcd_300to600_fixSig/*')
qcd_400to600 = glob.glob('/home/bjornb/Granularity_tests/ntuples/qcd_400to600_fixSig/*')
qcd_600to3000 = glob.glob('/home/bjornb/Granularity_tests/ntuples/qcd_600to3000_fixSig/*')

qcd_files = qcd_300to600 + qcd_400to600 + qcd_600to3000


#test = 'test_output/test_dz.root'
test = 'test_output/ttbar_new-production_test.root'
# Loop over decays
for d, decay in enumerate(decays):

    if d != args.label:
        #pass
        continue
   
    if d == 0:
        filelist = qcd_files
        tfile_idxs = range(1, len(qcd_files)+1)
        njets = 1
    elif d == 1:
        filelist = ttbar_files
        tfile_idxs = range(1, len(ttbar_files)+1)
    else:
        print 'decay must be equal to 0 or 1'
        break
    if nOutFiles == 0: nOutFiles = len(tfile_idxs)
    	
    print '>> Doing decay[%d]: %s'%(d, decay)

    tfile_idxs.sort()
    print '>> File idxs:', tfile_idxs
    idx_chunk = int(len(tfile_idxs)/nOutFiles)
    #break

    if args.chooseFile: #primarily for debug purposes
        tfile_idxs = [args.fileNum]
        idx_chunk = 1
        #tfile_idxs = [1,2]
        #idx_chunk = 1

    eventId_array = []
    lumiId_array = []
    runId_array = []
    X_ECAL_array = []
    X_ECAL_EEup_array = []
    X_ECAL_stacked_array = []
    X_EB_array = []
    X_EEm_array = []
    X_HBHE_array = []
    X_HBHE_EM_array = []
    X_HBHE_EB_up_array = []
    jetSeed_iphi_array = []
    jetSeed_ieta_array = []
    X_jets_array = []
    jetPt_array = []
    jetM_array = []
    y_array = []

    # Loop over root ntuples
    for idxs in [tfile_idxs[i:i+idx_chunk] for i in range(0, len(tfile_idxs), idx_chunk)]:
        for n in idxs:

            if n < args.file_idx_start:
                continue
            if n > args.file_idx_end:
                break

            1 == 1
            #tfile_str = glob.glob('%s/%s*_IMG/*/*/output_%d.root'%(eosDir,decay,n))
	        #tfile_str = glob.glob('%s/%s/*_%d.root'%(eosDir,decay,n))
	        #tfile_str = str('%s/%s/%s_%d.root'%(eosDir,decay,decay,n))
	        #tfile_str = str('%s/%s_%d.root'%(eosDir,decay,n))
            tfile_str = str(filelist[n])
	
            #tfile_str = 'output_dijet.root' # DEBUG mode: for single, local file
            print " >> For input file:", tfile_str
            #tfile = ROOT.TXNetFile(tfile_str)
            tfile = ROOT.TFile(tfile_str)
            tree = tfile.Get('fevt/RHTree')
            nevts = tree.GetEntries()

            #neff = (nevts//1000)*1000
            #neff = (nevts//100)*100
            #neff = 200
            neff = int(nevts)
            #tot_evts += int(nevts)
            if neff < chunk_size:
                chunk_size = neff
            if neff > nevts:
                neff = int(nevts)
            proc_range = list(range(0, neff, chunk_size))[:-1]
            print " >> Total events:", nevts
            print " >> Effective events:", neff
            print " >> Proccess list is:", list(proc_range)

            # eventId
            branches = ["eventId"]
            eventId = da.concatenate([\
                        da.from_delayed(\
                            load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                            shape=(get_chunk_size(i,neff,chunk_size),),\
                            dtype=np.int32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],eventId.shape)

            ## lumiId
            #branches = ["lumiId"]
            #lumiId = da.concatenate([\
            #            da.from_delayed(\
            #                load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
            #                shape=(get_chunk_size(i,neff,chunk_size),),\
            #                dtype=np.int32)\
            #            for i in proc_range])
            #print " >> %s: %s"%(branches[0],lumiId.shape)

            # runId
            branches = ["runId"]
            runId = da.concatenate([\
                        da.from_delayed(\
                            load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                            shape=(get_chunk_size(i,neff,chunk_size),),\
                            dtype=np.int32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],runId.shape)

            # ECAL
            readouts = [280,360]
            branches = ["ECAL_energy"]
            X_ECAL = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_ECAL.shape)

            # ECAL with resampled EE
            X_ECAL_EEup = X_ECAL.map_blocks(lambda x: block_resample_EE(x), dtype=np.float32)
            print " >> %s: %s"%('ECAL_EEup_energy',X_ECAL_EEup.shape)

             TODO get this working
             ECAL Upsample
            if args.granularity != 1:
                X_ECAL_EEup = tile_array(X_ECAL_EEup, args.granularity, args.granularity)

            # pT Tracks at ECAL
            if args.granularity != 1:
                readouts = [args.granularity*280, args.granularity*360]
                branches = ['ECALadj_tracksPt_%dx%d'%(args.granularity, args.granularity)]
            else:
                readouts = [280, 360]
                branches = ['ECAL_tracksPt_atECALfixIPfromPV']
            X_TracksAtECAL = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_TracksAtECAL.shape)

            # Average d0 Tracks at ECAL
            if args.granularity != 1:
                readouts = [args.granularity*280, args.granularity*360]
                branches = ['ECALadj_tracksD0Sig_%dx%d'%(args.granularity,args.granularity)]
            else:
                readouts = [280, 360]
                branches = ['ECAL_tracksD0Sig_atECALfixIP']
            X_D0TracksAtECAL = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_D0TracksAtECAL.shape)

            # Average dz Tracks at ECAL
            if args.granularity != 1:
                readouts = [args.granularity*280, args.granularity*360]
                branches = ['ECALadj_tracksDzSig_%dx%d'%(args.granularity,args.granularity)]
            else:
                readouts = [280, 360]
                branches = ['ECAL_tracksDzSig_atECALfixIP']
            X_DzTracksAtECAL = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_DzTracksAtECAL.shape)

            # Pixel L1 Rec Hits Fixed
            if args.granularity != 1:
                readouts = [args.granularity*280, args.granularity*360]
                branches = ['BPIX_layer1_ECALadj_%dx%d'%(args.granularity,args.granularity)]
            else:
                readouts = [280, 360]
                branches = ['BPIX_layer1_ECAL_atPV']
            X_PixelRecHitsL1 = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_PixelRecHitsL1.shape)

            # Pixel L2 Rec Hits Fixed
            if args.granularity != 1:
                readouts = [args.granularity*280, args.granularity*360]
                branches = ['BPIX_layer2_ECALadj_%dx%d'%(args.granularity,args.granularity)]
            else:
                readouts = [280, 360]
                branches = ['BPIX_layer2_ECAL_atPV']
            X_PixelRecHitsL2 = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_PixelRecHitsL2.shape)

            # Pixel L3 Rec Hits Fixed
            if args.granularity != 1:
                readouts = [args.granularity*280, args.granularity*360]
                branches = ['BPIX_layer3_ECALadj_%dx%d'%(args.granularity,args.granularity)]
            else:
                readouts = [280, 360]
                branches = ['BPIX_layer3_ECAL_atPV']
            X_PixelRecHitsL3 = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_PixelRecHitsL3.shape)

            # HBHE upsample
            readouts = [56,72]
            branches = ["HBHE_energy"]
            upscale = 5*args.granularity
            X_HBHE_up = da.concatenate([\
                        da.from_delayed(\
                            load_X_upsampled(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1], upscale),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s(upsampled): %s"%(branches[0],X_HBHE_up.shape)
        
            #X_MuonsAtECAL, 
            #X_ECAL_stacked = da.concatenate([X_ECAL_EEup, X_HBHE_up], axis=-1) #if no tracks
            #X_ECAL_stacked = da.concatenate([X_TracksAtECAL, X_ECAL_EEup, X_HBHE_up], axis=-1) # tracks, ECAL, HCAL
            #X_ECAL_stacked = da.concatenate([X_TracksAtECAL, X_ECAL_EEup, X_HBHE_up, X_D0TracksAtECAL, X_DzTracksAtECAL], axis=-1) #pt tracks, ECAL, HCAL, d0 tracks, dz tracks
            #X_ECAL_stacked = da.concatenate([X_TracksAtECAL, X_D0TracksAtECAL, X_DzTracksAtECAL, X_HBHE_up], axis=-1) #pt tracks, d0 tracks, dz tracks, HCAL
            #print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)
            #X_Pixel_stacked = da.concatenate([X_PixelRecHitsL1, X_PixelRecHitsL2, X_PixelRecHitsL3], axis=-1)
            #print " >> %s: %s"%('X_Pixel_stacked', X_Pixel_stacked.shape)
            #X_Jet_stacked = da.concatenate([X_ECAL_stacked, X_Pixel_stacked], axis=-1)
            #X_Jet_stacked = da.concatenate([X_TracksAtECAL, X_D0TracksAtECAL, X_DzTracksAtECAL, X_HBHE_up, X_PixelRecHitsL1, X_PixelRecHitsL2, X_PixelRecHitsL3], axis=-1)
            X_Jet_stacked = da.concatenate([X_TracksAtECAL, X_D0TracksAtECAL, X_DzTracksAtECAL, X_PixelRecHitsL1, X_PixelRecHitsL2, X_PixelRecHitsL3], axis=-1)
            print " >> %s: %s"%('X_Jet_stacked', X_Jet_stacked.shape)

            # EB
            readouts = [170,360]
            #branches = ["HBHE_energy_EB"]
            #branches = ["TracksQPt_EB","EB_energy"]
            #branches = ["TracksPt_EB","EB_energy"]
            branches = ["EndTracksPt_EB","EB_energy"]
            #branches = ["EB_energy"]
            #branches = ["EB_energy","HBHE_energy_EB","Tracks_EB"]
            X_EB = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_EB.shape)

            # EE-
            readouts = [100,100]
            branches = ["TracksPt_EEm","EEm_energy","HBHE_energy_EEm"]
            #branches = ["EndTracksPt_EEm","EEm_energy","HBHE_energy_EEm"] # for pt weighted at ECAL face
            #branches = ["TracksQPt_EEm","EEm_energy","HBHE_energy_EEm"] # for Qxpt weighted 
            #branches = ["EEm_energy","HBHE_energy_EEm","Tracks_EEm"]
            X_EEm = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_EEm.shape)

            # EE+
            readouts = [100,100]
            branches = ["TracksPt_EEp","EEp_energy","HBHE_energy_EEp"]
            #branches = ["EndTracksPt_EEp","EEp_energy","HBHE_energy_EEp"] # for pt weighted at ECAL face
            #branches = ["TracksQPt_EEp","EEp_energy","HBHE_energy_EEp"] # for Qxpt weighted
            #branches = ["EEp_energy","HBHE_energy_EEp","Tracks_EEp"]
            X_EEp = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_EEp.shape)

            # HBHE
            readouts = [56,72]
            branches = ["HBHE_energy"]
            X_HBHE = da.concatenate([\
                        da.from_delayed(\
                            load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s: %s"%(branches[0],X_HBHE.shape)

            # HBHE_EM
            #readouts = [56,72]
            #branches = ["HBHE_EMenergy"]
            #X_HBHE_EM = da.concatenate([\
            #            da.from_delayed(\
            #                load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
            #                shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
            #                dtype=np.float32)\
            #            for i in proc_range])
            #print " >> %s: %s"%(branches[0],X_HBHE_EM.shape)

            # HB_EB upsample
            readouts = [34,72]
            branches = ["HBHE_energy_EB"]
            upscale = 5
            X_HBHE_EB_up = da.concatenate([\
                        da.from_delayed(\
                            load_X_upsampled(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1], upscale),\
                            shape=(get_chunk_size(i,neff,chunk_size), readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                            dtype=np.float32)\
                        for i in proc_range])
            print " >> %s(upsampled): %s"%(branches[0],X_HBHE_EB_up.shape)

            X_EB = da.concatenate([X_EB, X_HBHE_EB_up], axis=-1)
            print " >> %s: %s"%('X_EB', X_EB.shape)
   
            # Loop over jets
            for ijet in range(njets):

                print ' >> jet index:',ijet

                # jet m0
                branches = ["jetM"]
                jetM = da.concatenate([\
                            da.from_delayed(\
                                load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                                shape=(get_chunk_size(i,neff,chunk_size),),\
                                dtype=np.float32)\
                            for i in proc_range])
                print "  >> jetM:", jetM.shape

                # jet pt
                branches = ["jetPt"]
                jetPt = da.concatenate([\
                            da.from_delayed(\
                                load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                                shape=(get_chunk_size(i,neff,chunk_size),),\
                                dtype=np.float32)\
                            for i in proc_range])
                print "  >> jetPt:", jetPt.shape

                # jet seed iphi
                branches = ["jetSeed_iphi"]
                jetSeed_iphi = da.concatenate([\
                            da.from_delayed(\
                                load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                                shape=(get_chunk_size(i,neff,chunk_size),),\
                                dtype=np.float32)\
                            for i in proc_range])
                print "  >> jetSeed_iphi:", jetSeed_iphi.shape

                # jet seed ieta
                branches = ["jetSeed_ieta"]
                jetSeed_ieta = da.concatenate([\
                            da.from_delayed(\
                                load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                                shape=(get_chunk_size(i,neff,chunk_size),),\
                                dtype=np.float32)\
                            for i in proc_range])
                print "  >> jetSeed_ieta:", jetSeed_ieta.shape

                # jet window
                X_jets = da.concatenate([\
                            da.from_delayed(\
                                 #crop_jet_block(X_ECAL_stacked[i:i+get_chunk_size(i,neff,chunk_size)],\
                                 crop_jet_block(X_Jet_stacked[i:i+get_chunk_size(i,neff,chunk_size)],\
                                                 jetSeed_iphi[i:i+get_chunk_size(i,neff,chunk_size)],\
                                                 jetSeed_ieta[i:i+get_chunk_size(i,neff,chunk_size)], jet_shape, args.granularity),\
                                shape=(get_chunk_size(i,neff,chunk_size), jet_shape, jet_shape, X_Jet_stacked.shape[-1]),\
                                #shape=(get_chunk_size(i,neff,chunk_size), jet_shape, jet_shape, X_ECAL_stacked.shape[-1]),\
                                dtype=np.float32)\
                            for i in proc_range])
                print "  >> X_jets:", X_jets.shape

                # Class label
                label = d
                print "  >> Class label:",label
                y = da.from_array(\
                        np.full(len(X_jets), label, dtype=np.float32),\
                        chunks=(get_chunk_size(i,neff,chunk_size),))
                print "  >> y shape:",y.shape

                print "blah blah"
                #eventId_array.append(eventId)
                #lumiId_array.append(lumiId)
                #runId_array.append(runId)
                #X_ECAL_array.append(X_ECAL)
                #X_ECAL_EEup_array.append(X_ECAL_EEup)
                #X_ECAL_stacked_array.append(X_ECAL_stacked)
                #X_EB_array.append(X_EB)
                #X_EEm_array.append(X_EEm)
                #X_HBHE_array.append(X_HBHE)
                #X_HBHE_EM_array.append(X_HBHE_EM)
                #X_HBHE_EB_up_array.append(X_HBHE_EB_up)
                #jetSeed_iphi_array.append(jetSeed_iphi)
                #jetSeed_ieta_array.append(jetSeed_ieta)
                #X_jets_array.append(X_jets)
                #jetPt_array.append(jetPt)
                #jetM_array.append(jetM)
                #y_array.append(y)

                #outPath = '%s/%s_IMGjet'%(outDir, decay)
                #outPath = outDir
                #if not os.path.isdir(outPath):
                #    os.makedirs(outPath)
                print y.dask
                #print da.concatenate(jetPt_array, axis=0).dask
                #file_out_str = "%s/%s_IMGjet_RH%d_n%d_label%d_jet%d_%d.hdf5"%(outPath,decay,int(scale[0]),neff,label,ijet,n)
                #file_out_str = "%s/%s_IMGjet_RH%d_n%d_label%d_%d-%d.hdf5"%(outPath,decay,int(scale[0]),sum([iy for iy in y_array]),label,idxs[0],idxs[-1])
                #file_out_str = "%s/%s_IMGjet_RH%d_n%d_label%d_jet%d-of-%d_file%d.hdf5"%(outPath,decay,int(scale[0]),y.shape[0],label,ijet+1,njets,idxs[0])  #use this one
                #file_out_str = "%s/%s_noTracks_IMGjet_RH%d_n%d_label%d_jet%d-of-%d_file%d.hdf5"%(outPath,decay,int(scale[0]),y.shape[0],label,ijet+1,njets,idxs[0])
                #file_out_str = 'test_output/test_new_ttbar.hdf5'
                #file_out_str = '%s/%s_f%d-%d_n%d.hdf5'%(outPath, outName, args.file_idx_start, args.file_idx_end,tot_evts)
                #file_out_str = '%s/%s_f%d-%d.hdf5'%(outPath, outName, args.file_idx_start, args.file_idx_end)

                outPath = '%s/%s_IMGjet_%d-Granularity'%(outDir, decay, args.granularity)
                if not os.path.isdir(outPath):
                    os.makedirs(outPath)
                file_out_str = '%s/%s_f%d_j%d_n%d.hdf5'%(outPath, decay, n, ijet+1, int(nevts))
                if os.path.isfile(file_out_str):
                    os.remove(file_out_str)

                print "  >> Writing to:", file_out_str
                da.to_hdf5(file_out_str, {
                                  #'eventId': da.concatenate(eventId_array, axis=0),
                                  #'lumiId': da.concatenate(lumiId_array, axis=0),
                                  #'runId': da.concatenate(runId_array, axis=0),
                                  #'X_ECAL': da.concatenate(X_ECAL_array, axis=0),
                                  #'X_ECAL_EEup': da.concatenate(X_ECAL_EEup_array, axis=0),
                                  #'X_ECAL_stacked': da.concatenate(X_ECAL_stacked_array, axis=0),
                                  #'X_EB': da.concatenate(X_EB_array, axis=0),
                                  #'X_EEm': da.concatenate(X_EEm_array, axis=0),
                                  #'X_EEp': da.concatenate(X_EEp_array, axis=0),
                                  #'X_HBHE': da.concatenate(X_HBHE_array, axis=0),
                                  #'X_HBHE_EM': da.concatenate(X_HBHE_EM_array, axis=0),
                                  #'X_HBHE_EB_up': da.concatenate(X_HBHE_EB_up_array, axis=0),
                                  #'jetSeed_iphi': da.concatenate(jetSeed_iphi_array, axis=0),
                                  #'jetSeed_ieta': da.concatenate(jetSeed_ieta_array, axis=0),
                                  #'X_jets': da.concatenate(X_jets_array, axis=0),
                                  #'jetPt': da.concatenate(jetPt_array, axis=0),
                                  #'jetM': da.concatenate(jetM_array, axis=0),
                                  #'/y': da.concatenate(y_array, axis=0)
                                  #'X_Pixel_stacked': X_Pixel_stacked,
                                  #'X_ECAL_stacked': X_ECAL_stacked,
                                  'X_jets': X_jets,
                                  'jetPt': jetPt,
                                  'jetM': jetM,
                                  '/y': y
                                  }, compression='lzf')

            tfile.Close()

    #print "Transfering to: %s/%s"%(xrootd,eosOut)
    #os.system('xrdcp %s %s/%s'%(file_out_str,xrootd,eosOut))
    print "  >> Done.\n"

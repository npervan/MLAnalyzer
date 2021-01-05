import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/step2_ttbar_p8_03/step2_qcd8_1109.root'
#inputFiles_='file:../aod_test.root'
#inputFiles_='/store/user/npervan/e2e/jmar_aodsim/ZprimeToTT_M-2000_W-20_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FA50EB2A-3EC8-E611-B9C5-02163E019C8B.root'
inputFiles_='/store/user/lpcml/npervan/qcd_aodsim/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/FE8523C4-95B6-E611-A75D-001C23BED6B6.root'
#70000/1817C734-42C8-E611-BB9F-02163E019DD8.root'
#inputFiles_='/store/user/npervan/e2e/jmar_aodsim/ZprimeToTT_M-500_W-5_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/022AC3C8-CABC-E611-9C50-0CC47A78A42C.root'
#inputFiles_='file:/uscms/home/bburkle/nobackup/working_area/CMSSW_5_3_32/src/opendatadnn/step2_test.root'
#inputFiles_='file:/uscms/home/bburkle/nobackup/working_area/CMSSW_5_3_32/src/MLAnalyzer/test/step2_ttbarOD.root'
#inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/CRAB_UserFiles/step2_QCD600to3000_01/190213_183439/0000/step2_QCDPt_15_3000_Flat_V27_961.root'
#inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/CRAB_UserFiles/step2_ttbarOD_EmBj_01/190308_200019/0000/step2_OpenData_10.root'

isTTbar_ = 1

#maxEvents_=100
#skipEvents_=0#
#outputFile_ = 'test.root'
outputFile_ = 'test/top_aod_test_10.root'

cmd="cmsRun %s inputFiles=%s outputFile=%s isTTbar=%d" %(cfg,inputFiles_,outputFile_,isTTbar_)
print '%s'%cmd
os.system(cmd)

from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
import random
config = config()


#mc='TTbar'
mc='QCD'
pT=200
m_num=3500
filelist_ = ''
files_ = []
the_name=''
splitting = 1 
nJets_=2
isTTbar_=0

#dataset_ = '/RSGluonToTT_M-4000_TuneCUETP8M1_14TeV-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
assert mc == 'TTbar' or mc == 'QCD'
if mc=='TTbar':
    filelist_ = 'aod_m-%d_filelist.txt' % (m_num)
    the_name = 'TTbar_pT%d_%dTjet_aod_m-%d' % (pT,nJets_,m_num)
    splitting = 1 #current ~50 events per file
    files_ = open(filelist_).readlines()
    isTTbar_=1
elif mc=='QCD':
    filelist_ = 'qcd_aod_filelist.txt'
    the_name = 'QCD_pT%d_%djet_aod_Dec20' % (pT,nJets_)
    splitting = 3 #current 300 events per file
    files_ = open(filelist_).readlines()
    random.shuffle(files_)
    #nJets_=1
    isTTbar_=0

outputFile_ = '%s.root' % the_name
maxEvents_ = -1
skipEvents_ = 0


config.General.requestName = the_name
config.General.workArea = '/uscms/home/npervan/nobackup/e2e/slc7_jmar/CMSSW_8_0_21/src/MLAnalyzer'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RecHitAnalyzer/python/ConfFile_cfg.py'
config.JobType.pyCfgParams = ['outputFile=%s'%outputFile_,'isTTbar=%d'%isTTbar_]
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
#config.JobType.pyCfgParams = ['maxEvents=%s'%maxEvents_,'outputFile=%s'%outputFile_,'skipEvents=%s'%skipEvents_]
#config.JobType.pluginName = 'RecHitAnalyzer/python/ConfFile_cfg.py outputFile=%s maxEvents=%s skipEvents=%s' % (outputFile_, maxEvents_, skipEvents_)

config.Data.userInputFiles = files_
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = splitting #if splitting by files
config.Data.publication = False
#config.Data.outLFNDirBase = '/store/user/%s/e2e/jmar_aod_ntuples/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/group/lpcml/npervan/QCD_ntuples_Dec20/'
config.Data.outputDatasetTag = the_name

config.Site.storageSite = 'T3_US_FNALLPC'

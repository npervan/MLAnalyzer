filelist_ = 'aod_m-2000_filelist.txt'
#the_name = 'TTbar_pT%d_%dTjet_aod_m-%d' % (pT,nJets_,m_num)
splitting = 20 #current ~50 events per file                                                                                                              
files_ = open(filelist_).readlines()
isTTbar_=1

for file in files_:
    print(file)

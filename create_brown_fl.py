import subprocess

file_nums = [1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 500, 750]
outfile = open('jmar_aod_ntuple_filelist.txt','w+')

for num in file_nums:
    folder='/store/user/npervan/e2e/top_jmar121019/CRAB_UserFiles/TTbar_pT200_1Tjet_aod_m-%s/' %(num)
    print(folder)
    #folder='/store/user/npervan/'                                                                                                                                                

    eosls_subdir=subprocess.check_output(['eos','root://cmseos.fnal.gov','ls',folder]).strip()
    print('eos','root://cmseos.fnal.gov','ls',folder+eosls_subdir+'/0000/')
    eosls=subprocess.check_output(['eos','root://cmseos.fnal.gov','ls',folder+eosls_subdir+'/0000/'])
    #eosls=subprocess.check_output(['ls',folder])                                                                                                                                 
        #print eosls                                                                                                                                                              
    names=eosls.splitlines()
    if 'log' in names:                                                                                                                                             
        names.remove('log')                                                                                                                                                      
    #if 'failed' in names:                                                                                                                                                        
    #    names.remove('failed')                                                                                                                                                   
    for k in names:
        outfile.write(folder+eosls_subdir+'/0000/'+k+'\n')
    #outfile.write(folder+'\n')

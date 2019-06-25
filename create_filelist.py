import subprocess

folder='/mnt/hadoop/store/mc/RunIISummer16DR80Premix/ZprimeToTT_M-2000_W-20_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/'
outfile = open('aod_m-2000_filelist.txt','w+')
eosls=subprocess.check_output(['ls',folder])
        #print eosls
names=eosls.splitlines()
if 'log' in names:
    names.remove('log')
if 'failed' in names:
    names.remove('failed')
for k in names:
    outfile.write(folder+k+'\n')

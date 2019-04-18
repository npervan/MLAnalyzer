#!/usr/bin/env python
import os, sys
import glob

import argparse
parser = argparse.ArgumentParser(description='Tell what label you are using')
parser.add_argument('-l', '--label', required=True, type=int, help='0 is QCD, 1 is TTbar')
parser.add_argument('-s', '--splitting', type=int, default=10, help='Number of input files per condor job')
parser.add_argument('-j', '--nJets', type=int, default=1, help='Number of jets in the root file, will always be reset to 1 for QCD')
parser.add_argument('-g', '--gran', type=int, default=1, help='Increased track and pixel granularity')
parser.add_argument('--start', type=int, default=0, help='Which file to start on')
parser.add_argument('--end', type=int, default=99999, help='Which file to end on')
args=parser.parse_args()

ttbar_files = len(glob.glob('/home/bjornb/Granularity_tests/ntuples/ttbar_fixSig/*'))
qcd_300to600 = glob.glob('/home/bjornb/Granularity_tests/ntuples/qcd_300to600_fixSig/*')
qcd_400to600 = glob.glob('/home/bjornb/Granularity_tests/ntuples/qcd_400to600_fixSig/*')
qcd_600to3000 = glob.glob('/home/bjornb/Granularity_tests/ntuples/qcd_600to3000_fixSig/*')

mc = ''
split = args.splitting
njets = args.nJets
qcd_files = len(qcd_300to600) + len(qcd_400to600) + len(qcd_600to3000)

assert args.label == 1 or args.label == 0
if args.label == 0:
        mc = 'QCD'
        jets = 1
        files = qcd_files
elif args.label == 1:
        mc = 'TTbar'
        files = ttbar_files
def condor(num):
        basedir = os.getcwd()
        if not os.path.exists(basedir+'/condor/%s'%mc):
                os.makedirs(basedir+'/condor/%s'%mc)
        for i in range(num):
                start = split*i + 1
                end = split*(i+1)
                if start < args.start or end > args.end:
                    continue
                dict={'start':start,'end':end,'dir':basedir,'mc':mc,'label':args.label, 'jets':njets, 'gran':args.gran}
                filename = '%(dir)s/condor/%(mc)s/condor_%(jets)sjets_%(start)s-%(end)s_x%(gran)s.job' % dict
                jdf = open(filename, 'w')
                jdf.write("""
universe = vanilla
executable = %(dir)s/doConvert.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Notification = Error

use_x509userproxy=true

Arguments = %(start)s %(end)s %(label)s %(jets)s %(gran)s

Output = %(dir)s/condor/%(mc)s/log_out/log_%(start)s-%(end)s_x%(gran)s.stdout
Error = %(dir)s/condor/%(mc)s/log_out/log_%(start)s-%(end)s_x%(gran)s.stderr
Log = %(dir)s/condor/%(mc)s/log_out/log_%(start)s-%(end)s_x%(gran)s.condorlog

Queue 1
                """%dict)

                jdf.close()
                os.system('condor_submit %s' % filename)

pass
condor(int(files/split)+1)


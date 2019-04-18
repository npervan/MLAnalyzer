import numpy as np
import os
import h5py
import glob

np.random.seed(1337)

import argparse
parser = argparse.ArgumentParser(add_help=True, description='Weight QCD files by pT')
parser.add_argument('-s', '--start', default=0, type=int, help='first file to weight')
parser.add_argument('-e', '--end', default=-1, type=int, help='last file to weight')
args = parser.parse_args()

#files = glob.glob('QCD_IMGjet/*.hdf5')
files = glob.glob('QCD-Granularity-Tests_IMGjet_1-Granularity/*')
print len(files)
outdir = 'QCD_weighted_x1'
keys = ['X_jets', 'jetPt', 'jetM', 'y']

weights = open('qcd_weighting.txt', 'r').read().splitlines()
weights = [ list(map(float, line.split('\t'))) for line in weights if float(line.split('\t')[0]) >= 400.0 ]
max_wt = max([wt[-1] for wt in weights])
wts = [ [wt[0], wt[1], wt[2]/max_wt] for wt in weights ]

def weight_pt(pts):
    print 'There are', len(pts), 'Jets'
    for i, pt in enumerate(pts):
        #print pt
        for low, high, wt in wts:
            #print low, high
            if pt > high:
                continue
            elif pt < low:
                break
            elif np.random.random() > wt:
                continue
            else:
                yield i

def main():
    nEvts = 0
    #fs = list(range(args.start, args.end))
    fs = files[args.start:args.end]
    for i, file in enumerate(fs):
        print i, file
        try:
            f = h5py.File(file, 'r')
        except:
            print 'Skipping File'
            continue
        keep = list(weight_pt(f['jetPt']))
        nevts = len(keep)
        print nevts, 'jets are being kept'
        if nevts == 0:
            continue
        outstr = '%s/qcd_f%d_n%d.hdf5'%(outdir, i+args.start, nevts)
        if os.path.isfile(outstr):
            continue
        fout = h5py.File(outstr, 'w')
        chunk_temp = 100
        if nevts < chunk_temp:
            chunk_temp = nevts
        fout.create_dataset(keys[0], (nevts, 125*1, 125*1, 6), chunks = (chunk_temp, 125*1, 125*1, 6), compression='lzf')
        for key in keys[1:]:
            fout.create_dataset(key, (nevts, ), chunks = (chunk_temp, ), compression='lzf')
        try:
            fout[keys[0]][...] = np.take(f[keys[0]],keep,axis=0)
        except:
            fout.close()
            os.system('rm %s' % outstr)
            continue
        for key in keys[1:]:
            fout[key][...] = np.take(f[key],keep,axis=0)
        fout.close()
        f.close()

if __name__ == '__main__':
    main()


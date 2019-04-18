import h5py
import numpy as np
import glob
import os
import random
random.seed(1337)

test = False

batch_sz = 32
chunk_sz = batch_sz*500
#if os.path.isfile(output_name):
#    os.remove(output_name)

output_name = 'BoostedJets-split_x1'

#TTbar_dir = 'TTbar_IMGjet'
#TTbar_dir = 'ttbar_decorrelated'
#TTbar_dir = 'TTbar-RecHits-noCal_IMGjet_fixedTracks_withPix'
TTbar_dir = 'TTbar-Granularity-Tests_IMGjet_1-Granularity'
#QCD_dir = 'qcd_weighted'
#QCD_dir = 'qcd_decorrelated'
#QCD_dir = 'qcd_weighted_rechit'
QCD_dir = 'QCD_weighted_x1'
outDir = 'merged_files'
TTbar_files = glob.glob('%s/*'%TTbar_dir)
QCD_files = glob.glob('%s/*'%QCD_dir)
random.shuffle(TTbar_files)
random.shuffle(QCD_files)

TTbar_jets = 0
QCD_jets = 0
nevts = 0
TTbar_split = [0]
QCD_split = [0]
for i, jet in enumerate(TTbar_files):
    TTbar_jets += int(jet.split('_')[-1][1:-5])
    nevts += int(jet.split('_')[-1][1:-5])
    if nevts >= chunk_sz:
        TTbar_split.append(i+1)
        nevts = 0
TTbar_split.append(len(TTbar_files)-1)
nevts = 0
for i, jet in enumerate(QCD_files):
    QCD_jets += int(jet.split('_')[-1][1:-5])
    nevts += int(jet.split('_')[-1][1:-5])
    if nevts >= chunk_sz:
        QCD_split.append(i+1)
        nevts = 0
QCD_split.append(len(QCD_files)-1)

print 'Number of TTbar Jets:', TTbar_jets
print 'Number of QCD Jets:', QCD_jets
print 'TTbar splits are:', TTbar_split
print 'QCD_splits are:', QCD_split


file_sz = int(2*min(TTbar_jets,QCD_jets) / chunk_sz) * chunk_sz
if test == True:
    file_sz = 2*chunk_sz

def WriteToFile(i, TTbar_start, TTbar_stop, QCD_start, QCD_stop, TTbar_files, QCD_files):
    TTbar_start = TTbar_split[i]
    TTbar_stop = TTbar_split[i+1]
    QCD_start = QCD_split[i]
    QCD_stop = QCD_split[i+1]

    files = TTbar_files[TTbar_start:TTbar_stop] + QCD_files[QCD_start:QCD_stop]
    random.shuffle(files)
    dsets = [h5py.File(file) for file in files]

    X = np.concatenate([dset[input_keys[0]] for dset in dsets])
    #try:
    #    X = np.concatenate([dset[input_keys[0]] for dset in dsets])
    #except:
    #    print 'Merge failed for this chunk, not sure why. These are the files:'
    #    print files
    #    return 0
    X[0,...][X[0,...] < 1.e-3] = 0 # zero-suppresion
    #X[X < 1.e-3] = 0
    #X[-1,...] = 25.*X[-1,...]
    #X[1,...] = X[1,...]*2.
    #X = X/100.
    X[0,...] = X[0,...]/200. #pt
    X[1,...] = X[1,...]/1000. #d0
    X[2,...] = X[2,...]/5000. #dz
    pt = np.concatenate([dset[input_keys[1]] for dset in dsets])
    m0 = np.concatenate([dset[input_keys[2]] for dset in dsets])
    y = np.concatenate([dset[input_keys[3]] for dset in dsets])

    for dset in dsets:
        dset.close()

    print '\nNumber of events in group %d is %d' % ((i+1), len(y))
    if len(y) < 2*chunk_sz:
        print "Number of events smaller than chunk size"
        return 0

    l = zip(X, pt, m0, y)
    random.shuffle(l)
    X, pt, m0, y = zip(*l)
    X = X[:2*chunk_sz]
    pt = pt[:2*chunk_sz]
    m0 = m0[:2*chunk_sz]
    y = y[:2*chunk_sz]

    print 'Checking first 10 y values are:'
    for j in range(10):
        print y[j]
    
    print 'Writing chunks to output file'
    #g = f.CreateGroup('group_%d'%i+1)
    #g.create_dataset(output_keys[0], data = X, chunks = batch_sz*10)
    #g.create_dataset(output_keys[1], data = pt, chunks = batch_sz*10)
    #g.create_dataset(output_keys[2], data = m0, chunks = batch_sz*10)
    #g.create_dataset(output_keys[3], data = y, chunks = batch_sz*10)
    f[output_keys[0]][i*2*chunk_sz:(i+1)*2*chunk_sz] = X
    f[output_keys[1]][i*2*chunk_sz:(i+1)*2*chunk_sz] = pt
    f[output_keys[2]][i*2*chunk_sz:(i+1)*2*chunk_sz] = m0
    f[output_keys[3]][i*2*chunk_sz:(i+1)*2*chunk_sz] = y

    return len(y)

max_per_file = 320000*10
total_jets = 0
file_number = 1

f = h5py.File('%s/%s_file-%d.hdf5'%(outDir,output_name,file_number), 'w')

input_keys = ['X_jets', 'jetPt', 'jetM', 'y']
output_keys = ['X_jets', 'pt', 'm0', 'y']

f.create_dataset(output_keys[0], (file_sz, 125*1, 125*1, 6), chunks = (10*batch_sz, 125*1, 125*1, 6), compression='lzf')
for key in output_keys[1:]:
    f.create_dataset(key, (file_sz, ), chunks = (10*batch_sz, ), compression='lzf')

for i in range( min(len(TTbar_split), len(QCD_split)) - 1 ):
    total_jets += WriteToFile(i, TTbar_split[i], TTbar_split[i+1], QCD_split[i], QCD_split[i+1], TTbar_files, QCD_files)
    print 'Succesfully wrote chunk %d to file' % (i+1)
    if total_jets >= max_per_file:
        f.close()
        total_jets = 0
        file_number += 1
        f = h5py.File('%s_file-%d.hdf5'%(output_name, file_number), 'w')
        f.create_dataset(output_keys[0], (file_sz, 125*1, 125*1, 1), chunks = (10*batch_sz, 125*1, 125*1, 6), compression='lzf')
        for key in output_keys[1:]:
            f.create_dataset(key, (file_sz, ), chunks = (10*batch_sz, ), compression='lzf')
    if test == True:
        break

f.close()

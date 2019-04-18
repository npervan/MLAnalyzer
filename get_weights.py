import h5py
import ROOT as r
import numpy as np
import rootpy
import glob

QCD_files = glob.glob('QCD-Granularity-Tests_IMGjet_1-Granularity/*.hdf5')
ttbar_files = glob.glob('TTbar-Granularity-Tests_IMGjet_1-Granularity/*.hdf5')

nbins = 100
xmin = 400
xmax = 2000
qcd_hist = r.TH1D('qcd_hist', 'qcd_hist', nbins, xmin, xmax)
ttbar_hist = r.TH1D('ttbar_hist', 'ttbar_hist', nbins, xmin, xmax)

qcd_list = []
ttbar_list = []

def makePlot(h1, h2, h3):
    c = r.TCanvas('c', 'c', 600, 600)
    pad1 = r.TPad('pad1', 'pad1', 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetGridx()
    pad1.Draw()
    pad1.cd()
    h1.SetStats(0)
    h1.Draw('hist')
    h2.Draw('same hist')

    c.cd()
    pad2 = r.TPad('pad2', 'pad2', 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()
    pad2.cd()

    h3.SetStats(0)
    h3.SetMarkerStyle(21)
    h3.Draw('eq')

    h1.SetLineColor(r.kBlue+1)
    h1.SetLineWidth(2)
    h1.GetYaxis().SetTitleSize(20)
    h1.SetTitle('TTbar and QCD pT Distribution')
    h1.GetYaxis().SetTitleFont(43)
    h1.GetYaxis().SetTitleOffset(1.55)
    h1.GetYaxis().SetTitle('Number of Jets')
    
    h2.SetLineColor(r.kRed)
    h2.SetLineWidth(2)
    h3.SetTitle('')
    h3.GetYaxis().SetTitle('ttbar/qcd')
    h3.GetYaxis().SetTitleSize(20)
    h3.GetYaxis().SetTitleFont(43)
    h3.GetYaxis().SetTitleOffset(1.55)

    c.SaveAs('TTbar_QCD_pT.png')

#print 'reading', len(QCD_files), 'qcd files'
for files in QCD_files:
    print 'Reading file', files
    try:
        file = h5py.File(files, 'r')
    except:
        print 'Skipping file'
        continue
    qcd_list.append(np.asarray(file['jetPt'], 'd'))
    file.close()
qcd = np.concatenate(qcd_list)
print 'There are', len(qcd), 'qcd jets'

f_bad = open('bad_files.txt', 'w')

#print 'reading', len(ttbar_files), 'ttbar files'
for files in ttbar_files:
    print 'Reading file', files
    try:
        file = h5py.File(files, 'r')
    except:
        print 'skipping file'
        f_bad.write('%s\n'%files)
        continue
    ttbar_list.append(np.asarray(file['jetPt'], 'd'))
    file.close()
ttbar = np.concatenate(ttbar_list)
f_bad.close()
print 'There are', len(ttbar), 'ttbar jets'
print 'There are', len(qcd), 'qcd jets'

for jet in qcd:
    qcd_hist.Fill(jet)
print 'qcd hist jets:', qcd_hist.GetEntries()
qcd_hist.Scale(1/qcd_hist.GetEntries())
#qcd_hist.Scale(1)
for jet in ttbar:
    ttbar_hist.Fill(jet)
print 'ttbar hist jets:', ttbar_hist.GetEntries()
ttbar_hist.Scale(1/ttbar_hist.GetEntries())
#ttbar_hist.Scale(1)

ratios = ttbar_hist.Clone()
ratios.Divide(qcd_hist)

bin_list = list(range(xmin, xmax, (xmax-xmin)/nbins))
bin_list.append(xmax)
max_ratio = ratios.GetMaximum()
f = open('qcd_weighting.txt', 'w')
for bin in range(nbins):
    #f.write('%d\t%d\t%.4f\n' % (bin_list[bin], bin_list[bin+1], ratios.GetBinContent(bin+1)/max_ratio))
    f.write('%d\t%d\t%.4f\n' % (bin_list[bin], bin_list[bin+1], ratios.GetBinContent(bin+1)))
f.close()

makePlot(ttbar_hist, qcd_hist, ratios)

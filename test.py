from ROOT import TFile

infile=TFile('output.root')
tree=infile.Get('fevt/RHTree')
index=100
length=10
for e in tree:
	tracksPt=[]
	tracksPtadj=[]
	for i in range(len(e.ECAL_tracksPt)):
		if(e.ECAL_tracksPt[i]>0.001):
			tracksPt.append((i,e.ECAL_tracksPt[i]))
	for i in range(len(e.ECALadj_tracksPt)):
		if(e.ECALadj_tracksPt[i]>0.001):
			tracksPtadj.append((i,e.ECALadj_tracksPt[i]))
	#tracksPtadj.append((200,0.27))
	print set(tracksPt).symmetric_difference(set(tracksPtadj))
	print 1,tracksPt[index:index+length]
	print 2,tracksPtadj[index:index+length]


	# print 1,e.ECAL_tracksPt.size(),list(e.ECAL_tracksPt[index:index+length])
	# print 2,e.ECALadj_tracksPt.size(),list(e.ECALadj_tracksPt[index:index+length])
	print '\n'

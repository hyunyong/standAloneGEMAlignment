from ROOT import *
from array import array

tf = TFile("gemRsidual.root","RECREATE")
tt = TTree("residual","dx")
res = array('f', [0])
dID = array('i', [0])
tt.Branch("dx", res, "dx")
tt.Branch("dID", dID, "dID/I")
r = TRandom()
for x in range(5000):
  res[0] = r.Gaus(0.5,1)
  dID[0] = 116 
  tt.Fill()
for x in range(5000):
  res[0] = r.Gaus(-0.3,1)
  dID[0] = -118 
  tt.Fill()


tf.Write()
tf.Close()

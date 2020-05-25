from ROOT import *
from array import array
import ctypes 
import csv

gROOT.SetBatch()

tf = TFile("initial.root")

def MuonResidualsFitter_logPureGaussian(residual, center, sigma):
  sigma = abs(sigma)
  cgaus = 0.5 * TMath.Log( 2.*TMath.Pi() )
  return (-pow(residual - center, 2) *0.5 / sigma / sigma) - cgaus - TMath.Log(sigma)

def MuonResiduals3DOFFitter_FCN(npar, gin, fval, par, iflag):
  fval = 0.
  dx = par[0]
  dy = par[1]
  dphiz = par[2]
  residsigma = par[3]
  for e in tt:
    fval += MuonResidualsFitter_logPureGaussian(e.RdPhi_CSC_GE11, getResidual3DOF(dx, dy, dphiz, e.prop_CSC_localx_GE11, e.prop_CSC_r_GE11, e.prop_CSC_localy_GE11), residsigma)
  return fval

def getResidual3DOF(delta_x, delta_y, delta_phiz, track_x, R, track_y): 
  return delta_x - (track_x/R - 3.*pow(track_x/R, 3)) * delta_y - track_y * delta_phiz 

def doFit3DOF(DoDx = True, DoDy = True, DoDphiz = True):
  tmpMean = 0.
  tmpMean2 = 0.
  for e in tt:
    tmpMean += e.RdPhi_CSC_GE11
    tmpMean2 += e.RdPhi_CSC_GE11*e.RdPhi_CSC_GE11
  n = tt.GetEntries()
  print n, tmpMean, tmpMean2
  par = [tmpMean/n, 0.0, 0.0, TMath.Sqrt(tmpMean2/n - pow(tmpMean/n,2))]
  mfit = TMinuit(4)
  mfit.SetFCN(MuonResiduals3DOFFitter_FCN)
  mfit.DefineParameter(0, "dx", par[0], 0.01, 0, 0)
  mfit.DefineParameter(1, "dy", par[1], 0.1, 0, 0)
  mfit.DefineParameter(2, "dphiz", par[2], 0.001, 0, 0)
  mfit.DefineParameter(3, "sig", par[3], 0.01*par[2], 0, 0)
  mfit.FixParameter(3)
  if not DoDx: mfit.FixParameter(0)
  if not DoDy: mfit.FixParameter(1)
  if not DoDphiz: mfit.FixParameter(2)

  arglist = array('d', 10*[0.])  
  arglist[0] = 0.5
  ierflg = ctypes.c_int(0)
  smierflg = ctypes.c_int(0)
  mfit.mnexcm("SET ERR", arglist, 1, ierflg)
  for i in range(10): arglist[i] = 0.0
  arglist[0] = 1
  mfit.mnexcm("SET STR", arglist, 1, ierflg)
  for i in range(10): arglist[i] = 0.0
  arglist[0] = 5000
  mfit.mnexcm("MIGRAD", arglist, 1, ierflg)
  try_again = False
  if ierflg != 0: try_again = True
  if try_again: 
    print "2nd fit"
    for i in range(10): arglist[i] = 0.0
    arglist[0] = 5000
    mfit.mnexcm("MIGRAD", arglist, 1, smierflg)

  result = [Double(0), Double(0), Double(0), Double(0)]
  error = [Double(0), Double(0), Double(0), Double(0)]
  mfit.GetParameter(0, result[0], error[0])
  mfit.GetParameter(1, result[1], error[1])
  mfit.GetParameter(2, result[2], error[2])
  return result
    
if __name__ == '__main__':
  DoDx = True
  DoDy = False
  DoDphiz = False
  f = open("fitter.csv", 'w')
  w = csv.writer(f)

  for j in [-1, 1]:
    for i in range(36):
      print "At chamber ", i
      detNum = j*(i+101)
      ttf = TFile("tmp","recreate")
      tt = tf.Get("analyser/MuonData").CopyTree("det_id=={}".format(detNum))
      fitResults = doFit3DOF(DoDx,DoDy,DoDphiz)
      dx = "%.3f" % fitResults[0]
      dy = "%.3f" % fitResults[1]
      dz = 0
      dphix = 0
      dphiy = 0
      dphiz = "%.3f" % fitResults[2]
      entries = tt.GetEntries()
      w.writerow([detNum, dx, dy, dz, dphix, dphiy, dphiz, entries])


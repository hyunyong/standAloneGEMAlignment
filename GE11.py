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

def MuonResiduals1DOFFitter_FCN(npar, gin, fval, par, iflag):
  fval = 0.
  residpeak = par[0]
  residsigma = par[1]
  for e in tt:
    fval += MuonResidualsFitter_logPureGaussian(e.RdPhi_CSC_GE11, residpeak, residsigma) 
  return fval

def MuonResiduals2DOFFitter_FCN(npar, gin, fval, par, iflag):
  fval = 0.
  dx = par[0]
  dphiz = par[1]
  residsigma = par[2]
  for e in tt:
    fval += MuonResidualsFitter_logPureGaussian(e.RdPhi_CSC_GE11, getResidual2DOF(dx, dphiz, e.prop_CSC_localy_GE11), residsigma)
  return fval

def getResidual2DOF(delta_x, delta_phiz, track_y): 
  return delta_x - track_y * delta_phiz 

def doFit():
  tmpMean = 0.
  tmpMean2 = 0.
  for e in tt:
    tmpMean += e.RdPhi_CSC_GE11
    tmpMean2 += e.RdPhi_CSC_GE11*e.RdPhi_CSC_GE11
  n = tt.GetEntries()
  print n, tmpMean, tmpMean2
  par = [tmpMean/n, TMath.Sqrt(tmpMean2/n - pow(tmpMean/n,2))]
  mfit = TMinuit(2)
  mfit.SetFCN(MuonResiduals1DOFFitter_FCN)
  arglist = array('d', 10*[0.])  
  arglist[0] = 0.5
  ierflg = ctypes.c_int(0)
  smierflg = ctypes.c_int(0)
  mfit.DefineParameter(0, "dx", par[0], 0.01*par[1], 0, 0)
  mfit.DefineParameter(1, "sig", par[1], 0.01*par[1], 0, 0)
  mfit.FixParameter(1)
  mfit.mnexcm("SET ERR", arglist, 1, ierflg)
  arglist[0] = 5000
  mfit.mnexcm("MIGRAD", arglist, 1, ierflg)
  result = [Double(0), Double(0)]
  mfit.GetParameter(0, result[0], result[1])
  return result

def doFit2DOF():
  tmpMean = 0.
  tmpMean2 = 0.
  for e in tt:
    tmpMean += e.RdPhi_CSC_GE11
    tmpMean2 += e.RdPhi_CSC_GE11*e.RdPhi_CSC_GE11
  n = tt.GetEntries()
  print n, tmpMean, tmpMean2
  par = [tmpMean/n, 0.0, TMath.Sqrt(tmpMean2/n - pow(tmpMean/n,2))]
  mfit = TMinuit(3)
  mfit.SetFCN(MuonResiduals2DOFFitter_FCN)
  arglist = array('d', 10*[0.])  
  arglist[0] = 0.5
  ierflg = ctypes.c_int(0)
  smierflg = ctypes.c_int(0)
  mfit.DefineParameter(0, "dx", par[0], 0.1, 0, 0)
  mfit.DefineParameter(1, "dphiz", par[1], 0.001, 0, 0)
  mfit.DefineParameter(2, "sig", par[2], 0.01*par[2], 0, 0)
  mfit.FixParameter(2)
  mfit.mnexcm("SET ERR", arglist, 1, ierflg)
  arglist[0] = 5000
  mfit.mnexcm("MIGRAD", arglist, 1, ierflg)
  result = [Double(0), Double(0), Double(0)]
  error = [Double(0), Double(0), Double(0)]
  mfit.GetParameter(0, result[0], error[0])
  mfit.GetParameter(1, result[1], error[1])
  return result
 
    
if __name__ == '__main__':

  f = open("fitter.csv", 'w')
  w = csv.writer(f)

  for j in [-1, 1]:

    for i in range(36):
      print "At chamber ", i
      detNum = j*(i+101)
      ttf = TFile("tmp","recreate")
      tt = tf.Get("analyser/MuonData").CopyTree("det_id=={}".format(detNum))
      #tt = tf.Get("analyser/MuonData").Clone("det{}".format(detNum))
      #tt.Draw(">>eventList", "det_id=={}".format(detNum), "")
      #el = tt.GetEntryNumber(gDirectory.Get("eventList"))
      #tt.SetEventList(el)
      fitResults = doFit2DOF()
      print "Going to doFit()"
      dx = "%.3f" % fitResults[0]
      print "Finished doFit()"
      dy = 0
      dz = 0
      dphix = 0
      dphiy = 0
      dphiz = "%.3f" % fitResults[1]
      entries = tt.GetEntries()
      w.writerow([detNum, dx, dy, dz, dphix, dphiy, dphiz, entries])



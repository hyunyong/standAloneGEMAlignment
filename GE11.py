from ROOT import *
from array import array
import ctypes 

tf = TFile("gemRsidual.root")


def MuonResidualsFitter_logPureGaussian(residual, center, sigma):
  sigma = abs(sigma)
  cgaus = 0.5 * TMath.Log( 2.*TMath.Pi() )
  return (-pow(residual - center, 2) *0.5 / sigma / sigma) - cgaus - TMath.Log(sigma)


def MuonResiduals1DOFFitter_FCN(npar, gin, fval, par, iflag):
  fval = 0.
  residpeak = par[0]
  residsigma = par[1]
  for e in tt:
    fval += MuonResidualsFitter_logPureGaussian(e.dx, residpeak, residsigma) 

  return fval


def doFit():
  tmpMean = 0.
  tmpMean2 = 0.
  for e in tt:
    tmpMean += e.dx
    tmpMean2 += e.dx*e.dx
  n = tt.GetEntries()
  print n
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
    
if __name__ == '__main__':

  tt = tf.Get("residual").CopyTree("dID==116")
  print doFit()
  tt = tf.Get("residual").CopyTree("dID==-118")
  print doFit()



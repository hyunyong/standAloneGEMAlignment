from ROOT import *
import csv
gROOT.SetBatch()
gStyle.SetOptStat(0)
c = TCanvas("","",800,600)
gemAl = csv.reader(open('initial.csv','r'))
fitter = csv.reader(open('iteration1_inner.csv','r'))
iter1 = csv.reader(open('iteration2_inner.csv','r'))
iter2 = csv.reader(open('iteration3_inner.csv','r'))

gemAlMap = {}
fitterMap = {}
iter1Map = {}
iter2Map = {}

for l in gemAl:
  gemAlMap[l[0]] = l[1:]

for l in fitter:
  fitterMap[l[0]] = l[1:-1]

for l in iter1:
  iter1Map[l[0]] = l[1:]

for l in iter2:
  iter2Map[l[0]] = l[1:]


h1 = TH1D("result1", "GEM Alignment (tracker)",72,0,72)
h2 = TH1D("result2", "GEM Alignment",72,0,72)
h3 = TH1D("result3", "GEM Alignment",72,0,72)
h4 = TH1D("result4", "GEM Alignment",72,0,72)

for i in range(1,73):
  if i < 37: h1.GetXaxis().SetBinLabel(i,"GE-1/1 ch{}".format(i))
  else: h1.GetXaxis().SetBinLabel(i,"GE1/1 ch{}".format(i-36))
av0 = 0.
av1 = 0.
av2 = 0.
av3 = 0.
for i in range(1,73):
  if i < 37: v = gemAlMap[str(-100-i)][0]
  else: v = gemAlMap[str(100+i-36)][0]
  av0 += abs(float(v))
  h1.SetBinContent(i,abs(float(v)))
for i in range(1,73):
  if i < 37: v = fitterMap[str(-100-i)][0]
  else: v = fitterMap[str(100+i-36)][0]
  av1 += abs(float(v))
  h2.SetBinContent(i,abs(float(v)))
for i in range(1,73):
  if i < 37: v = iter1Map[str(-100-i)][0]
  else: v = iter1Map[str(100+i-36)][0]
  av2 += abs(float(v))
  h3.SetBinContent(i,abs(float(v)))
for i in range(1,73):
  if i < 37: v = iter2Map[str(-100-i)][0]
  else: v = iter2Map[str(100+i-36)][0]
  av3 += abs(float(v))
  h4.SetBinContent(i,abs(float(v)))


h1.SetMarkerStyle(20)
h2.SetMarkerStyle(21)
h3.SetMarkerStyle(22)
h4.SetMarkerStyle(23)
h1.SetMarkerColor(kBlack)
h2.SetMarkerColor(kBlue)
h3.SetMarkerColor(kRed)
h4.SetMarkerColor(kViolet)

h1.Draw("p")
h2.Draw("p same")
h3.Draw("p same")
h4.Draw("p same")

le = TLegend(0.6, 0.8, 0.9, 0.89)
le.SetFillStyle(0)
le.SetBorderSize(0)
le.AddEntry(h1, "Input Geometry  mean: {:.2f} #mu m".format(av0/72.*10.*1000.))
le.AddEntry(h2, "1st Results mean: {:.2f} #mu m".format(av1/72.*10.*1000.))
le.AddEntry(h3, "2nd Results (iter1) mean: {:.2f} #mu m".format(av2/72.*10.*1000.))
le.AddEntry(h4, "3rd Results (iter2) mean: {:.2f} #mu m".format(av3/72.*10.*1000.))
le.Draw()
c.SaveAs("GEMAlTrk.png")

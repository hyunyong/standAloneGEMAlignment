import csv

def sumEle(key, map1, map2):
  tmp = map1[key]
  try: tmp2 = map2[key]
  except: tmp2 = [0,0,0,0,0,0]
  tmp3 = []
  for i,x in enumerate(tmp):
    tmp3.extend([float(x)+float(tmp2[i])])
  return tmp3

gemAl = csv.reader(open('gemAl.csv','r'))
fitter = csv.reader(open('fitter.csv','r'))

gemAlMap = {}
fitterMap = {}
for l in gemAl:
  gemAlMap[l[0]] = l[1:]

for l in fitter:
  fitterMap[l[0]] = l[1:-1]

out = csv.writer(open("iter.csv",'w'))
for i in gemAlMap.keys():
  out.writerow([int(i)]+sumEle(i, gemAlMap, fitterMap))
  

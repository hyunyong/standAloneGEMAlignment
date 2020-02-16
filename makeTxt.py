import csv, random

f = open("gemAl.csv", 'w')
w = csv.writer(f)

for r in [-1,1]:
  for s in [1,2]:
    if s == 1: nsec = 36
    if s == 2: nsec = 18
    for c in range(nsec):
      detNum = r*(s*100+c+1)
      dx = float("%1.3f"%random.gauss(0,1))
      dy = 0.0
      dz = 0.0
      dpx = 0.0
      dpy = 0.0
      dpz = 0.0
      w.writerow([detNum, dx, dy, dz, dpx, dpy, dpz])
      

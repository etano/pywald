import sys,os
import math,cmath
from numpy import dot, pi, zeros, array, linspace, log, logspace
from random import uniform

def mag(k):
  return math.sqrt(dot(k,k))

def Genps(species,N,D,b,r):
  p = zeros((len(species)))
  for s in range(len(species)):
    for i in range(N):
      C = zeros((d,2*mMax+1))
      for d in range(D):
        Cd = cmath.exp(complex(0.,dot(b[d],r[s,i])))
        for n in range(-nMax,nMax):
          C[d,n+nMax] = Cd**n
      for n in ns:
        prod = 1.
        for d in range(D):
          prod *= C[d,n[d]+nMax]**n[d]
        p[s] += prod
  return ps

def Genns(nMax,D):
  #ns = zeros(((2*nMax-1)**D,D))
  ns = []
  for n1 in range(-nMax,nMax+1):
    if D == 1:
      ns.append(array([n1]))
    else:
      for n2 in range(-nMax,nMax+1):
        if D == 2:
          ns.append(array([n1,n2]))
        else:
          for n3 in range(-nMax,nMax+1):
            ns.append(array([n1,n2,n3]))
  return ns

def ShortRanged(L,D,rs):
    V = lambda r: 1./r

    Vs = []
    for r in rs:
        Vs.append(V(r))

    return [Vs]

def GetUnique(ks, fs):
  uniqueKs, uniqueFs = [], []
  for i in range(len(ks)):
    if not (ks[i] in uniqueKs):
      uniqueKs.append(ks[i])
      uniqueFs.append(fs[i])
  return uniqueKs, uniqueFs

def TraditionalEwald(L,D,rs,nMax):
    alpha = 7./(L)

    VShort = lambda r,a: math.erfc(a*r)/r
    fVLong = lambda k2,a: 4*pi*math.exp(-k2/(4*a*a))/(k2*(L**D))

    # Short-ranged r-space
    Vss = []
    for r in rs:
      Vs = VShort(r,alpha)
      Vss.append(Vs)

    # Long-ranged k-space
    ns = Genns(nMax,D)
    ks = []
    fVls = []
    for n in ns:
      magn2 = dot(n,n)
      k2 = 4.*pi*pi*magn2/(L*L)
      ks.append(math.sqrt(k2))
      if magn2 != 0:
        fVl = fVLong(k2,alpha)
        fVls.append(fVl)
      else:
        fVls.append(0.)

    ks, fVls = GetUnique(ks,fVls)
    for (k,fVl) in sorted(zip(ks,fVls)):
      print k, fVl

    return [Vss,fVls]

def TestMadelung(L,D,nMax,args):
    ns = Genns(nMax,D)
    nMax2 = nMax*nMax

    # Fixed coordinates
    s1 = L/2
    N = 8
    x,Q = []*N,[]*N
    x.append(array([0, 0, 0]))
    x.append(array([s1, s1, 0]))
    x.append(array([s1, 0, s1]))
    x.append(array([0, s1, s1]))
    x.append(array([s1, 0, 0]))
    x.append(array([0, s1, 0]))
    x.append(array([0, 0, s1]))
    x.append(array([s1, s1, s1]))
    Q.append(1.0)
    Q.append(1.0)
    Q.append(1.0)
    Q.append(1.0)
    Q.append(-1.0)
    Q.append(-1.0)
    Q.append(-1.0)
    Q.append(-1.0)

    # Determine type of calculation
    V = 0.
    if len(args) == 1:
        # Naive potential in r-space
        Vr = lambda r: 1./r

        # Direct sum
        nMax = 10
        ns = Genns(nMax,D)
        V = 0.
        for i in range(len(x)):
          for j in range(len(x)):
            r = x[i] - x[j]
            for n in ns:
              if(sum(n*n) <= nMax2):
                rij = mag(r+n*L)
                if rij != 0:
                  V += Q[i]*Q[j]*Vr(rij)
        V *= 0.5
        print 'V', V
    else:
        dround = lambda x: math.floor(x+0.5)

        # Traditional Ewald
        rCut = 100.0
        alpha = 7./L
        vol = L**D

        VLong = lambda r,a: math.erf(a*r)/r
        VShort = lambda r,a: math.erfc(a*r)/r
        fVLong = lambda k2,a: 4*pi*math.exp(-k2/(4*a*a))/(k2*vol)

        # Short-ranged r-space
        Li = 1.0/L;
        Vs = 0.
        for i in range(0,len(x)-1):
          for j in range(i+1,len(x)):
              r = x[i] - x[j]
              for d in range(D):
                r[d] -= dround(r[d]*Li)*L # Put in box
              magr = mag(r)
              if (magr <= rCut):
                Vs += Q[i]*Q[j]*VShort(magr,alpha)
        print 'Vs', Vs

        # Long-ranged k-space
        fVl = []
        for n in ns:
          magn2 = dot(n,n)
          k2 = 4.*pi*pi*magn2/(L*L)
          if magn2 != 0:
            fVl.append(fVLong(k2,alpha))
          else:
            fVl.append(0.)

        Vl = 0.
        for i in range(len(ns)):
          n = ns[i]
          magn2 = dot(n,n)
          if(magn2 < nMax2):
            Re,Im = 0.,0.
            for j in range(len(x)):
              h = 2*pi*dot(n,x[j])/L
              Re += Q[j]*math.cos(h)
              Im -= Q[j]*math.sin(h)
            for j in range(len(x)):
              h = 2*pi*dot(n,x[j])/L
              Vl += Q[j] * (Re*math.cos(h) - Im*math.sin(h)) * fVl[i]
        Vl *= 0.5
        print 'Vl', Vl

        Vself = -sum([q*q for q in Q]) * alpha / math.sqrt(pi)
        print 'Vself', Vself

        print 'Vl+Vself', Vl+Vself
        V = Vs+Vl+Vself
        print 'V', V

    # Neutralizing Background
    rmax = rCut
    rmin = 0.00000001
    nr = 10000
    dr = (rmax-rmin)/nr
    s0 = 0.
    for r in linspace(rmin,rmax,nr):
      s0 += dr*4*pi*r*r*VShort(r,alpha)/vol

    Vb = 0.
    for i in range(N):
      Vb -= 0.5*N*N*s0
    print 'Vb', Vb

    # Madelung Constant
    estMad = L*V/N
    print 'Estimated madelung constant:', estMad
    extMad = -1.747564594633182190636212035544397403481
    print 'Exact madelung constant', extMad
    print 'Relative error:', abs(estMad-extMad)/extMad

def usage():
  print "Usage: python %s L D rMin rMax rGrid nPts nMax type" % os.path.basename(sys.argv[0])
  sys.exit(2)

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()

  try:
    L = float(sys.argv[1])
    D = int(sys.argv[2])
    rMin = float(sys.argv[3])
    rMax = float(sys.argv[4])
    rGrid = int(sys.argv[5])
    nPts = int(sys.argv[6])
    nMax = int(sys.argv[7])
    type = int(sys.argv[8])
  except:
    usage()

  if rGrid == 0: # Linear grid
    rs = linspace(rMin,rMax,nPts)
  elif rGrid == 1: # Log grid
    rs = logspace(log(rMin),log(rMax),nPts)
  else:
    sys.stderr.write('ERROR: Unknown grid type')
    usage()

  if type == 0: # No breakup
    Vs = ShortRanged(L,D,rs)
  elif type == 1: # Traditional Ewald
    Vs = TraditionalEwald(L,D,rs,nMax)
  elif type == 2: # Optimized Ewald
    print 'Optimized Ewald'
  else:
    sys.stderr.write('ERROR: Unknown break up type')
    usage()

  TestMadelung(L,D,nMax,Vs)

if __name__ == "__main__":
  sys.exit(main())



#!/usr/bin/python

### partially edited by Ravindra
import sys
import math

#Usage:
if (len(sys.argv)<3):
  print "Usage:"
  print "./molcas2champ.py (Molcas output) (Molcas orbital output) <threshold>"
  print "Remember to use >>EXPORT MOLCAS_PRINT=4 before the CASPT2 program."
  print "You can kill the CASPT2 right after it prints the CSF's (at the beginning)."
  print "Use SPHERICAL d-functions, i.e., -2, -1, 0, +1, +2, ..."
  print "f-functions are NOT implemented. To implement it, maybe add some more entries"
  print "to the dictionary ""dfuncs"", starting from 5 (to put them behind d-functions)."
  print "Good luck next time!"
  sys.exit()

# Open the input files
inputf=file(sys.argv[1])
inporbf=file(sys.argv[2])

#Third part of input: Threshold
if (len(sys.argv)>3):
  thresh=float(sys.argv[3])
else:
  thresh=0.0

# Functions definition

# Define orbitals
def gauss_s(x, alpha): return math.exp(-x*x*alpha)*(2.0*alpha)**0.75*2.0*math.pi**(-0.25)
def gauss_p(x, alpha): return math.exp(-x*x*alpha)*(2.0*alpha)**1.25*math.sqrt(8.0/3.0)*math.pi**(-0.25)
def gauss_d(x, alpha): return math.exp(-x*x*alpha)*(2.0*alpha)**1.75*math.sqrt(16.0/15.0)*math.pi**(-0.25)
def gauss_f(x, alpha): return math.exp(-x*x*alpha)*(2.0*alpha)**2.25*math.sqrt(32.0/105.0)*math.pi**(-0.25)

# Which type of run?
def run_type(inputf):
  line=inputf.readline()
  while ( (len(line.split())<4) or (line.split()[2]!='module') or (line.split()[3] not in ('RASSCF','SCF')) ):
    line=inputf.readline()
  run=line.split()[3]
  print "Detected a ", run, "run!"
  return run;

# Number of symmetry species:
def symmspec(inputf):
  line=inputf.readline()
  while ( (len(line.split())< 3) or (line.split()[1]+line.split()[2]!='Orbitalspecifications:') ):
    line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  nsymmspec = len(line.split())-2
  if (run=='RASSCF'):
     for l in range(4):
       line=inputf.readline(); l=l+1
     nactsymspec= 0
     for i in range(nsymmspec):
       if (int(line.split()[2+i])!=0):
         nactsymspec=nactsymspec+1
     for l in range(5):
       line=inputf.readline();l=l+1
     ndeleted=0
     symdel=[]
     for i in range(nsymmspec):
       symdel.append(int(line.split()[2+i]))
       ndeleted=ndeleted+int(line.split()[2+i])
  else:
     nactsymspec=0
  if (run=='SCF'):
    for l in range(7):
      line=inputf.readline();l=l+1
    ndeleted=0
    symdel=[]
    for i in range(nsymmspec):
      symdel.append(int(line.split()[2+i]))
      ndeleted=ndeleted+int(line.split()[2+i])
  print "Number of symmetry species:", nsymmspec
  print "Number of symmetry species with active orbitals:", nactsymspec
  print "Number of deleted orbitals:",ndeleted
  return nsymmspec, nactsymspec, ndeleted, symdel

#Set 1: Call functions
#run= run_type(inputf)
run = "SCF"
nsymmspec, nactsymspec, ndeleted, symdel=symmspec(inputf)

#Get the number of basis functions and number per irreducible representation:
inputf=file(sys.argv[1])
def numbas(inputf):
    for line in inputf:
        if line.startswith("++    Orbital specifications:"):
            while not line.startswith('--'):
                if line.strip().startswith('Number of basis functions'):
                    nbasis = line.split()[4:]
                    symmetry_count = len(nbasis)
                    nbas = sum(map(int, nbasis))
                    bas = []
                    for i in range(4,symmetry_count+4):
                        bas.append(line.split()[i])
                    bas = list(map(int, bas))

                # Get the total number of orbitals
                if line.strip().startswith('Total number of orbitals'):
                    nmos = line.split()[4:]
                    norb = sum(map(int, nmos))

                line = next(inputf)
    print (" Number of Basis functions     ::  ", nbas)
    print (" Basis functions per irrep     ::  ", bas)
    print (" Total Number of Orbitals      ::  ", norb)
    return nbas,bas, norb






#Set 2: Call functions
nbas, bas, norb =numbas(inputf)

#Create and read for geometry file
geo_file=open(sys.argv[1][:-4]+'.geo','w')
inputf=file(sys.argv[1])
def geometry(inputf):
  l=1
  line=inputf.readline()
  while ( (len(line.split())<3) or (line.split()[0]+line.split()[1]+line.split()[2]!='****CartesianCoordinates')):
    line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  atomweights={'O':6.0, 'N':5.0, 'C':4.0, 'H':1.0, 'S':6.0, 'LI':1.0}
  nctype=0
  natom=0
  types=0
  prevlabel=''
  labels=[]
  centers=[]
  while (len(line.split())>1):
    natom=natom+1
    label=str(line.split()[1][:-1])
    if (prevlabel!=label):
      prevlabel=label
      labels.append(label)
      nctype=nctype+1
    x=float(line.split()[2])
    y=float(line.split()[3])
    z=float(line.split()[4])
    centers.append((x,y,z,nctype))
    line=inputf.readline()
  return nctype, natom, labels, centers, atomweights;

#Set 3: Call function
nctype, natom, labels, centers, atomweights =geometry(inputf)

#Create basis function files and bfinfo file, Defining a function is probably not the easiest way here: Too many dependencies!
bfinfo_file=open(sys.argv[1][:-4]+'.bfinfo','w')
bfinfo_file.write('qmc_bf_info 1\n')
line=inputf.readline()
while ( (len(line.split())<5) or (line.split()[1]+line.split()[4]!='Primitive(Valence)')):
  line=inputf.readline()
line=inputf.readline()
line=inputf.readline()
line=inputf.readline()
atoms=[]
for label in labels:
  basisfcts=[]
  line=inputf.readline()
  if (label=='H'):
    nums=int(line.split()[1][19])#Modify
    nump=int(line.split()[1][21])#Modify
    numd=0
  else:
    nums=int(line.split()[1][23])#Modify depending on your basis set
    nump=int(line.split()[1][25])#-do-
    numd=int(line.split()[1][27])#-do-
  if (label!=line.split()[1][4]):
    if (len(label.split())>1):
      label=line.split()[1][4:5]
    else:
      print "mismatch!", label, line.split()[1]
  atoms.append([label, nums, nump, numd])
  line=inputf.readline()
  line=inputf.readline()
  while ((len(line.split())>0) and line.split()[0]=='Type'):
    line=inputf.readline()
    type=line.split()[0]
    line=inputf.readline()
    line=inputf.readline()
    primitives=[]
    bcoeffs=[]
    while (len(line.split())>0):
      primitives.append(float(line.split()[1].replace('D','E')))
      bcoeffs.append ([float(x.replace('D','E')) for x in line.split()[2:]])
      line=inputf.readline()
    for ifct in range(len(bcoeffs[0])):
      curr_fct=[type]
      for iprim in range(len(primitives)):
        curr_fct.append((primitives[iprim],bcoeffs[iprim][ifct]))
      basisfcts.append(curr_fct)
    line=inputf.readline()
  bas_file=open(sys.argv[1][:-4]+'.basis.'+label,'w')
  bas_file.write("%i"%len(basisfcts)+" 3 2000 1.003000 20.000000 0\n")
  npoints=2000
  gridarg=1.003
  gridr0=20.0
  gridr0=gridr0/(gridarg**(npoints-1)-1)
  for i in range(npoints):
    r=gridr0*gridarg**i-gridr0
    outstr="%20.12E"%r
    for basisfct in basisfcts:
      value=0.0
      type=basisfct[0]
      for element in basisfct[1:]:
        if (type=='s'):
          value=value+element[1]*gauss_s(r,element[0])
        elif (type=='p'):
          value=value+element[1]*gauss_p(r,element[0])*r
        elif (type=='d'):
          value=value+element[1]*gauss_d(r,element[0])*r**2
        elif (type=='f'):
          value=value+element[1]*gauss_f(r,element[0])*r**3
        else:
          print "I don't know what this type is! ",type
      if (value<1e-15):
        value=0.0
      outstr=outstr+"%20.12E"%value
    bas_file.write(outstr+"\n")
  s_str=""
  p_str=""
  d_str=""
  f_str=""
  ns=0
  np=0
  nd=0
  nf=0
  i=0
  for basisfct in basisfcts:
    i=i+1
    if (basisfct[0]=='s'):
      ns=ns+1
      s_str=s_str+"%3i"%i
    elif (basisfct[0]=='p'):
      np=np+1
      p_str=p_str+"%3i"%i
    elif (basisfct[0]=='d'):
      nd=nd+1
      d_str=d_str+"%3i"%i
    elif (basisfct[0]=='f'):
      nf=nf+1
      f_str=f_str+"%3i"%i
    else:
      print "I don't know what this type is! ",basisfct[0]
  bfinfostr="%3i"%(-ns)+"  0"
  for i in range(3):
    bfinfostr=bfinfostr+"%3i"%(-np)
  bfinfostr=bfinfostr+"  0  0  0  0"
  for i in range(5):
    bfinfostr=bfinfostr+"%3i"%(-nd)
  bfinfostr=bfinfostr+"  0  0  0  0"
  for i in range(7):
    bfinfostr=bfinfostr+"%3i"%(-nf)
  for i in range(12):
    bfinfostr=bfinfostr+"  0"
  bfinfo_file.write(bfinfostr+"\n")
  bfinfostr=s_str
  for i in range(3):
    bfinfostr=bfinfostr+p_str
  for i in range(5):
    bfinfostr=bfinfostr+d_str
  for i in range(7):
    bfinfostr=bfinfostr+f_str
  bfinfo_file.write(bfinfostr+"\n")
bfinfo_file.write("end\n")

#Create a list of tuples of equivalent centers
def equiv_centers(inputf):
  line=inputf.readline()
  while ( (len(line.split())<2) or ((line.split()[0]+line.split()[1])!='Irreduciblerepresentation')):
      line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  centers_list=[]
  while (len(line.split())>=5 and not 'Basis function(s)' in line):
    if (len(line.split())==5):
      if (line.split()[3] not in centers_list):
        centers_list.append((line.split()[3]))
      else:
        line=inputf.readline()
      line=inputf.readline()
    else:
      if (line.split()[3] not in centers_list):
        sublist=list()
        for i in range(3,len(line.split()),2):
          sublist.append((line.split()[i]))
        centers_list.append(sublist)
      else:
        line=inputf.readline()
      line=inputf.readline()
    line=inputf.readline()
  equiv_center=[]
  for i in centers_list:
    if i not in equiv_center:
      equiv_center.append(i)
  return equiv_center

#####
# Get orbital labels(kept in the original format): dictbas update
inputf=file(sys.argv[1])
def orb_label(inputf):
  l=0
  line=inputf.readline();l=l+1
  if (nsymmspec==1):
    while ( (len(line.split())<2) or (line.split()[1]!='Petite')):
      line=inputf.readline();l=l+1
    for n in range(7):
      line=inputf.readline();l=l+1
    #order of d-functions:
    dfuncs={'d0': 0, 'd2+': 1, 'd2-': 2, 'd1+': 3, 'd1-': 4}
    prevcent=int(line.split()[3])
    atombasis=list()
    i=0
    total_basis=0
    dictbas={}
    while ( len(line.split())==4  ):
      i=i+1
      type=line.split()[2][1:]
      center=int(line.split()[3])
      priority=-1
      if type in dfuncs.keys():
        priority=dfuncs[type]
      if (center!=prevcent):
        sortatombasis=sorted(atombasis,key=lambda x: x[1])
        j=0
        for element in sortatombasis:
          j=j+1
          dictbas.update({element[0]:total_basis+j})
        total_basis=total_basis+j
        atombasis=[(i,priority),]
      else:
        atombasis.append((i,priority))
      prevcent=center
      line=inputf.readline();l=l+1
    sortatombasis=sorted(atombasis,key=lambda x: x[1])#For the last atom center
    j=0
    for element in sortatombasis:
      j=j+1
      dictbas.update({element[0]:total_basis+j})
    total_basis=total_basis+j
    print "Built correct order of orbitals"
    return dictbas;
  else:
    while ( (len(line.split())<2) or ((line.split()[0]+line.split()[1])!='Irreduciblerepresentation')):
      line=inputf.readline();l=l+1
    symclass=line.split()[3]
    for n in range(4):
      line=inputf.readline();l=l+1
    pfuncs={'px': 0, 'py' : 1, 'pz' : 2}
    dfuncs={'d0': 3, 'd2+': 4, 'd2-': 5, 'd1+': 6, 'd1-': 7}
    prevcent=int(line.split()[3])
    prevlabel=line.split()[1][0]
    atombasis_full=list()
    i=0
    total_basis=0
    dictbas={}
    atomcenter={}
    atomcenter.update({prevcent:prevlabel})
    types=[]
    centers=[]
    centers.append(prevcent)
    while ( len(line.split())>=5  ):
      if (len(line.split())>2):
        atomlabel=line.split()[1][:-1]
      else:
        atomlabel=line.split()[1][0]
      type=line.split()[2][1:]
      types.append(type)
      priority=-1
      if type in pfuncs.keys():
        priority=pfuncs[type]
      if type in dfuncs.keys():
        priority=dfuncs[type]
      for k in range(3,len(line.split()),2):
        i=i+1
        center=int(line.split()[k])
        if (center not in centers):
          centers.append(center)
          atomcenter.update({center:atomlabel})
          atombasis_full.append((center, i,priority))
        else:
          atombasis_full.append((center,i,priority))
      line=inputf.readline();l=l+1
    line=inputf.readline
    line=inputf.readline
    sortatombasis_full=sorted(atombasis_full,key=lambda x: (x[0], x[2]))
    j=0
    for element in atombasis_full:
      j=j+1
      dictbas.update({element[1]:total_basis+j})
    total_basis=total_basis+j
    from itertools import groupby
    sortatombasis_list=[list(g) for k, g in groupby(sortatombasis_full, lambda s: s[0])]
    for k in range(len(sortatombasis_list)):
      sortatombasis=sortatombasis_list[k]
    print "Built correct order of orbitals"
    if (len(dictbas)< nbas):
      nzero= nbas-len(dictbas)
      print nzero,"AOs have 0 coefficients in symmetry class",symclass
      j=0
      for k in range(len(dictbas)+1,nbas+1):
         j=j+1
         dictbas.update({k:total_basis+j})
    return dictbas,sortatombasis_list,atomcenter;

#Count the number of states for the RASSCF runs
def state(inputf):
  if (run=='RASSCF'):
    line=inputf.readline()
    while ( (len(line.split())!=9) or (line.split()[2]!='CI-coefficients') ):
      line=inputf.readline()
    line=inputf.readline()
    line=inputf.readline()
    line=inputf.readline()
    csflist=[]
    states=0
    cicoeffs=[]
    while (len(line.split())==3+nactsymspec):
      states=states+1
      print "reading in state",states
      cicoeff_state={}
      while (len(line.split())==3+nactsymspec):
        strcsf=line.split()[1]
        for i in range(2,nactsymspec+1):
          strcsf=strcsf+line.split()[i]
        ccsf=float(line.split()[nactsymspec+1])
        cicoeff_state.update({strcsf:ccsf})
        if ((abs(ccsf)>=thresh) and (strcsf not in csflist)):
          csflist.append(strcsf)
        line=inputf.readline()
      line=inputf.readline()
      line=inputf.readline()
      line=inputf.readline()
      line=inputf.readline()
      cicoeffs.append(cicoeff_state)
    for strcsf in csflist:
      istate=0
      for statedict in cicoeffs:
        istate=istate+1
        if (strcsf not in statedict):
          print "unknown CI coefficient due to a threshold that is too brutal"
          print "state:",istate
          print "CSF:",strcsf
          print "Will assume 0.00 and carry on."
          statedict.update({strcsf:0.00})
    nstates=states
    print "States:",nstates
    print "CSF's above the threshold:",len(csflist)
  return nstates, csflist, cicoeffs;

#CSFs
def csf_rasscf(inputf):
  line=inputf.readline()
  while ( (len(line.split())!=6) or (line.split()[2]+line.split()[3]+line.split()[4]!='closedshellelectrons')):
    line=inputf.readline()
  line=inputf.readline()
  actel=int(line.split()[6])#got it
  line=inputf.readline()
  line=inputf.readline()
  line=inputf.readline()
  inactorb=int(line.split()[4])
  line=inputf.readline()
  actorb=int(line.split()[4])
  line=inputf.readline()
  extorb=int(line.split()[4])
  while ( (len(line.split())!=4) or (line.split()[2]!='CSFs')):
    line=inputf.readline()
  ncsf=int(line.split()[3])
  print "It looks like we have a CAS(",actel,",",actorb,") with ",inactorb,"inactive and ",extorb," external orbitals."
  frozenslate=""
  for i in range(0,inactorb):
    frozenslate=frozenslate+"%4s"%(i+1)
#frozenslate is a string that has the integers from core orbitals, e.g. "  1  2  3  4"
  while ( (len(line.split())!=4) or ((line.split()[0]+line.split()[1]+line.split()[2])!='SGUGAinfois') ):
    line=inputf.readline()
  inputf.readline()
  slaters=()
  nslat=0
  csfmap=[]
  line=inputf.readline()
  print "starting going through CSF's in CASPT2"
  i=0
  for strcsf in csflist:
    i=i+1
    csfmap.append(())
  while (len(line.split())>0):
    strcsf=line.split()[-3]
    for j in range(1,nactsymspec):
      strcsf=line.split()[-3-j]+strcsf
    if (strcsf in csflist):
      imatch=csflist.index(strcsf)
      inputf.readline()
      line=inputf.readline()
      map=()
      while ( len(line.split())>0 ):
        sign=line.split()[0]
        coeff=line.split()[1]
        slat=line.split()[2][1:-1]
        fcoeff=float(coeff.split('/')[0])/float(coeff.split('/')[1])
        if(sign=="+"):
          fsign=1.0
        elif(sign=="-"):
          fsign=-1.0
        else:
          print "problem! no +/- in line:"
          print line
          sys.exit()
        strslata=frozenslate
        strslatb=frozenslate
        for i in range(0,actorb):
          if (slat[i]=="2"):
            strslata=strslata+"%4s"%(i+inactorb+1)
            strslatb=strslatb+"%4s"%(i+inactorb+1)
          elif (slat[i]=="a"):
            strslata=strslata+"%4s"%(i+inactorb+1)
          elif (slat[i]=="b"):
            strslatb=strslatb+"%4s"%(i+inactorb+1)
        match=-1
        for islat in range(0,nslat):
          if (slat==slaters[islat][0]):
            match=islat
            break
        if (match==-1):
          match=nslat
 #remember, it is 0-indexed
          nslat=nslat+1
          slatsign=1
          for i in range(0,actorb):
            if (slat[i]=="b"):
              for j in range(i+1,actorb):
                if (slat[j]=="a"):
                  slatsign=-slatsign
          slaters=slaters+((slat,str(strslata+"  "+strslatb),slatsign),)
        map=map+((match+1,fcoeff,fsign*slaters[match][2]),)
#still 0-indexed, so you need match+1...
        line=inputf.readline()
      csfmap[imatch]=(strcsf,map)
      line=inputf.readline()
    else:
      inputf.readline()
      line=inputf.readline()
      while (len(line.split())>0):
        line=inputf.readline()
      line=inputf.readline()
#####################################
# slaters contains a list of tuples identifying Slater determinants, like:
  # (slater_string_molcas, slater_string_champ, sign)
  # the sign can be +1 or -1 to include to the Molcas coefficient to get to CHAMP.
  ############
  # csfmap contains a list of tuples for each CSF.
  # each tuple has a list of Slater determinants and their respective coefficients.
  #####################################
  detcoeffs=[]
  i=0
  for map in csfmap:
    i=i+1
  i=0
  for state in cicoeffs:
    row=[]
    for i in range(nslat):
      row=row+[0.0]
    detcoeffs.append(list(row))
  icsf=0
  for map in csfmap:
    for slat in map[1]:
      x=math.sqrt(slat[1])*slat[2]
      islat=slat[0]-1
#back to 0-indexed stuff...
      istate=0
      for state in cicoeffs:
        detcoeffs[istate][islat]=detcoeffs[istate][islat]+x*state[map[0]]
        istate=istate+1
    icsf=icsf+1
  return actel, inactorb, nslat, detcoeffs, csfmap, slaters;

# def csf_scf(inputf):
#   line=inputf.readline()
#   while ( (len(line.split())<4) or (line.split()[2]!='module') or (line.split()[3]!=('SCF')) ):
#     line=inputf.readline()
#   while ( (len(line.split())<3) or ((line.split()[0]+line.split()[1]+line.split()[2])!='Symmetryspecies1') ):
#     line=inputf.readline()
#   lnlen=len(line.split())
#   norb=0
#   line=inputf.readline()
#   line=inputf.readline()
#   for i in range(2,lnlen):
#     norb=norb+int(line.split()[i])# Frozen orbitals
#   line=inputf.readline()#"aufbau" orbitals
#   if (line.split()[1]!= 'orbitals'):
#      norb=norb+int(line.split()[1])
#   else:
#      norb=norb+int(line.split()[2])
#   for i in range(5):
#     line=inputf.readline()
#   print "SCF calculation with ",norb," orbitals. Creating a determinant file..."
#   return norb

#Find the number of active orbitals in each symmetry species and sort them as per occupation number and index
def active_orb(inputf):
  l=0
  line=inputf.readline()
  orbital=[]
  energy=[]
  occup_no=[]
  orbtype=[]
  comb_list=[]
  counter1 = 0; counter2 = 0
  if (run=='SCF'):
    for line in inputf:
        if line.strip().startswith("Title: SCF orbitals"):
            while not line.startswith('--'):
            # At this line Orbital index are listed like "Orbital            1         2         3         4         5         6         7         8         9        10"
              tokens = line.split()
              # print ("tokens ", tokens )

              if line.strip().startswith('Molecular orbitals for symmetry species'):
                counter2 += 1

              if line.strip().startswith('Orbital'):
                orbital.extend(tokens[1:])
                counter1 += 1
                orbtype.extend(str(counter2)*len(tokens[1:]))
                line = next(inputf)

                if line.strip().startswith('Energy'):
                  tokens = line.split()
                  energy.extend(tokens[1:])
                  line = next(inputf)

                if 'Occ. No.' in line:
                  tokens = line.split()
                  occup_no.extend(tokens[2:])
                  line = next(inputf)
              line = next(inputf)

    orbital = list(map(int, orbital))
    energy  = list(map(float, energy) )
    occup_no = list(map(float, occup_no))
    orbtype = list(map(int, orbtype))

    for i in range(len(orbital)):
      comb_list.append([orbtype[i],orbital[i],energy[i],occup_no[i]])
    comb_list=sorted(comb_list,key=lambda x : (x[1],-x[2]))
    occup_list=[]
    for i in range(len(comb_list)):
      if (comb_list[i][2]>0.0):
        occup_list.append(comb_list[i])
  return occup_list,comb_list

#Create a new orbital file to read from for the final orb file; instead of ScfOrb or RasOrb
inputf=file(sys.argv[1])
inporbf=file(sys.argv[2])
def new_orb(inputf,inporbf):
  l=0; m=0
  line1=inputf.readline()
  line2=inporbf.readline()
  s_list={'s':0}
  p_list={'px':0,'py':1,'pz':2}
  d_list={'d0':0,'d2+':1,'d2-':2,'d1+':3,'d1-':4}
  ref_list=list()
  for i in range(len(comb_list)):
    print ("nbas - ndel ", i)
    ref_list.append(comb_list[i][0])
  final_list=list()
  for i in range(nsymmspec):
    while ((len(line1.split())<2) or ((line1.split()[0]+line1.split()[1])!='Irreduciblerepresentation')):
      line1=inputf.readline(); l=l+1
    for c in range(4):
      line1=inputf.readline(); l=l+1
    orb_list=[]
    for j in range(len(equiv_center)):
        s_type=[]
        p_type=[]
        d_type=[]
        print ("line1 ", line1)
        center=line1.split()[3]
        label=line1.split()[1][:-1]
        while((len(line1.split())>=5) and (line1.split()[4]=='1') and (line1.split()[3]==center)):
          type=line1.split()[2][1:]
          phase=[]
          for iphase in range(4,len(line1.split()),2):
            phase.append(line1.split()[iphase])
          if (type[0]=='s'):
            s_type.append([type,phase])
          elif (type[0]=='p'):
            p_type.append([type,phase])
          elif (type[0]=='d'):
            d_type.append([type,phase])
          else:
            print "Do not recognize orbital type!"
            sys.exit()
          line1=inputf.readline();m=m+1
        np_type=[]
        nd_type=[]
        for ptype in range(len(p_type)):
          if p_type[ptype][0] not in np_type:
             np_type.append(p_type[ptype][0])
        for dtype in range(len(d_type)):
          if d_type[dtype][0] not in nd_type:
           nd_type.append(d_type[dtype][0])
        nptype=len(np_type)
        ndtype=len(nd_type)
        orb_list.append([center,label,s_type,p_type,d_type,nptype,ndtype,np_type,nd_type])
    new_coeffs=[]
    while ( (len(line2.split())<3) or (line2.split()[1]!='ORBITAL') or (int(line2.split()[2])!=i+1) ):
        line2=inporbf.readline(); m=m+1
    for k in range(bas[i]-symdel[i]):
      symclass=int(line2.split()[2])
      orbnum=int(line2.split()[3])
      line2=inporbf.readline()
      testorb=[]
      for t in range(int(math.ceil(bas[i]/4.0))):
        extlist=[ float(line2.split()[0][n:n+18]) for n in range(0,len(line2.split()[0]),18)]
        testorb.extend(extlist)
        line2=inporbf.readline()
      lc1=0
      lc2=0
      tabline_new=[]
      for count in range(nbas):
        tabline_new.append(0.0)
      for n in range(len(orb_list)):
        center=orb_list[n][0]
        for element in equiv_center:
          if center in element[0]:
            div=int(len(element))
            copy=len(element)-1
        label=orb_list[n][1]
        for num in range(len(atoms)):
          if (atoms[num][0]==label):
            nums=atoms[num][1]
            nump=3*atoms[num][2]
            numd=5*atoms[num][3]
        for ns in range(len(orb_list[n][2])):
          tabline_new[lc2+ns]=(testorb[lc1+ns]*int(orb_list[n][2][ns][1][0]))/math.sqrt(div)
        if (orb_list[n][5]==1):
          if (orb_list[n][7]==['px']):
            for np in range(len(orb_list[n][3])):
              tabline_new[lc2+nums+np]=(testorb[lc1+len(orb_list[n][2])+np]*int(orb_list[n][3][np][1][0]))/math.sqrt(div)
          elif (orb_list[n][7]==['py']):
            for np in range(len(orb_list[n][3])):
              tabline_new[lc2+nums+(nump/3)+np]=(testorb[lc1+len(orb_list[n][2])+np]*int(orb_list[n][3][np][1][0]))/math.sqrt(div)
          elif (orb_list[n][7]==['pz']):
            for np in range(len(orb_list[n][3])):
              tabline_new[lc2+nums+(2*nump/3)+np]=(testorb[lc1+len(orb_list[n][2])+np]*int(orb_list[n][3][np][1][0]))/math.sqrt(div)
        if (orb_list[n][5]==2):
          if (orb_list[n][7]==['px','py']):
            for np in range(len(orb_list[n][3])):
              tabline_new[lc2+nums+np]=(testorb[lc1+len(orb_list[n][2])+np]*int(orb_list[n][3][np][1][0]))/math.sqrt(div)
	  elif (orb_list[n][7]==['py','pz']):
            for np in range(len(orb_list[n][3])):
              tabline_new[lc2+nums+(nump/3)+np]=(testorb[lc1+len(orb_list[n][2])+np]*int(orb_list[n][3][np][1][0]))/math.sqrt(div)
          elif (orb_list[n][7]==['px','pz']):
            for np1 in range(len(orb_list[n][3])/2):
              tabline_new[lc2+nums+np1]=(testorb[lc1+len(orb_list[n][2])+np1]*int(orb_list[n][3][np1][1][0]))/math.sqrt(div)
            for np2 in range(len(orb_list[n][3])/2):
              tabline_new[lc2+nums+(2*nump/3)+np2]=(testorb[lc1+len(orb_list[n][2])+(len(orb_list[n][3])/2)+np2]*int(orb_list[n][3][(len(orb_list[n][3])/2)+np2]))/math.sqrt(div)
        if (orb_list[n][5]==3):
          for np in range(len(orb_list[n][3])):
            tabline_new[lc2+nums+np]=(testorb[lc1+len(orb_list[n][2])+np]*int(orb_list[n][3][np][1][0]))/math.sqrt(div)
        if (orb_list[n][6]==1):
          if (orb_list[n][8]==['d0']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
        if (orb_list[n][6]==2):
          if (orb_list[n][8]==['d0','d2+']):
            for nd1 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/2)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/2)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d0']):
            for nd1 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/2)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/2)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d0','d1+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-','d0']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d2+']):
            for nd1 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/2)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/2)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1+','d2+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-','d2+']):
            for nd1 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(numd/5)+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/2)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/2)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1+']):
            for nd1 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/2)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/2)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1-']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-','d1+']):
            for nd1 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(numd/5)+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/2)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/2)+nd2][1][0]))/math.sqrt(div)
        if (orb_list[n][6]==3):
          if (orb_list[n][8]==['d2-','d0','d2+']):
            for nd1 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
            for nd3 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(2*len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(2*len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d0','d1+','d2+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-','d0','d2+']):
            for nd1 in range(2*len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(numd/5)+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(2*len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(2*len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d0','d1+']):
            for nd1 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(2*len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1-','d0']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-','d0','d1+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1+','d2+']):
            for nd1 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(2*len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1-','d2+']):
            for nd1 in range(2*len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(2*len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(2*len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-','d1+','d2+']):
            for nd1 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(numd/5)+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(2*len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1-','d1+']):
            for nd1 in range(2*len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/3):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(2*len(orb_list[n][4])/3)+nd2]*int(orb_list[n][4][(2*len(orb_list[n][4])/3)+nd2][1][0]))/math.sqrt(div)
        if (orb_list[n][6]==4):
          if (orb_list[n][8]==['d2-','d0','d1+','d2+']):
            for nd1 in range(len(orb_list[n][4])/4):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(3*len(orb_list[n][4])/4):
              tabline_new[lc2+nums+nump+(2*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/4)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/4)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1-','d0','d2+']):
            for nd1 in range(3*len(orb_list[n][4])/4):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/4):
              tabline_new[lc2+nums+nump+(4*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(3*len(orb_list[n][4])/4)+nd2]*int(orb_list[n][4][(3*len(orb_list[n][4])/4)+nd2][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d1-','d0','d1+','d2+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+(numd/5)+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1-','d0','d1+']):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
          elif (orb_list[n][8]==['d2-','d1-','d1+','d2+']):
            for nd1 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+nd1]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd1]*int(orb_list[n][4][nd1][1][0]))/math.sqrt(div)
            for nd2 in range(len(orb_list[n][4])/2):
              tabline_new[lc2+nums+nump+(3*numd/5)+nd2]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+(len(orb_list[n][4])/2)+nd2]*int(orb_list[n][4][(len(orb_list[n][4])/2)+nd2][1][0]))/math.sqrt(div)
        if (orb_list[n][6]==5):
            for nd in range(len(orb_list[n][4])):
              tabline_new[lc2+nums+nump+nd]=(testorb[lc1+len(orb_list[n][2])+len(orb_list[n][3])+nd]*int(orb_list[n][4][nd][1][0]))/math.sqrt(div)
        lc2=lc2+nums+nump+numd
        lc2c=lc2
        lc1=lc1+len(orb_list[n][2])+len(orb_list[n][3])+len(orb_list[n][4])
        lc3=0
        lc4=lc2c-nums-nump-numd
        if (copy>0):
          for icopy in range(0,copy):
            sc=0
            pc=0
            dc=0
            for ncops in range(nums):
              if (tabline_new[lc4+ncops]!=0):
                tabline_new[lc3+lc2+ncops]=tabline_new[lc4+ncops]*int(orb_list[n][2][sc][1][1+icopy])
                sc=sc+1
            for ncopp in range(nump):
              if (tabline_new[lc4+nums+ncopp]!=0):
                tabline_new[lc3+lc2+nums+ncopp]=tabline_new[lc4+nums+ncopp]*int(orb_list[n][3][pc][1][1+icopy])
                pc=pc+1
            if (numd>0):
              for ncopd in range(numd):
                if (tabline_new[lc4+nums+nump+ncopd]!=0):
                  tabline_new[lc3+lc2+nums+nump+ncopd]=tabline_new[lc4+nums+nump+ncopd]*int(orb_list[n][4][dc][1][1+icopy])
                  dc=dc+1
            lc3=lc3+nums+nump+numd
            lc2=lc2+nums+nump+numd
      final_list.append([(symclass,orbnum),tabline_new])
  final_list= sorted(final_list, key=lambda x: ref_list.index(x[0]))
  return final_list;

#open and read the orbital file:
inporbf=file(sys.argv[2])
def orbf(inporbf):
  print "I opened the RasOrb file now... fingers crossed"
  print "nbas=",nbas
  coeffs=[]
  tabline=[]
  for i in range(nbas):
    tabline.append(0.0)
  for i in range(nbas):
    coeffs.append(list(tabline))
  line=inporbf.readline()
  if (nsymmspec==1):
    while ( (len(line.split())<2) or (line.split()[1]!='ORBITAL')):
      line=inporbf.readline()
    for iorb in range(nbas):
      line=inporbf.readline()
      orbcoeffs=[]
      for j in range(int(math.ceil(1.0*nbas/4.0))):
        extlist=[ float(line.split()[0][i:i+18]) for i in range(0,len(line.split()[0]),18)]
        orbcoeffs.extend(extlist)
        line=inporbf.readline()
####
      l_edit=coeffs[iorb]
      for ibas in range(nbas):
        l_edit[dictbas[ibas+1]-1]=orbcoeffs[ibas]
      coeffs[iorb]=l_edit
  else:
    for iorb in range(nbas-ndeleted):
      orbcoeffs=final_list[iorb][1]
      l_edit=coeffs[iorb]
      for ibas in range(nbas):
        l_edit[dictbas[ibas+1]-1]=orbcoeffs[ibas]
      coeffs[iorb]=l_edit
  while ( (len(line.split())!=2) or (line.split()[1]!='1234567890')):
    line=inporbf.readline()
  deleted=0
  while ( (len(line.split())!=0)):
    for char in line.split()[1]:
      if (char=='d'):
        deleted=deleted+1
    line=inporbf.readline()
  print "deleted orbitals:",deleted
  return coeffs,deleted;

# Set 4: Call functions
dictbas_sep=[]
sortatombasis_listsep=[]
if (nsymmspec>1):
  for i in range(nsymmspec):
    dictbas,sortatombasis_list,atomcenter=orb_label(inputf)
    dictbas_sep.append(dictbas)
    sortatombasis_listsep.append(sortatombasis_list);i=i+1
else:
  dictbas=orb_label(inputf)
if (run=='RASSCF'):
  nstates, csflist, cicoeffs=state(inputf)
#If nstates>1, create mstates files
  if (nstates>1):
      mstates_file=open(sys.argv[1][:-4]+'.mstates','w')
      mstates_file.write('multiple_cistates %i\n'%nstates)
  actel, inactorb, nslat, detcoeffs, csfmap, slaters=csf_rasscf(inputf)
else:
  _, _, norb=numbas(inputf)

#Write output files
#Write geometry file
geo_file.write("&atoms nctype%3i natom%4i\n"%(nctype,natom))
geo_file.write("&atom_types")
i=0
for label in labels:
  i=i+1
  geo_file.write("%3i %s"%(i,label))
geo_file.write("\ngeometry\n")
for center in centers:
  geo_file.write("%10.6f %10.6f %10.6f %3i\n"%center)
geo_file.write("end\nznuc\n")
for label in labels:
  if (label in atomweights):
    geo_file.write("%6.2f"%atomweights[label])
  else:
    print "WARNING! Unknown center: ",label
    print "Will write 0.00"
geo_file.write("\nend\n")

#Write orbital files
outputf=open(sys.argv[1][:-4]+'.orb','w')
if (nsymmspec>1):
  inputf=file(sys.argv[1])
  inporbf=file(sys.argv[2])
  equiv_center=equiv_centers(inputf)
  print ("equivalent centers ", equiv_center)
  for element in equiv_center:
    if (len(element) > 1):
      print "Orbital coefficients for center ", element[0], "will be copied", len(element)-1,"times."
  inputf=file(sys.argv[1])
  occup_list,comb_list=active_orb(inputf)
  inputf=file(sys.argv[1])
  inporbf=file(sys.argv[2])
  final_list=new_orb(inputf,inporbf)
coeffs,deleted=orbf(inporbf)
norb=nbas-deleted
outputf.write("lcao"+"%5i"%norb+"%5i"%nbas+"   1"+"\n")
if (deleted>0):
  del(coeffs[-deleted:])
for coeff in coeffs:
  printstr=""
  for lcao in coeff:
    printstr=printstr+"%20E"%lcao
  outputf.write(printstr+"\n")
outputf.write("end\n")

#Write CSF/determinant files
csf_file=open(sys.argv[1][:-4]+'.csf','w')
if (run=='RASSCF'):
  csf_file.write("&electrons nelec %5i nup %5i\n"%(actel+2*inactorb,actel/2+inactorb))
  istate=0
  for state in cicoeffs:
    csf_file.write("# State %4i\n"%(istate+1))
    csf_file.write("determinants %4i %4i\n"%(nslat,1))
    outstring=""
    icsf=0
    for islat in range(nslat):
       outstring=outstring+"%15.8f"%detcoeffs[istate][islat]
    csf_file.write(outstring+'\n')
    if (nstates>1):
      mstates_file.write(outstring+'\n')
    istate=istate+1
    for slater in slaters:
      i=i+1
      csf_file.write(slater[1]+"\n")
  csf_file.write("end\n")
  if (nstates>1):
    mstates_file.write("end\n")
  csf_file.write( "csf %4i %4i \n"%(len(csflist),len(cicoeffs)) )
  for statedict in cicoeffs:
    outstring=""
    icsf=0
    for strcsf in csflist:
      outstring=outstring+"%15.8f"%statedict[strcsf]
    csf_file.write(outstring+'\n')
  csf_file.write("end\n")
  csf_file.write("csfmap\n")
  nslat_redundant=0
  for map in csfmap:
    nslat_redundant=nslat_redundant+len(map[1])
  csf_file.write("%5i%5i%5i\n"%(len(csflist),nslat,nslat_redundant))
  i=0
  for map in csfmap:
    i=map[0]
  #  print "Map for CSF number:",i
  #  csf_file.write("%5i\n"%i)
    j=0
    csf_file.write("%4i\n"%len(map[1]))
    for slat in map[1]:
      j=j+1
      x=math.sqrt(slat[1])*slat[2]
      csf_file.write("%4i"%slat[0]+"%10.6f"%x+'\n')
  csf_file.write("end\n")
else:
  csfstr=""
  csf_file.write("&electrons nelec %5i nup %5i\n"%(norb*2,norb))
  csf_file.write("determinants 1 1\n")
  csf_file.write("1.00000\n")
  for i in range(norb):
    csfstr=csfstr+"%4i"%(i+1)
  csfstr=csfstr+" "+csfstr+"\n"
  csf_file.write(csfstr)
  csf_file.write("end\n")

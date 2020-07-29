#!/usr/bin/python

import io
import readline
import sys
import argparse
from numpy import *

def search_array(array,value):
# returns array index for value
    for index,item in enumerate(array):
       if item.strip() == value:
         return index

def atnr2res(nr,atoms):
# search array atoms for atom number nr, returns atom residue name

    atom_res=atoms[int(nr)-1].split()[3]
    return atom_res

def atnr2resNr(nr,atoms):
# search array atoms for atom number nr, returns atom residue number as integer

    atom_resnr=atoms[int(nr)-1].split()[2]
    atom_resnr=int(atom_resnr)
    return atom_resnr

def inter(pair1,pair2,atoms,lastresNr,lastatom,array):

    vW6 = {
    "H":0.,
    "CH0R":6.873E-2,
    "CH1R":6.873E-2,
    "CH2R":6.873E-2,
    "CH0":4.838E-2,
    "CH1":5.396E-2,
    "CH2":6.873E-2,
    "CH3":8.278E-2,
    "HO":0.,
    "OR":4.756E-2,
    "OA":4.756E-2,
    "OE":4.756E-2,
    "COO":.00489, #avery added 7/29/20 needs checked
    "OOC":.00236 #same
    }
    vW12 = {
    "H":0.,
    "CH0R":1.36E-3,
    "CH1R":1.36E-3,
    "CH2R":1.36E-3,
    "CH0":1.837E-3,
    "CH1":1.933E-3,
    "CH2":2.1775964E-3,
    "CH3":2.455782E-3,
    "HO":0.,
    "OR":0.685E-3,
    "OA":1.125E-3,
    "OE":1.125E-3,
    "COO":1.359E-5, #avery added 7/29/20
    "OOC":1.59E-6 #added 7/29/20
    }

    nvW6 = {
    "HO1O5":0.,
    "HO1C3":0.,
    "HO1C5":0.,
    "HO2C4":0.,
    "HO2O5":0.,
    "HO3C1":0.,
    "HO3C5":0.,
    "HO4C2":0.,
    "HO4O5":0.,
    "C6C1":5.689469E-3,
    "C6C3":5.689469E-3,
    "O1C3":3.268799E-3,
    "O1C5":3.268799E-3,
    "O2C4":3.268799E-3,
    "O3C1":3.268799E-3,
    "O3C5":3.268799E-3,
    "O4C2":3.268799E-3,
    "O1C6":4.110135E-3,
    "O3C6":4.110135E-3,
    "O1C2":3.268799E-3,
    "O1C1":3.268799E-3,
    "O1C4":3.268799E-3
    }
    nvW12 = {
    "HO1O5":0.7E-6,
    "HO1C3":0.35E-6,
    "HO1C5":0.35E-6,
    "HO2C4":0.35E-6,
    "HO2O5":0.35E-6,
    "HO3C1":0.35E-6,
    "HO3C5":0.35E-6,
    "HO4C2":0.35E-6,
    "HO4O5":0.35E-6,
    "C6C1":3.33986E-6,
    "C6C3":3.33986E-6,
    "O1C3":0.7E-6,
    "O1C5":0.7E-6,
    "O2C4":2.5E-6,
    "O3C1":2.5E-6,
    "O3C5":2.5E-6,
    "O4C2":2.5E-6,
    "O1C6":8.4108E-6,
    "O3C6":8.4108E-6,
    "O1C2":2.5E-6,
    "O1C1":2.5E-6,
    "O1C4":2.5E-6
    }
    atom1_name = atoms[int(pair1)-1].split()[4]
    atom2_name = atoms[int(pair2)-1].split()[4]
    atom1_type = atoms[int(pair1)-1].split()[1]
    atom2_type = atoms[int(pair2)-1].split()[1]


    atom1atom2 = atom1_name+atom2_name
    atom2atom1 = atom2_name+atom1_name
    at1resNr = int(atnr2resNr(pair1,atoms))
    at2resNr = int(atnr2resNr(pair2,atoms))
    lastresNr = int(lastresNr)
    if( (at1resNr == lastresNr or at2resNr == lastresNr) and lastatom == 'CH3' ):
      nvW12re={
       "O1C5":0.7E-6,
       "O1C3":0.7E-6
      }
    else:
      nvW12re={
       "O1C5":2.5E-6,
       "O1C3":2.5E-6
      }


    link = args.c
    if ( link == '1-2'):
      condition = ( atom1atom2=='O1C4' or atom1atom2=='O1C1' or atom1atom2=='O1C5' ) and at1resNr+1 == at2resNr
    if ( link == '1-3'):
      condition = ( atom1atom2=='O1C1' or atom1atom2=='O1C5' or atom1atom2=='O1C6' ) and at1resNr+1 == at2resNr
    if ( link == '1-4'):
      condition = atom1atom2=='O1C2'  and at1resNr+1 == at2resNr
    if ( link == '1-6'):
      condition = atom1atom2=='O1C4'  and at1resNr+1 == at2resNr
      nvW12['O1C4']=1.53e-06


# 1. O1C2 between rings
    if ( condition ):
       if (link == '1-3'):
         nvW12['O1C1'] = 2.5E-6
         nvW12['O1C6'] = 8.4108E-6
         nvW12['O1C5'] = 2.5E-6
       array.append('%+5s %+5s   1   %4e %4e; %4s %4s NRINGS' % (str(pair1),str(pair2),nvW6[atom1_name+atom2_name],(nvW12[atom1_name+atom2_name]), str(atom1_name), str(atom2_name) ))
    else:

# 2. N interactions in ring - a1a2
      if ( nvW6.has_key(atom1_name+atom2_name) and (at1resNr == at2resNr)  ):
        if (at2resNr == lastresNr and (atom1atom2 == 'O1C3' or atom1atom2 == 'O1C5' )):
          array.append('%+5s %+5s   1   %4e %4e; %4s %4s NRE' % (str(pair1),str(pair2),nvW6[atom1_name+atom2_name],(nvW12re[atom1_name+atom2_name]), str(atom1_name), str(atom2_name) ))
        else:
          array.append('%+5s %+5s   1   %4e %4e; %4s %4s N' % (str(pair1),str(pair2),nvW6[atom1_name+atom2_name],(nvW12[atom1_name+atom2_name]), str(atom1_name), str(atom2_name) ))
      else:
# 3. N interactions in ring - a2a1
       if (nvW6.has_key(atom2_name+atom1_name) and at1resNr == at2resNr ):
        if (at2resNr == lastresNr and (atom2atom1 == 'O1C3' or atom2atom1 == 'O1C5' )):
            array.append('%+5s %+5s   1   %4e %4e; %4s %4s NRE' % (str(pair1),str(pair2),nvW6[atom2_name+atom1_name],(nvW12re[atom2_name+atom1_name]), str(atom2_name), str(atom1_name) ))
        else:
            array.append('%+5s %+5s   1   %4e %4e; %4s %4s N' % (str(pair1),str(pair2),nvW6[atom2_name+atom1_name],(nvW12[atom2_name+atom1_name]), str(atom2_name), str(atom1_name) ))
       else:
# 4. standard interactions in ring
          array.append('%+5s %+5s   1   %4e %4e; %4s %4s' % (str(pair1),str(pair2),vW6[atom1_type]*vW6[atom2_type],vW12[atom1_type]*vW12[atom2_type], str(atom1_name), str(atom2_name)))

    return array

def atnr2name(nr,atoms):
#returns atom name and atom type for nr in atoms array
    atom_name=atoms[int(nr)-1].split()[4]
    atom_type=atoms[int(nr)-1].split()[1]
    return atom_name, atom_type


def array_remove(array,exclude):
    newlist = [item for item in array if not item.strip().startswith(exclude)]
    return newlist


parser = argparse.ArgumentParser()
parser.add_argument('-c', help='Connection: 1-2, 1-3, 1-4, 1-6, either option will work for monomers ',metavar='1-X',required=True)
parser.add_argument('-i', help='Input topology', metavar='input.top', required=True)
parser.add_argument('-o',help='Output topology', metavar='output.top',required=True)
args = parser.parse_args()

print ("Connection: %s" % args.c )
print ("Input file: %s" % args.i )
print ("Output file: %s" % args.o )

input_top = args.i
out_top = args.o
lines = []
f = open(args.i,'r')

for line in f:
    lines.append(line)
f.close

atoms_index = search_array(lines,'[ atoms ]')
bonds_index = search_array(lines,'[ bonds ]')
pairs_index = search_array(lines,'[ pairs ]')
exclus_index = search_array(lines,'[ exclusions ]')
angles_index = search_array(lines,'[ angles ]')


head = lines[:(pairs_index+1)]
atoms0 = lines[(atoms_index+1):(bonds_index-1)]
pairs0 = lines[(pairs_index+1):(exclus_index-1)]
exclus0 = lines[(exclus_index+1):(angles_index-1)]
tail = lines[exclus_index:]

atoms = array_remove(atoms0,";")
pairs = array_remove(pairs0,";")
exclus = array_remove(exclus0,";")

lastresnr = atoms[-1].split()[2]
lastatom = atoms[-1].split()[1]


pairs_new = []
for index,it in enumerate(pairs):
# LJ interactions for [ pairs ] section
    item = it.strip()
    pair1 = item.split()[0]
    pair2 = item.split()[1]
    inter(pair1,pair2,atoms,lastresnr,lastatom,pairs_new)

for index,it in enumerate(exclus):
# LJ interactions for [ exclusions ] section
    item = it.strip()
    if (len(item.split()) == 2):
      pair1 = item.split()[0]
      pair2 = item.split()[1]
      inter(pair1,pair2,atoms,lastresnr,lastatom,pairs_new)
    if (len(item.split()) == 3):
     pair1=item.split()[0]
     pair2=item.split()[1]
     pair3=item.split()[2]
     inter(pair1,pair2,atoms,lastresnr,lastatom,pairs_new)
     inter(pair1,pair3,atoms,lastresnr,lastatom,pairs_new)

#for index,it in enumerate(atoms):
#    item=it.strip()
#    print item

top_new = open(args.o,'w')

for i,item in enumerate(head):
     top_new.write((item))

for item in pairs_new:
     top_new.write((item+"\n"))

top_new.write("\n")


for item in tail:
     top_new.write((item))

#!/usr/bin/env python3

## This script will take the coordinates of alginate monomer units located in
## ./structures.csv's, concatenate them, and convert the polymer to .pdb format


"""
    These are the steps for this script

    1) take in the input that defines the monomer organization
        example: 'MMMGGMMGMGMGMMMM'
        for now, has to be even number of monomers

    2) Call the function that defines the header to the output.pdb

    3) Call the function that takes monomer (M/G) and spits out .pdb for it
        this depends on the location of the monomer

    2) open .pdb file for writing output

    3) write headers of the .pdb file

    4) for each monomer
        open proper .csv file to view
        write xyz coordinates in .csv file to .pdb file
        note that I do not need to write the partial charges
"""
## imports
import sys
import os
import csv
from csv import DictReader

## Grabbing input
polymer = sys.argv[1]

## Defining output
output_dir = polymer + '/'
try:
    os.makedirs(output_dir)
    os.makedirs(output_dir + '00_setup/')
    os.makedirs(output_dir + '01_min/')
    os.makedirs(output_dir + '02_eqm/')
    os.makedirs(output_dir + '03_prod/')
    print("Directory ", output_dir, " created")
except FileExistsError:
    print("Directory ", output_dir, " already exists")

output_name = output_dir + '00_setup/' + polymer + '.pdb' #where the output will be written
output_header = "This is what goes at the top of the .pdb file"
output_remark = "Add any remarks here"

# Size of simulation box, assuming periodic
#box_size = [10, 10, 10]


## Writing the Headers to the pdb output file
with open(output_name, 'w') as output:
# Headers
    output.write("TITLE     " + output_header + "\n")
    output.write("REMARK    " + output_remark + "\n")
    #output.write("CRYST1    " + str(box_size[0]) + "   " + str(box_size[2]) + \
    #"   " + str(box_size[2]) + "   " + "90.00  90.00  90.00 P 1           1\n")
    #if you want to go back to the beginning, f.seek(0)

## initialize atom and res number
atomNum = 0
resNum = 0

## initialize spatial translation of monomer to fit into polymer
new_move = [0,0,0]

## importing and reading structure csv file
path2structures = '/Users/averyagles/google/research/projects/gromacs/alginate/alginator/'
mon_pos = 1 #initializing position of monomer in the polymer
for mon in polymer:
    #beginning of polymer
    if mon_pos == 1:
        if mon == 'G':
            structure = path2structures + 'structures/GNRD.csv'
        if mon == 'M':
            structure = path2structures + 'structures/MNRD.csv'
    #end of polymer
    elif mon_pos == len(polymer):
        if mon == 'G':
            if mon_pos % 2 == 0:
                structure = path2structures + 'structures/GRED_evn.csv'
            else:
                structure = path2structures + 'structures/GRED_odd.csv'
        if mon == 'M':
            if mon_pos % 2 == 0:
                structure = path2structures + 'structures/MRED_evn.csv'
            else:
                structure = path2structures + 'structures/MRED_odd.csv'
    #middle, even
    elif mon_pos % 2 == 0:
        if mon == 'G':
            structure = path2structures + 'structures/GEVN.csv'
        if mon == 'M':
            structure = path2structures + 'structures/MEVN.csv'
    #middle, odd
    else:
        if mon == 'G':
            structure = path2structures + 'structures/GODD.csv'
        if mon == 'M':
            structure = path2structures + 'structures/MODD.csv'

    ## Grabbing info in structure and appending to .csv
    with open(structure, 'r') as f:
        dict_reader = DictReader(f)
        resNum += 1

        # loopin through lines of data in csv
        move = new_move #resetting the translation for this monomer
        for row in dict_reader:
            atomName = row['name']
            atomNum += 1 #index for atom number
            resName = row['res']
            elem = row['element']
            x = float('%.3f'%(float(row['x']) + move[0])) #have to convert to number before addition
            y = float('%.3f'%(float(row['y']) + move[1]))
            z = float('%.3f'%(float(row['z']) + move[2]))
            # This step redefines the move for the next monomer
            if atomName == 'O1':
                new_move = [x, y, z]
            occupancy = '1.00'
            tempFact = '0.00'
            none = "" #for columns in pdb that will be empty

            atomName = atomName #left justified, 4 columns
            # atomNum
            altLoc = none.rjust(0)
            resName = resName.rjust(4) #right adjusted, 4 columns (ASSUMING 4 characters)
            chainID = none.rjust(0)
            # resNum
            code = none.rjust(1)
            x = str(x).rjust(8) #convert back to string
            y = str(y).rjust(8)
            z = str(z).rjust(8)
            occupancy = occupancy.rjust(6)
            tempFact = tempFact.rjust(6)
            segID = none.ljust(4)
            elem = elem.rjust(2)


            with open(output_name, 'a') as output: #a for append!
                output.write("ATOM  " + str(atomNum).rjust(5) + none.rjust(2) + atomName.ljust(4) + altLoc + resName + " " + chainID + str(resNum).rjust(4) + code + none.rjust(3) + x + y + z + occupancy + tempFact + none.rjust(6) + segID + elem + '\n')
    mon_pos += 1

with open(output_name, 'a') as output:
    output.write("TER\nENDMDL")

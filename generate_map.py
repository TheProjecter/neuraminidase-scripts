"""
Generate map for mutations in neuraminidase (H1N1) sequences.
@author: Carlos Olivares
@email: alfaceor@gmail.com
"""

import domain_data
import sys
import os
import Mutation
# lee el archivo de input solo el basename
alignBasename = sys.argv[1].split(".")[0].split("/")[-1]
if alignBasename=="":
    quit()

#########################
# Read the aligment files 
#########################
from Bio import SeqIO

#handle = open(alignBasename+".aln","rU")
handle = open(sys.argv[1],"rU")
alignRecords = list(SeqIO.parse(handle,"clustal"))
handle.close()
print "clustalw time"
#print str(tocks-ticks)
print "sequence reference"
print alignRecords[0].id
# print alignRecords[0].seq
#########################
# Use alligment to obtain the mutations
#########################
# wildtype Mexico: ACQ73395
wildtype=alignRecords[0]
# print '-'*10
# print wildtype==records[3]

#open files for the output
outputfile_1d_structure=open(alignBasename+"mapa_mutaciones-1d.csv",'w')
outputfile_3d_structure=open(alignBasename+"mapa_mutaciones-3d.csv",'w')

outputfile_2d_structure=open(alignBasename+"mapa_mutaciones-2d.csv",'w')
outputfile_as_structure=open(alignBasename+"mapa_mutaciones-as.csv",'w')

csv_token=","
# Print header for csv files.
head="description"+csv_token+"id"
for i in range(len(wildtype.seq)):
    head=head+csv_token+str(i+1)
outputfile_1d_structure.write(head+'\n')
outputfile_2d_structure.write(head+'\n')
for i in range(len(wildtype.seq)):
    head=head+csv_token+str(i+1)

outputfile_3d_structure.write(head+'\n')


# region list for the 29 region (active site and others).
numRegion_AS=29
mutRegion_AS=[[] for i in range(numRegion_AS)]
numRegion_2D=61
mutRegion_2D=[[] for i in range(numRegion_2D)]

# Regions for Secondary Structure
regions_SS=[81,95,103,104,109,113,118,120,125,127,136,137,140,155,161,172,178,179,186,188,198,200,209,210,218,231,234,236,244,251,259,261,269,279,284,286,292,300,307,310,317,350,353,355,361,367,377,387,393,397,400,401,408,409,413,418,430,436,448,469]

# Regions for Active Site
regions_AS=[82,117,118,119,150,151,152,155,156,176,177,222,223,224,225,227,228,276,277,278,292,293,294,295,367,368,401,402]


# print csv file with the sequences
for i in range(len(alignRecords)):
    mut=""
    mut3d=""
    for j in range(len(wildtype.seq)):
        # if a mutation appears then show that
        if wildtype.seq[j] != alignRecords[i].seq[j]:
			# record the Accession Number, the mutation and country
            j3d=domain_data.notation_3D(j)
            notation3d = str(wildtype.seq[j])+j3d+str(alignRecords[i].seq[j]) # notation for 3D structure.
            # check what ist the region in the sequence by active site.
            if not('-' in notation3d):
                # to notation issues j+1 to print.
                mut = mut+str(wildtype.seq[j])+str(j+1)+str(alignRecords[i].seq[j])
                # Active Site
                mutRegion_AS[domain_data.whatRegion_AS(j+1)].append(notation3d)
                # Secondary structure
                mutRegion_2D[domain_data.whatRegion_2D(j+1)].append(notation3d)
                mut3d = mut3d+notation3d

        mut = mut+csv_token
        mut3d = mut3d+csv_token
    line=str(alignRecords[i].description)+csv_token+str(alignRecords[i].id)+csv_token+mut
    line3d=str(alignRecords[i].description)+csv_token+str(alignRecords[i].id)+csv_token+mut3d
    outputfile_1d_structure.write(line+'\n')
    outputfile_3d_structure.write(line3d+'\n')

# For the active site map.
for i in range(numRegion_2D):
    # outputfile_2d
    outputfile_2d_structure.write(str(i)+","+str(set(mutRegion_2D[i]))+'\n')

for i in range(numRegion_AS):
    # outputfile_as
    outputfile_as_structure.write(str(i)+","+str(set(mutRegion_AS[i]))+'\n')

outputfile_1d_structure.close()
outputfile_2d_structure.close()
outputfile_3d_structure.close()
outputfile_as_structure.close()


# TODO:
"""
# imprime el diccionario
# ordenar mutaciones por posicion
# N3Y, R10D, T70U
# genera el mapa
 # leer las posiciones claves (sitio activo)
 # recorrer diccionario
 # si es una posicion clave anadir a lista

"""

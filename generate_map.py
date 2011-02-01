"""
Generate map for mutations in neuraminidase (H1N1) sequences.
@author: Carlos Olivares
@email: alfaceor@gmail.com
"""



import domain_data
import sys
import os
# lee el archivo de input solo el basename
alignBasename = sys.argv[1].split(".")[0].split("/")[-1]
if alignBasename=="":
    quit()
"""
#########################
# Run a clustalw
#########################
# FIXME: Why cannot use the package Bio.Align.Applications command not run
#from Bio.Align.Applications import ClustalwCommandline
#commandString = ClustalwCommandline(infile=alignBasename+".fasta")
commandString = "mpirun -n 7 /home/alfaceor/instaladores/clustalw-mpi-0.13/clustalw-mpi -OUTORDER=INPUT -infile="+alignBasename+".fasta"
# FIXME: This is a temporal solution for the fail with Bio.Align.Applications package.
# Reference: http://www.ibm.com/developerworks/aix/library/au-python/
import commands
import time
try:
    #run a 'clustalw' command and assign results to a variable
    #commandString = cline
    ticks = time.time()
    commandOutput = commands.getoutput(commandString)
    tocks = time.time()
#    print commandOutput
    
except:
    print "Clustalw can not run, maybe you must install first"
"""

#########################
# Read the aligment files 
#########################
from Bio import SeqIO

handle = open(alignBasename+".aln","rU")
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

# region list for the 29 region (active site and others).
numRegion_AS=29
mutRegion_AS=[[] for i in range(numRegion_AS)]
numRegion_2D=61
mutRegion_2D=[[] for i in range(numRegion_2D)]


outputfile_3d_structure.write(head+'\n')

# print csv file with the sequences
for i in range(len(alignRecords)):
    mut=""
    mut3d=""
    for j in range(len(wildtype.seq)):
        # if a mutation appears then show that
        if wildtype.seq[j] != alignRecords[i].seq[j]:
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

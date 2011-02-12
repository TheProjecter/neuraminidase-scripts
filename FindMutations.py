#import domain_data
import sys
import os
import Mutation

# Aligment sequences in clustal format
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
print "sequence reference"
print alignRecords[0].id
# print alignRecords[0].seq

#########################
# Use alligment to obtain the mutations
#########################
# wildtype Mexico: ACQ73395
wildtype=alignRecords[0]


#########################
# Local Database Neuraminidase
#########################
from BioSQL import BioSeqDatabase
server = BioSeqDatabase.open_database(driver="MySQLdb", user="root", passwd = "rootpass", host = "localhost", db="bioseqdb")
db = server["neuraminidase"]

#########################
# Google country name
#########################
import string
from geopy import geocoders
""" Get the country name from the genbank description field for the sequence

args:
     description: record.description field in genbank format
return:
        countryLocation: string name for the country in a uniform format, based on google API.
"""

def getCountry_fromGBDescription(description)
   # separete the fields by '/' to get the ubication
   recordGeoLocation=string.split(description,'/')[1] # A/Mexico City/s
   # localhost key for the google api
   g=geocoders.Google('ABQIAAAAjB-oGyPYZ_fAg6eMOY4uoxT2yXp_ZAY8_ufC3CFXhHIE1NvwkxRM-xD84KmxCLors0B9Ii1PetMJwQ')
   # place in google api format, latitude and longitud
   place, (lat,lng)=g.geocode(recordGeoLocation)
   print place
   # last field is the country so splited and then get the last one
   splitedLocation=string.split(place,',')
   print splitedLocation
   # get the last one
   countryLocation=splitedLocation[len(splitedLocation)-1]
   return countryLocation

#########################
# getMutations
#########################
from time import sleep
from Bio import Entrez

mutations_list=[]
print len(alignRecords)
for i in range(len(alignRecords)):
   mut=""
   mut3d=""
   # TODO: make a module for this process
   # Insert sequence in database to actualice database.
   if False:
      try:
         handle = Entrez.efetch(db="protein", id=alignRecords[i].id, rettype="gb")
         db.load(SeqIO.parse(handle, "genbank"))
         server.commit() #On Biopython 1.49 or older, server.adaptor.commit()
         sleep(3) # wait for 3 seconds
         print ''+str(i)+'The accession number '+accessionNumber+'has been inserted'
      except Exception:
         print 'Error can\'t insert in database - '+str(i)+' - '+str(alignRecords[i].id)
   
   # Search in the sequence for mutations
   for j in range(len(wildtype.seq)):
      # if a mutation appears then show that
      seq_record=db.lookup(accession=alignRecords[i].id)
      alignRecords[i].description=seq_record.description
      
      if (wildtype.seq[j] != alignRecords[i].seq[j]) and (alignRecords[i].seq[j]!='-'):
      # record the Accession Number, the mutation and country
         mutations_list.append(Mutation.Mutation(wildtype.id,alignRecords[i].id,wildtype.seq[j],j,alignRecords[i].seq[j],getCountry_fromGBDescription(alignRecords[i].description)))
         
         
         
         





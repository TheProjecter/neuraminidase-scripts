import domain_data
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
server = BioSeqDatabase.open_database(driver="MySQLdb", user="root", passwd = "123456", host = "localhost", db="bioseqdb")
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

def getCountry_fromGBDescription(description):
   # separete the fields by '/' to get the ubication
   recordGeoLocation=string.split(description,'/')[1] # A/Mexico City/s
   # localhost key for the google api
   g=geocoders.Google('ABQIAAAAjB-oGyPYZ_fAg6eMOY4uoxT2yXp_ZAY8_ufC3CFXhHIE1NvwkxRM-xD84KmxCLors0B9Ii1PetMJwQ')
   # place in google api format, latitude and longitud
   place, (lat,lng)=g.geocode(recordGeoLocation)
   #print place
   # last field is the country so splited and then get the last one
   splitedLocation=string.split(place,',')
   #print splitedLocation
   # get the last one
   countryLocation=splitedLocation[len(splitedLocation)-1]
   return countryLocation.strip()

#########################
# Regions SS and AS
#########################

# Regions for Secondary Structure
regions_SS=[81,95,103,104,109,113,118,120,125,127,136,137,140,155,161,172,178,179,186,188,198,200,209,210,218,231,234,236,244,251,259,261,269,279,284,286,292,300,307,310,317,350,353,355,361,367,377,387,393,397,400,401,408,409,413,418,430,436,448,469]
# Regions for Active Site
regions_AS=[82,117,118,119,150,151,152,155,156,176,177,222,223,224,225,227,228,276,277,278,292,293,294,295,367,368,401,402]


#########################
# getMutations
#########################
from time import sleep
from Bio import Entrez
import domain_data
import csv
dataWriter = csv.writer(open('mutationData.csv', 'wb'), delimiter=',')
dataWriter.writerow(["Accession Number","old_aa","position3D","new_aa","Country","Secondary Structure","Active Site"])

mutations_list=[]
print len(alignRecords)
for i in range(len(alignRecords)):
#for i in range(3):
   mut=""
   mut3d=""
   # TODO: make a module for this process
   # Insert sequence in database to actualice database.
   if False:
      try:
         handle = Entrez.efetch(db="protein", id=alignRecords[i].id, rettype="gb")
         db.load(SeqIO.parse(handle, "genbank"))
         server.commit() #On Biopython 1.49 or older, server.adaptor.commit()
         sleep(1) # wait for 1 seconds
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
         country=getCountry_fromGBDescription(alignRecords[i].description)
         mutations_list.append(Mutation.Mutation(wildtype.id,alignRecords[i].id,wildtype.seq[j],j,alignRecords[i].seq[j],getCountry_fromGBDescription(alignRecords[i].description)))
         # write data in csv file in 3D NOTATION
         mut_act_pos=len(mutations_list)-1
         # Define regions for SS and AS
         ss_reg=domain_data.whatRegion(mutations_list[mut_act_pos].position_1D,regions_SS)
         as_reg=domain_data.whatRegion(mutations_list[mut_act_pos].position_1D,regions_AS)
         csvrow=[mutations_list[mut_act_pos].sequence_AN,mutations_list[mut_act_pos].old_aa,domain_data.notation_3D(mutations_list[mut_act_pos].position_1D),mutations_list[mut_act_pos].new_aa,str(mutations_list[mut_act_pos].country),ss_reg,as_reg]
         dataWriter.writerow(csvrow)
         
# put country_list for accession number mutation
country_list={}
country_csv_files={}
for m in mutations_list:
   # check if country in list
   if m.country not in country_list:
      country_list.update({m.country:1})
      # Create csv for file
      country_csv_files.update({m.country:csv.writer(open('mutation_'+m.country+'.csv', 'wb'), delimiter=',')})
      ss_reg=domain_data.whatRegion(m.position_1D,regions_SS)
      as_reg=domain_data.whatRegion(m.position_1D,regions_AS)
      country_csv_files[m.country].writerow(["Accession Number","mutation3D","ss_region","as_region"])
      csvrow=[m.sequence_AN,m.old_aa+domain_data.notation_3D(m.position_1D)+m.new_aa,ss_reg,as_reg]
      country_csv_files[m.country].writerow(csvrow)
   else:
      country_list[m.country]=country_list[m.country]+1
      ss_reg=domain_data.whatRegion(m.position_1D,regions_SS)
      as_reg=domain_data.whatRegion(m.position_1D,regions_AS)
      csvrow=[m.sequence_AN,m.old_aa+domain_data.notation_3D(m.position_1D)+m.new_aa,ss_reg,as_reg]
      country_csv_files[m.country].writerow(csvrow)
  

# Open csv with the 

print country_list

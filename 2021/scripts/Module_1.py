#!/usr/bin/env python

##1st Module of the Interpro wrapper. Main idea is to be able to quick get an owerview inside all possible domain rearrangiaments of a 
##protein family or a proteoms (!To Be Tested!).


import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='First module of Interproscan output wrapper.')
#TSV input file
parser.add_argument('--in_tsv', required=True, help='The input tsv file.')
#Prefix for output file
parser.add_argument('--out', required=True, help='prefix for output files.')
#Enable Multi-Specie mode, the workflow is the same but storing also Species informations. Useful in comparative genomics
parser.add_argument('--Multiple_Species', action='store_true')
#Domain db to use
parser.add_argument('--domain_db',choices=['Pfam', 'CDD', 'SMART', 'ProSiteProfiles'],help='Domain db to use', required=True, type=str)


args = parser.parse_args()

#Setting up arguments to pass to the "table reformatting" function.

DOMAIN_DB = args.domain_db
IN_FILE = args.in_tsv
OUT = args.out
if args.Multiple_Species :
    MODE = True
else :
    MODE = False

#Setting up output directories
path = os.getcwd()
TABLE_OUT = "Tables"
SEQUENCES_OUT = "Sequences"
TABLE_path = "/".join((path,TABLE_OUT))
SEQUENCES_path = "/".join((path,SEQUENCES_OUT))
os.mkdir(TABLE_path)
os.mkdir(SEQUENCES_path)

#--------------------MAIN----------------------------------#

#Function for correcting parsing of Interproscan tsv file:
def table_reformatting(IN_FILE, MODE, DOMAIN_DB) :
    
#Creation of list were interested values will be stored
#Query ID
    query = []
#Pfam ID
    Domain_Hit = []
#Starting coordinates
    Start = []
#Ending coordinates
    End = []
#If Multi-specie option enable, we need one more list to store species informations    
    if MODE :
        specie = []
#Open tsv file and keep only interested columns and rows
    with open(IN_FILE, 'r') as file:
        l = file.readline()
        for l in file :
#Split file by tab and keep rows in which the 3rd field match the "Pfam" pattern
            sl = l.split('\t')
            if sl[3] == DOMAIN_DB :
#Append all interested values (Query ID, Pfam ID, Starting and Ending coordinates):
                query.append(str(sl[0]))
                Domain_Hit.append(str(sl[4]))
                Start.append(int(sl[6]))             
                End.append(int(sl[7]))
                if MODE :
                    sl_2 = l.split("|")
                    specie.append(str(sl_2[0]))
    if MODE :
        #Creation of a sorted dataframe with species Name ( == multi-specie mode)
        df = pd.DataFrame(list(zip(specie, query, Domain_Hit,Start,End)),columns = ['Specie','Query', 'Domain','Start','End'])
    else :
        #Creation of a sorted dataframe without any species name ( == Single Specie mode)
        df = pd.DataFrame(list(zip(query, Domain_Hit,Start,End)),columns = ['Query', 'Domain','Start','End']) #Creation of a sorted dataframe
    return df

#-----------------Results summary and printing---------------#

df = table_reformatting(IN_FILE, MODE,DOMAIN_DB)
df = df.sort_values(['Query','Start'])

#Summary table of number of domains
item_counts = pd.DataFrame(df["Domain"].value_counts())
item_counts.to_csv(path_or_buf=open('Tables/%s_Domain.Summary.tsv' % OUT, 'w'),sep="\t",header = False)
print(" Summary number of total %s domains :" %DOMAIN_DB, item_counts, sep="\n")
print("-----> Total sum of", item_counts.sum())
print("-----> Accounting for :", len(item_counts.index), "Unique domains")

#Summary table with occurence of each domain in each specie if MODE = Multi-Specie
if MODE :
    Domain_specie = df.groupby(['Specie','Domain']).size().reset_index(name='count')
    Domain_specie.to_csv(path_or_buf=open('Tables/%s_Domain.Occurrence.tsv' % OUT, 'w'),sep="\t",header = True, index = False)
    print("Printing Species Domain occurrence to %s_Domain.Occurrence.tsv" % OUT)

print("\n")
print("Finiding unique domains compositions...")
print("\n")

#Summary table with domain arrangeament for each Query
sep = '-'
if MODE :
    Domain_arr = pd.DataFrame(df.groupby(['Specie','Query'], sort=False)['Domain'].apply(lambda x: sep.join(x))).reset_index()
else :
    Domain_arr = pd.DataFrame(df.groupby('Query', sort=False)['Domain'].apply(lambda x: sep.join(x))).reset_index()

Domain_arr.to_csv(path_or_buf=open('Tables/%s_Domain.Arrangiament.tsv' % OUT, 'w'),sep="\t", header = True, index = False)
print("Printing Domain arrangeaments to %s_Domain.Arrangiament.tsv" % OUT)
Domain_arr_Count = Domain_arr['Domain'].value_counts()
Domain_arr_Count.to_csv(path_or_buf=open('Tables/%s_Domain.Arrangiament_Count.tsv' % OUT, 'w'),sep="\t", header = False)
print("Brief Summary :", pd.DataFrame(Domain_arr['Domain'].value_counts()))
print("----> Founded :", len(Domain_arr.Domain.unique()),"unique domains arrangeaments accounting for :", len(Domain_arr.index), "total queries.")
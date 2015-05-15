#!/usr/bin/python
import _mysql as mysql
import _mysql_exceptions as mysql_exceptions
import MySQLdb.cursors
import os
import pickle
import numpy as np
import sys
import StringIO
from scipy.misc import comb
import getopt
import csv
import re

def usage():
    print "\n-----------------------------------------------------------------"
    print "Usage: "
    print "   sifter_find_families.py [options] <output_file>"
    print "-----------------------------------------------------------------\n"
    print "Examples:"
    print "   sifter_find_families.py -p C0JYY2_HUMAN ../examples/family_list.txt\n"
    print "   sifter_find_families.py -s 9823 ../examples/family_list.txt\n"
    print "   sifter_find_families.py -A ../examples/family_list.txt\n"
    print "   sifter_find_families.py --ip ../example/protein_list.txt ./examples/family_list.txt\n"
    print "   sifter_find_families.py -p C0JYY2_HUMAN --dbaddr www.example.org --dbuser jack --dbpass 1234 ./examples/family_list.txt\n"    
    print "This function reports list of Pfam families for your query protein or species."
    print "@author Sayed Mohammad Ebrahim Sahraeian (mohammad@compbio.berkeley.edu)"
    print "Please cite new paper:"
    print "-Sahraeian SME, Luo KR, Brenner SE (2015)"
    print "\nThe SIFTER algorithm presented in the following paper:"
    print "- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. Genome-scale phylogenetic function annotation of large and diverse protein families. Genome Research 21:1969-1980. \n"
    print "inputs:"
    print "        <output_file>            the output file where the list of"
    print "                                 of Pfam families will be placed"
    print "options: (you should only use one of '-p -s --ip -A' options.)"
    print "           -p          STRING    The query protein (use Uniprot ID"
    print "                                 or Accession)."
    print "           -s          STRING    NCBI taxonomy ID for input species."
    print "           --ip        STRING    Path to the input file where the list"
    print "                                 of proteins are placed."
    print "           -A          STRING    Use to find all pfam families."
    print "           --dbaddr    STRING    Address of the MySQL database that"
    print "                                 has neccessary data for SIFTER"
    print "                                 [Default: localhost]"
    print "           --dbname    STRING    Name of the MySQL database that"
    print "                                 has neccessary data for SIFTER"
    print "                                 [Default: sifter_db]"
    print "           --dbuser    STRING    Name of the MySQL database that"
    print "                                 has neccessary data for SIFTER"
    print "                                 [Default: root]"
    print "           --dbpass    STRING    Password of the user for the MySQL"
    print "                                 database that has neccessary data"
    print "                                 for SIFTER [Default: '']"
    print "           -h                    Help. Print Usage."


def msql(query, db):
    c = db.cursor()
    c.execute(query)
    results = c.fetchall()
    c.close()
    return results


def find_pfam_for_genes(genes):
    sql="""SELECT
      pfamseq.pfamseq_acc,
      pfamseq.pfamseq_id,
      pfamA.pfamA_acc
     FROM   pfamseq
      INNER JOIN pfamA_reg_full_significant on (pfamA_reg_full_significant.auto_pfamseq=pfamseq.auto_pfamseq)
      INNER JOIN pfamA on (pfamA.auto_pfamA=pfamA_reg_full_significant.auto_pfamA)
     WHERE
       pfamseq.pfamseq_acc in ('%s')
       OR pfamseq.pfamseq_id in ('%s')
    """%("','".join(genes),"','".join(genes))
    seq_anns = msql(sql, db_mysql)
    return seq_anns

def find_pfam_for_taxid(taxid):
    sql="""SELECT
      pfamseq.pfamseq_acc,
      pfamseq.pfamseq_id,
      pfamA.pfamA_acc
     FROM   pfamseq
      INNER JOIN pfamA_reg_full_significant on (pfamA_reg_full_significant.auto_pfamseq=pfamseq.auto_pfamseq)
      INNER JOIN pfamA on (pfamA.auto_pfamA=pfamA_reg_full_significant.auto_pfamA)
     WHERE
       pfamseq.ncbi_taxid = '%s'
    """%(taxid)
    seq_anns = msql(sql, db_mysql)
    return seq_anns

def find_all_pfams():
    sql="""SELECT
      DISTINCT(pfamA_acc)
     FROM   pfamA
    """
    seq_anns = msql(sql, db_mysql)
    res=list(set([w['pfamA_acc']for w in seq_anns]))
    return res    
    
def find_pfams(res):
    pfams=[w['pfamA_acc'] for w in res]
    return pfams
    

# See how many sequences are in each family, add it to pfds
def get_pfds(pfams):
    pfds={}
    for p in pfams:
        sql_q="""select pfamA.num_full, pfamA.number_species, pfamA.pfamA_id, pfamA.description, group_concat(go_id) as go_ids, group_concat(term) as go_terms from pfamA left join gene_ontology on gene_ontology.auto_pfamA = pfamA.auto_pfamA where pfamA.pfamA_acc='%s' group by pfamA_acc
        """%(p)
        #AND Locus.type=1
        #AND Synonym.type=2
        r = msql(sql_q, db_mysql)
        if r:
            pfds[p]={}
            for w in r[0].keys():
                pfds[p][w]=r[0][w]
    return pfds

if __name__=="__main__":
    
    # Initialization
    params_mysql = {\
    'db_address': 'localhost',
    'db_username': 'root',
    'db_password': '',
    'db_name': 'sifter_db'
    }

    taxid=''
    query_proteins=[]
    protein=''
    input_file=''
    all_fams=0
    # Check for options
    
    opts, args = getopt.getopt(sys.argv[1:], "hAp:s:",['ip=','dbname=','dbpass=','dbuser=','dbaddr=']) 
    if len(args) != 1:
        usage()
        sys.exit()
    choices=[]
    if len(opts)>0:
        for o, a in opts:
            if o == "-p":
                protein=a
                choices.append('p')
            elif o == "-s":
                taxid = a
                choices.append('s')                
            elif o == "--ip":
                input_file = a
                choices.append('ip')                
            elif o == "-A":
                all_fams = 1
                choices.append('A')                
            elif o == "--dbname":
                params_mysql['db_name']= a
            elif o == "--dbaddr":
                params_mysql['db_address']= a
            elif o == "--dbpass":
                params_mysql['db_password']= a
            elif o == "--dbuser":
                params_mysql['db_username']= a
            else:
                usage()
                sys.exit()

    if len(choices)==0:
        print "\nERROR: No queries are entered."
        print "Please use one of the '-p -s --ip -A' options to enter your query.\n"
        sys.exit()
    elif len(choices)>1:
        print "\nERROR: Please use ONLY one of the '-p -s --ip -A' options to enter your query.\n"
        sys.exit()

   
    output_file=args[0]
    
    ###
    db_mysql = MySQLdb.connect(host=params_mysql['db_address'],
                   user=params_mysql['db_username'],
                   passwd=params_mysql['db_password'],
                   db=params_mysql['db_name'],
                   cursorclass=MySQLdb.cursors.DictCursor)
    print "\n\n--------------Searching for Pfams families------------"
    if protein:
        query_proteins=[protein]
        res=find_pfam_for_genes(query_proteins)
        pfams=find_pfams(res)
        print "Found %s Pfam families for query protein %s"%(len(pfams),protein)
    elif input_file:
        if not os.path.exists(input_file):
            print "\nERROR: No file exists at %s\n"%input_file
            sys.exit()
        f = open(input_file, 'r')
        a=f.read()
        splited =re.split(' |,|;|\n',a.strip())
        query_proteins=list(set([w for w in splited if w]))
        res=find_pfam_for_genes(query_proteins)
        pfams=find_pfams(res)
        print "Found %s Pfam families for %s query proteins"%(len(pfams),len(query_proteins))
    elif taxid:
        res=find_pfam_for_taxid(taxid)
        pfams=find_pfams(res)
        print "Found %s Pfam families for query species (taxid=%s)"%(len(pfams),taxid)
    elif all_fams==1:
        pfams=find_all_pfams()
        print "Found %s Pfam families"%(len(pfams))

    if not pfams:
        print "There are no pfam families for your input query."

    pfds=get_pfds(pfams)

    tree_sizes = {}
    for p in pfds.keys():
        tree_sizes[p] = pfds[p]['num_full']
    sorted_fams = sorted(pfds.keys(), key=lambda k:pfds[k]['num_full'])
    
    print "Number of families:" ,len(sorted_fams)
    for p in sorted_fams:
        print p, "(Family size: %s)"%tree_sizes[p]
    print 

    with open(output_file, 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\n',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(sorted_fams)

    print "\nResults are written to %s"%output_file
    print "\n----------------Finding Families is Done------------------"    
        


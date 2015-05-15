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
    print "   sifter_extract.py [options] <results_folder> <output_file>"
    print "-----------------------------------------------------------------\n"
    print "Examples:"
    print "   sifter_extract.py -p C0JYY2_HUMAN ../example/results ../examples/preds.txt\n"
    print "   sifter_extract.py -s 9823 ../example/results ../examples/preds.txt\n"
    print "   sifter_extract.py -f PF12491 ../example/results ../examples/preds.txt\n"
    print "   sifter_extract.py --ip input_file -r ../example/results ../examples/preds.txt\n"
    print "   sifter_extract.py -s 9823 --dbaddr www.example.org --dbuser jack --dbpass 1234 ../example/results ../examples/preds.txt\n"    
    print "   sifter_extract.py -A --hit_file ../example/pfam_res.txt ../example/results ../examples/preds.txt\n"    
    print "This function extracts SIFTER predictions for your query protein/species/family using the outputs of 'sifter_run.py'"
    print "@author Sayed Mohammad Ebrahim Sahraeian (mohammad@compbio.berkeley.edu)"
    print "Please cite new paper:"
    print "-Sahraeian SME, Luo KR, Brenner SE (2015)"
    print "\nThe SIFTER algorithm presented in the following paper:"
    print "- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. Genome-scale phylogenetic function annotation of large and diverse protein families. Genome Research 21:1969-1980. \n"
    print "inputs:"
    print "        <results_folder>         Path to the results folder where"
    print "                                 results have been written to." 
    print "                                 Use the same directory used as" 
    print "                                 results_folder in 'sifter_run.py' script."
    print "        <output_file>            Output file where the extracted"
    print "                                 SIFTER predictions will be written to"
    print "options: (you should only use one of '-p -s -f --ip -A' options.)"
    print "           -p          STRING    List of query proteins (use Uniprot ID"
    print "                                 or Accession) in comma seperated format."
    print "           -s          STRING    NCBI taxonomy ID for input species."
    print "           -f          STRING    List of Pfam families for which you"
    print "                                 want to prepare data."
    print "                                 (in comma seperated format)"
    print "           --ip        STRING    Path to the input file where the list"
    print "                                 of proteins are placed."
    print "           --hit_file  STRING    Output of pfam_scan.pl file on the "
    print "                                 novel genome. This file consists of"
    print "                                 the list of pfam hits for the genome."
    print "                                 If this option is uded, we will"
    print "                                 look in this file to find Pfams"
    print "                                 instead of the SQL database."
    print "           -A                    Prepare for all Pfam families of queried"
    print "                                 novel genome. (hit_file should be provided)"
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

def find_genes_for_pfams(pfam_ids):
    sql="""SELECT
      pfamseq.pfamseq_id,
      pfamA.pfamA_acc
     FROM   pfamseq
      INNER JOIN pfamA_reg_full_significant on (pfamA_reg_full_significant.auto_pfamseq=pfamseq.auto_pfamseq)
      INNER JOIN pfamA on (pfamA.auto_pfamA=pfamA_reg_full_significant.auto_pfamA)
     WHERE
       pfamA.pfamA_acc in ('%s')
    """%("','".join(pfam_ids))
    seq_anns = msql(sql, db_mysql)
    return seq_anns


def find_pfam_2_gene(res):
    pfam_2_gene={}
    gene_2_pfam={}
    for w in res:
        my_pfam=w['pfamA_acc']
        my_gene=w['pfamseq_id']
        if my_pfam not in pfam_2_gene:
           pfam_2_gene[my_pfam]=set([])
        pfam_2_gene[my_pfam].add(my_gene)
        if my_gene not in gene_2_pfam:
           gene_2_pfam[my_gene]=set([])
        gene_2_pfam[my_gene].add(my_pfam)
    return pfam_2_gene,gene_2_pfam



def find_pfam_2_gene_from_file(hit_file):
    pfam_2_gene={}
    gene_2_pfam={}
    with open(hit_file, 'rb') as infile:
        for line in infile:
            line=line.strip()
            if not line:
                continue
            if len(line)<3:
                continue
            if line[0]=="#" and not line[2]=="<":
                continue
            if line[0]=="#" and line[2]=="<":
                keys=line.split('> <')
                keys[0]=keys[0].split('<')[1]
                keys[-1]=keys[-1].split('>')[0]
                continue
            row=line.split()
            if not len(row)==15:
                print "ERR"
                break
            r={k:row[i] for i,k in enumerate(keys)}
            if r['significance']=='1':
                pfam_id=r['hmm acc'][0:r['hmm acc'].find('.')]                
                my_gene=r['seq id']
                if pfam_id not in pfam_2_gene.keys():
                    pfam_2_gene[pfam_id]=set([])
                pfam_2_gene[pfam_id].add(my_gene)
                if my_gene not in gene_2_pfam.keys():
                    gene_2_pfam[my_gene]=set([])
                gene_2_pfam[my_gene].add(pfam_id)

    print "Your queried novel genome has:"
    print len(gene_2_pfam), "genes in pfam"
    print len(pfam_2_gene), "pfam families\n"
    return  pfam_2_gene,gene_2_pfam

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



# ##Process Evidence
def parse_GO_OBO(obo_file):
    go_dict={}
    new_is_comming=-1
    with open(obo_file, "r") as infile:
        currentGOTerm = None
        for line in infile:
            line = line.strip()
            if not line: continue #Skip empty
            if new_is_comming==1:
                key, sep, val = line.partition(":")  
                key=key.strip()
                val=val.strip()
                currentGOTerm=val
                go_dict[currentGOTerm]={}
                new_is_comming=0
                continue
            if line == "[Term]":
                new_is_comming=1
            elif line == "[Typedef]":
                #Skip [Typedef sections]
                new_is_comming=-1
            elif  new_is_comming==0:
                #Only process if we're inside a [Term] environment
                key, sep, val = line.partition(":")
                key=key.strip()
                val=val.strip()            
                if key not in go_dict[currentGOTerm]:
                    go_dict[currentGOTerm][key]=[]
                go_dict[currentGOTerm][key].append(val.strip())
        #Add last term


    #remove obsoletes
    obseletes=[]
    for term in go_dict:
        if 'is_obsolete' in go_dict[term]:
            if go_dict[term]['is_obsolete'][0]== 'true':
                obseletes.append(term)
                continue
    for term in obseletes:
        del go_dict[term]


    ontologies=['biological_process','molecular_function','cellular_component']
    DAGs={w:{} for w in ontologies}
    DAGs_r={w:{} for w in ontologies}
    roots={w:{} for w in ontologies}
    for term in go_dict.keys():
        ont=go_dict[term]['namespace'][0]
        DAGs[ont][term]=[]
        DAGs_r[ont][term]=[]
    for term in go_dict.keys():
        ont=go_dict[term]['namespace'][0]
        if 'is_a' in go_dict[term]:
            for pa in go_dict[term]['is_a']:
                term_2=pa.split(' ! ')[0]
                DAGs[ont][term].append(term_2)
                DAGs_r[ont][term_2].append(term)
        else:
           roots[ont]=term 
    
    
    return go_dict,DAGs,DAGs_r,roots

def trace_to_ontology_root(cur_node):
    """
    Generator to recursively visit all nodes on each path
    from a node up to the root node.
    """
    #print "Graph node:", cur_node
    yield cur_node
    for pa in DAGs[ont][cur_node]:
        for n in trace_to_ontology_root(pa):
            yield n

def get_ontology_subdag(annotated_term_nodes):
    """
    Given evidence_set, returns a filtered subgraph of evidence_ontology
    that only contains those nodes or their ancestors.
    """
    # For each annotated node, traverse to the root node of the ontology
    # to include all its less-specific terms
    all_term_nodes = set([])
    for go_term in annotated_term_nodes:
        traced=trace_to_ontology_root(go_term)
        all_term_nodes.update(set(traced))
    sub_dag = all_term_nodes
    return sub_dag

def get_leaves_from_node(sub_dag, top_node):
    descendant_leaves = set()

    #print "Top node is: %s"%str(top_node)
    #print "Successors: %s"%str(godag.successors(top_node))
    for ch in set(DAGs_r[ont][top_node])&set(sub_dag):
        if not set(DAGs_r[ont][ch])&set(sub_dag):
            descendant_leaves.add(ch)
        else:
            descendant_leaves.update(get_leaves_from_node(sub_dag,ch))
    return descendant_leaves

def find_candidate_fcns(unique_terms):
    '''os.devnull
    Using the parsed evidence, this places the evidence set
    and modifies the gene ontology graph in the SIFTER 2.0 way.
    '''
    # For each protein in the evidence set, store the annotation
    # into the evidence graph            
    annotated_term_nodes = []
    for go_term in unique_terms:
        if go_term not in DAGs[ont]:
            print "GO term, %s doesn't seem to be named in your ontology."%go_term            
            continue
        annotated_term_nodes.append(go_term)

    go_subdag = get_ontology_subdag(annotated_term_nodes=annotated_term_nodes)

    root_node = roots[ont]
    candidate_fcns=get_leaves_from_node(go_subdag, root_node)

    return candidate_fcns

def extract_for_each_family(pfam_id):
    r_flag=0
    res_pickle_file=results_folder+'/%s_result.pickle'%pfam_id
    results={}
    if os.path.exists(res_pickle_file):    
        results=pickle.load(open(res_pickle_file,'r')) 
        r_flag=1
    else:
        processing_file=results_folder+'/%s.sifterj.processing'%pfam_id
        if os.path.exists(processing_file):    
            print "SIFTER is currently processing the family %s "%(pfam_id)
            r_flag=2    
        else:
            print "SIFTER has not yet run on family %s "%(pfam_id)    
            r_flag=0    
    if results:
        results['results']={g:v for g,v in results['results'].iteritems() if g in pfam_2_gene[pfam_id]}
    return results,r_flag
if __name__=="__main__":
    
    # Initialization
    params_mysql = {\
    'db_address': 'localhost',
    'db_username': 'root',
    'db_password': '',
    'db_name': 'sifter_db'
    }

    main_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    obo_file=main_dir+'/data/go.obo'

   
    taxid=''
    query_families=[]
    query_proteins=[]
    input_file=''
    pfams=[]
    hit_file=''
    all_genes=0
    # Check for options
    
    opts, args = getopt.getopt(sys.argv[1:], "hAp:s:f:",['ip=','dbname=','dbpass=','dbuser=','dbaddr=','hit_file=']) 
    if len(args) != 2:
        usage()
        sys.exit()
    choices=[]
    if len(opts)>0:
        for o, a in opts:
            if o == "-p":
                splited =a.strip().split(',')
                query_proteins=list(set([w for w in splited if w]))  
                choices.append('p')
            elif o == "-s":
                taxid = a
                choices.append('s')                
            elif o == "-f":
                splited =a.strip().split(',')
                query_families=list(set([w for w in splited if w]))  
                choices.append('f')                
            elif o == "--ip":
                input_file = a
                choices.append('ip')                
            elif o == "-A":
                all_genes = 1
                choices.append('A')                
            elif o == "--hit_file":
                hit_file = a
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
        print "Please use one of the '-p -s -f --ip -A' options to enter your query.\n"
        sys.exit()
    elif len(choices)>1:
        print "\nERROR: Please use ONLY one of the '-p -s -f --ip -A' options to enter your query.\n"
        sys.exit()
   
    ###
    results_folder=args[0]
    output_file=args[1]
   
    
    db_mysql = MySQLdb.connect(host=params_mysql['db_address'],
                   user=params_mysql['db_username'],
                   passwd=params_mysql['db_password'],
                   db=params_mysql['db_name'],
                   cursorclass=MySQLdb.cursors.DictCursor)
    print "\n\n--------------Reading the query information------------"
    if hit_file:
        if not os.path.exists(hit_file):
            print "\nERROR: No Pfam hit file at %s.\n"%hit_file
            sys.exit()
        else:
            pfam_2_gene_hit,gene_2_pfam_hit=find_pfam_2_gene_from_file(hit_file)        
            
    
    if query_families:
        res=find_genes_for_pfams(query_families)
        pfam_2_gene,gene_2_pfam=find_pfam_2_gene(res)        
        pfams=pfam_2_gene.keys()            
        print "Extract SIFTER results for Pfam family %s"%','.join(query_families)
    elif query_proteins or input_file:
        if input_file:
            if not os.path.exists(input_file):
                print "\nERROR: No file exists at %s\n"%input_file
                sys.exit()
            f = open(input_file, 'r')
            a=f.read()
            splited =re.split(' |,|;|\n',a.strip())
            query_proteins=list(set([w for w in splited if w]))
        if not hit_file:
            res=find_pfam_for_genes(query_proteins)
            pfam_2_gene,gene_2_pfam=find_pfam_2_gene(res)
            pfams=pfam_2_gene.keys()
        else:
            gene_2_pfam={p:gene_2_pfam_hit[p] for p in query_proteins if p in gene_2_pfam_hit}
            pfam_2_gene={}
            for g,fs in gene_2_pfam.iteritems():
                for f in fs:
                    if not f in pfam_2_gene:                    
                        pfam_2_gene[f]=set([])
                    pfam_2_gene[f].add(g)
            pfams=pfam_2_gene.keys()
        print "Extract SIFTER results for %s Pfam families for  %s query proteins"%(len(pfams),query_proteins)
    elif taxid:
        if not hit_file:
            res=find_pfam_for_taxid(taxid)
            pfam_2_gene,gene_2_pfam=find_pfam_2_gene(res)
            pfams=pfam_2_gene.keys()
            print "Extract SIFTER results for %s Pfam families for query species (taxid=%s) with %s proteins"%(len(pfams),taxid,len(gene_2_pfam))
        else:
            gene_2_pfam=gene_2_pfam_hit;
            pfam_2_gene=pfam_2_gene_hit;
            pfams=pfam_2_gene.keys()
            print "-s will be ignored. We will extract SIFTER results for all %s genes in the hit-file"%(len(gene_2_pfam))
    elif all_genes==1:
        if not hit_file:
            print "\nERROR: -A option can only used for novel genomes (hit_file should be provided)\n"
            sys.exit()
        else:
            gene_2_pfam=gene_2_pfam_hit;
            pfam_2_gene=pfam_2_gene_hit;
            pfams=pfam_2_gene.keys()
            print "We extract SIFTER results for all %s genes in the hit-file"%(len(gene_2_pfam))

    if not pfams:
        print "\nERROR: There are no pfam families for your input query."
        print "Please use one of the '-p -s -f --ip -A --hit-file' options to enter your query.\n"
        sys.exit()

    pfds=get_pfds(pfams)

    tree_sizes = {}
    for p in pfds.keys():
        tree_sizes[p] = pfds[p]['num_full']
    sorted_fams = sorted(pfds.keys(), key=lambda k:pfds[k]['num_full'])
    print "Number of families:" ,len(sorted_fams)

    print "\n-----------------Reading the ontology file----------------"
    ont='molecular_function'
    go_dict,DAGs,DAGs_r,roots=parse_GO_OBO(obo_file)


    print "\n------------Prepare the necessary query files-------------"
    pfams_to_process = []    
    results={}
    r_flags={}
    for i,pfam_id in enumerate(sorted_fams):
        res,r_flag=extract_for_each_family(pfam_id)
        if r_flag==1:
            results[pfam_id]=res
            print "Results for %s domains extracted from %s "%(len(results[pfam_id]),pfam_id) 
        r_flags[pfam_id]=r_flag   
    for pfam_id,res in results.iteritems():
        candids=res['candids']
        print '\t\t'+'%s'%candids
        for g in res['results']:
            print g, '%s'%res['results'][g]
    
    genes_with_pred=set([])
    with open(output_file, 'w') as myfile:
        for g,fs in gene_2_pfam.iteritems():
            g_printed=0
            for f in fs:
                if f in results:
                    res=results[f]
                    if g not in res['results']:
                        continue
                    if g_printed==0:
                        g_printed=1
                        line="Protein: %s"%g
                        myfile.write(line+'\n')
                        print line
                    genes_with_pred.add(g)
                    candids=res['candids']
                    for d in res['results'][g]:
                        line='Domain: %s/%s\t%s'%(g,d[0],f)
                        myfile.write(line+'\n')
                        print line
                        line='Predictions:\tGO ID\t\tTerm name\tConfidence Score'
                        myfile.write(line+'\n')
                        print line
                        
                        preds=[[i,score] for i,score in enumerate(d[1])]
                        preds=sorted(preds,key=lambda x:x[1],reverse=True)
                        for i,r in preds:
                            name=''
                            if candids[i] in go_dict:
                                if 'name' in go_dict[candids[i]]:
                                    if go_dict[candids[i]]['name']:
                                        name=go_dict[candids[i]]['name'][0]
                            line='\t\t%s\t%s\t%0.2f'%(candids[i],name,max(min(r,1),0))
                            myfile.write(line+'\n')
                            print line
            if g_printed==1:
                myfile.write('\n')                
                print
            
    print "Predictions for %s genes are written to the output file %s"%(len(genes_with_pred),output_file)
    print "-------------------Extracting the results is Done----------------------"
    
        


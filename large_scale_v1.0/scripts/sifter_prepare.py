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
import subprocess
import random
import glob

def usage():
    print "\n-----------------------------------------------------------------"
    print "Usage: "
    print "   sifter_prepare.py [options] <families_data_folder> <output_folder>"
    print "-----------------------------------------------------------------\n"
    print "Examples:"
    print "   sifter_prepare.py -p C0JYY2_HUMAN ../example/fam_data ../example/queries\n"
    print "   sifter_prepare.py -s 9823  ../example/fam_data ../example/queries\n"
    print "   sifter_prepare.py -f PF03818 -a ../example/fam_data ../example/queries\n"
    print "   sifter_prepare.py --ip ../example/protein_list.txt -r ../example/fam_data ../example/queries\n"
    print "   sifter_prepare.py --if ../example/family_list.txt  ../example/fam_data ../example/queries\n"
    print "   sifter_prepare.py -p C0JYY2_HUMAN -x 1e5  ../example/fam_data ../example/queries\n"
    print "   sifter_prepare.py -p C0JYY2_HUMAN -t 2  ../example/fam_data ../example/queries\n"    
    print "   sifter_prepare.py -s 9823 --dbaddr www.example.org --dbuser jack --dbpass 1234  ../example/fam_data ../example/queries\n"    
    print "   sifter_prepare.py -A --hit_file ../example/pfam_res.txt  ../example/fam_data ../example/queries\n"    
    print "This function prepares necessary files for your query to run SIFTER on."
    print "@author Sayed Mohammad Ebrahim Sahraeian (mohammad@compbio.berkeley.edu)"
    print "Please cite new paper:"
    print "-Sahraeian SME, Luo KR, Brenner SE (2015)"
    print "\nThe SIFTER algorithm presented in the following paper:"
    print "- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. Genome-scale phylogenetic function annotation of large and diverse protein families. Genome Research 21:1969-1980. \n"
    print "inputs:"
    print "        <families_data_folder>   Path to the folder where the"
    print "                                 families data are placed. You can"
    print "                                 download the precomputed data"
    print "                                 or build it using the"
    print "                                 'sifter_gather_family_data.py' script." 
    print "        <output_folder>          Path to the output folder where"
    print "                                 the necessary query files and"
    print "                                 results will be written to." 
    print "options: (you should only use one of '-p -s -f --ip -A' options.)"
    print "           -p          STRING    List of query proteins (use Uniprot ID"
    print "                                 or Accession) in comma seperated format."
    print "           -s          STRING    NCBI taxonomy ID for input species."
    print "           -f          STRING    List of Pfam families for which you"
    print "                                 want to prepare data."
    print "                                 (in comma seperated format)"
    print "           --ip        STRING    Path to the input file where the list"
    print "                                 of proteins are placed."
    print "           --if        STRING    Path to the input file where the list"
    print "                                 of families are placed."
    print "           --hit_file  STRING    Output of pfam_scan.pl file on the "
    print "                                 novel genome. This file consists of"
    print "                                 the list of pfam hits for the genome."
    print "                                 If this option is uded, we will"
    print "                                 look in this file to find Pfams"
    print "                                 instead of the SQL database."
    print "           -A                    Prepare for all Pfam families of queried"
    print "                                 novel genome. (hit_file should be provided)"
    print "           -a                    Include all experimental and"
    print "                                 non-experimental evidence" 
    print "                                 in the inference. (Defualt [if"
    print "                                 this option is not used]: only"
    print "                                 experimental evidence will be used)."
    print "           -r                    Remove all query files already prepared"
    print "                                 and rebuild the queries."
    print "           -x          INT       Maximum number of nonzero elements"
    print "                                 in the transition matrix. Should be" 
    print "                                 a number in [1e5,1e7] for reasonable" 
    print "                                 time and accuracy balance (Default=2250000)"
    print "                                 Smaller value leads to faster running time."
    print "           -t          INT       Number of functions to truncate" 
    print "                                 to in approximation [Default:"
    print "                                 adaptive based on -x option]"
    print "                                 Smaller value leads to faster running time."
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
      pfamseq.pfamseq_acc,
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
        my_gene=w['pfamseq_acc']
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

def max_fun_possible(i,thr):
    max_f=0
    for j in range(1,i+1):
        max_f_temp=max_f+comb(i,j,exact=0)
        if max_f_temp>thr:
            return [j-1,max_f]
        else:
            max_f=max_f_temp
    return [i,max_f]

def calc_numel(numTerms, maxFun):
    return pow(sum([float(comb(numTerms, i, exact=0)) for i in range(1, maxFun + 1)]), 2)

def get_criteria(numTerms, maxFun):
    if numTerms > 8:
        if maxFun > 1:
            return 1
        else:
            return 2
    else:
        return 3

def get_category(numel, famSize):
    # List of dividers based upon the number of elements NUMEL in transition matrix
    numelDivs = [65025.0, 330625.0, 1046529.0]

    # List of dividers based upon the family size FAMSIZE
    famSizeDivs = [567.0, 1637.0, 4989.0]
    
    n = sum(map(lambda x: numel > x, numelDivs))
    s = sum(map(lambda x: famSize > x, famSizeDivs))
    return (n, s)        

def est_processing_time(numTerms, famSize, maxFun,numel):
    paramsDict = {1: [-6.6940979152046394, 1.2175437752942884,   0.61437156459022535],
              2: [-3.6107074614976109, 0.91343454244972999,  0.45521131812635984],
              3: [-2.7026843343076519, 0.052132418536663394, 0.93755721899494526]}
    crit = get_criteria(numTerms, maxFun)
    line = paramsDict[crit]
    return pow(10, line[0]) * pow(numel, line[1]) * pow(famSize, line[2])

def get_upper_bound(eTime, cat, per):
    percentileDict={(0, 0): {'95': 8.3435056315411593, '99.9': 10.953643510480756},
      (0, 1): {'95': 9.4040189875556379, '99.9': 10.175590194144538},
      (0, 2): {'95': 7.0857310513064657, '99.9': 10.031292126553355},
      (0, 3): {'95': 4.3471755740354761, '99.9': 8.7766092407283836},
      (1, 0): {'95': 4.0445760101251587, '99.9': 9.5270816900332136},
      (1, 1): {'95': 2.3310236959309329, '99.9': 3.4547033474036422},
      (1, 2): {'95': 1.8195072570575042, '99.9': 2.9109043732685018},
      (1, 3): {'95': 2.0892177205927638, '99.9': 7.8978069638688924},
      (2, 0): {'95': 2.2542718513558571, '99.9': 2.9746194223225029},
      (2, 1): {'95': 2.6775509810516125, '99.9': 4.4976310858312294},
      (2, 2): {'95': 2.9809620961392786, '99.9': 4.8087748272548554},
      (2, 3): {'95': 4.4914777165287258, '99.9': 6.7709753345612205},
      (3, 0): {'95': 2.6439743599924892, '99.9': 3.3485478896514702},
      (3, 1): {'95': 2.883955861280195, '99.9': 3.9323761482164077},
      (3, 2): {'95': 3.156846158873563, '99.9': 3.904755873693849},
      (3, 3): {'95': 3.898056279279821, '99.9': 4.4261063907623219}}
    percentiles = percentileDict[(cat[0],cat[1])][per]
    return eTime*percentiles

def format_times(times):
    if not times:
        return times
    t = times[0]
    if t < 1:
        return ['%.1f seconds' % (60 * t) for t in times]
    elif t < 60:
        return ['%.1f minutes' % t for t in times]
    elif t < 60 * 24:
        return ['%.1f hours' % (t / 60) for t in times]
    elif t < 60 * 24 * 365:
        return ['%.1f days' % (t / 60 / 24) for t in times]
    else:
        return ['%.1f years' % (t / 60 / 24 / 365) for t in times]

def estimate_time(numTerms, famSize,t_lev):
    tableBody = []
    pers = ['95','99.9']
    maxFun=min(t_lev,numTerms)
    numel = calc_numel(numTerms, maxFun)
    eTime = est_processing_time(numTerms, famSize, maxFun,numel)
    eTime = max(eTime, 1.0) # set minimum estimated time to 1 minute
    cat = get_category(numel, famSize)
    row = [maxFun]
    times = [eTime]
    for j in range(len(pers)):
        upper = get_upper_bound(eTime, cat, pers[j])
        times.append(upper)
    row.extend(times)
    row.extend(format_times(times))
    return row


def store_run_data(pfam_id):   
    data={}
    data['pfam_id']=pfam_id
    data['query_proteins']=[]#pplacer_queries[pfam_id]
    data['query_protein_accs']={}#{k['id']:pfamseq_acc_for_id[k['id']] for k in pplacer_queries[pfam_id]}
    data['tree_size']=tree_sizes[pfam_id]
    data['evidence_constraints']=evidence_allowed
    data['tree_loc']=reconciled_folder+'/%s'%pfam_id+"_reconciled.xml"
    data['tree_format']='phyloxml',
    data['annotation_loc']=evidence_folder+'/%s.pli'%pfam_id
    data['annotation_loc_pickle']=evidence_folder+'/%s.pickle'%pfam_id
    data['annotation_format']='pli'
        
    print "Loading goa annotations for %s..."%pfam_id   
    evidence_pickle_file = evidence_folder+'/%s.pickle'%pfam_id  # file with annotations 
    rand_id_1=random.randint(1000000,9999999)                    
    if os.path.exists(evidence_pickle_file+'.gz'):
        if os.path.exists('%s.%d'%(evidence_pickle_file,rand_id_1)):
            subprocess.check_call("rm %s"%(evidence_pickle_file),shell=True)                    
        subprocess.check_call("gunzip -c %s.gz > %s.%d"%(evidence_pickle_file,evidence_pickle_file,rand_id_1),shell=True)

    [evidence_file2,pfam_anns, pp, seq_lookup] = pickle.load(open('%s.%d'%(evidence_pickle_file,rand_id_1), 'rb'))
    if os.path.exists('%s.%d'%(evidence_pickle_file,rand_id_1)):
        subprocess.check_call("rm %s.%d"%(evidence_pickle_file,rand_id_1),shell=True)                    

    # Filter for only experimental annotations.
    unique_terms=set([])
    num_ev = 0    
    for prot_id, anns in pfam_anns.iteritems():
        alwd_ev = [a['acc'] for a in anns if a['code'] in evidence_allowed]
        unique_terms=unique_terms.union(set(alwd_ev))
        # a is an annotation with a function and a code that says where the function came from
        # keep this annotation if it was gotten through experiments
        if len(alwd_ev) > 0:
            num_ev += 1
    print pfam_id,'has' ,num_ev, "annotated proteins with allowed evidence type"
    data['num_ev_prots']=num_ev
    data['num_any_ev_prots']=len(pfam_anns)
    
    if len(unique_terms)==1 and ('GO:0003674' in unique_terms):
        num_ev=0
    if num_ev>0:
        # Input evidence
        candidate_fcns=find_candidate_fcns(unique_terms)  
        evidence_format = 'pli'
        data['n_terms'] = len(candidate_fcns)
        data['candids'] = candidate_fcns
        thr=max_fun_possible(data['n_terms'],np.sqrt((mx_numel)))[0]
        if truncation_level:
            thr=min(thr,truncation_level)
        row=estimate_time(data['n_terms'],data['tree_size'],thr)
        data['e_time']=row
        print "Number of functions:",data['n_terms']
        print "We will use truncation level = %s"%row[0]
        print "Estimated running time for family %s = %s (95%% confidence upper bound = %s)"%(pfam_id,row[4],row[5])            
        pickle.dump(data, open(queries_folder+'/%s_query.pickle'%pfam_id, 'wb'))
        print "Processed evidence from:", pfam_id
    else:
        print "No candidate functions: SIFTER will not be run on this family."    
        data['n_terms'] = 0
        pickle.dump(data, open(queries_folder+'/NQ/%s_query.pickle'%pfam_id, 'wb')) 
    
def prepare_for_each_family(pfam_id):
    reconciled_fname = reconciled_folder+'/%s'%pfam_id            
    evidence_file = evidence_folder+'/%s.pli'%pfam_id
    evidence_pickle_file = evidence_folder+'/%s.pickle'%pfam_id              
    queries_to_process=[]
   
   
    skip_flag=0
    if not(os.path.isfile(reconciled_fname+"_reconciled.xml.gz")):   
        print "\nERROR: No tree file %s. Skip this family.\n"%(reconciled_fname+"_reconciled.xml.gz")
        skip_flag=1
    if not(os.path.isfile(evidence_file+'.gz')):
        print "\nERROR: No evidence file %s.gz. Skip this family.\n"%(evidence_file)
        skip_flag=1

    if not(os.path.isfile(evidence_pickle_file+'.gz')):
        print "\nERROR: No evidence file %s.gz. Skip this family.\n"%(evidence_pickle_file)
        skip_flag=1

    q_flag=0
    if (skip_flag==0):
        if not(os.path.isfile(queries_folder+'/%s_query.pickle'%pfam_id)) and  not(os.path.isfile(queries_folder+'/NQ/%s_query.pickle'%pfam_id)):
            store_run_data(pfam_id)
        else:
            print "Family %s already prepared."%(pfam_id)
        if (os.path.isfile(queries_folder+'/%s_query.pickle'%pfam_id)):
            q_flag=1
         
    return q_flag
     

if __name__=="__main__":
    
    # Initialization
    params_mysql = {\
    'db_address': 'localhost',
    'db_username': 'root',
    'db_password': '',
    'db_name': 'sifter_db'
    }

    evidence_constraints_exp = [
        # Experimental
        'EXP',  # Experiment
        'IDA',  # Direct Assay
        'IPI',  # Physical Interaction
        'IMP',  # Mutant Phenotype
        'IGI',  # Genetic Interaction
        'IEP',  # Expression Pattern
        # Author Statements
        'TAS',  # Traceable Author Statement
        'NAS',  # Non-traceable Author Statement
        ]
    
        
    evidence_constraints_all = [
        # Experimental
        'EXP',  # Experiment
        'IDA',  # Direct Assay
        'IPI',  # Physical Interaction
        'IMP',  # Mutant Phenotype
        'IGI',  # Genetic Interaction
        'IEP',  # Expression Pattern
        # Author Statements
        'TAS',  # Traceable Author Statement
        'NAS',  # Non-traceable Author Statement
        # Computational Analysis Evidence Codes
        'ISS',  # Sequence/Structural Similarity
        'ISO', # Sequence Orthology
        'ISA', # Sequence Alignment
        'ISM', # Sequence Model
        'IGC', # Genomic Context
        'IBA', # Biological aspect of ancestor
        'IBD', # Biological aspect of descendant
        'IKR', # Key Residues
        'IRD', # Rapid Divergence
        'RCA',  # Reviews Computational Analysis
        # Curator Statement
        'IC',  # Curator
        'ND',  # No biological data available
        # Automatically assigned
        'IEA',  # Electronic Annotation
        # Obsolete
        'NR'  # Not recorded
        ]
    main_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    obo_file=main_dir+'/data/go.obo'
    evidence_allowed = evidence_constraints_exp


    taxid=''
    query_families=[]
    query_proteins=[]
    p_input_file=''
    f_input_file=''
    truncation_level=0
    pfams=[]
    remove_query_files=0
    mx_numel=2250000
    hit_file=''
    all_fams=0
    # Check for options
    
    opts, args = getopt.getopt(sys.argv[1:], "hraAp:s:f:t:x:",['ip=','if=','dbname=','dbpass=','dbuser=','dbaddr=','hit_file=']) 
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
                p_input_file = a
                choices.append('ip')                
            elif o == "--if":
                f_input_file = a
                choices.append('if')                
            elif o == "-A":
                all_fams = 1
                choices.append('A')                
            elif o == "--hit_file":
                hit_file = a
            elif o == "-a":
                evidence_allowed = evidence_constraints_all
            elif o == "-r":
                remove_query_files=1
            elif o == "-x":
                mx_numel=int(float(a))
            elif o == "-t":
                truncation_level=int(a)
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
        print "Please use one of the '-p -s -f --ip --if -A' options to enter your query.\n"
        sys.exit()
    elif len(choices)>1:
        print "\nERROR: Please use ONLY one of the '-p -s -f --ip --if -A' options to enter your query.\n"
        sys.exit()

   
    families_data_path=args[0]
    
    if not os.path.exists(families_data_path):
        print "\nERROR: families_data directory ( %s ) does not exist\n"%families_data_path
        sys.exit()
    evidence_folder=families_data_path+'/annotations'
    if not os.path.exists(evidence_folder):
        print "\nERROR: annotations directory ( %s ) not exists\n"%evidence_folder
        sys.exit()
    reconciled_folder=families_data_path+'/reconciled_trees'    
    if not os.path.exists(reconciled_folder):
        print "\nERROR: reconciled_trees directory ( %s ) not exists\n"%reconciled_folder
        sys.exit()
    alignment_folder=families_data_path+'/alignments'
    if not os.path.exists(alignment_folder):
        print "\nERROR: alignment directory( %s ) not exists\n"%alignment_folder
        sys.exit()
    
    ###
    output_path=args[1]
    queries_folder=output_path
    if remove_query_files==1:
        os.system('rm -rf %s'%queries_folder)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    
    db_mysql = MySQLdb.connect(host=params_mysql['db_address'],
                   user=params_mysql['db_username'],
                   passwd=params_mysql['db_password'],
                   db=params_mysql['db_name'],
                   cursorclass=MySQLdb.cursors.DictCursor)
    queries_folder_NQ=output_path+'/NQ'
    if not os.path.exists(queries_folder_NQ):
        os.mkdir(queries_folder_NQ)
    
    prepared_queries=glob.glob(output_path+'/*.pickle')+glob.glob(queries_folder_NQ+'/*.pickle')
    prepared_queries=[(w.split('/')[-1]).split('_query')[0] for w in prepared_queries]
    already_prepared_fams=[]
    
    print "\n\n--------------Reading the query information------------"
    if hit_file:
        if not os.path.exists(hit_file):
            print "\nERROR: No Pfam hit file at %s.\n"%hit_file
            sys.exit()
        else:
            pfam_2_gene_hit,gene_2_pfam_hit=find_pfam_2_gene_from_file(hit_file)        
            
    
    if query_families or f_input_file:
        if f_input_file:
            if not os.path.exists(f_input_file):
                print "\nERROR: No file exists at %s\n"%f_input_file
                sys.exit()
            f = open(f_input_file, 'r')
            a=f.read()
            splited =re.split(' |,|;|\n',a.strip())
            query_families=list(set([w for w in splited if w]))            
        
        already_prepared_fams=list(set(query_families)&set(prepared_queries))
        toprep_families=list(set(query_families)-set(prepared_queries))
        print "%s out of %s Families have already prepared. We will Check %s others."%(len(already_prepared_fams),len(query_families),len(toprep_families))               
        query_families=toprep_families
        res=find_genes_for_pfams(query_families)
        pfam_2_gene,gene_2_pfam=find_pfam_2_gene(res)        
        pfams=pfam_2_gene.keys()
        for f in set(query_families)-set(pfams):
            print "Family %s is not in the SQL database."%f       
        if query_families:
            print "Run SIFTER for Pfam families: %s"%','.join(query_families)
    elif query_proteins or p_input_file:
        if p_input_file:
            if not os.path.exists(p_input_file):
                print "\nERROR: No file exists at %s\n"%p_input_file
                sys.exit()
            f = open(p_input_file, 'r')
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
        print "Run SIFTER for %s Pfam families for %s query proteins"%(len(pfams),len(query_proteins))
    elif taxid:
        if not hit_file:
            res=find_pfam_for_taxid(taxid)
            pfam_2_gene,gene_2_pfam=find_pfam_2_gene(res)
            pfams=pfam_2_gene.keys()
            print "Run SIFTER for %s Pfam families for query species (taxid=%s) with %s proteins"%(len(pfams),taxid,len(gene_2_pfam))
        else:
            gene_2_pfam=gene_2_pfam_hit;
            pfam_2_gene=pfam_2_gene_hit;
            pfams=pfam_2_gene.keys()
            print "-s will be ignored. We will run on all %s Pfam families in the hit-file"%(len(pfams))
    elif all_fams==1:
        if not hit_file:
            print "\nERROR: -A option can only used for novel genomes (hit_file should be provided)\n"
            sys.exit()
        else:
            gene_2_pfam=gene_2_pfam_hit;
            pfam_2_gene=pfam_2_gene_hit;
            pfams=pfam_2_gene.keys()
            print "We will run on all %s Pfam families in the hit-file"%(len(pfams))

    
    if (not pfams):
        if (not already_prepared_fams):
            print "\nERROR: There are no pfam families for your input query."
            print "Please use one of the '-p -s -f --ip -A --hit-file' options to enter your query.\n"
            sys.exit()
        else:
            print "-------------------Preperation is Done----------------------"
            print "All of your %s query families have been already prepared."%(len(already_prepared_fams))
            print "\nNext step is to run 'sifter_run.py'."
            print "You may exclude some of the more complex families there.\n"
    else:
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
        for i,pfam_id in enumerate(sorted_fams):
            q_flag=prepare_for_each_family(pfam_id)
            if q_flag==1:
                pfams_to_process.append(pfam_id)
                print "Input file prepared for %s (%d out of %d families)"%(pfam_id,i+1,len(sorted_fams))
        nqs=0
        for pfam_id in sorted_fams:
            nqf=queries_folder_NQ+'/%s_query.pickle'%pfam_id
            if (os.path.isfile(nqf)):
                nqs+=1
        errors=len(sorted_fams)-len(pfams_to_process)-nqs

        

        if len(pfams_to_process)>0:  
            
            e_times = []
            total_e=0
            total_95=0
            with open(output_path+'/running_estimation.csv', 'w') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
                spamwriter.writerow(['Family','number of candidate functions','Family size','Truncation level','Estimated running time','95%% confidence upper bound','99.9%% confidence upper bound'])
                for pfam_id in pfams_to_process:
                    qfile=queries_folder+'/%s_query.pickle'%pfam_id
                    query_data = pickle.load(open(qfile, "rb" ))
                    row=query_data['e_time']
                    spamwriter.writerow([pfam_id,query_data['n_terms'],query_data['tree_size'],row[0],row[4],row[5],row[6]])
                    total_e +=(row[1])
                    total_95 +=(row[2])
            
            fe=format_times([total_e,total_95])
            if already_prepared_fams:
                print "%s of your query families have been already prepared."%(len(already_prepared_fams))
                print "Here is the statistics for the rest of queries."

                
            print "\nFiles are prepared for %d out of %d families. (%s missed due to errors, %s are skipped duo to no candidate functions)"%(len(pfams_to_process),len(sorted_fams),errors,nqs)
            print "-------------------Preperation is Done----------------------"
            print "There are %s families to run SIFTER on."%(len(pfams_to_process))
            print "\nTotal estimated time for your query is %s (95%% confidence upper bound = %s)."%(fe[0],fe[1])
            print "Details for individual families are written in '%s/running_estimation.csv'"%output_path
            print "\nNext step is to run 'sifter_run.py'."
            print "You may exclude some of the more complex families there.\n"
            
        else:
            print "\nFiles are prepared for %d out of %d families. (%s missed due to errors, %s are skipped duo to no candidate functions)"%(len(pfams_to_process),len(sorted_fams),errors,nqs)

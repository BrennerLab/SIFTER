#!/usr/bin/python
import os
import pickle
from Bio import AlignIO,Phylo,SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo import PhyloXMLIO,PhyloXML
import ete2
from ete2 import Phyloxml
import sys
import StringIO
import gzip
import subprocess
import getopt
import copy
import _mysql as mysql
import _mysql_exceptions as mysql_exceptions
import MySQLdb.cursors
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+'/lib/pplacer_python_scripts')
from PplacerWrapper import GuppyWrapper
from PplacerWrapper import PplacerWrapper
from PfamHMMERAlign import PfamHMMERAlign
import Queue
import threading
import re
import random

def usage():
    print "\n-----------------------------------------------------------------"
    print "Usage: "
    print "   sifter_gather_family_data.py [options] <families_data_folder>"
    print "-----------------------------------------------------------------\n"
    print "Examples:"
    print "   sifter_gather_family_data.py -f PF12491,PF13820 ../example/fam_data\n"
    print "   sifter_gather_family_data.py -i ../example/family_list.txt ../example/fam_data\n"
    print "   sifter_gather_family_data.py -f PF09172 --dbaddr www.example.org --dbuser jack --dbpass 1234 ../example/fam_data\n"    
    print "This function gather necessary 'alignment', 'tree', and 'evidence' files needed to run SIFTER for each query family."
    print "\nTo gather data for a set of sequences (of a novel genome) which is not already in Pfam, you may enter the sequences file, the pfam hit file (found by pfam_scan.pl), and the NCBI taxonomy ID of your geneome of interest:"    
    print "Examples:"
    print "   sifter_gather_family_data.py -f PF07083,PF14343,PF03818 --seq_file ../example/myseq.fasta --hit_file ../example/pfam_res.txt --taxid 1192197 ../example/fam_data\n"
    print "   sifter_gather_family_data.py -A --seq_file myseq.fasta --hit_file ../example/pfam_res.txt --taxid 1192197 ../example/fam_data\n"
    print "@author Sayed Mohammad Ebrahim Sahraeian (mohammad@compbio.berkeley.edu)"
    print "Please cite new paper:"
    print "-Sahraeian SME, Luo KR, Brenner SE (2015)"
    print "The SIFTER algorithm presented rin the following paper:"
    print "- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. Genome-scale phylogenetic function annotation of large and diverse protein families. Genome Research 21:1969-1980. \n"
    print "inputs:"
    print "        <families_data_folder>   Path to the folder where the"
    print "                                 families data will be placed." 
    print "options: (you should only use one of -i or -f  or -A options.)"
    print "           -f          STRING    List of Pfam families for which you"
    print "                                 want to gather necessary data."
    print "                                 (in comma seperated format)"
    print "           -i          STRING    Path to the input file where the lis"
    print "                                 of families are placed."
    print "           --seq_file  STRING    The fasta format input sequences file"
    print "                                 of the novel genome."
    print "           --hit_file  STRING    Output of pfam_scan.pl file on the "
    print "                                 novel genome. This file consists of"
    print "                                 the list of pfam hits for the genome."
    print "           --taxid     INT       The NCBI taxonomy ID of the genome. If"
    print "                                 the tax-ID is not in the species"
    print "                                 tree, you may enter the NCBI taxonomy"
    print "                                 ID of a close species."
    print "           -A                    Run on all Pfam families of queried"
    print "                                 novel genome."
    print "           -n          INT       Number of threads (Default=4)"
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


# See how many sequences are in each family, add it to pfds
def get_pfds(pfams,db):
    pfds={}
    for p in pfams:
        sql_q="""select pfamA.num_full, pfamA.number_species, pfamA.pfamA_id, pfamA.description, group_concat(go_id) as go_ids, group_concat(term) as go_terms from pfamA left join gene_ontology on gene_ontology.auto_pfamA = pfamA.auto_pfamA where pfamA.pfamA_acc='%s' group by pfamA_acc
        """%(p)
        #AND Locus.type=1
        #AND Synonym.type=2
        r = msql(sql_q, db)
        if r:
            pfds[p]={}
            for w in r[0].keys():
                pfds[p][w]=r[0][w]
        else:
            print "%s is a wrong Pfam ID and will be excluded."%p
    return pfds

# ##2-Extract alignment information for each Pfam family
def get_pfam_alignment_by_id(pfam_id, outpt_fname,db):
    # Get file from MySQL query, gunzip, and save into query directory.
    mysql_aq =  "(select auto_pfamA from pfamA where pfamA_acc='" + pfam_id + "')"
    mysql_q = "select alignment from alignments_and_trees "             + "where auto_pfamA = "+mysql_aq+" "             + "and type='full';"
    gzipped = True# # this field is gzipped in the table
    print mysql_q
    res = msql(mysql_q, db)
    print len(res)
    if len(res) == 0:
        return
    
    t_data = ''
    f_output = open(outpt_fname, "w")
    if gzipped:
        f_gzipped = StringIO.StringIO(res[0]['alignment'])
        f = gzip.GzipFile(fileobj=f_gzipped, mode='rbU')
        t_data = f.read()
        f_output.write(t_data)
        f.close()
        f_gzipped.close()
    else:
        t_data = res[0]['alignment']
        f_output.write(t_data)
    f_output.close()


def get_alignment(pfam_id,my_db):
    outpt_fname = alignment_folder+'/%s'%pfam_id
    if not(os.path.isfile(outpt_fname+".fasta.gz")):          
        print "Saving alignment for", pfam_id
        print ""
        get_pfam_alignment_by_id(pfam_id=pfam_id, outpt_fname=outpt_fname+".sth",db=my_db)
        AlignIO.convert(outpt_fname+".sth","stockholm",outpt_fname+".fasta","fasta")
        if os.path.exists('%s.fasta.gz'%(outpt_fname)):
            subprocess.check_call("rm  %s.fasta.gz"%(outpt_fname),shell=True)
        subprocess.check_call("gzip %s.fasta"%(outpt_fname),shell=True)
        



def find_pfam_hits_from_file(hit_file,my_sequence_file):
    pplacer_queries={}
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
                print "Number of columns are not writ in this line, will skip the line:\n %s"%line
                continue
            r={k:row[i] for i,k in enumerate(keys)}
            if r['significance']=='1':
                pfam_id=r['hmm acc'][0:r['hmm acc'].find('.')]                
                gene_pplacer_id=r['seq id']
                gene_pplacer_id+='/'+r['envelope start']+'-'+r['envelope end']
                if pfam_id not in pplacer_queries.keys():
                    pplacer_queries[pfam_id]=[{'pplacer_id':gene_pplacer_id,'id':r['seq id'],'seq_start':int(r['envelope start']),'seq_end':int(r['envelope end'])}]
                else:
                   pplacer_queries[pfam_id].append({'pplacer_id':gene_pplacer_id,'id':r['seq id'],'seq_start':int(r['envelope start']),'seq_end':int(r['envelope end'])}) 
    gene_seq={}
    handle = open(my_sequence_file, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        gene_seq[record.id]=record.seq

    print "Your queried novel genome has:"
    print len(gene_seq), "sequences"
    print len(set([v['id'] for w in pplacer_queries.values() for v in w])), "genes in pfam"
    print len(pplacer_queries), "pfam families\n"
    return     gene_seq,    pplacer_queries





def find_tax_name(ncbi_taxid,db):
    sql_q="""select species from ncbi_taxonomy where ncbi_taxid='%s' limit 1;
    """%(ncbi_taxid)
    my_id = msql(sql_q, db)
    if my_id:
        my_id=my_id[0]['species']
    return my_id


def find_ncbi_taxid(pfamseq_id,db):
    sql_q="""select ncbi_taxid,species from pfamseq where pfamseq_id='%s';
    """%(pfamseq_id)
    my_id = msql(sql_q, db)
    if not my_id:
        if pfamseq_id in gene_seq:
            my_id=({'ncbi_taxid': my_taxid, 'species': my_taxid},)
    return max([int(w['ncbi_taxid']) for w in my_id])

def find_taxid(unip_id,db):
    sql_q="""select tax_id from goa_db where uniprot_id='%s';
    """%(unip_id)
    my_id = msql(sql_q, db)
    if my_id:
        return max([int(w['tax_id']) for w in my_id])
    else:
        return my_id

def find_tax_id_unip(unip_id,db):
    sql_q="""SELECT pfamseq_acc
                FROM pfamseq
                WHERE pfamseq_id='%s';
    """%(unip_id)
    my_id = msql(sql_q, db)
    if not(my_id):
        tax_id=find_ncbi_taxid(unip_id,db)
        if tax_id>0:
            sp=find_tax_name(tax_id,db)            
            return(tax_id,sp)
    if not(my_id):      
        if unip_id in gene_seq:
            tax_id=my_taxid
            sp=my_taxid
            return(tax_id,sp)
        else:
            my_id=({'pfamseq_acc':unip_id.split('_')[0]},)
    
    tax_id=find_taxid(my_id[0]['pfamseq_acc'],db)
    if tax_id:
        sp=find_tax_name(tax_id,db)
        if sp:
            return (tax_id,sp)
        else:
            return(tax_id,'')
        
    sql_q="""SELECT ncbi_id 
                FROM uniprot_2_NCBI 
                WHERE uniprot_id='%s';
    """%(my_id[0]['pfamseq_acc'])
    my_id2 = msql(sql_q, db)
    if my_id2:
        tax_id=max([w['ncbi_id'] for w in my_id2])
        sp=find_tax_name(tax_id,db)
        return (tax_id,sp)
    else:
        tax_id=find_ncbi_taxid(unip_id,db)
        if tax_id>0:
            sp=find_tax_name(tax_id,db)  
            return(tax_id,sp)
        
        if unip_id in gene_seq:
            tax_id=my_taxid
            sp=my_taxid
            return(tax_id,sp)            
        else:
            p_code=unip_id.split('_')[1]
            if p_code in zero_taxids:
                tax_id=zero_taxids[p_code]
                sp=find_tax_name(tax_id,db)
                return (tax_id,sp)
            else:
                tax_id=0
                sp=''       
                return (tax_id,sp)

def find_best_taxid(my_id,db):
    my_id=int(my_id)
    mynode=orig_sp_tree_0.search_nodes(name='%d'%my_id)
    if mynode:
        des=list(set([int(w.name) for w in mynode[0].iter_descendants()]) & set(all_species_txids))
        if des:
            return des[0],find_tax_name(des[0],db)
        else:
            for node in mynode[0].iter_ancestors():
                if node.name=="NoName":
                    continue
                des=list(set([int(w.name) for w in node.iter_descendants()]) & set(all_species_txids))
                if des:
                    return des[0],find_tax_name(des[0],db)                
    return -my_id,''
    


#reconcile the tree
def reconcile_tree(gene_tree_file,reconciled_file,rec_tag,pfam_id,db):
    if (os.path.isfile(rec_tag+'ids.pickle')) and  (pplacer_flag==1): 
        id_information = pickle.load(open(rec_tag+'ids.pickle', 'rb'))      
        existing_genes=id_information['existing_genes']
        Sequnces=[]
        p_ids=[]
        new_genes=set([w['id'] for w in pplacer_queries[pfam_id]])
        if not (new_genes-set(existing_genes)):
            print "All %s Genes for family %s have already been placed in the reconciled tree."%(len(new_genes),pfam_id)
            print "Skip Reconciliation for %s"%pfam_id
            return

    txid_file=rec_tag+'txid.xml'       
    if not(os.path.isfile(rec_tag+'ids.pickle')) or not(os.path.isfile(reconciled_file+'.gz')) or  (pplacer_flag==1): 
        print "Running Reconciliation for: %s"%pfam_id
        
        rand_id=random.randint(1000000,9999999)        
        subprocess.check_call("gunzip -c %s/%s.nw.gz > %s.%d"%(tree_folder,pfam_id,gene_tree_file,rand_id),shell=True)
        tree = ete2.PhyloTree('%s.%d'%(gene_tree_file,rand_id), format=0)
        tree.resolve_polytomy()
        tree.write(format=0, outfile=txid_file+'.tmp.nw')
        if os.path.exists('%s.%d'%(gene_tree_file,rand_id)):
            subprocess.check_call("rm  %s.%d"%(gene_tree_file,rand_id),shell=True)

        Phylo.convert(txid_file+'.tmp.nw', 'newick', txid_file+'.tmp.xml', 'phyloxml')
        treexml = PhyloXMLIO.read(open(txid_file+'.tmp.xml','r'))
        tree = treexml[0]
        treexml.attributes.pop('schemaLocation', None)  # not supported by Forester
        tree.rooted = True
        my_ids=set([])
        my_query_by_taxid={}
        for leaf in tree.clade.find_clades(terminal=True):
            up_name = leaf.name.split('/')[0]
            tax_id,tax_name=find_tax_id_unip(up_name,db)
            if tax_id not in all_species_txids:
                if tax_id in merged_taxid.keys():
                    tax_id=merged_taxid[tax_id]
                    tax_name=find_tax_name(tax_id,db)
                if tax_id in best_taxid_map.keys():
                    tax_id=best_taxid_map[tax_id]
                    tax_name=find_tax_name(tax_id,db)
                else:
                    tax_id0=tax_id
                    tax_id,tax_name=find_best_taxid(tax_id,db)
                    if tax_id>0:
                        best_taxid_map[tax_id0]=tax_id
            if tax_id<0:
                if (-tax_id) in merged_taxid.keys():
                    tax_id=merged_taxid[-tax_id]
                    tax_name=find_tax_name(tax_id,db)
            if tax_id in my_query_by_taxid:
               my_query_by_taxid[tax_id].append(up_name)
            else:
               my_query_by_taxid[tax_id]=[up_name]
            my_ids.add(tax_id)
            my_tax_id = PhyloXML.Id(tax_id, provider='ncbi_taxonomy')
            taxon=PhyloXML.Taxonomy(id=my_tax_id)
            taxon.scientific_name = tax_name
            leaf._set_taxonomy(taxon)
        PhyloXMLIO.write(treexml, open(txid_file,'w'))    
        os.system('rm '+txid_file+'.tmp.nw')
        os.system('rm '+txid_file+'.tmp.xml')
        print "Taxid file done for: %s"%pfam_id
        existing_ids=list(set(my_ids)&set(all_species_txids))
        existing_genes=[g for txid in my_query_by_taxid.keys() for g in my_query_by_taxid[txid] if txid in existing_ids]        
        pickle.dump({'pfam_id':pfam_id,'existing_ids':existing_ids,'existing_genes':existing_genes}, open(rec_tag+'ids.pickle', 'wb'))      
        print "Pickle file done for: %s"%pfam_id
        
       
    if os.path.exists(reconciled_file):
        os.system('rm '+reconciled_file)
    os.system("java -Xmx4g -cp %s/forester_1038.jar org.forester.application.gsdi -g %s %s/ncbi_2_fixed.xml %s"%(lib_path, txid_file, species_tree_data_path, reconciled_file))
    if os.path.exists(reconciled_file):
        if os.path.exists(reconciled_file+'.gz'):
            subprocess.check_call("rm  %s.gz"%(reconciled_file),shell=True)
        subprocess.check_call("gzip %s"%(reconciled_file),shell=True)
    os.system('rm '+rec_tag+'reconciled_species_tree_used.xml')
    os.system('rm '+rec_tag+'reconciled_gsdi_log.txt')
    os.system('rm '+txid_file)
    print "Reconciliation file done for: %s"%pfam_id




def update_hmmer_alignment(sequences, orig_alignment, hmm,tmpf):
    #tmpf = tempfile.NamedTemporaryFile(delete=False,
    #                                   suffix='.sto')
    tmpf_o=open(tmpf,'w')
    SeqIO.write(sequences, tmpf_o, "fasta")
    tmpf_o.close()
    
    pfam_hmm_align = PfamHMMERAlign()
    # Make call and do initial parsing of results into format for ProteinInformation retrievers.
    pfam_hmm_align.setup_caller(executable_locations = {
                                    'hmmpress': path_to_hmmpress,
                                    'hmmalign': path_to_hmmalign
                                },
                                params={
                                    'hmm_file': hmm,
                                    'orig_alignment': orig_alignment,
                                    'query_sequences_fasta_file': tmpf
                                })
    pfam_hmm_align.call()
    
    pfam_hmm_align.parse_results()
    return pfam_hmm_align.parsed_results

def pplacer_call(pplacer_package, aln_file, jplace_output_file):
    pplacer = PplacerWrapper()
    
    # Make call and do initial parsing of results into format for ProteinInformation retrievers.
    pplacer.setup_caller(executable_locations = {
                            'pplacer': path_to_pplacer
                        },
                        params={
                            'package_location': pplacer_package,
                            'orig_alignment': aln_file,
                            'jplace_output_file': jplace_output_file
                        })
    pplacer.call()
    pplacer.parse_results()
    return pplacer.parsed_results

def guppy_tree(jplace_file, tree_output_file):
    guppy = GuppyWrapper()
    guppy.setup_caller(executable_locations = {
                            'guppy': path_to_guppy
                        },
                        params={
                            'jplace_file': jplace_file,
                            'tree_output_file': tree_output_file
                        })
    guppy.call()
    guppy.parse_results()
    return guppy.parsed_results

def hmm_build(hmmbuild_executable_loc,
              sequence_file,
              output_file):
    '''
    Calls hmmbuild
    '''
    cmd = hmmbuild_executable_loc \
        + " " + output_file \
        + " " + sequence_file
    raw_data = subprocess.check_call(cmd, shell=True)
    
    
def taxit_create(taxit_executable_loc,
                aln_fasta,
                hmm_file,
                tree_file,
                tree_stats,
                pfam_acc,
                output_location,
                aln_stockholm):
    '''
    Calls taxit
    '''
    #taxit create --clobber --aln-fasta ./PF14424.dedup.fasta --profile ./PF14424.wholefam.hmm --tree-file ./PF14424.dedup.nh  --locus PF14424 --package-name PF14424.pplacer
    cmd = taxit_executable_loc \
        + " create --clobber" \
        + " --aln-fasta " + aln_fasta \
        + " --profile " + hmm_file \
        + " --tree-file " + tree_file \
        + " --tree-stats " + tree_stats \
        + " --locus " + pfam_acc \
        + " --package-name " + output_location
    raw_data = subprocess.check_call(cmd, shell=True)   
    input_handle = open(aln_fasta, "rU")
    output_handle = open(aln_stockholm, "w")
    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "stockholm")
    output_handle.close()
    input_handle.close()
    

def add_pplaced(pfam_id):
    if pfam_id in pplacer_queries.keys():
        print "Running PPlacer for: %s"%pfam_id

        pplace_log=pplacer_folder+'/%s_pplace_log.txt'%pfam_id
        Already_placed=[]
        if os.path.exists(pplace_log):
            with open(pplace_log, "r") as myfile:
                for line in myfile:
                    line=line.strip()
                    if not line:
                        continue
                    line=line.split('\t')
                    if not len(line)==2:
                        continue
                    Already_placed.extend(line[1].split(','))
        Sequnces=[]
        p_ids=[]
        for new_gene in pplacer_queries[pfam_id]:    
            p_id = new_gene['pplacer_id']
            if p_id in Already_placed:
                continue
            p_ids.append(p_id)
            p_seq = gene_seq[new_gene['id']][(new_gene['seq_start']-1):new_gene['seq_end']]
            Sequnces.append(SeqRecord(p_seq, id=p_id))
        if not p_ids:
            print "All %s domains for family %s have already been pplaced."%(len(Already_placed),pfam_id)
            return

		

        rand_id_1=random.randint(1000000,9999999)        
        rand_id_2=random.randint(1000000,9999999)        
        rand_id_3=random.randint(1000000,9999999)        

        subprocess.check_call("gunzip -c %s/%s.log.gz > %s/%s.log.%d"%(tree_folder,pfam_id,tree_folder,pfam_id,rand_id_1),shell=True)
        subprocess.check_call("gunzip -c %s/%s.nw.gz > %s/%s.nw.%d"%(tree_folder,pfam_id,tree_folder,pfam_id,rand_id_2),shell=True)
        subprocess.check_call("gunzip -c %s/%s.fasta.gz > %s/%s.fasta.%d"%(alignment_folder,pfam_id,alignment_folder,pfam_id,rand_id_3),shell=True)

        AlignIO.convert("%s/%s.fasta.%d"%(alignment_folder,pfam_id,rand_id_3),"fasta","%s/%s.sth.%d"%(alignment_folder,pfam_id,rand_id_3),"stockholm")

        hmm_build(hmmbuild_executable_loc=path_to_hmmbuild,
              sequence_file='%s/%s.sth.%d'%(alignment_folder,pfam_id,rand_id_3),
              output_file='%s/%s.hmm'%(pplacer_folder,pfam_id))


        taxit_create(taxit_executable_loc=path_to_taxit,
            aln_fasta='%s/%s.fasta.%d'%(alignment_folder,pfam_id,rand_id_3),
            hmm_file='%s/%s.hmm'%(pplacer_folder,pfam_id),
            tree_file='%s/%s.nw.%d'%(tree_folder,pfam_id,rand_id_2),
            tree_stats='%s/%s.log.%d'%(tree_folder,pfam_id,rand_id_1),
            pfam_acc=pfam_id,
            output_location='%s/%s_pplacer'%(pplacer_folder,pfam_id),
            aln_stockholm='%s/%s_pplacer/%s.sto.%d'%(pplacer_folder,pfam_id,pfam_id,rand_id_3),                        
            )


        if os.path.exists("%s/%s.log.%d"%(tree_folder,pfam_id,rand_id_1)):
            subprocess.check_call("rm  %s/%s.log.%d"%(tree_folder,pfam_id,rand_id_1),shell=True)
        if os.path.exists("%s/%s.nw.%d"%(tree_folder,pfam_id,rand_id_2)):
            subprocess.check_call("rm  %s/%s.nw.%d"%(tree_folder,pfam_id,rand_id_2),shell=True)
        if os.path.exists("%s/%s.fasta.%d"%(alignment_folder,pfam_id,rand_id_3)):
            subprocess.check_call("rm  %s/%s.fasta.%d"%(alignment_folder,pfam_id,rand_id_3),shell=True)
        if os.path.exists("%s/%s.sth.%d"%(alignment_folder,pfam_id,rand_id_3)):
            subprocess.check_call("rm  %s/%s.sth.%d"%(alignment_folder,pfam_id,rand_id_3),shell=True)

        output_prefix = '%s/%s_pplaced'%(pplacer_folder,pfam_id)
        updated_aln = output_prefix + '.sto'
        jplace_output_file = output_prefix + '.jplace'
        tree_output_file = output_prefix + '.tre'
        sequence_file='%s/%s.sth'%(alignment_folder,pfam_id)
        aln_fasta='%s/%s.fasta'%(alignment_folder,pfam_id)
        tree_file='%s/%s.nw'%(tree_folder,pfam_id)

        pplacer_pkg_dir ='%s/%s_pplacer'%(pplacer_folder,pfam_id)
        pplacer_pkg_hmm = '%s/%s.hmm'%(pplacer_pkg_dir,pfam_id)
        pplacer_pkg_aln = '%s/%s.sto.%d'%(pplacer_pkg_dir,pfam_id,rand_id_3)
        tmpf='%s/%s.tmpf'%(pplacer_pkg_dir,pfam_id)
        # Update alignment to include the query sequence for the hypothetical domain.
        aln_res = update_hmmer_alignment(Sequnces,
                                         orig_alignment=pplacer_pkg_aln,
                                         hmm=pplacer_pkg_hmm,tmpf=tmpf)
        aln_out = open(updated_aln,'w')
        AlignIO.write(aln_res[0], aln_out, 'stockholm')
        aln_out.close()

        # Call pplacer to generate placements onto the tree.
        pplaced = pplacer_call(pplacer_package=pplacer_pkg_dir,
                               aln_file=updated_aln,
                               jplace_output_file=jplace_output_file)

        # Use the "guppy" tool to generate the best-placement tree with query as a leaf.
        gt = guppy_tree(jplace_file=jplace_output_file,
                        tree_output_file=tree_output_file)
        #Phylo.convert(tree_output_file, 'newick', tree_output_file_xml, 'phyloxml')


        os.system('rm -rf %s'%(pplacer_pkg_dir))
        os.system('rm %s/%s.hmm'%(pplacer_folder,pfam_id))
        os.system('rm %s/%s_pplaced.jplace'%(pplacer_folder,pfam_id))
        os.system('mv %s %s'%(updated_aln,sequence_file))
        AlignIO.convert(sequence_file,"stockholm",aln_fasta,"fasta")
        if os.path.exists(aln_fasta+'.gz'):
            subprocess.check_call("rm  %s.gz"%(aln_fasta),shell=True)
        subprocess.check_call("gzip %s"%(aln_fasta),shell=True)

        
        
        cmd='mv %s %s'%(tree_output_file,tree_file)
        os.system(cmd)
        if os.path.exists(tree_file+'.gz'):
            subprocess.check_call("rm  %s.gz"%(tree_file),shell=True)
        
        subprocess.check_call("gzip %s"%(tree_file),shell=True)
        with open(pplace_log, "a") as myfile:
            myfile.write("%s\t%s\n"%(my_sequence_file,','.join(p_ids)))
                


def process_tree(pfam_id,db):
    align_file = alignment_folder+'/%s'%pfam_id
    reconciled_fname = reconciled_folder+'/%s'%pfam_id            
    if (os.path.isfile(align_file+".fasta.gz")): #if you have an alignment file...
        if not(os.path.isfile(reconciled_fname+"_reconciled.xml.gz")) or (pplacer_flag==1):   
            print "Process Tree for", pfam_id, "with", tree_sizes[pfam_id],"leaves"
            # make a gene tree
            if not(os.path.isfile("%s/%s.nw.gz"%(tree_folder,pfam_id))): #make a tree based on the alignment
                print "Running FastTree for: %s"%pfam_id
                rand_id_1=random.randint(1000000,9999999)        
                subprocess.check_call("gunzip -c %s.fasta.gz > %s.fasta.%d"%(align_file,align_file,rand_id_1),shell=True)
                subprocess.check_call(lib_path + "/FastTree -log %s/%s.log %s/%s.fasta.%d > %s/%s.nw"%(tree_folder,pfam_id,alignment_folder,pfam_id,rand_id_1,tree_folder,pfam_id),shell=True);
                if os.path.exists("%s.fasta.%d"%(align_file,rand_id_1)):
                    subprocess.check_call("rm  %s.fasta.%d"%(align_file,rand_id_1),shell=True)                
                if os.path.exists('%s/%s.log.gz'%(tree_folder,pfam_id)):
                    subprocess.check_call("rm  %s/%s.log.gz"%(tree_folder,pfam_id),shell=True)
                if os.path.exists('%s/%s.nw.gz'%(tree_folder,pfam_id)):
                    subprocess.check_call("rm  %s/%s.nw.gz"%(tree_folder,pfam_id),shell=True)
                subprocess.check_call("gzip %s/%s.log"%(tree_folder,pfam_id),shell=True)
                subprocess.check_call("gzip %s/%s.nw"%(tree_folder,pfam_id),shell=True)
            if pplacer_flag:
                pplaced_file = add_pplaced(pfam_id)
            gene_tree_file = tree_folder+ '/%s.nw'%(pfam_id)
            reconciled_file='%s/%s_reconciled.xml'%(reconciled_folder,pfam_id)    
            rec_tag='%s/%s_'%(reconciled_folder,pfam_id)
            if not(os.path.isfile(reconciled_file+'.gz')) or (pplacer_flag==1):  
                reconcile_tree(gene_tree_file, reconciled_file, rec_tag, pfam_id,db)
            return pfam_id



def get_goa_annotations_for_pfam_acc_new(domain_acc, evidence_constraints_all,db):
    sql="""SELECT
      goa_db.term_id as acc,
      goa_db.evidence_code as code,
      pfamseq.pfamseq_acc,
      pfamseq.pfamseq_id
     FROM   goa_db
      INNER JOIN pfamseq on (goa_db.uniprot_id=pfamseq.pfamseq_acc)
      INNER JOIN pfamA_reg_full_significant on (pfamA_reg_full_significant.auto_pfamseq=pfamseq.auto_pfamseq)
      INNER JOIN pfamA on (pfamA.auto_pfamA=pfamA_reg_full_significant.auto_pfamA)
     WHERE
       pfamA.pfamA_acc='%s'
       AND goa_db.term_type = 'F'
       AND goa_db.evidence_code in ('%s')
    """%(domain_acc, "','".join(evidence_constraints_all))
    seq_anns = msql(sql, db)
    return seq_anns


def parse_goa_pfam_annots(seq_anns):
    anns = {}
    seq_lookup = {}
    for a in seq_anns:
        pid = a['pfamseq_acc']
        seq_lookup[a['pfamseq_acc']] = a['pfamseq_id']
        if pid not in anns:
            anns[pid] = []
        anns[pid].append(a)
    return anns, seq_lookup

def write_goa_anns_to_pli(evidence_file,goa_anns, fam_id, seq_lookup):
    '''
    This converts the database rows to B. Engelhardt's arbitrary evidence XML foramt.
    Input looks like:
            {'A2VE79': [{'acc': 'GO:0000287',
             'code': 'ISS',
             'full_name': 'Diphosphoinositol polyphosphate phosphohydrolase 1',
             'genus': 'Bos',
             'is_not': 0L,
             'name': 'magnesium ion binding',
             'species': 'taurus',
             'symbol': 'NUDT3',
             'xref_dbname': 'UniProtKB',
             'xref_key': 'A2VE79'},
            {'acc': 'GO:0008486',
             'code': 'ISS',
             'full_name': 'Diphosphoinositol polyphosphate phosphohydrolase 1',
             'genus': 'Bos',
             'is_not': 0L,
             'name': 'diphosphoinositol-polyphosphate diphosphatase activity',
             'species': 'taurus',
             'symbol': 'NUDT3',
             'xref_dbname': 'UniProtKB',
             'xref_key': 'A2VE79'},
            ...
    '''
    f = open(evidence_file, 'w')
    f.write("<?xml version=\"1.0\"?>\n<Family>\n")
    f.write("  <FamilyID>%s</FamilyID>\n"%fam_id)
    
    for p_id, anns in goa_anns.iteritems():
        f.write("  <Protein>\n")
        f.write("    <ProteinName>%s</ProteinName>\n"%seq_lookup[p_id])
        f.write("    <ProteinNumber>%s</ProteinNumber>\n"%p_id)
        go_str = ''
        moc_str = ''
        for i,a in enumerate(anns):
            go_str += a['acc'][3:]
            moc_str += a['code']
            if i < len(anns)-1:
                go_str += ', '
                moc_str += ', '
        f.write("    <GONumber>%s</GONumber>\n"%('['+go_str+']'))
        f.write("    <MOC>%s</MOC>\n"%('['+moc_str+']'))
        f.write("  </Protein>\n")
    f.write("</Family>\n")
    f.close()



# Get all evidence and then write each to file.
def get_evidence(p,my_db):    
    evidence_file = evidence_folder+'/%s.pli'%p
    if not(os.path.isfile(evidence_file+'.gz')):
        evidence_pickle_file = evidence_folder+'/%s.pickle'%p
        print "Retrieving goa annotations for %s..."%p
        seq_anns = get_goa_annotations_for_pfam_acc_new(domain_acc=p, evidence_constraints_all=evidence_constraints_all,db=my_db)
        anns, seq_lookup = parse_goa_pfam_annots(seq_anns=seq_anns)
        print "got %i results."%len(anns)
 
        id_information = pickle.load(open(reconciled_folder+'/%s_ids.pickle'%p, 'rb'))      
        miss_flag=0
        del_anns=[]
        for gene,ev in anns.iteritems():
            if ev[0]['pfamseq_id'] not in id_information['existing_genes']:
                miss_flag=1
                del_anns.append(gene)
                continue
        if miss_flag==1:
            for gene in  del_anns:
               anns.pop(gene) 
        write_goa_anns_to_pli(evidence_file,anns, p, seq_lookup)
        pickle.dump([evidence_file,anns, p, seq_lookup], open(evidence_pickle_file, 'wb'))
        
        if os.path.exists(evidence_file):
            if os.path.exists(evidence_file+'.gz'):
                subprocess.check_call("rm  %s.gz"%(evidence_file),shell=True)
            subprocess.check_call("gzip %s"%(evidence_file),shell=True)
        if os.path.exists(evidence_pickle_file):
            if os.path.exists(evidence_pickle_file+'.gz'):
                subprocess.check_call("rm  %s.gz"%(evidence_pickle_file),shell=True)
            subprocess.check_call("gzip %s"%(evidence_pickle_file),shell=True)
        print "Wrote evidence to %s"%evidence_file


def gather_for_each_family(pfam_id,my_db):
    align_file = alignment_folder+'/%s'%pfam_id
    reconciled_fname = reconciled_folder+'/%s'%pfam_id            
    evidence_file = evidence_folder+'/%s.pli'%pfam_id
    evidence_pickle_file = evidence_folder+'/%s.pickle'%pfam_id              
    queries_to_process=[]
    if not(os.path.isfile(align_file+".fasta.gz")):          
        get_alignment(pfam_id,my_db)
    
    if (os.path.isfile(align_file+".fasta.gz")): #if you have an alignment file...
        if not(os.path.isfile(reconciled_fname+"_reconciled.xml.gz")) or (pplacer_flag==1):   
            process_tree(pfam_id,my_db)

    if not(os.path.isfile(evidence_pickle_file+'.gz')):
        get_evidence(pfam_id,my_db)



class ProcessingThread_gather(threading.Thread):
    """Thread for running sequence alignments on a given input homolog cluster."""
    def __init__(self, thread_queue,db):
        threading.Thread.__init__(self)
        self.thread_queue = thread_queue
        self.db=db
    def thread_operation(self, thread_data):
        pfam_id = thread_data
        my_db=self.db

        try:
            print "--------------------------------------------------"
            print "Gathering family data for %s"%pfam_id
            
            # Input evidence
            gather_for_each_family(pfam_id,my_db)

            print "Family data files gathered for %s"%pfam_id
            print "---------------"
        
        except Exception as e:
            print >> sys.stderr, "Error gathering family data for %s"%pfam_id
            print >> sys.stderr, "Error: ", e
            exit(1)
        
    
    def run(self):
        while True:
            # Spawn a thread with data from the queue
            thread_data = self.thread_queue.get()
            # Run thread's function on the data
            try:
                self.thread_operation(thread_data)
            except:
                print "Unexpected thread error:", sys.exc_info()[0]
                print "Thread data:", thread_data
            # Send signal that this task finished
            self.thread_queue.task_done()


if __name__=="__main__":
    
    # Initialization
    params_mysql = {\
    'db_address': 'localhost',
    'db_username': 'root',
    'db_password': '',
    'db_name': 'sifter_db'
    }

      
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
    species_tree_data_path = main_dir+'/data/species_tree_data'
    lib_path=main_dir+'/lib'

    path_to_hmmbuild='hmmbuild'
    path_to_hmmpress='hmmpress'
    path_to_hmmalign='hmmalign'
    path_to_taxit='taxit'
    path_to_pplacer=lib_path+'/pplacer-v1.1/pplacer'
    path_to_guppy=lib_path+'/pplacer-v1.1/guppy'


    #???
    best_taxid_map_file=species_tree_data_path+'/best_taxid_map.pickle'
    #???
    merged_taxid_file = species_tree_data_path+'/merged.dmp'

    num_threads=4
    pplacer_flag = 0
    my_sequence_file=''
    hit_file=''
    my_taxid=''
    all_fams=0
    pfams=[]
    input_file=''
    # Check for options
    opts, args = getopt.getopt(sys.argv[1:], "hi:f:n:A",['dbname=','dbpass=','dbuser=','dbaddr=','seq_file=','hit_file=','taxid=']) 
    if len(args) != 1:
        usage()
        sys.exit()
    
    choices=[]
    new_genome_choices=[]
    if len(opts)>0:
        for o, a in opts:
            if o == "-f":
                splited =a.strip().split(',')
                pfams=list(set([w for w in splited if w]))  
                choices.append('f')
            elif o == "-i":
                input_file = a
                choices.append('i')
            elif o == "-A":
                all_fams = 1
                choices.append('A')
            elif o == "--seq_file":
                my_sequence_file = a
                new_genome_choices.append('seq_file')
            elif o == "--hit_file":
                hit_file = a
                new_genome_choices.append('hit_file')
            elif o == "--taxid":
                my_taxid = a
                new_genome_choices.append('taxid')
            elif o == "-n":
                num_threads=int(a)
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
   
    if len(new_genome_choices)>0 and len(new_genome_choices)<3:
        print "\nERROR: To gather data for a new genome, please enter information for all '--seq_file', '--hit_file', and '--taxid' options.\n"
        sys.exit()

    if my_sequence_file:
        if not os.path.exists(my_sequence_file):
            print "\nERROR: No sequence file at %s.\n"%my_sequence_file
            sys.exit()

    if hit_file:
        if not os.path.exists(hit_file):
            print "\nERROR: No Pfam hit file at %s.\n"%hit_file
            sys.exit()

    if len(choices)==0:
        print "\nERROR: No pfam families are entered."
        print "Please use one of the -f or -i or -A options to enter your query.\n"
        sys.exit()
    elif len(choices)>1:
        print "\nERROR: Please use ONLY one of the -f or -i or -A options to enter your query.\n"
        sys.exit()
    
    if (len(new_genome_choices)==0) and ('A' in choices):
        print "\nERROR: Option -A can only be used when gather data for a novel species.\n"
        sys.exit()

    output_families_data_path=args[0]
    if not os.path.exists(output_families_data_path):
        os.mkdir(output_families_data_path) 
        

    evidence_folder=output_families_data_path+'/annotations'
    if not os.path.exists(evidence_folder):
        os.mkdir(evidence_folder)
    alignment_folder=output_families_data_path+'/alignments'
    if not os.path.exists(alignment_folder):
        os.mkdir(alignment_folder)
    reconciled_folder=output_families_data_path+'/reconciled_trees'
    if not os.path.exists(reconciled_folder):
        os.mkdir(reconciled_folder)
    tree_folder=reconciled_folder+'/trees'    
    if not os.path.exists(tree_folder):
        os.mkdir(tree_folder)

    db_mysql = MySQLdb.connect(host=params_mysql['db_address'],
                   user=params_mysql['db_username'],
                   passwd=params_mysql['db_password'],
                   db=params_mysql['db_name'],
                   cursorclass=MySQLdb.cursors.DictCursor)


    print "\n\n--------------Reading the input famiy information------------"
    if pfams:
        print "Run SIFTER for %s Pfam Families"%len(pfams)
    elif input_file:
        if not os.path.exists(input_file):
            print "\nERROR: No file exists at %s\n"%input_file
            sys.exit()
        f = open(input_file, 'r')
        a=f.read()
        splited =re.split(' |,|;|\n',a.strip())
        pfams=list(set([w for w in splited if w]))
        print "Run SIFTER for %s Pfam Families"%len(pfams)

    if (not pfams) and (len(new_genome_choices)==0) :
        print "\nERROR: No pfam families are entered.\n"
        sys.exit()

    if len(new_genome_choices)==3:
        print "\n\n--------------Reading the Pfam hit file------------"
        gene_seq,pplacer_queries=find_pfam_hits_from_file(hit_file,my_sequence_file)
        if not pfams:
            pfams=pplacer_queries.keys()
            print "We will run on All %s Pfam Families of your novel species."%len(pfams)

    pfds=get_pfds(pfams,db_mysql)

    tree_sizes = {}
    for p in pfds.keys():
        tree_sizes[p] = pfds[p]['num_full']
    sorted_fams = sorted(pfds.keys(), key=lambda k:pfds[k]['num_full'])
    print "Number of families to process:" ,len(sorted_fams)


    print "\n\n--------------Reading the species tree data------------"
    # ##3-Extract tree information for each Pfam family
    all_species_txids_pickled=species_tree_data_path + '/all_species_txids.pickle'
    all_species_txids=pickle.load(open(all_species_txids_pickled))
    orig_sp_tree_0 = ete2.PhyloTree(species_tree_data_path + '/ncbi.nw', format=0)
        
        
    zero_taxids={'HUMAN':9606, '9CAUD':70702, 'BABHY':37552, 'ARATH':3702, '9STAP':1077965, 'SALTM':99287,  '9MYCO':512402, '9RETR':31697,
                 'BEABA':176275, '9EURO':1194637, '9BACE':871324, '9CAEN':1181719 }

    best_taxid_map=pickle.load(open(best_taxid_map_file,'r'))
    #???
    merged_f=open(merged_taxid_file,'r')
    merged_taxid={}
    for line in merged_f.readlines():
        my_line=line.split('\t')
        merged_taxid[int(my_line[0])]=int(my_line[2])
    merged_f.close()

    print "\n------------Gather the necessary data for families-------------"

    my_taxid0=my_taxid
    if my_taxid0:
        success=1
        if my_taxid0 not in all_species_txids:
            success=0
            if my_taxid0 in merged_taxid.keys():
                my_taxid0=merged_taxid[my_taxid0]
                tax_name=find_tax_name(my_taxid0,db_mysql)
            if my_taxid0 in best_taxid_map.keys():
                my_taxid0=best_taxid_map[my_taxid0]
                tax_name=find_tax_name(my_taxid0,db_mysql)
                success=1
            else:
                tax_id0=my_taxid0
                my_taxid0,tax_name=find_best_taxid(my_taxid0,db_mysql)
                if my_taxid0>0:
                    best_taxid_map[tax_id0]=my_taxid0
                    success=1
                    
        if success==0:
            print "\nThe taxonomy ID you entered does not exist in our database, please enter the correct NCBI taxonomy ID for your species. You may also enter the NCBI taxonomy ID for a close species to your query that exist in our dataset.%s\n"
            sys.exit()
        else:
            pplacer_flag=1
            pplacer_folder=reconciled_folder+'/pplaced'
            if not os.path.exists(pplacer_folder):
                os.mkdir(pplacer_folder)
            
                
    

    pfams_to_process = []
    for i,pfam_id in enumerate(sorted_fams):
        if pplacer_flag==1:
            if pfam_id not in pplacer_queries:
                print "Your queried species does not have a gene in family %s)"%(pfam_id)
                print "We will Skip this family"
                print "---------------"
                continue
        pfams_to_process.append(pfam_id)
 
 
    thread_queue = Queue.Queue()
    for i in range(num_threads):
        my_db=MySQLdb.connect(host=params_mysql['db_address'],
                   user=params_mysql['db_username'],
                   passwd=params_mysql['db_password'],
                   db=params_mysql['db_name'],
                   cursorclass=MySQLdb.cursors.DictCursor)
        t = ProcessingThread_gather(thread_queue,my_db)
        t.setDaemon(True)
        t.start()

     
    for pfam_id in pfams_to_process:
        thread_queue.put(item=pfam_id, block=False)
    
    # Wait on the queue until everything has been processed         
    thread_queue.join()


    errors=0
    for pfam_id in pfams_to_process:
        evidence_pickle_file = output_families_data_path+'/annotations/%s.pickle'%pfam_id
        if not (os.path.isfile(evidence_pickle_file+'.gz')):
            errors+=1



    if errors==0:  
        print "-------------------Data gadering is Done----------------------"
        print "\nNext step is to run 'sifter_prepare.py' to prepares necessary files for your query to run SIFTER on."
    else:
        print "\nData files are gatherd for %d out of %d families. (%s missed due to errors)"%(len(pfams_to_process)-errors,len(pfams_to_process),errors)
        
  

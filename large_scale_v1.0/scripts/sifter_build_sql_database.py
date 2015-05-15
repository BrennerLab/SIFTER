#!/usr/bin/python
import os
import getopt
import _mysql as mysql
import _mysql_exceptions as mysql_exceptions
import MySQLdb.cursors
import sys
import csv
import subprocess
from Bio import AlignIO,Phylo,SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo import PhyloXMLIO,PhyloXML
import ete2
from ete2 import Phyloxml
import pickle

def usage():
    print "\n-----------------------------------------------------------------"
    print "Usage: "
    print "   sifter_build_SQLdb.py [options] temp_output_dir"
    print "-----------------------------------------------------------------\n"
    print "Examples:"
    print "   sifter_build_sql_database.py ../example/tmp_output\n"
    print "   sifter_build_sql_database.py --goa --ont ../example/tmp_output\n"
    print "   sifter_build_sql_database.py --id ../example/tmp_output\n"
    print "   sifter_build_sql_database.py --dbaddr www.example.org --dbuser jack --dbpass 1234  ../example/tmp_output\n"    
    print "This function build the MySQL database that is needed to run all other scripts."
    print "NOTE: Only use this script if you don't wish to use the precomputed SQL database built based on latest releases of Pfam and GOA databases."
    print "@author Sayed Mohammad Ebrahim Sahraeian (mohammad@compbio.berkeley.edu)"
    print "Please cite new paper:"
    print "-Sahraeian SME, Luo KR, Brenner SE (2015)"
    print "The SIFTER algorithm presented in the following paper:"
    print "- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. Genome-scale phylogenetic function annotation of large and diverse protein families. Genome Research 21:1969-1980. \n"
    print "inputs:"
    print "        <temp_output_dir>        Path to the folder where the"
    print "                                 temporary files will be written." 
    print "options: (you should only use one of -i or -f  or -A options.)"
    print "           --pfam                Only update the Pfam tables."
    print "           --goa                 Only update the Gene ontology"
    print "                                 annotation (GOA) tables."
    print "           --id                  Only update the 'uniprot ID" 
    print "                                 map to NCBI' table."
    print "           --ont                 Only update the go_term.sqlite file."
    print "           --sp                  Only update the species tree files."
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


def find_tax_name(ncbi_taxid):
    sql_q="""select species from ncbi_taxonomy where ncbi_taxid='%s' limit 1;
    """%(ncbi_taxid)
    my_id = msql(sql_q, db_mysql)
    if my_id:
        my_id=my_id[0]['species']
    return my_id

def prepare_species_tree(FILE_TREE_IN,FILE_TREE_OUT):
    clan_taxa = {}
    treexml = PhyloXMLIO.read(open(FILE_TREE_IN, 'r'))
    tree = treexml[0]
    treexml.attributes.pop('schemaLocation', None)  # not supported by Forester
    tree.rooted = True
    leaf_dict = {}
    for node in tree.clade.find_clades():
        if node.name:
            tax_id = node.name
            if tax_id.startswith('INT'):
                tax_id = tax_id[3:]
            taxon = PhyloXML.Taxonomy(id=PhyloXML.Id(tax_id, provider='ncbi_taxonomy'))
            try:
                taxon.scientific_name = find_tax_name(tax_id)
            except KeyError:
                taxon.scientific_name = '(NA)'
            node._set_taxonomy(taxon)
            node.name = None
        else:
            pass
    PhyloXMLIO.write(treexml, FILE_TREE_OUT)
        
if __name__=="__main__":
    
    # Initialization
    params_mysql = {\
    'db_address': 'localhost',
    'db_username': 'root',
    'db_password': '',
    'db_name': 'sifter_db'
    }

    # Check for options
    opts, args = getopt.getopt(sys.argv[1:], "h",['dbname=','dbpass=','dbuser=','dbaddr=','pfam','goa','id','ont','sp']) 
    if len(args) != 1:
        usage()
        sys.exit()
    
    choices=[]
    if len(opts)>0:
        for o, a in opts:
            if o == "--pfam":
                choices.append('pfam')
            elif o == "--goa":
                choices.append('goa')
            elif o == "--id":
                choices.append('id')
            elif o == "--ont":
                choices.append('ont')
            elif o == "--sp":
                choices.append('sp')
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
   
    output_dir=args[0]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir) 
        
    main_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    lib_path=main_dir+'/lib'
    data_path=main_dir+'/data'
    species_tree_data_folder = main_dir+'/data/species_tree_data'
   
    if not choices:
        choices=['pfam','goa','id','ont','sp']
    
    if ('pfam' in choices) or ('goa' in choices) or 'id' in choices:
        print "\n\n--------------Creating the database------------"
        db_mysql = MySQLdb.connect(host=params_mysql['db_address'],
                           user=params_mysql['db_username'],
                           passwd=params_mysql['db_password'])
        cursor = db_mysql.cursor()
        sql = 'CREATE DATABASE IF NOT EXISTS %s'%params_mysql['db_name']
        cursor.execute(sql)


        db_mysql = MySQLdb.connect(host=params_mysql['db_address'],
                       user=params_mysql['db_username'],
                       passwd=params_mysql['db_password'],
                       db=params_mysql['db_name'],
                       cursorclass=MySQLdb.cursors.DictCursor)

        if ('pfam' in choices):        
            print "\n\n----Download Pfam SQL tables and put into the database-----"
            pfam_ftp_path='ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/'
            pfam_tables=['pfamA','gene_ontology','ncbi_taxonomy','pfamseq','pfamA_reg_full_significant','alignments_and_trees'];#['pfamA','alignments_and_trees']
            for table in pfam_tables:
                if not os.path.exists("%s/%s.processed"%(output_dir,table)):
                    msql("drop table if exists %s"%table,db_mysql)    
                    if os.path.exists("%s/%s.sql"%(output_dir,table)):
                        print "Removing previously downloaded files for: " + table
                        cmd='rm %s/%s.sql'%(output_dir,table)
                        os.system(cmd)
                    if os.path.exists("%s/%s.sql.gz"%(output_dir,table)):
                        print "Removing previously downloaded files for: " + table
                        cmd='rm %s/%s.sql.gz'%(output_dir,table)
                        os.system(cmd)
                    if os.path.exists("%s/%s.txt"%(output_dir,table)):
                        print "Removing previously downloaded files for: " + table
                        cmd='rm %s/%s.txt'%(output_dir,table)
                        os.system(cmd)
                    if os.path.exists("%s/%s.txt.gz"%(output_dir,table)):
                        print "Removing previously downloaded files for: " + table
                        cmd='rm %s/%s.txt.gz'%(output_dir,table)
                        os.system(cmd)
                        
                    try:
                        print "Downloading schema for: " + table
                        cmd='wget %s/%s.sql.gz --directory-prefix %s'%(pfam_ftp_path,table,output_dir)
                        retcode=subprocess.call(cmd,shell=True)
                        if not retcode == 0:
                            raise Exception(cmd)
                        print "Unzipping schema for: " + table
                        cmd='gunzip %s/%s.sql.gz'%(output_dir,table)
                        retcode=subprocess.call(cmd,shell=True)
                        if not retcode == 0:
                            raise Exception(cmd)
                        print "Loading schema for: " + table
                        cmd = "mysql -u " + params_mysql['db_username'] \
                                + " --password=" + params_mysql['db_password'] \
                                + " --database=" + params_mysql['db_name'] \
                                + " < "+"%s/%s.sql"%(output_dir,table)
                        retcode=subprocess.call(cmd,shell=True)
                        if not retcode == 0:
                            raise Exception(cmd)

                        print "Downloading data for: " + table
                        cmd='wget %s/%s.txt.gz --directory-prefix %s'%(pfam_ftp_path,table,output_dir)
                        retcode=subprocess.call(cmd,shell=True)
                        if not retcode == 0:
                            raise Exception(cmd)
                        print "Unzipping data for: " + table
                        cmd='gunzip %s/%s.txt.gz'%(output_dir,table)
                        retcode=subprocess.call(cmd,shell=True)
                        if not retcode == 0:
                            raise Exception(cmd)
                        print "Loading data for: " + table
                        cmd = "mysql --local-infile=1 -u " + params_mysql['db_username'] \
                                + " --password=" + params_mysql['db_password'] \
                                + " --database=" + params_mysql['db_name'] \
                                + ' -e "LOAD DATA LOCAL INFILE \''+"%s/%s.txt"%(output_dir,table)+'\' into table '+table+'"' 
                        retcode=subprocess.call(cmd,shell=True)
                        if not retcode == 0:
                            raise Exception(cmd)
                    except Exception as e:
                        print >> sys.stderr, "Error: ", e
                        exit(1)
                        

                    print "Removing downloaded files for: " + table
                    cmd='rm %s/%s.sql'%(output_dir,table)
                    os.system(cmd)
                    cmd='rm %s/%s.txt'%(output_dir,table)
                    os.system(cmd)

                    print "Wrting the .processed file for: " + table
                    f = open("%s/%s.processed"%(output_dir,table), 'w')
                    f.close() 
                else:
                    print table+"%s.processed file exist --> Already processed:"+table  
           
        if ('goa' in choices):        
            print "\n\n----Download GOA data and put into the database-----"
            goa_ftp_path="ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/"
            table="goa_db"
            if not os.path.exists("%s/goa_db.processed"%(output_dir)):
                msql("drop table if exists goa_db ",db_mysql)    
                msql("create table goa_db (uniprot_id VARCHAR(10) not null, symbol VARCHAR(10) not null, term_id VARCHAR(10) not null, evidence_code VARCHAR(6) not null, term_type CHAR(1) not null, tax_id int(10) unsigned not null,   key (uniprot_id,term_id,evidence_code,term_type,tax_id));",db_mysql)    
                if os.path.exists("%s/gene_association.goa_uniprot.filtered"%(output_dir)):
                    print "Removing previously downloaded files for: " + table
                    cmd='rm %s/gene_association.goa_uniprot.filtered'%(output_dir)
                    os.system(cmd)
                if os.path.exists("%s/gene_association.goa_uniprot"%(output_dir)):
                    print "Removing previously downloaded files for: " + table
                    cmd='rm %s/gene_association.goa_uniprot'%(output_dir)
                    os.system(cmd)
                if os.path.exists("%s/gene_association.goa_uniprot.gz"%(output_dir)):
                    print "Removing previously downloaded files for: " + table
                    cmd='rm %s/gene_association.goa_uniprot.gz'%(output_dir)
                    os.system(cmd)
                            
                try:
                    print "Downloading data for: " + table
                    cmd='wget %s/gene_association.goa_uniprot.gz --directory-prefix %s'%(goa_ftp_path,output_dir)
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)
                    print "Unzipping data for: " + table    
                    cmd='gunzip %s/gene_association.goa_uniprot.gz'%(output_dir)
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)
                    print "Filtering out non-necessary columns for: " + table        
                    with open('%s/gene_association.goa_uniprot'%output_dir, 'rb') as csvfile_r:
                        spamreader = csv.reader(csvfile_r, delimiter='\t', quotechar='|')
                        cnt=0
                        with open('%s/gene_association.goa_uniprot.filtered'%output_dir, 'wb') as csvfile_w:
                            spamwriter = csv.writer(csvfile_w, delimiter='\t',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                            for row in spamreader:
                                if row[0][0]=='!':
                                    continue
                                cnt+=1
                                if 'taxon:' in row[12]:
                                    new_row=[row[1],row[2].split('"')[0],row[4],row[6],row[8],str(int((row[12].split('|')[0]).split(':')[1]))] 
                                if 'taxon:' in row[11]:
                                    new_row=[row[1],row[2].split('"')[0],row[4],row[6],row[8],str(int((row[11].split('|')[0]).split(':')[1]))]                 
                                spamwriter.writerow(new_row)
                    
                    print "Removing gene_association.goa_uniprot" 
                    cmd='rm %s/gene_association.goa_uniprot'%(output_dir)
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)

                    print "Loading data for: " + table        
                    cmd = "mysql --local-infile=1 -u " + params_mysql['db_username'] \
                            + " --password=" + params_mysql['db_password'] \
                            + " --database=" + params_mysql['db_name'] \
                            + ' -e "LOAD DATA LOCAL INFILE \''+"%s/gene_association.goa_uniprot.filtered"%(output_dir)+'\' into table goa_db"' 
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)
                except Exception as e:
                    print >> sys.stderr, "Error: ", e
                    exit(1)
                    

                print "Removing downloaded files for: " + table
                cmd='rm %s/gene_association.goa_uniprot.filtered'%(output_dir)
                os.system(cmd)
                
                print "Wrting the .processed file for: " + table
                f = open("%s/goa_db.processed"%(output_dir), 'w')
                f.close()    
            else:
                print "goa_db.processed file exist --> Already processed:"+table  

        if ('id' in choices):        
            print "\n\n----Download Uniprot 2 NCBI ID  mapping data and put into the database-----"
            idmap_ftp_path="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/"
            table="uniprot_2_NCBI"
            if not os.path.exists("%s/uniprot_2_NCBI.processed"%(output_dir)):
                msql("drop table if exists uniprot_2_NCBI",db_mysql)    
                msql("create table uniprot_2_NCBI (uniprot_id VARCHAR(10) not null, ncbi_id VARCHAR(15) not null, primary key (uniprot_id));",db_mysql)
                
                if os.path.exists("%s/idmapping.dat"%(output_dir)):
                    print "Removing previously downloaded files for: " + table
                    cmd='rm %s/idmapping.dat'%(output_dir)
                    os.system(cmd)
                if os.path.exists("%s/idmapping.dat.gz"%(output_dir)):
                    print "Removing previously downloaded files for: " + table
                    cmd='rm %s/idmapping.dat.gz'%(output_dir)
                    os.system(cmd)
                if os.path.exists("%s/idmapping_NCBI_filtered.tab"%(output_dir)):
                    print "Removing previously downloaded files for: " + table
                    cmd='rm %s/idmapping_NCBI_filtered.tab'%(output_dir)
                    os.system(cmd)
                
                try:
                    print "Downloading data for: " + table
                    cmd='wget %s/idmapping.dat.gz --directory-prefix %s'%(idmap_ftp_path,output_dir)
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)
                    
                    print "Unzipping data for: " + table    
                    cmd='gunzip %s/idmapping.dat.gz'%(output_dir)
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)
                    
                    print "Filtering out non-necessary columns for: " + table        
                    with open('%s/idmapping.dat'%output_dir, 'rb') as csvfile_r:
                        spamreader = csv.reader(csvfile_r, delimiter='\t', quotechar='|')
                        id_map={}
                        cnt=0
                        with open('%s/idmapping_NCBI_filtered.tab'%output_dir, 'wb') as csvfile_w:
                            spamwriter = csv.writer(csvfile_w, delimiter='\t',quotechar='|', quoting=csv.QUOTE_MINIMAL)
                            for row in spamreader:
                                cnt+=1
                                if row[1]=='NCBI_TaxID':
                                    spamwriter.writerow([row[0],row[2]])


                    print "Removing idmapping.dat"
                    cmd='rm %s/idmapping.dat'%(output_dir)
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)


                    cmd = "mysql --local-infile=1 -u " + params_mysql['db_username'] \
                            + " --password=" + params_mysql['db_password'] \
                            + " --database=" + params_mysql['db_name'] \
                            + ' -e "LOAD DATA LOCAL INFILE \''+"%s/idmapping_NCBI_filtered.tab"%(output_dir)+'\' into table uniprot_2_NCBI"' 
                    
                    retcode=subprocess.call(cmd,shell=True)
                    if not retcode == 0:
                        raise Exception(cmd)
                except Exception as e:
                    print >> sys.stderr, "Error: ", e
                    exit(1)
                    

                print "Removing downloaded files for: " + table
                cmd='rm %s/idmapping_NCBI_filtered.tab'%(output_dir)
                os.system(cmd)
                
                print "Wrting the .processed file for: " + table
                f = open("%s/uniprot_2_NCBI.processed"%(output_dir), 'w')
                f.close()    
            else:
                print "uniprot_2_NCBI.processed file exist --> Already processed:"+table  


    if ('ont' in choices):        
        print "\n\n----Updata the go_term.sqlite file-----"
        #check here for more info: 
        #https://code.google.com/p/variationtoolkit/wiki/GeneOntologyDbManager
        #https://code.google.com/p/variationtoolkit/wiki/HowToInstall
        if not os.path.exists("%s/goterms_sqlite.processed"%(output_dir)):
            go_rdf_download_path="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz"
            if os.path.exists("%s/go_daily-termdb.rdf-xml.gz"%(output_dir)):
                print "Removing previously downloaded files for: go_daily-termdb.rdf-xml.gz"
                cmd='rm %s/go_daily-termdb.rdf-xml.gz'%(output_dir)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
            try:
                print "Downloading data for go_term.sqlite"
                cmd='wget %s --directory-prefix %s'%(go_rdf_download_path,output_dir)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                
                print "Unzipping data for go_term.sqlite"    
                cmd=' gunzip -c %s/go_daily-termdb.rdf-xml.gz |  %s/variationtoolkit-read-only/bin/godbmgr loadrdf -f %s/goterms.sqlite'%(output_dir,lib_path,data_path)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                
            except Exception as e:
                print >> sys.stderr, "Error: ", e
                exit(1)

            print "Removing downloaded files go_daily-termdb.rdf-xml.gz"
            cmd='rm %s/go_daily-termdb.rdf-xml.gz'%(output_dir)
            os.system(cmd)

            print "Wrting the .processed file for go_term.sqlite"
            f = open("%s/goterms_sqlite.processed"%(output_dir), 'w')
            f.close() 
        else:
            print "Already processed: go_term.sqlite"  
            print "goterms_sqlite.processed file exist --> Already processed: go_term.sqlite"  

    if ('sp' in choices): 
        print "\n\n----Updata the species tree files-----"
        #check here for more info: 
        #https://github.com/jhcepas/ncbi_taxonomy
        
        
        db_mysql = MySQLdb.connect(host=params_mysql['db_address'],
               user=params_mysql['db_username'],
               passwd=params_mysql['db_password'],
               db=params_mysql['db_name'],
               cursorclass=MySQLdb.cursors.DictCursor)


        if not os.path.exists("%s/ncbi_tree.processed"%(output_dir)):
            taxdump_path="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
            if os.path.exists("%s/taxdump.tar.gz"%(output_dir)):
                print "Removing previously downloaded files for ncbi_tree"
                cmd='rm %s/taxdump.tar.gz'%(output_dir)
                os.system(cmd)            
            if os.path.exists("%s/ncbi_taxonomy/taxdump"%(lib_path)):
                print "Removing previously downloaded files for ncbi_tree"
                cmd='rm -rf %s/ncbi_taxonomy/taxdump'%(lib_path)
                os.system(cmd)            
            if os.path.exists("%s/ncbi_taxonomy/tmp"%(lib_path)):
                print "Removing previously downloaded files for ncbi_tree"
                cmd='rm -rf %s/ncbi_taxonomy/tmp'%(lib_path)
                os.system(cmd)            
            try:
                print "Downloading data for ncbi_tree"
                cmd='wget %s --directory-prefix %s'%(taxdump_path,output_dir)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                
                
                print "Unzipping data for ncbi_tree"    
                cmd="mkdir %s/ncbi_taxonomy/taxdump"%(lib_path)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='tar zxf %s/taxdump.tar.gz --directory %s/ncbi_taxonomy/taxdump'%(output_dir,lib_path)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                print "Obtain the latest taxadb"    
                cmd='python %s/ncbi_taxonomy/update_taxadb_modified.py'%(lib_path)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                    
                print "Convet tree to phyloxml"
                tmp_folder='%s/ncbi_taxonomy/tmp'%(lib_path)
                cmd="mkdir %s"%tmp_folder
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='rm  %s/ncbi_taxonomy/taxa.tab'%(lib_path)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='rm  %s/ncbi_taxonomy/syn.tab'%(lib_path)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='rm  %s/ncbi_taxonomy/taxa.sqlite'%(lib_path)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='mv  %s/ncbi_taxonomy/ncbi.nw %s/ncbi.nw'%(lib_path,tmp_folder)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)


                    
                sp_tree = ete2.PhyloTree(tmp_folder+"/ncbi.nw", format=0)
                sp_tree.write(format=0, outfile=tmp_folder+'/ncbi_2.nw')
                    
                os.system("java -Xmx4g -cp %s/forester_1038.jar org.forester.application.phyloxml_converter -f=nn %s/ncbi_2.nw %s/ncbi_2.xml"%(lib_path,tmp_folder, tmp_folder))

                print "Add taxid information"
                prepare_species_tree(FILE_TREE_IN = tmp_folder+'/ncbi_2.xml', FILE_TREE_OUT = tmp_folder+'/ncbi_2_fixed.xml')

                
                print "Pickle all_species_txids"
                sp_tree_org = Phylo.read(tmp_folder + '/ncbi_2_fixed.xml', 'phyloxml')
                a1=[int(w) for w in set([(w.taxonomy.id.value) for w in sp_tree_org.get_nonterminals() if w.taxonomy])-set(['1'])]
                a2=[int(w.taxonomy.id.value) for w in sp_tree_org.get_terminals()]
                all_species_txids=list(a1+a2)
                pickle.dump(all_species_txids,open(tmp_folder+'/all_species_txids.pickle','w'))

                print "Moving generated files to data/species_tree_data "
                cmd='mv  %s/ncbi.nw %s/ncbi.nw'%(tmp_folder,species_tree_data_folder)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='mv  %s/ncbi_2_fixed.xml %s/ncbi_2_fixed.xml'%(tmp_folder,species_tree_data_folder)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='mv  %s/all_species_txids.pickle %s/all_species_txids.pickle'%(tmp_folder,species_tree_data_folder)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)
                cmd='mv  %s/ncbi_taxonomy/taxdump/merged.dmp %s/merged.dmp'%(lib_path,species_tree_data_folder)
                retcode=subprocess.call(cmd,shell=True)
                if not retcode == 0:
                    raise Exception(cmd)

                
            except Exception as e:
                print >> sys.stderr, "Error in command '%s'"%e
                exit(1)

            print "Removing downloaded files taxdump.tar.gz"
            cmd='rm %s/taxdump.tar.gz'%(output_dir)
            os.system(cmd)
            print "Removing directory taxdump"
            cmd='rm -rf %s/ncbi_taxonomy/taxdump'%(lib_path)
            os.system(cmd)

            print "Removing directory tmp"
            cmd='rm -rf %s/ncbi_taxonomy/tmp'%(lib_path)
            os.system(cmd)

            print "Wrting the .processed file for ncbi_tree"
            f = open("%s/ncbi_tree.processed"%(output_dir), 'w')
            f.close() 
        else:
             print "ncbi_tree.processed file exist --> Already processed: ncbi_tree"  


    print "-------------------Data gadering is Done----------------------"
    print "\nNext step is to run 'sifter_prepare.py' to prepares necessary files for your query to run SIFTER on."
   
  

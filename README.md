# SIFTER
A pipeline for large-scale phylogeny-based protein function prediction using SIFTER (Statistical Inference of Function Through Evolutionary Relationships)

========================================================================

SIFTER is a statistical approach to predicting protein function 
that uses a protein family's phylogenetic tree, as the natural structure 
for representing protein relationships, overlaid with all known protein 
functions in the family.

This package provides a pipeline for large scale protein function 
prediction using SIFTER algorithm. Thus, it can be used for genome-wide
protein function prediction.

The large Scale implementation of SIFTER is developed by Sayed Mohammad 
Ebrahim Sahraeian at Department of Plant and Microbial Biology, 
University of California, Berkeley.

Please cite new paper:
       Sahraeian SME, Luo KR, Brenner SE (2015)

The orignial SIFTER algorithm is developed by Barbara E Engelhardt.
Original paper:
- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. 
Genome-scale phylogenetic function annotation of large and diverse 
protein families.Genome Research 21:1969-1980. doi:10.1101/gr.104687.109 

You can also use the SIFTER webserver at http://sifter.berkeley.edu 
to access online the predictions on 16,863,537 proteins across 
232,403 species. More information at:
- Sahraeian SME, Luo KR, Brenner SE. 2015. SIFTER Search: A web server 
for accurate phylogeny-based protein function prediction. Nucleic Acids 
Research, to appear


Other previous developers:
  Philip Johnson, Steven R. Chan, Micheal Souza


========================================================================
##Download Package
========================================================================

    mkdir sifter_large_scale
    cd sifter_large_scale
    git clone https://github.com/BrennerLab/SIFTER.git
    
    Extract lib and data folders:
    
    cd SIFTER/large_scale_v1.0
    tar -xzvf data.tar.gz
    tar -xzvf lib.tar.gz

========================================================================
##Scripts:
========================================================================

    sifter_find_families.py             Finds Pfam families for your 
                                        query protein or species.
    
    sifter_gather_family_data.py        Gathers necessary 'alignment', 
                                        'tree', and 'evidence' files 
                                        needed to run SIFTER for each 
                                        query family.
                                        [NOTE: Only use this script if you 
                                        don't wish to use the precomputed 
                                        data files built based on latest 
                                        releases of Pfam and GOA 
                                        databases, OR if you have a set 
                                        of sequences (of a novel genome) 
                                        which is not already in Pfam]

    sifter_prepare.py                    Prepares necessary files for 
                                        your query to run SIFTER on.
    
    sifter_run.py                        Runs SIFTER on the prepared 
                                        files generated by 
                                        'sifter_prepare.py'.

    sifter_extract.py                    Extracts SIFTER predictions for 
                                        your query proteins, species, or
                                        families. (using the outputs of 
                                        'sifter_run.py'.)

    sifter_build_sql_database.py        Gathers necessary 'alignment', 
                                        'tree', and 'evidence' files 
                                        needed to run SIFTER for each 
                                        query family.
                                        [NOTE: Only use this script if 
                                        you don't wish to use the 
                                        precomputed SQL database built 
                                        based on latest releases of 
                                        Pfam and GOA databases.

========================================================================
##Requirments:
========================================================================
1-Install necessary Prerequisite packages:

    sudo apt-get update
    sudo apt-get install mysql-server build-essential python-dev pip \
        libmysqlclient-dev liblapack-dev libatlas-dev gfortran \
        default-jdk git 

    Python packages:
        sudo pip install numpy 
        sudo pip install scipy
        sudo pip install Biopython
        sudo pip install ete2
        sudo pip install MySQL-python
        sudo pip install sqlite

TO BE ABLE TO RUN ON NEW GENOMES (which don't exist in Pfam current release):

2-Install PfamScan source code and data
    #Documentation at: https://wiki.gacrc.uga.edu/wiki/Pfam_scan

    wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz
    tar -zxvf PfamScan.tar.gz
    cd PfamScan
    wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
    gunzip Pfam-A.hmm.dat.gz
    wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
    gunzip active_site.dat.gz
    hmmpress Pfam-A.hmm

3-Install Perl packages
    sudo apt-get install libmoose-perl
    sudo cpan BioPerl
    sudo cpan IPC:RUN

4-Install hmmer:
    wget ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-3.1b1.tar.gz
    tar -zxvf hmmer-3.1b1.tar.gz
    cd hmmer-3.1b1
    mkdir build
    ./configure
    make
    sudo make install

5-Install TAXTASTIC
    sudo pip install taxtastic

6-Fix the paths to hmmbuild, hmmpress, hmmalign, and taxit if not already in your PATH.



========================================================================
##Usage
========================================================================
Before running SIFTER you should performs the one time data setup (STEPs
1 and 2 below).

NOTE: To run the scripts first change directory to the scripts folder
		
		cd scripts/large_scale_v1.0/


STEP 1-Build SIFTER MySQL databases:
    To Run the Python scripts in this package you need to have a MySQL 
    database that encompasses the necessary information into the tables 
    (The .gz file is ~34GB).
    
    --You can download use the precomputed databases as follows:

        cd data
        wget sifter.berkeley.edu/media/sifter_db.gz    
        [You may need to set max_allowed_packet=1024M in the [mysqld] 
        section of your MySQL my.cnf and then restart the MySQL]
        mysqladmin -u user -p pass create sifter_db
        gunzip < sifter_db.gz | mysql sifter_db  

    --Alternatively you can build the MySQL from scratch using latest
    available data using the 'sifter_build_sql_database.py' script 
    (this may take 1-2 days)

        python sifter_build_sql_database.py ../example/tmp_output
    

STEP 2- Obtain family data
    To run SIFTER you need to prepare necessary tree and evidence files.
    We have prepared necessary data files based on latest releases of 
    Pfam and Gene Ontology Annotation databases. You can download those 
    files as follows (The .tar.gz file is ~25GB):

        wget sifter.berkeley.edu/media/families_data.tar.gz
        tar -xzvf families_data.tar.gz
    
    Alternatively, you may gather necessary family data using the 
    'sifter_gather_family_data.py' script. You can use this script to
    gather data for ALL Pfam families ((this may take 1-2 days) or ONLY 
    the Pfam families you may need for running your query.
    
    Example:
    -To create family data for two Pfam families PF12491 and PF13820:
        python sifter_gather_family_data.py -f PF12491,PF13820 ../example/fam_data

    -To create family data for families which a spesific gene has domain in:
        python sifter_find_families.py -p C0JYY2_HUMAN ../examples/family_list.txt
        python sifter_gather_family_data.py -i ../example/family_list.txt ../example/fam_data
        
    -To create family for All Pfam families:
        python sifter_find_families.py -A ../examples/family_list.txt
        python sifter_gather_family_data.py -i ../example/family_list.txt ../example/fam_data
        
    -If you already have the precomputed damily data and you want to keep 
    the alignments and trees data (since Pfam database does not change 
    so fast), but update the evidence annotations to most recent version,
    you may use the following:
        mv families_data/annotations families_data/old_annotations
        python sifter_build_sql_database.py --goa --ont --id ../example/tmp_output
        python sifter_find_families.py -A ../examples/family_list.txt
        python sifter_gather_family_data.py -i ../example/family_list.txt path/to/families_data


STEP 3-Run SIFTER for your query
    Given that the mysql (step 1) and family data (step 2) are ready, you can
    start running SIFTER on you queries.

    Example:
    -To run SIFTER on some families:
        python sifter_prepare.py -f PF12491,PF13820 path/to/families_data ../example/queries
        python sifter_run.py ../example/queries ../example/results
        python sifter_extract.py -f PF12491,PF13820 ../example/results ../examples/preds.txt

    -To run SIFTER on a gene:
        python sifter_prepare.py -p C0JYY2_HUMAN path/to/families_data ../example/queries
        python sifter_run.py ../example/queries ../example/results
        python sifter_extract.py -p C0JYY2_HUMAN ../example/results ../examples/preds.txt

    -To run SIFTER on a species (NOTE: This usually will take a long time as the species has 
    domains in many families).
        python sifter_prepare.py -s 9823 path/to/families_data ../example/queries
        python sifter_run.py ../example/queries ../example/results
        python sifter_extract.py -s 9823 ../example/results ../examples/preds.txt


========================================================================
##Run SIFTER a new genome (or set of genes) not already in Pfam database
========================================================================
    If your query genes are not already in the Pfam database, you can not
    use the precomputed families data to build necessary phylogenetic tree
    and evodence files for SIFTER to run.
    
    So, you need to find the families where the query genes have domains in
    and add those genes    to the current trees. The following steps help you
    to do this task: (NOTE:You should have  followd Steps 2-6 of the 
    Installation)
        a) First make sure your query data is not in Pfam. Let say you have
           a list of genes (UniProt ids) in the protein_list.txt. Run:
           python sifter_find_families.py --ip ../examples/protein_list.txt ../examples/family_list.txt
           
           If it returned zero hits in family_list.txt then you have a 
           set of new genes.
           
        b) Find Pfam domains of your query genes (We assume the fasta file 
           of your query genes is at 'myseq.fasta'). 
           
               perl /path/to/PfamScan/pfamscan.pl --fasta ../example/myseq.fasta --dir /path/to/PfamScan/ -e_dom 1 -e_seq 1 -outfile  ../example/pfam_res.txt
               python sifter_gather_family_data.py -A --seq_file ../example/myseq.fasta --hit_file ../example/pfam_res.txt --taxid 1192197 path/to/families_data\n"
               python sifter_prepare.py -A --hit_file ../example/pfam_res.txt path/to/families_data ../example/queries
               python sifter_run.py ../example/queries ../example/results
               python sifter_extract.py -A --hit_file ../example/pfam_res.txt ../example/results ../examples/preds.txt
           
           If you want to run only on specific set of families you can use:
               perl /path/to/PfamScan/pfamscan.pl --fasta ../example/myseq.fasta --dir /path/to/PfamScan/ -e_dom 1 -e_seq 1 -outfile  ../example/pfam_res.txt
               python scripts/sifter_gather_family_data.py -f PF12491,PF13820 --seq_file ../example/myseq.fasta --hit_file ../example/pfam_res.txt --taxid 1192197  path/to/families_data\n"
               python sifter_prepare.py -f PF12491,PF13820 --hit_file ../example/pfam_res.txt path/to/families_data ../example/queries
               python sifter_run.py ../example/queries ../example/results
               python sifter_extract.py -f PF12491,PF13820 --hit_file ../example/pfam_res.txt ../example/results ../examples/preds.txt
        

========================================================================
##Estimating the running time
========================================================================

    Once you run the 'sifter_prepare.py' script to prepare you query, a 
    file ('running_estimation.csv') will be generated in the output folder
    that shows the information for individual families the SIFTER will 
    be run on along with the estimated running time on each family.
    You can also use the following webpage to get an estimate of SIFTER 
    running time based on different family features.
    http://sifter.berkeley.edu/complexity/
    

========================================================================
##Control on Speed
========================================================================
   You can control the speed of running SIFTER using the -x and -t options
   in the 'sifter_prepare.py' script:
   
        -x   INT     Maximum number of nonzero elements
                     in the transition matrix. Should be
                     a number in [1e5,1e7] for reasonable
                     time and accuracy balance (Default=1e6)
                     Smaller value leads to faster runningtime.

        -t   INT     Number of functions to truncate
                     to in approximation [Default: adaptive based
                     on -x option]
                     Smaller value leads to faster runningtime.
   
  
========================================================================
##Multi-threading:
========================================================================
You can run the 'sifter_gather_family_data.py' and 'sifter_run.py' 
scripts on multi-threads using -n option. (Default is 4)
  
========================================================================
##Datasets used:
========================================================================
Used in Version 1.0:
	Pfam 27.0 
	(March 2013, 14831 families)

	GOA 
	ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz 
	(Updated on March 31, 2015)

	Uniprot ID Mapping 
	(ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/)
	(Updated on Jan 7th, 2015)

	Gene Ontology
	http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz
	(Updated on March 31, 2015)

	NCBI Species tree
	ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	(Updated on March 3rd, 2015)

========================================================================
##Version History
========================================================================
5/15/2015	large_scale_v1.0 source code and data released


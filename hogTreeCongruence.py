from __future__ import division
import re
import os
import sys
sys.path.append('/home/gfhtk/workspace/projects/UCLrepos/family-analyzer/familyanalyzer')
import familyanalyzer as fa
import random
import dendropy
import argparse
from Bio import SeqIO
from Bio import AlignIO
import imp
#imp.load_source('fa', '/home/gfhtk/workspace/projects/UCLrepos/family-analyzer/familyanalyzer/familyanalyzer.py')
from dendropy.calculate import treecompare
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import PhymlCommandline

## GET ARGUMENTS FROM COMMAND LINE ##
parser = argparse.ArgumentParser(description='Calculate Gene Tree Congruence on HOGs.')
parser.add_argument('--input', help='orthoxml file', required=True)
parser.add_argument('--file', help='file with list of HOGs')
parser.add_argument('--num', help='number of HOGs to randomly select (default=1)', type=int, default=1)
parser.add_argument('--bootstrap', help='bootstrap number', type=int, default=100)
parser.add_argument('--collapse', help='if specified, collapse threshold', type=float, default=50)
parser.add_argument('--seqs', help='fasta file containing sequences (include complete path)', required=True)
parser.add_argument('--taxlevel', help='taxonomic level of interest', default='root')
parser.add_argument('--family', help='specify specific family')
args = parser.parse_args()

## FAMILY ANALYZER STUFF ##
orthoxml = open(args.input, 'r')
op = fa.OrthoXMLParser(orthoxml)
tax = fa.TaxonomyFactory.newTaxonomy(op)
op.augmentTaxonomyInfo(tax)
tax.get_histories(op) #Generates FamHistory object
tax.get_comparisons(op)
alltax = tax.histories[tax.root] #Finds the root of the taxonomy
hoglist = alltax.geneFamList
taxlevel = args.taxlevel 

## IF USER SPECIFIES A PARTICULAR LIST OF HOGS ##
if args.file:
    hoglist_sample = []
    with open(args.file, 'r') as infile:
        #print('The following list of GeneFamilies (HOGs) will be analyzed:\n'
        #for line in infile:
        famid = infile.readline()
        famid = famid.strip()
        hoglist_sample.append(alltax[famid])
        print(hoglist_sample)

if args.family:
    hoglist_sample = []
    famid = args.family
    hoglist_sample.append(alltax[famid])
    print(hoglist_sample)

## OTHERWISE TAKE 100 RANDOM HOG (OR NUMBER SPECIFIED IN OPTIONS) ##
else:
    random.seed()
    hoglist_sample = random.sample(hoglist,args.num)
    print('The following list of GeneFamilies (HOGs) will be analyzed:\n'+str(hoglist_sample))
    match = re.match(r'.*id#=(.*)\)', str(hoglist_sample[0]))
    if match:
        famid = str(match.group(1))      

## FUNCTIONS ##

def create_hog_tree(hog, taxlevel, hogtree):
   "Infers newick tree from HOG"
   print(hog, str(hog))
   gtt = fa.GeneTreeTracer(op,tax)
   tax_node = tax[tax.root] #taxNode object
   gene_tree = gtt.trace_gene_family(tax_node, hog, None, False)
   newick_long = gene_tree.write()
   with open(FAtree, 'w') as outfile1:
       outfile1.write(str(newick_long))
   #newick_short = re.sub(r'\".*?\" \[.*?\],','', newick_long)
   newick_short = re.sub(r'\[.*?\]|\".*?\"', '',newick_long)
   newick_short =newick_short+';'
   newick_short = re.sub(r', \)', ')', newick_short)
   newick_short = re.sub(r'\(, ', '(', newick_short)
   newick_short = re.sub(r' , ', ' ', newick_short)
   with open(hogtree, 'w') as outfile:
       outfile.write(str(newick_short))
   return (newick_short)

def create_fasta(hog, seqs, fastaoutfile):
    "Creates fasta file for gene tree with IDs from hog"
    omaseqs_dict = SeqIO.index(seqs, 'fasta')
    line = ''
    for geneId in hog.getMemberGenes():
        line += ">"+omaseqs_dict[op.mapGeneToXRef(geneId,"protId")].id+"\n"
        line += str(omaseqs_dict[op.mapGeneToXRef(geneId,"protId")].seq)+"\n"
    with open(fastaoutfile, 'w') as outfile:
        outfile.write(line)

def create_msa(fasta_infile, msa_fasta,msa_phy):
    "Creates a multiple sequence alignment with mafft in phylip format"
    mafft_cline = MafftCommandline(input=fasta_infile) #Create mafft command line
    stdout,stderr = mafft_cline() #save mafft output into variable
    with open(msa_fasta, 'w') as handle:
         handle.write(stdout) #write mafft output in fasta format
    AlignIO.convert(msa_fasta,"fasta", msa_phy, "phylip-relaxed") #convert mafft output from fasta to phylip

def create_gene_tree(msa_phy, bootstrap):
    'Runs PhyML to make gene trees'
    print("Creating ML gene tree...")
    cmdline = PhymlCommandline("phyml", input=msa_phy, datatype='aa', bootstrap=bootstrap)
    out_log, err_log = cmdline()

def collapse(phyml_tree, collapse_threshold, collapsed_tree):
    print("Collapsing poorly supported nodes...")
    cmd_collapse = "~/bin/collapse.py --threshold "+ str(collapse_threshold) + " " + str(phyml_tree)+" |sed \"s/.*/&;/\" > " + collapsed_tree 
    os.system(cmd_collapse)
    
def tree_comparison(collapsed_tree, hogtree, phyml_tree, gene_tree_congruence, report):
    print("Calculating gene tree congruence metric...")
    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(file=open(collapsed_tree, 'r'), schema='newick', taxon_namespace=tns)
    tree2 = dendropy.Tree.get(file=open(hogtree, 'r'), schema='newick', taxon_namespace=tns)
    tree3 = dendropy.Tree.get(file=open(phyml_tree, 'r'), schema='newick', taxon_namespace=tns)
    diff = treecompare.symmetric_difference(tree1,tree2, is_bipartitions_updated=False) #same as unweighted RF distance
    tree1_node_num = (len(tree1.internal_nodes())) #gets number of splits (nodes) in tree
    tree2_node_num = (len(tree2.internal_nodes()))
    tree3_node_num = (len(tree3.internal_nodes()))
    hog_metric = 1 - diff/(tree1_node_num + tree2_node_num)
    num_nodes_collapsed = tree3_node_num - tree1_node_num

    with open(gene_tree_congruence, 'w') as outfile1:
        outfile1.write(str(hog_metric) + "\n")
    with open(report, 'w') as outfile2:
        with open(collapsed_tree, 'r') as collapsedtreefile:
            with open(hogtree, 'r') as hogtreefile:
                with open(phyml_tree, 'r') as phymltreefile:
                    phyml_genetree = phymltreefile.read()
                    hogtree        = hogtreefile.read()
                    collapsed_genetree = collapsedtreefile.read()
                    collapsedplot = tree1.as_ascii_plot()
                    hogplot = tree2.as_ascii_plot()
                    outfile2.write("HOG tree:\n" + str(hogtree) + "\n"+ str(hogplot)+ "\n\nUncollapsed Gene tree:\n" + str(phyml_genetree) + "\n\nCollapsed Gene tree:\n" + str(collapsed_genetree) + "\n"+ str(collapsedplot)+ "\n\n# nodes collapsed: "+str(num_nodes_collapsed)+ "\n\nSymmetric distance: "+str(diff)+ "\n\n# Nodes (Hogtree): "+str(tree2_node_num) + "\n\n# Nodes (Uncollaped Gene tree): "+ str(tree3_node_num)+ "\n\n# Nodes (Collapsed Gene tree): " + str(tree1_node_num)+"\n\nHOG Tree Congruence:" + str(hog_metric)+ "\n\n")

### MAIN ###
for hog in hoglist_sample:    
    os.mkdir(str(famid))
    cwd                   = os.getcwd()+ "/" + str(famid)
    FAtree                = cwd + "/" + str(famid) + '.FAtree' 
    hogtree               = cwd + "/" + str(famid) + '.hogtree'
    fastaoutfile          = cwd + "/" + str(famid) + '.fasta'
    msa_fasta             = cwd + "/" + str(famid) + '.alignment'
    msa_phy               = cwd + "/" + str(famid) + '.phy'
    phyml_tree            = cwd + "/" + str(famid) + '.phy_phyml_tree'
    gene_tree_congruence  = cwd + "/" + str(famid) + '.gtc'
    report                = cwd + "/" + str(famid) + '.txt'
    collapsed_tree        = cwd + "/" + str(famid) + '.phy_phyml_tree.collapsed'
    
    create_hog_tree(hog, taxlevel, hogtree)
    create_fasta(hog, args.seqs, fastaoutfile)
    create_msa(fastaoutfile,msa_fasta,msa_phy)
    genetree= create_gene_tree(msa_phy, args.bootstrap)
    collapse(phyml_tree, args.collapse, collapsed_tree)
    tree_comparison(collapsed_tree, hogtree, phyml_tree, gene_tree_congruence, report)
 
exit()


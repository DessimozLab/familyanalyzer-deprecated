import sys, os
import familyanalyzer as fa

# Checks for right input parameters
if len(sys.argv) < 3:
    print("Use: python compute_species_coverage <OrthoXML> TaxLevel [threshold]")
    sys.exit(0)
    
XML = sys.argv[1]
if not os.path.exists(XML):
    print("Couldn't find OrthoXML file: " + XML)
    sys.exit(0)
if not XML.endswith(".orthoxml"):
    print("OrthoXML has wrong format or file ending (expected: file.orthoxml)! ")
    sys.exit(0)

level = sys.argv[2]

if len(sys.argv) < 4:
    print("WARNING: No threshold given! Using default 0.5!")
    threshold = 0.5
else:
    threshold = float(sys.argv[3])

    
# Parse OrthoXML and compute taxonomic information
op = fa.OrthoXMLParser(XML)
if level not in op.getLevels():
    print("Couldn't find level: " + level + "! Available levels:")
    for l in op.getLevels():
        print(l)
    sys.exit(0)
tax = fa.TaxonomyFactory.newTaxonomy(op)
op.augmentTaxonomyInfo(tax)
all_levels = tax.get_histories(op)
all_branches = tax.get_comparisons(op)

# Compute species coverages for each HOG
HOGs = tax.histories[level]
out = "HOG ID, species coverage\n"
out_below = ""
out_above = ""

for HOG in HOGs:
    if level == "LUCA":
        num_all_species = len(op.getSpeciesSet())
    else:
        num_all_species = len(level.split("/"))
    HOG_species = set()
    for gene in HOG.getMemberGenes():
        HOG_species.add(op.mapGeneToSpecies(gene))
    coverage = len(HOG_species)/num_all_species
    out += HOG.getFamId() + ", " + str(coverage) + "\n"
    if (coverage < threshold):
        out_below += HOG.getFamId() + "\n"
    else:
        out_above += HOG.getFamId() + "\n"

# Write output
with open(XML.replace(".orthoxml", ".coverages.above.txt"), "w") as f:
    f.write(out_above)
with open(XML.replace(".orthoxml", ".coverages.below.txt"), "w") as f:
    f.write(out_below)
with open(XML.replace(".orthoxml", ".coverages.csv"), "w") as f:
    f.write(out)

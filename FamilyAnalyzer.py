#!/usr/bin/env python
#
#  This scripts has the purpose of analyzing an orthoXML file.
#  It does so by providing several methods operating on the file.
#  
#                            Adrian Altenhoff, June 2013
#  
try:
    import xml.etree.cElementTree as etree
except ImportError:
    import xml.etree.ElementTree as etree
import collections
import sys

class OrthoXMLParser(object):
    ns={"ns0":"http://orthoXML.org/2011/"}
    def __init__(self, filename):
        """creates a OrthoXMLParser object. the parameter filename needs to
        be a path pointing to the orthoxml file to be analyzed."""
        doc = etree.parse(filename)
        self.root = doc.getroot()

        self.__buildMappings()

    def mapGeneToXRef(self, id, typ='geneId'):
        if typ is None:
            res = id
        else:
            tup = (id,typ)
            res = self.__xrefs.get(tup,None)
            if res is None:
                # fallback, if not in dict
                gene = self.root.find('.*//{{{}}}gene[@id="{}"]'.format(
                    self.ns['ns0'], id))
                # return desired ID typ, otherwise whole element
                if typ is not None:
                    res = gene.get(typ, gene) 
                else:
                    res = gene
        return res
        
    def getSpeciesSet(self):
        return self.__species

    def getGeneIds(self):
        return self.__gene2species.keys()
    
    def getToplevelGroups(self):
        """A function yielding the toplevel orthologGroups from the file.
        This corresponds to gene families for the purpose of this project."""
        
        return self.root.findall(".//{{{ns0}}}groups/{{{ns0}}}orthologGroup".
           format(**self.ns))

    def getSubFamilies(self, level, root=None):
        """return a forest of orthologGroup nodes with roots at the given 
        taxonomic level. This function requires that the orthologGroup nodes
        are annotated with a 'property' element with 'TaxRange' as a name and
        the actual level in a 'value' attribute. This is not required by the
        orthoxml schema."""
        if root is None:
            root = self.root
        return root.findall(
            ".//{{{0}}}property[@name='TaxRange'][@value='{1}']/.."
            .format(self.ns['ns0'], level))

    def mapGeneToSpecies(self, id, typ='name'):
        if self.__gene2species is None:
            self.__buildMappings()
        return self.__gene2species[id].get(typ)

    def __buildMappings(self):
        mapping=dict()
        xref = dict()
        for species in self.root.findall(".//{{{ns0}}}species".format(**self.ns)):
            genes = species.findall(".//{{{ns0}}}gene".format(**self.ns))
            for gene in genes:
                id=gene.get('id')
                mapping[gene.get('id')] = species
                for tag in gene.keys():
                    if tag!="id":
                        xref[(id,tag)]=gene.get(tag) 
        self.__gene2species=mapping
        self.__xrefs = xref
        self.__species = frozenset({z.get('name') for z in mapping.values()})

    def getUbiquitusFamilies(self, minCoverage=.5):
        families = self.getToplevelGroups();
        return filter(
            lambda x:len(self.getGenesPerSpeciesInFam(x))>=minCoverage*len(self.getSpeciesSet()), 
            families)

    def getGenesPerSpeciesInFam(self, fam):
        genes=collections.defaultdict(set)
        geneRefs = fam.findall(".//{{{ns0}}}geneRef".format(**self.ns));
        for gref in geneRefs:
            gid = gref.get('id')
            sp = self.mapGeneToSpecies(gid)
            genes[sp].add(gid)
        return genes


    def getFamHistory(self, species=None, level=None):
        gene2copies = collections.defaultdict(list)
        famWhereLost = list()
        # assure that orthologGroup xml elements annotated with an 'og' attr
        if self.root.find(".//*[@og]") is None:
            GroupAnnotator(self).annotateDoc()

        famHist = FamHistory(self, species, level)
        for fam in self.getSubFamilies(level):
            genesInFam = self.getGenesPerSpeciesInFam(fam)
            famHist.addFamily(fam)
        return famHist


class FamHistory(object):
    XRefTag = None;
    def __init__(self, parser, species, level):
        self.parser = parser
        self.species = set(species)
        self.level = level
        self.__gene2copies = dict()
        for g in self.species:
            self.__gene2copies[g] = collections.defaultdict(list)
        self.__gene2fam = dict()
        self.__famWhereLost = collections.defaultdict(list)
    
    def setXRefTag(self, tag):
        self.XRefTag = tag

    def addFamily(self, fam):
        genesPerSpec = self.parser.getGenesPerSpeciesInFam(fam)
        famId = fam.get('og')
        for species in self.species:
            gids = genesPerSpec.get(species)
            if gids is None:
                self.__famWhereLost[species].append(famId)
            else:
                for gid in gids:
                    self.__gene2copies[species][gid] = gids
                    self.__gene2fam[gid] = famId
 
    def __getFamOfGeneIds(self, gids):
        coveredFams = set(map(lambda x: self.__gene2fam.get(x,None), gids))
        return coveredFams
    
    def write(self):
        print('\nFamily Analysis:')
        for species in self.species:
            gids = filter(lambda gid:self.parser.mapGeneToSpecies(gid)==species,
                self.parser.getGeneIds())
            gids.sort(cmp=lambda x,y:len(self.__gene2copies[species][x]) - len(self.__gene2copies[species][y]) )
            coveredFams = set(map(lambda x: self.__gene2fam.get(x,None), gids))
            print("{} - {} of {} sub-families covered".
                format(species, len(coveredFams), 
                    len(coveredFams)+ len(self.__famWhereLost[species])))
            for gid in gids:
                if len(self.__gene2copies[species][gid])<=0:
                    print(" {}: n/a (singleton not in any family)".format(
                        self.parser.mapGeneToXRef(gid, self.XRefTag)))
                else:
                    args = dict(gXref=self.parser.mapGeneToXRef(gid,self.XRefTag),
                        famId=self.__gene2fam[gid],
                        cnt=len(self.__gene2copies[species][gid]),
                        sib=";".join([self.parser.mapGeneToXRef(z,self.XRefTag) 
                            for z in self.__gene2copies[species][gid]]))
                    print(" {gXref}: {famId} ({cnt}): {sib}".format(**args))
            for fam in self.__famWhereLost[species]:
                print(" n/a: {} (0) no member in subfamily".format(fam))
                
        
class GroupAnnotator(object):
    def __init__(self, parser):
        self.parser = parser
        self.ns=parser.ns

    def __getNextSubId(self, idx):
        while len(self.dupCnt)<idx:
            self.dupCnt.append(0)
        self.dupCnt[idx-1] += 1
        return self.dupCnt[idx-1]

    def __encodeParalogClusterId(self, prefix, nr):
        while nr/25>0:
            prefix += chr(97+nr%25)
            nr /= 25
        prefix += chr(97+nr%25)
        return prefix

    def __annotateGroupR(self, node, og, idx=0):
        if node.tag=="{{{ns0}}}orthologGroup".format(**self.ns):
            node.set('og',og)
            for child in node.getchildren():
                self.__annotateGroupR(child, og, idx)
        elif node.tag=="{{{ns0}}}paralogGroup".format(**self.ns):
            idx += 1
            nextOG = "{}.{}".format(og, self.__getNextSubId(idx))
            for i, child in enumerate(node.getchildren()):
                self.__annotateGroupR(child, self.__encodeParalogClusterId(nextOG,i), idx);
    
    def annotateDoc(self):
        for i, fam in enumerate(self.parser.getToplevelGroups()):
            self.dupCnt=list()
            self.__annotateGroupR(fam, fam.get('id',str(i)))

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Analyze Hierarchical OrthoXML families.')
    parser.add_argument('-t','--xreftag', default=None, help='xref tag of genes to report')
    parser.add_argument('path', help='path to orthoxml file')
    parser.add_argument('level', help='taxonomic level at which analysis should be done')
    parser.add_argument('species', nargs="+", help='(list of) species to be analyzed')
    args = parser.parse_args()

    op = OrthoXMLParser( args.path )
    print("Analyzing {} on taxlevel {}".format(args.path, args.level))
    print("Species found:")
    print("; ".join(op.getSpeciesSet()))
    print("--> analyzing " + "; ".join(args.species))

    #print op.getSubFamilies("mouse2_mouse")
    #for fam in op.getUbiquitusFamilies(minCoverage=.75):
    #    print fam.get('id');

    hist=op.getFamHistory(args.species, args.level)
    hist.setXRefTag(args.xreftag)
    hist.write()

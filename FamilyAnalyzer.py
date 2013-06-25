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

class ElementError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return str(self.msg)

class OrthoXMLParser(object):
    ns={"ns0":"http://orthoXML.org/2011/"} # xml namespace
    def __init__(self, filename):
        """creates a OrthoXMLParser object. the parameter filename needs to
        be a path pointing to the orthoxml file to be analyzed."""
        doc = etree.parse(filename)
        self.root = doc.getroot()

        self.__buildMappings() # builds three dictionaries - see def below

    def mapGeneToXRef(self, id, typ='geneId'):
        """
        Looks up id (integer id number as a string) and type in self.__xrefs
        dictionary (aka self._OrthoXMLParser__xrefs)
        Returns lookup value, which is a gene name
        """
        if typ is None:
            res = id
        else:
            tup = (id,typ)
            res = self.__xrefs.get(tup,None)
            if res is None:
                # fallback, if not in dict
                gene = self.root.find('.*//{{{}}}gene[@id="{}"]'.format(self.ns['ns0'], id))
                # return desired ID typ, otherwise whole element
                if typ is not None:
                    res = gene.get(typ, gene) 
                else:
                    res = gene
        return res
        
    def getSpeciesSet(self):
        return self.__species # all species in the xml tree

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
            ".//{{{0}}}property[@name='TaxRange'][@value='{1}']/..".format(self.ns['ns0'], level))

    def getSubFamiliesByRecursion(self, level, root=None, subfamilies=None):
        """return a forest of orthologGroup nodes with roots at the given 
        taxonomic level or lower. Each node is the topmost matching node - any
        matching children are not returned, to avoid multiple counting.

        This function requires that the orthologGroup nodes
        are annotated with a 'property' element with 'TaxRange' as a name and
        the actual level in a 'value' attribute, and that the value is a 
        '/' joined list of species names. This is not required by the
        orthoxml schema."""

        if subfamilies is None:
            subfamilies = []

        if root is None:
            root = self.root

        # Stopping criterion for matching top-level node
        if (self.__is_ortholog_group(root)
            and self.__is_at_desired_level(root, level) ):
            subfamilies.append(root)
            return

        # Stopping criterion no children
        if len(root) == 0:
            return

        # Recurse on child nodes
        for child in root.getchildren():
            self.getSubFamiliesByRecursion(level, child, subfamilies)

        return subfamilies

    def __is_ortholog_group(self, element):
        return element.tag == '{http://orthoXML.org/2011/}orthologGroup'

    def __is_at_desired_level(self, element, querylevel):
        """
        Tests if the element is at a taxonomic level equal to or below
        `querylevel`.
        If querylevel == "LUCA" -> return True, because necessarily all
        nodes in the XML are equal to or below LUCA.
        If element level == "LUCA", but query level is not "LUCA" (which it 
        isn't, because it would have returned True) return false, because the 
        element's level is necessarily above the querylevel.
        Otherwise, test that the element's level is a subset of querylevel, by
        splitting the level string on '/' and making a set. NB issubset returns
        True when sets are equal.
        """
        if not self.__is_ortholog_group(element): 
            raise ElementError('Not an orthologGroup node')
            # return False
        if querylevel == 'LUCA':
            return True
        if not isinstance(querylevel, set):
            querylevel = set(querylevel.split('/'))
        prop = element.find('{http://orthoXML.org/2011/}property')
        level = prop.get('value')
        if level == "LUCA":
            return False
        level = set(level.split('/'))
        return level.issubset(querylevel)

    def mapGeneToSpecies(self, id, typ='name'):
        """
        Does a lookup in the self.__gene2species dict:
        key = idnum, return = species name
        """
        if self.__gene2species is None:
            self.__buildMappings()
        return self.__gene2species[id].get(typ)

    def __buildMappings(self):
        """
        Builds two dictionaries:
          self.__gene2species - keys are ID numbers, values are species
          self.__xrefs - keys are tuples (idnum, idtype ['geneId','protId']), values are gene names
        Also builds the set:
          self.__species - All species names in the xml tree
        """
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
        self.__levels = set(n.get('value') for n in self.root.findall(
            ".//{{{0}}}property".format(self.ns['ns0'])))

    def getUbiquitusFamilies(self, minCoverage=.5):
        families = self.getToplevelGroups();
        return filter(
            lambda x:len(self.getGenesPerSpeciesInFam(x))>=minCoverage*len(self.getSpeciesSet()), 
            families)

    def getLevels(self):
        return self.__levels

    def getGenesPerSpeciesInFam(self, fam):
        """
        Takes a gene family, returns a dictionary:
          keys   = species names
          values = set of geneIds belonging to that species at a level descended
                from the family
        """
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

    def getFamHistoryByRecursion(self, species=None, level=None):
        gene2copies = collections.defaultdict(list)
        famWhereLost = list()
        # assure that orthologGroup xml elements annotated with an 'og' attr
        if self.root.find(".//*[@og]") is None:
            GroupAnnotator(self).annotateDoc()

        famHist = FamHistory(self, species, level)
        for fam in self.getSubFamiliesByRecursion(level):
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
        letters = []
        while nr/26 > 0:
            letters.append(chr(97+nr%26))
            nr = nr/26 - 1
        letters.append(chr(97+nr%26))
        return prefix+''.join(letters[::-1]) # letters were in reverse order

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
    parser.add_argument('--show_levels', action='store_true', help='show available levels and species and quit')
    parser.add_argument('-r', '--use-recursion', action='store_true', help='Use recursion to sample families that are a subset of the query')
    parser.add_argument('path', help='path to orthoxml file')
    parser.add_argument('level', help='taxonomic level at which analysis should be done')
    parser.add_argument('species', nargs="+", help='(list of) species to be analyzed')
    args = parser.parse_args()

    op = OrthoXMLParser( args.path )
    if args.show_levels:
        print("Species:\n{0}\n\nLevels:\n{1}".format('\n'.join(sorted(list(op.getSpeciesSet()))),
                '\n'.join(sorted(op.getLevels()))))
        sys.exit()    
    print("Analyzing {} on taxlevel {}".format(args.path, args.level))
    print("Species found:")
    print("; ".join(op.getSpeciesSet()))
    print("--> analyzing " + "; ".join(args.species))

    #print op.getSubFamilies("mouse2_mouse")
    #for fam in op.getUbiquitusFamilies(minCoverage=.75):
    #    print fam.get('id');

    if args.use_recursion:
        hist=op.getFamHistoryByRecursion(args.species, args.level)
    else:
        hist=op.getFamHistory(args.species, args.level)
    hist.setXRefTag(args.xreftag)
    hist.write()

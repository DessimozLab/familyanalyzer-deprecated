#!/usr/bin/env python
#
#  This scripts has the purpose of analyzing an orthoXML file.
#  It does so by providing several methods operating on the file.
#
#                            Adrian Altenhoff, June 2013
#
import lxml.etree as etree
import collections
import itertools
import io

class ElementError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return str(self.msg)


class OrthoXMLQuery(object):
    """Helper class with predefined queries on an orthoxml tree."""

    ns = {"ns0": "http://orthoXML.org/2011/"}   # xml namespace

    @classmethod
    def getToplevelOrthologGroups(cls, root):
        """returns a list with the toplevel orthologGroup elements
        of the given root element."""
        xquery = ".//{{{ns0}}}groups/{{{ns0}}}orthologGroup".format(**cls.ns)
        return root.findall(xquery)

    @classmethod
    def getGeneFromId(cls, id_, root):
        xquery = ".*//{{{}}}gene[@id='{}']".format(cls.ns['ns0'], id_)
        return root.findall(xquery)

    @classmethod
    def getGroupsAtLevel(cls, level, root):
        """returns a list with the orthologGroup elements which have a
        TaxRange property equals to the requested level."""
        xquery = (".//{{{0}}}property[@name='TaxRange'][@value='{1}']/..".
                  format(cls.ns['ns0'], level))
        return root.findall(xquery)

    @classmethod
    def getSubNodes(cls, targetNode, root, recursivly=True):
        """method which returns a list of all (if recursively
        is set to true) or only the direct children nodes
        having 'targetNode' as their tagname.
        The namespace is automatically added to the tagname."""
        xPrefix = ".//" if recursivly else "./"
        xquery = "{}{{{}}}{}".format(xPrefix, cls.ns['ns0'], targetNode)
        return root.findall(xquery)

    @classmethod
    def is_geneRef_node(cls, element):
        """check whether a given element is an instance of a geneRef
        element."""
        return element.tag == '{{{ns0}}}geneRef'.format(**cls.ns)

    @classmethod
    def getLevels(cls, element):
        """returns a list of the TaxRange levels associated to the
        passed orthologGroup element. If the element does not have
        any TaxRange property tags associated, an empty list is
        returned."""
        propTags = cls.getSubNodes("property", element, recursivly=False)
        res = [t.get('value') for t in propTags if t.get('name')=='TaxRange']
        return res


class OrthoXMLParser(object):
    ns = {"ns0": "http://orthoXML.org/2011/"}   # xml namespace

    def __init__(self, filename):
        """creates a OrthoXMLParser object. the parameter filename
        needs to be a path pointing to the orthoxml file to be
        analyzed."""
        self.doc = etree.parse(filename)
        self.root = self.doc.getroot()
        self.tax = None

        self._buildMappings()   # builds three dictionaries - see def below

    def write(self, filename):
        """Write out the (modified) orthoxml file into a new file."""
        self.doc.write(filename)

    def mapGeneToXRef(self, id_, typ='geneId'):
        """
        Looks up id_ (integer id_ number as a string) and type in self._xrefs
        dictionary (aka self._xrefs)
        Returns lookup value, which is a gene name
        """
        if typ is None:
            res = id_
        else:
            tup = (id_, typ)
            res = self._xrefs.get(tup, None)
            if res is None:
                # fallback, if not in dict
                gene = OrthoXMLQuery.getGeneFromId(id_, self.root)
                # return desired ID typ, otherwise whole element
                if typ is not None:
                    res = gene.get(typ, gene)
                else:
                    res = gene
        return res

    def getSpeciesSet(self):
        return self._species  # all species in the xml tree

    def getGeneIds(self, speciesFilter=None, tag="name"):
        genes = list(self._gene2species.keys())
        if speciesFilter is not None:
            genes = [g for g in genes if self._gene2species[g].get(tag, None) in speciesFilter]
        return genes

    def getToplevelGroups(self):
        """A function yielding the toplevel orthologGroups from the file.
        This corresponds to gene families for the purpose of this project."""
        return OrthoXMLQuery.getToplevelOrthologGroups(self.root)

    def getSubFamilies(self, level, root=None):
        """return a forest of orthologGroup nodes with roots at the
        given taxonomic level. This function requires that the
        orthologGroup nodes are annotated with a 'property' element
        with 'TaxRange' as a name and the actual level in a 'value'
        attribute. This is not required by the orthoxml schema."""
        if root is None:
            root = self.root

        return OrthoXMLQuery.getGroupsAtLevel(level, root)

    def getSubFamiliesByRecursion(self, level, root=None, subfamilies=None):
        """return a forest of orthologGroup nodes with roots at the
        given taxonomic level or lower. Each node is the topmost
        matching node - any matching children are not returned, to
        avoid multiple counting.

        This function requires that the orthologGroup nodes
        are annotated with a 'property' element with 'TaxRange' as a
        name and the actual level in a 'value' attribute, and that the
        value is a '/' joined list of species names. This is not
        required by the orthoxml schema."""

        if subfamilies is None:
            subfamilies = []

        if root is None:
            root = self.root

        # Stopping criterion for matching top-level node
        if (self.is_ortholog_group(root) and
                self._is_at_desired_level(root, level)):
            subfamilies.append(root)
            return

        # Stopping criterion no children
        if len(root) == 0:
            return

        # Recurse on child nodes
        for child in root.getchildren():
            self.getSubFamiliesByRecursion(level, child, subfamilies)

        return subfamilies

    @classmethod
    def is_ortholog_group(cls, element):
        """
        Returns true if the passed element is an orthologGroup xml node
        """
        return element.tag == '{{{ns0}}}orthologGroup'.format(**cls.ns)

    @classmethod
    def is_paralog_group(cls, element):
        """
        Returns true if the passed element is an paralogGroup xml node
        """
        return element.tag == '{{{ns0}}}paralogGroup'.format(**cls.ns)

    @classmethod
    def is_evolutionary_node(cls, element):
        """Returns true if the passed element is an evolutionary event
        xml node, i.e. if it is either an orthologGroup or a
        paralogGroup element."""
        return (cls.is_ortholog_group(element) or
                cls.is_paralog_group(element))

    def _is_at_desired_level(self, element, querylevel):
        """Tests if the element is at a taxonomic level equal to or
        below `querylevel`.
        If querylevel == "LUCA" -> return True, because necessarily
        all nodes in the XML are equal to or below LUCA.
        If element level == "LUCA", but query level is not "LUCA"
        (which it isn't, because it would have returned True) return
        false, because the element's level is necessarily above the
        querylevel.
        Otherwise, test that the element's level is a subset of
        querylevel, by splitting the level string on '/' and making
        a set. NB issubset returns True when sets are equal."""

        if not self.is_ortholog_group(element):
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

    def mapGeneToSpecies(self, id_, typ='name'):
        """
        Does a lookup in the self._gene2species dict:
        key = idnum, return = species name
        """
        if self._gene2species is None:
            self._buildMappings()
        return self._gene2species[id_].get(typ)

    def _findSubNodes(self, targetNode, root=None):
        """return all (recursively) found elements with tagname
        'targetNode' below the root element. If no root element
        is provided, search starts at the document root."""
        rootNode = root if root is not None else self.root
        return rootNode.findall(".//{{{0}}}{1}".
                                format(self.ns['ns0'], targetNode))

    def _buildMappings(self):
        """Builds two dictionaries:
          self._gene2species - keys are ID numbers, values are species
          self._xrefs - keys are tuples
              (idnum, idtype ['geneId','protId']),
          values are gene names
        Also builds the set:
          self._species - All species names in the xml tree"""
        mapping = dict()
        xref = dict()
        for species in self._findSubNodes("species"):
            genes = self._findSubNodes("gene", root=species)
            for gene in genes:
                id_ = gene.get('id')
                mapping[id_] = species
                for tag in gene.keys():
                    if tag != "id":
                        xref[(id_, tag)] = gene.get(tag)
        self._gene2species = mapping
        self._xrefs = xref
        self._species = frozenset({z.get('name') for z in mapping.values()})
        self._levels = frozenset({n.get('value')
            for n in self._findSubNodes("property")
            if n.get('name') == "TaxRange"})

    def getUbiquitusFamilies(self, minCoverage=.5):
        families = self.getToplevelGroups()
        return [x for x in families if len(self.getGenesPerSpeciesInFam(x)) >=
                minCoverage*len(self.getSpeciesSet())]

    def getLevels(self):
        return self._levels

    def getGenesPerSpeciesInFam(self, fam):
        """
        Takes a gene family, returns a dictionary:
          keys   = species names
          values = set of geneIds belonging to that species at a level
                descended from the family
        """
        genes = collections.defaultdict(set)
        geneRefs = self._findSubNodes("geneRef", fam)
        for gref in geneRefs:
            gid = gref.get('id')
            sp = self.mapGeneToSpecies(gid)
            genes[sp].add(gid)
        return genes

    def getFamHistory(self):
        """This method returns a FamHistory object initialized with
        the most powerful LevelAnalysis handler. The handler depends
        on whether the parser contains a tax attribute, set by the
        augmentTaxonomyInfo method if called."""

        # assure that orthologGroup xml elements annotated with an 'og' attr
        if self.root.find(".//*[@og]") is None:
            GroupAnnotator(self).annotateDoc()

        analyzer = LevelAnalysisFactory().newLevelAnalysis(self)
        return FamHistory(self, analyzer)

    def getFamHistoryByRecursion(self, species=None, level=None):
        # assure that orthologGroup xml elements annotated with an 'og' attr
        if self.root.find(".//*[@og]") is None:
            GroupAnnotator(self).annotateDoc()

        famHist = FamHistory(self, species, level)
        for fam in self.getSubFamiliesByRecursion(level):
            famHist.addFamily(fam)
        return famHist

    def augmentTaxonomyInfo(self, tax, propagate_top):
        """Assign a taxonomy to the orthoxml file. this taxonomy
        is used to augment the xml with the relevant level infos
        as 'TaxRange' property tags in orthologGroup elements.

        The propagate_top parameter can be should be used to enable
        or disable the propagation of all the levels which are
        older than the families topmost level. In other words, if
        enabled, all families arouse at LUCA, otherwise families
        can be invented later on in evolution."""
        if self.tax is not None:
            raise Exception("a taxonomy can be assigned only once")
        self.tax = tax
        GroupAnnotator(self).annotateMissingTaxRanges(tax, propagate_top)


class TaxonomyInconsistencyError(Exception):
    pass


class TaxonomyFactory(object):
    @classmethod
    def newTaxonomy(cls, arg):
        if isinstance(arg, str):
            if arg.endswith('.xml'):
                return XMLTaxonomy(arg)
            else:
                suffix = arg[arg.rindex('.'):]
                if suffix in ['.nwk', '.tree', '.newick']:
                    return NewickTaxonomy(arg) 
        elif isinstance(arg, OrthoXMLParser):
            return TaxRangeOrthoXMLTaxonomy(arg)
        else:
            raise NotImplementedError("unknown type of Taxonomy")


class Taxonomy(object):
    def __init__(self):
        raise NotImplementedError("abstract class")

    def iterParents(self, node, stopBefor=None):
        """iterates over all the taxonomy nodes towards the root
        which are above 'node' and below 'stopBefor'."""

        if node == stopBefor:
            return
        tn = self.hierarchy[node]
        while tn.up is not None and tn.up.name != stopBefor:
            tn = tn.up
            yield tn.name

    def _countParentAmongLevelSet(self, levels):
        """helper method to count for each level how many levels
        are parent levels. e.g. (arrow: is-partent-of)
          A->B->C
              \>D->E
        will return A=0,B=1,C=2,D=2,E=3"""
        levelSet = set(levels)
        cnts = dict()
        for lev in levelSet:
            t = set(self.iterParents(lev)).intersection(levelSet)
            cnts[lev] = len(t)
        return cnts

    def mostSpecific(self, levels):
        """returns the most specific (youngest) level among a set of
        levels. it is required that all levels are on one monophyletic
        lineage, otherwise an Exception is raised."""
        # count who often each element is a child of any other one.
        # the one with len(levels)-1 is the most specific level
        cnts = self._countParentAmongLevelSet(levels)
        for lev, cnt in cnts.items():
            if cnt==len(levels)-1:
                return lev
        raise Exception("Non of the element is subelement of all others")

    def mostGeneralLevel(self, levels):
        """returns the most general (odest) level among a set of levels."""
        # count who often each element is a child of any other one.
        # the one with len(levels)-1 is the most specific level
        cnts = self._countParentAmongLevelSet(levels)
        for lev, cnt in cnts.items():
            if cnt==0:
                return lev
        raise Exception("Non of the element is the root of all others")

    def printSubTreeR(self, fd, lev=None, indent=0):
        if lev is None:
            lev = self.root
        fd.write("{}{}\n".format(" "*2*indent, lev))
        for child in self.hierarchy[lev].down:
            self.printSubTreeR(fd, child.name, indent+1)

    def __str__(self):
        fd = io.StringIO()
        self.printSubTreeR(fd)
        res = fd.getvalue()
        fd.close()
        return res


class XMLTaxonomy(Taxonomy):
    def __init__(self, filename):
        raise NotImplementedError("XML Taxonomies have not yet been implemented")


class TaxRangeOrthoXMLTaxonomy(Taxonomy):
    def __init__(self, parser):
        self.parser = parser
        self.extractAdjacencies()
        self.bloat_all()
        self.extractHierarchy()

    def _parseParentChildRelsR(self, grp):
        levels = None
        if OrthoXMLParser.is_ortholog_group(grp):
            levels = [l.get('value') for l in grp.findall(
                './{{{ns0}}}property[@name="TaxRange"]'
                .format(**OrthoXMLParser.ns))]
        directChildNodes = list(grp)
        children = [child for child in directChildNodes
                    if OrthoXMLParser.is_evolutionary_node(child)]
        geneRefs = [node for node in directChildNodes if OrthoXMLQuery.is_geneRef_node(node)]
        speciesOfGenes = {self.parser.mapGeneToSpecies(x.get('id')) for x in geneRefs}


        # recursively process childreen nodes
        subLevs = speciesOfGenes
        for child in children:
            subLevs.update(self._parseParentChildRelsR(child))

        if levels is not None:
            for parent in levels:
                for child in subLevs:
                    self.adj.add((parent, child))
            subLevs = set(levels)
        return subLevs

    def extractAdjacencies(self):
        self.adj = set()
        for grp in self.parser.getToplevelGroups():
            self._parseParentChildRelsR(grp)

        self.nodes = set(itertools.chain(*self.adj))

    def bloat_all(self):
        """build transitive closure of all parent - child relations"""
        while(self.bloat()):
            pass

    def bloat(self):
        found = False
        for pair in itertools.product(self.nodes, repeat=2):
            first, second = pair
            for node in self.nodes:
                if ((pair not in self.adj) and ((first, node) in self.adj) and
                        ((node, second) in self.adj)):
                    found = True
                    self.adj.add(pair)
        return found

    def extractHierarchy(self):
        self.hierarchy = dict(zip(self.nodes, map(TaxNode, self.nodes)))
        for pair in itertools.product(self.nodes, repeat=2):
            if pair in self.adj:
                if self.good(pair):
                    first, second = pair
                    #print "%s,%s is good" % pair
                    self.hierarchy[first].addChild(self.hierarchy[second])
                    self.hierarchy[second].addParent(self.hierarchy[first])
        noParentNodes = [z for z in self.nodes if self.hierarchy[z].up is None]
        if len(noParentNodes) != 1:
            raise TaxonomyInconsistencyError(
                "Warning: several/none TaxonomyNodes are roots: {}"
                .format(noParentNodes))
        self.root = noParentNodes[0]

    def good(self, pair):
        first, second = pair
        for node in self.nodes:
            if (first, node) in self.adj and (node, second) in self.adj:
                return False
        return True


class TaxNode(object):
    def __init__(self, name):
        self.name = name
        self.up = None
        self.down = list()

    def addChild(self, c):
        if not c in self.down:
            self.down.append(c)

    def addParent(self, p):
        if self.up is not None and self.up != p:
            raise TaxonomyInconsistencyError(
                "Level {} has several parents, at least two: {}, {}"
                .format(self.name, self.up.name, p.name))
        self.up = p

    def iterLeaves(self):
        if len(self.down) == 0:
            yield self
        else:
            for child in self.down:
                for elem in child.iterLeaves():
                    yield elem


class GeneFamily(object):
    """GeneFamily(root_element)

    Represents one gene family rooted at an orthologous group. """
    def __init__(self, root_element):
        if not OrthoXMLParser.is_ortholog_group(root_element):
            raise ElementError('Not an orthologGroup node')
        self.root = root_element

    def getMemberGenes(self):
        members = self.root.findall('.//{{{ns0}}}geneRef'.
                format(**OrthoXMLParser.ns))
        return [x.get('id') for x in members]

    def getFamId(self):
        return self.root.get('og')

    def getLevels(self):
        return OrthoXMLQuery.getLevels(self.root)

    def analyzeLevel(self, level):
        """analyze the structure of the family at a given taxonomic
        level.
        returns a list of LevelAnalysis object, one per sub-family"""

        subFamNodes = OrthoXMLQuery.getGroupsAtLevel(level, self.root)
        subFams = [GeneFamily(fam) for fam in subFamNodes]
        return subFams

    def analyze(self, strategy):
        """analyze the history of the GeneFamily using the strategy
        passed to the method. The strategy arguement must be an
        object providing a analyzeGeneFam method,
        e.g. a LevelAnalysis object."""
        self.summary = strategy.analyzeGeneFam(self)

    def write(self, fd, speciesFilter=None, idFormatter=lambda x: x):
        species = list(self.summary.keys())
        if speciesFilter is not None:
            species = [g for g in species if g in speciesFilter]
        for spec in self.summary.keys():
            if not spec in species:
                continue
            for sumElem in self.summary[spec]:
                refs = "; ".join([idFormatter(gid) for gid in sumElem.genes])
                fd.write("{}\t{}\t{}\t{}:{}\n".format(
                    self.getFamId(), spec, len(sumElem.genes),
                    sumElem.typ, refs))


class Singletons(GeneFamily):
    def __init__(self, element):
        memb = set()
        if isinstance(element, str):
            memb.add(element)
        else:
            memb.update(element)
        self.members = memb

    def getMemberGenes(self):
        return self.members

    def getFamId(self):
        return "n/a"

    def getLevels(self):
        return None


    def analyzeLevel(self, level, parser):
        return self

    def analyze(self, strategy):
        super().analyze(strategy)
        for specSum in self.summary.values():
            for sumElement in specSum:
                sumElement.typ = "SINGLETON"


def enum(*sequential, **named):
    """creates an Enum type with given values"""
    enums = dict(zip(sequential, range(len(sequential))), **named)
    enums['reverse'] = dict((value, key) for key, value in enums.items())
    return type('Enum', (object, ), enums)


class SummaryOfSpecies(object):
    def __init__(self, typ, genes):
        self.typ = typ
        self.genes = genes


class LevelAnalysisFactory(object):
    def __init__(self):
        pass

    def newLevelAnalysis(self, parser):
        """return a the appropriate LevelAnalysis instance based on
        the presents/absence of a taxonomy in the parser"""
        if parser.tax is None:
            return BasicLevelAnalysis(parser)
        else:
            return TaxAwareLevelAnalysis(parser, parser.tax)


class BasicLevelAnalysis(object):
    GeneClasses = enum("MULTICOPY", "SINGLECOPY", "ANCIENT_BUT_LOST", "LATER_GAINED", "SINGLETON")

    def __init__(self, parser):
        self.parser = parser

    def analyzeGeneFam(self, fam):
        """analyzes a single gene family and returns a summary dict.

        This method classifies all genes in the family depending on
        the number of copies per genome into MULTICOPY or SINGLECOPY
        genes."""

        spec2genes = collections.defaultdict(set)
        for geneId in fam.getMemberGenes():
            spec = self.parser.mapGeneToSpecies(geneId)
            spec2genes[spec].add(geneId)
        summary = dict()
        for spec in iter(spec2genes.keys()):
            nrMemb = len(spec2genes[spec])
            gclass = self.GeneClasses.MULTICOPY if nrMemb > 1 else self.GeneClasses.SINGLECOPY
            summary[spec] = [SummaryOfSpecies(self.GeneClasses.reverse[gclass],
                                              spec2genes[spec])]
        return summary


class TaxAwareLevelAnalysis(BasicLevelAnalysis):
    def __init__(self, parser, tax):
        super().__init__(parser)
        self.tax = tax

    def analyzeGeneFam(self, fam):
        """analyzes a singel gene family in the context of a known
        taxonomic tree.

        in addition to the method defined in the base class, this
        method adds information of lost genes. It does this by
        checking whether a species within the taxonomic range of
        the family contains a copy of the gene. if not, it had
        been lost."""
        summary = super().analyzeGeneFam(fam)
        lev = fam.getLevels()
        if lev is not None:
            # if several levels exist at this node, use oldest one
            mostGeneralLevel = self.tax.mostGeneralLevel(lev)
            speciesCoveredByLevel = {l for l in self.tax.hierarchy[mostGeneralLevel].iterLeaves()}
            lostSpecies = speciesCoveredByLevel.difference(summary.keys())
            for lost in lostSpecies:
                summary[lost] = [(SummaryOfSpecies("ANCIENT_BUT_LOST",[]))]
        return summary


class FamHistory(object):
    XRefTag = None

    def __init__(self, parser, analyzer):
        self.parser = parser
        self.analyzer = analyzer

    def setXRefTag(self, tag):
        """set the attribute name of the 'gene' elements which should
        be used for report. defined by orthoxml are 'geneId' and
        'protId'. If not defined, the (numerical) internal ids are used."""
        self.XRefTag = tag

    def analyzeLevel(self, level):
        gfamList = list()
        for fam in self.parser.getToplevelGroups():
            gfamList.extend(GeneFamily(fam).analyzeLevel(level))

        gene2FamIdx = dict()
        specInTaxRange = set()
        for idx, gfam in enumerate(gfamList):
            gfam.analyze(self.analyzer)
            specInTaxRange.update(gfam.summary.keys())
            for mem in gfam.getMemberGenes():
                gene2FamIdx[mem] = idx
        # get genes in taxrange not belonging to any family. those are
        # considered to be singletons and added to a special GeneFamily.
        allGenesInRange = self.parser.getGeneIds(speciesFilter=specInTaxRange)
        singletonsSet = set(allGenesInRange).difference(gene2FamIdx.keys())
        singletons = Singletons(singletonsSet)
        singletons.analyze(self.analyzer)

        gfamList.append(singletons)
        for singleton in singletonsSet:
            gene2FamIdx[singleton] = len(gfamList)-1

        self.geneFamList = gfamList
        self.gene2FamIdx = gene2FamIdx
        self.analyzedLevel = level

    def write(self, fd, speciesFilter=None):
        """writes the FamHistory object to a given stream object
        in a human readable format.
        The stream object needs to have a write(str) method defined.
        The optional speciesFilter argument accepts a set of
        species names for which the genes in the families are
        reported."""

        fd.write("FamilyAnalysis at {}\n".format(self.analyzedLevel))
        for fam in self.geneFamList:
            fam.write(fd, speciesFilter, idFormatter=lambda gid:
                    self.parser.mapGeneToXRef(gid, self.XRefTag))

    def __str__(self):
        fd = io.StringIO()
        self.write(fd)
        res = fd.getvalue()
        fd.close()
        return res

    def _find_subfamilies(self, query, targetlist):
        """ Used in compare method - fixes s.startswith bug where
        queries such as '20' would return subfamilies '200', '2000', etc. 
        Now '20' will only return subfamilies of the form '20.1a', 20.1b.2c.3e', 
        etc. """
        result = []
        q = query.split('.') # e.g. 8.1b.2a -> ['8', '1b', '2a']

        for family in targetlist:
            name_elements = family.split('.')
            if name_elements[:len(q)] == q: # exact prefix match to query
                result.append(family)

        return result

    def compare(self, other, fd):
        """compares two FamilyHistory objects.

        The two FamilyHistory objects are ment to operate on the same
        dataset, but on different taxonomic levels. It compares all
        families with respect to these two levels and clusters them
        according to whether they are
         - identical (no duplication/loss on the lineage between the
           two levels
         - duplications (old family duplicated along the lineage into
           two or more subfamilies)
         - novel families (did not exist at older level)

        It is expected that the 'other' family object corresponds to
        the younger taxonomic level, but this is not enforced. """

        famIds = [gfam.getFamId() for gfam in self.geneFamList]
        otherfamIds = [gfam.getFamId() for gfam in other.geneFamList]
        for f in famIds:
            if f == "n/a":
                continue
            if f in otherfamIds:
                fd.write("{} identical\n".format(f))
            else:
                subfam = self._find_subfamilies(f, otherfamIds)
                if len(subfam) == 0:
                    fd.write("{} -> LOST\n".format(f))
                else:
                    fd.write("{} -> {}\n".format(f, "; ".join(subfam)))
        for f in otherfamIds:
            topId = f[:f.find('.')]
            if not any(map(lambda x:x.startswith(topId), famIds)):
                fd.write("n/a -> {}\n".format(f))


class GroupAnnotator(object):
    """this class annotates orthologGroup elements with the numbering
    schema presented in the LOFT paper:
    van der Heijden, Snel, van Noort, Huynen
    Orthology prediction at scalable resolution by phylogenetic tree analysis.
    BMC Bioinformatics, 2007, 8, 83

    and adding additional property tags with skipped TaxRange levels."""

    def __init__(self, parser):
        self.parser = parser
        self.ns = parser.ns

    def _getNextSubId(self, idx):
        """helper method to return the next number at a given depth of
        duplication (idx)"""
        while len(self.dupCnt) < idx:
            self.dupCnt.append(0)
        self.dupCnt[idx-1] += 1
        return self.dupCnt[idx - 1]

    def _encodeParalogClusterId(self, prefix, nr):
        """encode the paralogGroups at the same level, e.g. 1a, 1b, 1c
        for 3 paralogGroups next to each other. the nr argument
        identifies the individual indices of those 3 paralogGroups."""
        letters = []
        while nr//26 > 0:
            letters.append(chr(97 + (nr % 26)))
            nr = nr//26 - 1
        letters.append(chr(97 + (nr % 26)))
        return prefix+''.join(letters[::-1])  # letters were in reverse order

    def _annotateGroupR(self, node, og, idx=0):
        """create the og attributes at the orthologGroup elements
        according to the naming schema of LOFT. ParalogGroup elements
        do not get own attributes (not possible in the xml schema),
        but propagate their sub-names for the subsequent orthologGroup
        elements."""
        if self.parser.is_ortholog_group(node):
            node.set('og', og)
            for child in list(node):
                self._annotateGroupR(child, og, idx)
        elif self.parser.is_paralog_group(node):
            idx += 1
            nextOG = "{}.{}".format(og, self._getNextSubId(idx))
            for i, child in enumerate(list(node)):
                self._annotateGroupR(child,
                                     self._encodeParalogClusterId(nextOG, i),
                                     idx)

    def _addTaxRangeR(self, node, last=None, noUpwardLevels=False):
        """recursive method to ad TaxRange property tags."""
        if self.parser.is_ortholog_group(node):
            levels = {z.get('value') for z in node.findall(
                './{{{ns0}}}property[@name="TaxRange"]'
                .format(**self.ns))}
            mostSpecificLevel = self.tax.mostSpecific(levels)
            if noUpwardLevels:
                levelsToParent = set()
            else:
                levelsToParent = {l for l in self.tax.iterParents(mostSpecificLevel, last)}
                levelsToParent.add(mostSpecificLevel)
                if not levels.issubset(levelsToParent):
                    raise Exception("taxonomy not in correspondance with found hierarchy: {} vs {}"
                                    .format(levels, levelsToParent))
            addLevels = levelsToParent - levels
            for lev in addLevels:
                node.append(self._createTaxRangeTag(lev))
            for child in list(node):
                self._addTaxRangeR(child, mostSpecificLevel)

        elif self.parser.is_paralog_group(node):
            for child in list(node):
                self._addTaxRangeR(child, last)
        elif OrthoXMLQuery.is_geneRef_node(node):
            # we check whether the parent node is a direct ancester in the
            # tax or not. if not, we creaete a fake orthologGroup.
            spec = self.parser.mapGeneToSpecies(node.get('id'))
            expRange = self.tax.hierarchy[spec].up.name
            directParent = parent = node.getparent()
            while not self.parser.is_ortholog_group(parent):
                parent = parent.getparent()
            levOfParent = OrthoXMLQuery.getLevels(parent)
            mostSpecific = self.tax.mostSpecific(levOfParent)
            if expRange != mostSpecific:
                self._insertOG(directParent, node, expRange, mostSpecific)

    def _insertOG(self, parent, child, specificLev, beforeLev):
        pos = parent.index(child)
        el = etree.Element('{{{ns0}}}orthologGroup'.format(**self.parser.ns))
        el.append(self._createTaxRangeTag(specificLev))
        for lev in self.tax.iterParents(specificLev, stopBefor=beforeLev):
            el.append(self._createTaxRangeTag(lev))
        el.append(child)
        parent.insert(pos, el)

    def _createTaxRangeTag(self, lev):
        return etree.Element('{{{ns0}}}property'.format(**self.parser.ns),
                attrib=dict(name='TaxRange', value=lev))

    def annotateMissingTaxRanges(self, tax, propagate_top=False):
        """This function adds left-out taxrange property elements to
        the orthologGroup elements in the xml. It will add all the levels
        defined in the 'tax'-Taxonomy between the parents most specific
        level and the current nodes level. If no parent exists, all
        tax-levels above the current one are used."""
        self.tax = tax
        for fam in self.parser.getToplevelGroups():
            self._addTaxRangeR(fam, noUpwardLevels=not propagate_top)
        del self.tax

    def annotateDoc(self):
        """apply the LOFT naming schema to all the orthologGroups."""
        for i, fam in enumerate(self.parser.getToplevelGroups()):
            self.dupCnt = list()
            self._annotateGroupR(fam, fam.get('id', str(i)))


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Analyze Hierarchical OrthoXML families.')
    parser.add_argument('--xreftag', default=None, help='xref tag of genes to report')
    parser.add_argument('--show_levels', action='store_true', help='show available levels and species and quit')
    parser.add_argument('-r', '--use-recursion', action='store_true', help='Use recursion to sample families that are a subset of the query')
    parser.add_argument('--taxonomy', default='implicit', help='Taxonomy used to reconstruct intermediate levels. Has to be either "implicit" (default) or a path to a file. If set to "implicit", the taxonomy is extracted from the input OrthoXML file')
    parser.add_argument('--propagate_top', action='store_true', help='propagate taxonomy levels up to the toplevel. If not set, only intermediate levels are propagated.')
    parser.add_argument('--show_taxonomy', action='store_true', help='show taxonomy used to infer missing levels')
    parser.add_argument('--store_augmented_xml', default=None, help='if set to a filename, the input orthoxml file with augmented annotations is written')
    parser.add_argument('--compare_second_level',default=None, help='compare secondary level with primary one')
    parser.add_argument('path', help='path to orthoxml file')
    parser.add_argument('level', help='taxonomic level at which analysis should be done')
    parser.add_argument('species', nargs="+", help='(list of) species to be analyzed')
    args = parser.parse_args()

    op = OrthoXMLParser(args.path)
    if args.show_levels:
        print("Species:\n{0}\n\nLevels:\n{1}".format(
                '\n'.join(sorted(list(op.getSpeciesSet()))),
                '\n'.join(sorted(op.getLevels()))))
        sys.exit()
    print("Analyzing {} on taxlevel {}".format(args.path, args.level))
    print("Species found:")
    print("; ".join(op.getSpeciesSet()))
    print("--> analyzing " + "; ".join(args.species))

    if args.taxonomy == "implicit":
        tax = TaxonomyFactory.newTaxonomy(op)
    else:
        from newick import NewickTaxonomy
        tax = TaxonomyFactory.newTaxonomy(args.taxonomy)
        if isinstance(tax, NewickTaxonomy):
            tax.annotate_from_orthoxml(op)

    if args.show_taxonomy:
        print("Use following taxonomy")
        print(tax)

    # add taxonomy to parser
    op.augmentTaxonomyInfo(tax, args.propagate_top)

    if args.use_recursion:
        hist = op.getFamHistoryByRecursion(args.species, args.level)
    else:
        hist = op.getFamHistory()
        hist.analyzeLevel(args.level)
    if args.compare_second_level is None:
        hist.setXRefTag(args.xreftag)
        hist.write(sys.stdout, speciesFilter=args.species)
    else:
        hist2 = op.getFamHistory()
        hist2.analyzeLevel(args.compare_second_level)
        print("Comparing taxlevel {}\n to taxlevel {}".format(
            args.level, args.compare_second_level))
        hist.compare(hist2, sys.stdout)

    if args.store_augmented_xml is not None:
        op.write(args.store_augmented_xml)

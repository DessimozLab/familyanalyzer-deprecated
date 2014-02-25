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
import re
import sys
try:
    from progressbar import ProgressBar, Percentage, Timer, ETA, Bar
    PROGRESSBAR = True
except ImportError:
    PROGRESSBAR = False


MAXINT = sys.maxsize


def setup_progressbar(msg, size):
    if not msg.endswith(': '):
        msg += ': '

    widgets = [msg,
               Percentage(), ' ',
               Bar(), ' ',
               Timer(), ' ',
               ETA()]

    pbar = ProgressBar(widgets=widgets, maxval=size)
    return pbar


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
        genes = root.findall(xquery)
        if len(genes) > 1:
            raise ElementError('several gene nodes with id {} '
                               'exist'.format(id_))
        gene = genes[0] if len(genes)>0 else None
        return gene

    @classmethod
    def getGroupsAtLevel(cls, level, root):
        """returns a list with the orthologGroup elements which have a
        TaxRange property equals to the requested level."""
        xquery = (".//{{{0}}}property[@name='TaxRange'][@value='{1}']/..".
                  format(cls.ns['ns0'], level))
        return root.findall(xquery)

    @classmethod
    def getSubNodes(cls, targetNode, root, recursively=True):
        """method which returns a list of all (if recursively
        is set to true) or only the direct children nodes
        having 'targetNode' as their tagname.
        The namespace is automatically added to the tagname."""
        xPrefix = ".//" if recursively else "./"
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
        propTags = cls.getSubNodes("property", element, recursively=False)
        res = [t.get('value') for t in propTags if t.get('name') == 'TaxRange']
        return res

    @classmethod
    def getInputGenes(cls, root, species=None):
        """returns a list of all gene elements in the orthoxml inside
        <species><database> tags, i.e. the list of genes prior to running
        OMA-HOGS. Optionally filtered by species."""
        filter_ = ('[@name="{}"]'.format(species)
                   if species is not None else '')
        if filter_ > '':
            xquery = ('/ns:orthoXML/ns:species{}/ns:database/'
                      'ns:genes//ns:gene'.format(filter_))
        else:
            xquery = '//ns:gene'
        return root.xpath(xquery, namespaces={'ns': cls.ns['ns0']})

    @classmethod
    def getGroupedGenes(cls, root, species=None):
        """ returns a list of all geneRef elements inside <group> tags, i.e.
        the list of genes clustered into families after running OMA-HOGS.
        Optionally filtered by species."""
        filter_ = ('[@name="TaxRange"and@value="{}"]'.format(species)
                   if species is not None else '')
        if filter_ > '':
            xquery = ('/ns:orthoXML/ns:groups/ns:orthologGroup//ns:property{}/'
                      'following-sibling::ns:geneRef'.format(filter_))
        else:
            xquery = '//ns:geneRef'
        return root.xpath(xquery, namespaces={'ns': cls.ns['ns0']})


class OrthoXMLParser(object):
    ns = {"ns0": "http://orthoXML.org/2011/"}   # xml namespace

    def __init__(self, filename):
        """creates a OrthoXMLParser object. the parameter filename
        needs to be a path pointing to the orthoxml file to be
        analyzed."""
        self.doc = etree.parse(filename)
        self.root = self.doc.getroot()
        self.tax = None
        self.singletons = None

        self._buildMappings()   # builds three dictionaries - see def below

    def write(self, filename, **kwargs):
        """Write out the (modified) orthoxml file into a new file.
        kwargs include:
            pretty_print=[True/False],
            xml_declaration=[True/False]
            encoding=[e.g. 'UTF-8']"""
        if 'pretty_print' in kwargs:
            self._remove_whitespace()
        self.doc.write(filename, **kwargs)

    def _remove_whitespace(self):
        """Remove whitespace from text and tail fields so that lxml can
        pretty_print the tree"""
        for element in self.root.iter():
            if element.text is not None and element.text.isspace():
                element.text = None
            element.tail = None

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
            genes = [g for g in genes
                     if self._gene2species[g].get(tag, None) in speciesFilter]
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

    def augmentTaxonomyInfo(self, tax, propagate_top=False):
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

    def augmentSingletons(self):
        GroupAnnotator(self).annotateSingletons()


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

    def __iter__(self):
        for n in self.hierarchy[self.root].iterDescendents():
            yield n

    def iterParents(self, node, stopBefore=None):
        """iterates over all the taxonomy nodes towards the root
        which are above 'node' and below 'stopBefore'."""

        if node == stopBefore:
            return
        tn = self.hierarchy[node]
        while tn.up is not None and tn.up.name != stopBefore:
            tn = tn.up
            yield tn.name

    def _countParentAmongLevelSet(self, levels):
        """helper method to count for each level how many levels
        are parent levels. e.g. (arrow: is-parent-of)
          A->B->C
              \>D->E
        will return A=0,B=1,C=2,D=2,E=3"""
        levelSet = set(levels)
        counts = dict()
        for lev in levelSet:
            t = set(self.iterParents(lev)).intersection(levelSet)
            counts[lev] = len(t)
        return counts

    def mostSpecific(self, levels):
        """returns the most specific (youngest) level among a set of
        levels. it is required that all levels are on one monophyletic
        lineage, otherwise an Exception is raised."""
        # count how often each element is a child of any other one.
        # the one with len(levels)-1 is the most specific level
        counts = self._countParentAmongLevelSet(levels)
        for lev, count in counts.items():
            if count == len(levels)-1:
                return lev
        raise Exception("Non of the element is subelement of all others")

    def mostGeneralLevel(self, levels):
        """returns the most general (oldest) level among a set of levels."""
        # count who often each element is a child of any other one.
        # the one with len(levels)-1 is the most specific level
        counts = self._countParentAmongLevelSet(levels)
        for lev, count in counts.items():
            if count == 0:
                return lev
        raise Exception("Non of the element is the root of all others")

    def younger_than_filter(self, levels, oldest_permitted):
        """
        Filters a set of levels, removing any that are older than the
        oldest_permitted (oldest_permitted=string: node_name)
        """
        try:
            oldest_permitted_node = self.hierarchy[oldest_permitted]
        except KeyError:
            raise Exception('No node with name {} found in '
                            'Taxonomy'.format(oldest_permitted))
        permitted = [node.name
                     for node in oldest_permitted_node.iterDescendents()]

        return [lev for lev in levels if lev in permitted]

    def printSubTreeR(self, fd, lev=None, indent=0):
        if lev is None:
            lev = self.root
        fd.write("{}{}\n".format(" "*2*indent, lev))
        for child in self.hierarchy[lev].down:
            self.printSubTreeR(fd, child.name, indent+1)

    def get_histories(self, parser):

        histories = {}

        if PROGRESSBAR:
            pbar = setup_progressbar('Getting histories', len(self.hierarchy))
            pbar.start()

        for i, level in enumerate(self.hierarchy, start=1):
            history = parser.getFamHistory()
            history.analyzeLevel(level)
            histories[level] = history
            self.hierarchy[level].attachFamHistory(history)
            if PROGRESSBAR:
                pbar.update(i)

        if PROGRESSBAR:
            pbar.finish()

        self.histories = histories
        return histories

    def get_comparisons(self, parser):

        if getattr(self, 'histories', None) is None:
            self.get_histories(parser)

        comparisons = {}
        to_compare = [(node, child) for node in self
                                    for child in node.down]

        if PROGRESSBAR:
            pbar = setup_progressbar('Comparing', len(to_compare))
            pbar.start()

        for i, (parent, child) in enumerate(to_compare, start=1):
            parent_history = self.histories[parent.name]
            child_history = self.histories[child.name]
            parent_child_comparison = parent_history.compare(child_history)
            comparisons[(parent.name, child.name)] = parent_child_comparison
            child.attachLevelComparisonResult(parent_child_comparison)
            if PROGRESSBAR:
                pbar.update(i)

        if PROGRESSBAR:
            pbar.finish()

        return comparisons

    def newick(self):
        return str(self.hierarchy[self.root]) + ';'

    def __str__(self):
        fd = io.StringIO()
        self.printSubTreeR(fd)
        res = fd.getvalue()
        fd.close()
        return res


class XMLTaxonomy(Taxonomy):
    def __init__(self, filename):
        raise NotImplementedError("XML Taxonomies have not "
                                  "yet been implemented")


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
        geneRefs = [node for node in directChildNodes
                    if OrthoXMLQuery.is_geneRef_node(node)]
        speciesOfGenes = {self.parser.mapGeneToSpecies(x.get('id'))
                          for x in geneRefs}

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

    reg = re.compile(r'\W') # matches anything that's NOT a-z, A-Z, 0-9 or _

    def __init__(self, name):
        self.name = name
        self.up = None
        self.down = list()
        self.history = None
        self.comparison = None

    def __str__(self):
        if self.reg.search(self.name):
            label = '"{}"'.format(self.name)
        else:
            label = self.name
        if self.isLeaf():
            return label
        subtree = ', '.join(str(ch) for ch in self.down)
        return '({0}){1}'.format(subtree, label)

    def addChild(self, c):
        if not c in self.down:
            self.down.append(c)

    def addParent(self, p):
        if self.up is not None and self.up != p:
            raise TaxonomyInconsistencyError(
                "Level {} has several parents, at least two: {}, {}"
                .format(self.name, self.up.name, p.name))
        self.up = p

    def isLeaf(self):
        return len(self.down) == 0

    def isInner(self):
        return not self.isLeaf()

    def isRoot(self):
        return self.up is None

    def iterDescendents(self):
        yield self
        for child in self.down:
            for elem in child.iterDescendents():
                yield elem

    def iterLeaves(self):
        for elem in self.iterDescendents():
            if elem.isLeaf():
                yield elem

    def iterInnerNodes(self):
        for elem in self.iterDescendents():
            if elem.isInner():
                yield elem

    def attachFamHistory(self, history):
        self.history = history

    def attachLevelComparisonResult(self, comparison):
        self.comparison = comparison


class GeneFamily(object):
    """GeneFamily(root_element)

    Represents one gene family rooted at an orthologous group. """
    def __init__(self, root_element):
        if not OrthoXMLParser.is_ortholog_group(root_element):
            raise ElementError('Not an orthologGroup node')
        self.root = root_element

    def __repr__(self):
        return '{} (id#={})'.format(self.__class__.__name__, self.getFamId())

    def __lt__(self, other):
        return self._cmp() < other._cmp()

    def __gt__(self, other):
        return self._cmp() > other._cmp()

    def __le__(self, other):
        return self._cmp() <= other._cmp()

    def __ge__(self, other):
        return self._cmp() >= other._cmp()

    def __eq__(self, other):
        return self._cmp() == other._cmp()

    def prefix_match(self, other):
        query = self._cmp()
        target = other._cmp()
        return target[:len(query)] == query

    def getMemberGenes(self):
        members = self.root.findall('.//{{{ns0}}}geneRef'.
                                    format(**OrthoXMLParser.ns))
        return [x.get('id') for x in members]

    def getFamId(self):
        return self.root.get('og')

    def _cmp(self):
        """
        Enable sorting based on LOFT numbering scheme,
        such that (e.g.)
            1.1a < 1.1b,
            1.1b < 1.1c,
            1.1z < 1.1aa,
            2 < 10,
            10 < n/a,
        """
        fam = self.getFamId()
        comp = tuple((0, int(num),) if num else (len(alpha.strip('.')),
            alpha.strip('.'),) for (num, alpha) in re.findall(r'(\d+)|(\D+)',
            fam))
        return comp

    def getLevels(self):
        return OrthoXMLQuery.getLevels(self.root)

    def analyzeLevel(self, level):
        """analyze the structure of the family at a given taxonomic
        level.
        returns a list of GeneFamily object, one per sub-family"""

        subFamNodes = OrthoXMLQuery.getGroupsAtLevel(level, self.root)
        subFams = [GeneFamily(fam) for fam in subFamNodes]
        return subFams

    def analyze(self, strategy, level):
        """analyze the history of the GeneFamily using the strategy
        passed to the method. The strategy arguement must be an
        object providing a analyzeGeneFam method,
        e.g. a LevelAnalysis object."""
        self.summary = strategy.analyzeGeneFam(self, level)

    def write(self, fd, speciesFilter=None, idFormatter=lambda x: x):
        species = list(self.summary.keys())
        if speciesFilter is not None:
            species = [g for g in species if g in speciesFilter]
        for spec in self.summary.keys():
            if not spec in species:
                continue
            sumElem = self.summary[spec]
            refs = "; ".join([idFormatter(gid) for gid in sumElem.genes])
            fd.write("{}\t{}\t{}\t{}:{}\n".format(
                self.getFamId(), spec, len(sumElem.genes),
                sumElem.typ, refs))

    def is_singleton(self):
        if not hasattr(self, 'summary'):
            return False
        return 'SINGLETON' in {x.typ for x in self.summary.values()}


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
        for sumElement in self.summary.values():
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
        elif parser.singletons is None:
            return TaxAwareLevelAnalysis(parser, parser.tax)
        else:
            return SingletonAwareLevelAnalysis(parser,
                                               parser.tax,
                                               parser.singletons)


class BasicLevelAnalysis(object):
    GeneClasses = enum("MULTICOPY",
                       "SINGLECOPY",
                       "ANCIENT_BUT_LOST",
                       "LATER_GAINED",
                       "SINGLETON")

    def __init__(self, parser):
        self.parser = parser

    def analyzeGeneFam(self, fam, level=None):
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
            gclass = (self.GeneClasses.MULTICOPY
                        if nrMemb > 1
                        else self.GeneClasses.SINGLECOPY)

            summary[spec] = SummaryOfSpecies(self.GeneClasses.reverse[gclass],
                                             spec2genes[spec])
        return summary


class TaxAwareLevelAnalysis(BasicLevelAnalysis):
    def __init__(self, parser, tax):
        super().__init__(parser)
        self.tax = tax

    def addLosses(self, fam, summary, level):
        lev = fam.getLevels()
        if lev is not None:
            # if several levels exist at this node, use oldest one
            # that is not older than `level'
            lev = self.tax.younger_than_filter(lev, level)
            mostGeneralLevel = self.tax.mostGeneralLevel(lev)
            speciesCoveredByLevel = {
                l.name for l in
                self.tax.hierarchy[mostGeneralLevel].iterLeaves()
            }

            lostSpecies = speciesCoveredByLevel.difference(summary.keys())
            for lost in lostSpecies:
                summary[lost] = SummaryOfSpecies("ANCIENT_BUT_LOST", [])

    def analyzeGeneFam(self, fam, level):
        """analyzes a single gene family in the context of a known
        taxonomic tree.

        in addition to the method defined in the base class, this
        method adds information of lost genes. It does this by
        checking whether a species within the taxonomic range of
        the family contains a copy of the gene. if not, it had
        been lost."""
        summary = super().analyzeGeneFam(fam, level)
        self.addLosses(fam, summary, level)
        return summary


class SingletonAwareLevelAnalysis(TaxAwareLevelAnalysis):
    def __init__(self, parser, tax, singletons):
        super().__init__(parser, tax)
        self.singletons = singletons

    def analyzeGeneFam(self, fam, level):
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
            if nrMemb > 1:
                gclass = self.GeneClasses.MULTICOPY
            else:
                if fam.getFamId() in self.parser.singletons:
                    gclass = self.GeneClasses.SINGLETON
                else:
                    gclass = self.GeneClasses.SINGLECOPY

            summary[spec] = SummaryOfSpecies(self.GeneClasses.reverse[gclass],
                                             spec2genes[spec])

        self.addLosses(fam, summary, level)
        return summary


class FamHistory(object):
    XRefTag = None

    def __init__(self, parser, analyzer):
        self.parser = parser
        self.analyzer = analyzer

    def __len__(self):
        return self.get_number_of_fams(singletons=True)

    def setXRefTag(self, tag):
        """set the attribute name of the 'gene' elements which should
        be used for report. defined by orthoxml are 'geneId' and
        'protId'. If not defined, the (numerical) internal ids are used."""
        self.XRefTag = tag

    def analyzeLevel(self, level):
        subFamNodes = OrthoXMLQuery.getGroupsAtLevel(level, self.parser.root)
        gfamList = [GeneFamily(fam) for fam in subFamNodes]

        for gfam in gfamList:
            gfam.analyze(self.analyzer, level)

        self.geneFamList = gfamList
        self.analyzedLevel = level

    def write(self, fd, speciesFilter=None):
        """writes the FamHistory object to a given stream object
        in a human readable format.
        The stream object needs to have a write(str) method defined.
        The optional speciesFilter argument accepts a set of
        species names for which the genes in the families are
        reported."""

        formatter = lambda gid: self.parser.mapGeneToXRef(gid, self.XRefTag)
        fd.write("FamilyAnalysis at {}\n".format(self.analyzedLevel))
        for fam in self.geneFamList:
            fam.write(fd, speciesFilter, idFormatter=formatter)

    def __str__(self):
        fd = io.StringIO()
        self.write(fd)
        res = fd.getvalue()
        fd.close()
        return res

    def compare(self, other):
        """ compares two FamilyHistory objects in linear time
            algorithm implemented in Comparer class
        """
        c = Comparer(self, other)
        c.run()
        return c.comp

    def get_number_of_fams(self, singletons=False):
        if not hasattr(self, 'geneFamList'):
            return 0
        if singletons:
            return len(self.geneFamList)
        return len([x for x in self.geneFamList if not x.is_singleton()])


class Comparer(object):
    """
    Compares two FamilyHistory objects in linear time.

    Algorithm:

    Init:
    preprocess input by sorting FamilyHistory geneFamLists
    maintain pointers to the lists that only move forwards

    Run:
    If families exactly match:
        Mark as FamIdent
        advance both lists one element
    Else if family2 is a prefix match of family 1:
        While family1.prefix_match(family2):
            Mark as FamDupl
            advance list 2
        advance list 1
    Else if family 1 < family 2 (according to sort comparison)
        While family 1 < family 2:
            Mark as FamLost
            advance list 1
    Else if family 1 > family 2
        While family 1 > family 2:
            Mark as FamNovel
            advance list 2
    When list 1 is exhausted:
        mark all remaining members of list 2 as FamNovel
    When list 2 is exhausted:
        mark all remaining members of list 1 as FamLost

    End:
    When both lists are exhausted
    """

    def __init__(self, fam_history_1, fam_history_2):
        self.i1 = (x for x in sorted(fam_history_1.geneFamList)
                    if not x.getFamId() == 'n/a')
        self.i2 = (x for x in sorted(fam_history_2.geneFamList)
                    if not x.getFamId() == 'n/a')
        self.f1 = None
        self.f2 = None
        self.advance_i1()
        self.advance_i2()
        self.comp = LevelComparisonResult(fam_history_1.analyzedLevel,
            fam_history_2.analyzedLevel)

    def run(self):
        while self.f1 is not None and self.f2 is not None:
            if self.f1 == self.f2:
                self.ident()

            elif self.f1.prefix_match(self.f2):
                self.dupl()

            elif self.f1 < self.f2:
                self.lost()

            elif self.f1 > self.f2:
                self.novel()

            else:
                raise Exception('Unexpected state')

        if self.f1 is None:
            self.l1_exhausted()

        if self.f2 is None:
            self.l2_exhausted()

    def ident(self):
        self.comp.addFamily(FamIdent(self.f1.getFamId()))
        self.advance_i1()
        self.advance_i2()

    def dupl(self):
        m = list()
        while self.f1.prefix_match(self.f2):
            m.append(self.f2)
            self.advance_i2()
            if self.f2 is None:
                break
        self.comp.addFamily(FamDupl(self.f1.getFamId(),
            '; '.join(gf.getFamId() for gf in m)))
        self.advance_i1()

    def lost(self):
        while self.f1 < self.f2 and not self.f1.prefix_match(self.f2):
            self.comp.addFamily(FamLost(self.f1.getFamId()))
            self.advance_i1()
            if self.f1 is None:
                break

    def novel(self):
        while self.f1 > self.f2 and not self.f1.prefix_match(self.f2):
            _Event = (FamSingleton if self.f2.is_singleton()
                      else FamNovel) # this check is probably redundant
                                     # because there shouldn't be any
                                     # singletons if l1 is not exhausted
                                     # (singletons are always annotated last)
            self.comp.addFamily(_Event(self.f2.getFamId()))
            self.advance_i2()
            if self.f2 is None:
                break

    def advance_i1(self):
        try:
            val = next(self.i1)
        except StopIteration:
            val = None
        self.f1 = val

    def advance_i2(self):
        try:
            val = next(self.i2)
        except StopIteration:
            val = None
        self.f2 = val

    def l1_exhausted(self):
        while self.f2 is not None:
            _Event = (FamSingleton if self.f2.is_singleton()
                      else FamNovel)
            self.comp.addFamily(_Event(self.f2.getFamId()))
            self.advance_i2()

    def l2_exhausted(self):
        while self.f1 is not None:
            self.comp.addFamily(FamLost(self.f1.getFamId()))
            self.advance_i1()


class FamEvent(object):
    event = None

    def __init__(self, fam):
        self.fam = fam

    def __str__(self):
        return "{}: {}\n".format(self.fam, self.event)

    def __eq__(self, other):
        return self.fam == other.fam and self.event == other.event


class FamIdent(FamEvent):
    event = "identical"


class FamNovel(FamEvent):
    event = "novel"


class FamLost(FamEvent):
    event = "lost"


class FamSingleton(FamEvent):
    """A single-member 'family' consisting of a gene that doesn't
    match any other family. Only occurs in a leaf."""
    event = "singleton"


class FamDupl(FamEvent):
    event = "duplicated"

    def __init__(self, fam, subfam):
        super().__init__(fam)
        if isinstance(subfam, list):
            subfam = "; ".join(subfam)
        self.into = subfam

    def __str__(self):
        return "{} --> {}\n".format(self.fam, self.into)

    def __eq__(self, other):
        return super().__eq__(other) and self.into == other.into


class LevelComparisonResult(object):

    """
    Result of comparing two 'FamHistory' objects at different levels
    on the tree.

    """

    @staticmethod
    def sort_key(item):
        if item.fam == 'n/a':
            return (MAXINT,)
        return tuple((int(num) if num else alpha) for
                    (num, alpha) in re.findall(r'(\d+)|(\D+)', item.fam))

    def __init__(self, lev1, lev2):
        self.fams = list()
        self.lev1 = lev1
        self.lev2 = lev2

    def __str__(self):
        fd = io.StringIO()
        self.write(fd)
        res = fd.getvalue()
        fd.close()
        return res

    def addFamily(self, famEvent):
        self.fams.append(famEvent)

    def write(self, fd):
        fd.write("\nLevelComparisonResult between taxlevel {} and {}\n".
                 format(self.lev1, self.lev2))
        self.fams.sort(key=self.sort_key)
        for fam in self.fams:
            fd.writelines(str(fam))


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
                levelsToParent = {l for l in
                                  self.tax.iterParents(mostSpecificLevel,
                                                       last)}
                levelsToParent.add(mostSpecificLevel)
                if not levels.issubset(levelsToParent):
                    raise Exception("taxonomy not in correspondance with found"
                                    " hierarchy: {} vs {}".format(
                                                                levels,
                                                                levelsToParent))
            addLevels = levelsToParent - levels
            for lev in addLevels:
                node.append(self._createTaxRangeTag(lev))
            for child in list(node):
                self._addTaxRangeR(child, mostSpecificLevel)

        elif self.parser.is_paralog_group(node):
            for child in list(node):
                self._addTaxRangeR(child, last)
        elif OrthoXMLQuery.is_geneRef_node(node):
            # to simplify analyses down to the taxlevel of a single species
            # we add an aditional fake orthologGroup just above each geneRef
            # element with all the taxRanges between the most specific level
            # of the parent orthologGroup node and the species itself.
            spec = self.parser.mapGeneToSpecies(node.get('id'))
            expRange = self.tax.hierarchy[spec].name
            directParent = parent = node.getparent()
            while not self.parser.is_ortholog_group(parent):
                parent = parent.getparent()
            levOfParent = OrthoXMLQuery.getLevels(parent)
            mostSpecific = self.tax.mostSpecific(levOfParent)
            self._insertOG(directParent, node, expRange, mostSpecific)

    def _insertOG(self, parent, child, specificLev, beforeLev):
        pos = parent.index(child)
        el = etree.Element('{{{ns0}}}orthologGroup'.format(**self.parser.ns))
        el.append(self._createTaxRangeTag(specificLev))
        for lev in self.tax.iterParents(specificLev, stopBefore=beforeLev):
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

        top_level_groups = self.parser.getToplevelGroups()

        if PROGRESSBAR:
            pbar = setup_progressbar(
                'Adding missing taxonomy annotation: ',
                len(top_level_groups)
                )
            pbar.start()

        for i, fam in enumerate(top_level_groups, start=1):
            self._addTaxRangeR(fam, noUpwardLevels=not propagate_top)
            if PROGRESSBAR:
                pbar.update(i)

        if PROGRESSBAR:
            pbar.finish()

        del self.tax

    def annotateDoc(self):
        """apply the LOFT naming schema to all the orthologGroups."""
        for i, fam in enumerate(self.parser.getToplevelGroups()):
            self.dupCnt = list()
            self._annotateGroupR(fam, fam.get('id', str(i)))

    def annotateSingletons(self):
        """Any input genes that aren't assigned to ortholog groups are
        singletons, which are added to the xml as extra ortholog groups"""
        if PROGRESSBAR:
            pbar = setup_progressbar('Adding singletons: ', 1)
            pbar.start()

        highest_group = max(self.parser.getToplevelGroups(),
                            key=lambda x: int(x.get('id')))
        input_genes = set(n.get('id') for n in
                          OrthoXMLQuery.getInputGenes(self.parser.root))
        grouped_genes = set(n.get('id') for n in
                            OrthoXMLQuery.getGroupedGenes(self.parser.root))
        singletons = input_genes - grouped_genes
        groups_node = OrthoXMLQuery.getSubNodes('groups', self.parser.root)[0]

        fam_num = int(highest_group.get('id')) + 1
        singleton_families = set()

        if PROGRESSBAR:
            pbar.maxval = len(singletons)

        for i, gene in enumerate(singletons, start=1):
            singleton_families.add(str(fam_num))
            species = self.parser.mapGeneToSpecies(gene)
            new_node = etree.Element('{{{ns0}}}orthologGroup'.format(
                                                            **self.parser.ns),
                                                            id=str(fam_num))
            new_node.append(self._createTaxRangeTag(species))
            new_node.append(etree.Element('{{{ns0}}}geneRef'.format(
                                                            **self.parser.ns),
                                                            id=gene))
            groups_node.append(new_node)
            fam_num += 1
            if PROGRESSBAR:
                pbar.update(i)

        self.parser.singletons = singleton_families
        if PROGRESSBAR:
            pbar.finish()


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Analyze Hierarchical OrthoXML families.')
    parser.add_argument('--xreftag', default=None,
                        help=("xref tag of genes to report. OrthoXML allows to "
                              "store multiple ids and xref annotations per gene "
                              "as attributes in the species section. If not set, "
                              "the internal (purely numerical) ids are reported."))
    parser.add_argument('--show_levels', action='store_true',
                        help='print the levels and species found in the orthoXML file and quit')
    parser.add_argument('-r', '--use-recursion', action='store_true',
                        help=("DEPRECATED: Use recursion to sample families that are a "
                              "subset of the query"))
    parser.add_argument('--taxonomy', default='implicit',
                        help=("Taxonomy used to reconstruct intermediate levels. "
                              "Has to be either 'implicit' (default) or a path to "
                              "a file in Newick format. The taxonomy might be "
                              "multifurcating. If set to 'implicit', the "
                              "taxonomy is extracted from the input OrthoXML file. "
                              "The orthoXML level do not have to cover all the "
                              "levels for all families. In order to infer gene losses "
                              "Family-Analyzer needs to infer these skipped levels "
                              "and reconcile each family with the complete taxonomy."))
    parser.add_argument('--propagate_top', action='store_true',
                        help=("propagate taxonomy levels up to the toplevel. As an "
                              "illustration, consider a gene family in an eukaryotic "
                              "analysis that has only mammalian genes. Its topmost "
                              "taxonomic level will therefor be 'Mammalia' and an "
                              "ancestral gene was gained at that level. However, if "
                              "'--propagete-top' is set, the family is assumed to have "
                              "already be present in the topmost taxonomic level, i.e. "
                              "Eukaryota in this example, and non-mammalian species "
                              "have all lost this gene."))
    parser.add_argument('--show_taxonomy', action='store_true',
                        help='write the taxonomy used to standard out. ')
    parser.add_argument('--store_augmented_xml', default=None,
                        help=("filename to which the input orthoxml file with "
                              "augmented annotations is written. The augmented "
                              "annotations include for example the additional "
                              "taxonomic levels of orthologGroup and unique HOG "
                              "IDs."))
    parser.add_argument('--add_singletons', action='store_true',
                        help=("Take singletons - genes from the input set that "
                              "weren't assigned to orthologous groups and "
                              "appear as novel gains in the taxonomy leaves - "
                              "and add them to the xml as single-member "
                              "orthologGroups"))
    parser.add_argument('--compare_second_level', default=None,
                        help=("Compare secondary level with primary one, i.e. "
                              "report what happend between the secondary and primary "
                              "level to the individual histories. Note that the "
                              "Second level needs to be younger than the primary."))
    parser.add_argument('orthoxml', help='path to orthoxml file to be analyzed')
    parser.add_argument('level', help='taxonomic level at which analysis should be done')
    parser.add_argument('species', nargs="+", help=("(list of) species to be analyzed. "
                              "Note that only genes of the selected species are "
                              "reported. In order for the output to make sense, "
                              "the selected species all must be part of the "
                              "linages specified in 'level' (and --compare_second_level)."))
    args = parser.parse_args()

    op = OrthoXMLParser(args.orthoxml)
    if args.show_levels:
        print("Species:\n{0}\n\nLevels:\n{1}".format(
              '\n'.join(sorted(list(op.getSpeciesSet()))),
              '\n'.join(sorted(op.getLevels()))))
        sys.exit()
    print("Analyzing {} on taxlevel {}".format(args.orthoxml, args.level))
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

    if args.add_singletons:
        op.augmentSingletons()

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
        comp = hist.compare(hist2)
        comp.write(sys.stdout)

    if args.store_augmented_xml is not None:
        op.write(args.store_augmented_xml)

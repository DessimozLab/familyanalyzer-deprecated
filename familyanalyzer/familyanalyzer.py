from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future.builtins import super
from future.builtins import next
from future.builtins import chr
from future.builtins import dict
from future.builtins import int
from future.builtins import str
from future import standard_library

standard_library.install_hooks()

# !/usr/bin/env python
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
from .tools import enum, PROGRESSBAR, setup_progressbar
from .orthoxmlquery import ElementError, OrthoXMLQuery
from .taxonomy import NewickTaxonomy, TaxRangeOrthoXMLTaxonomy, XMLTaxonomy

MAXINT = sys.maxsize


class OrthoXMLParser(object):
    ns = {"ns0": "http://orthoXML.org/2011/"}  # xml namespace

    def __init__(self, filename):
        """creates a OrthoXMLParser object. the parameter filename
        needs to be a path pointing to the orthoxml file to be
        analyzed."""
        self.doc = etree.parse(filename)
        self.root = self.doc.getroot()
        self.tax = None
        self.singletons = None

        self._buildMappings()  # builds three dictionaries - see def below

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

    def mapGeneToXRef(self, id_, typ='protId'):
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

    def get_species_below_node(self, node, return_gene_total_count=False):
        """ return a set of all species that have a geneRef present beneath
        the specified node

        :param node: Node of interest
        :param return_gene_total_count: boolean flag whether or not to return
               the number of genes below this node.
        :return:
        """
        generef_nodes = OrthoXMLQuery.getGeneRefNodes(node)
        species_covered = {self.mapGeneToSpecies(gr.get('id'))
                           for gr in generef_nodes}
        if return_gene_total_count:
            return species_covered, len(generef_nodes)
        else:
            return species_covered

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
                minCoverage * len(self.getSpeciesSet())]

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

    def __str__(self):
        return self.root.get('og')

    def prefix_match(self, other):
        query = self._cmp()
        target = other._cmp()
        return target[:len(query)] == query

    def getMemberGenes(self):
        """
        Get all genes belonging to this family. Cached to speed up repeated calls.
        """
        if hasattr(self, '_member_genes'):
            return self._member_genes
        else:
            members = self.root.findall('.//{{{ns0}}}geneRef'.
                                        format(**OrthoXMLParser.ns))
            self._member_genes = [x.get('id') for x in members]
            return self._member_genes

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

    def analyzeLevel(self, level):
        return self

    def analyze(self, strategy):
        super().analyze(strategy)
        for sumElement in self.summary.values():
            sumElement.typ = "SINGLETON"


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
            speciesCoveredByLevel = self.tax.descendents[mostGeneralLevel]
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
        famId = fam.getFamId()
        for spec, set_ in spec2genes.items():
            if len(set_) > 1:
                gclass = self.GeneClasses.MULTICOPY
            else:
                if famId in self.parser.singletons:
                    gclass = self.GeneClasses.SINGLETON
                else:
                    gclass = self.GeneClasses.SINGLECOPY

            summary[spec] = SummaryOfSpecies(self.GeneClasses.reverse[gclass],
                                             set_)

        self.addLosses(fam, summary, level)
        return summary


class FamHistory(object):
    XRefTag = None

    def __init__(self, parser, analyzer):
        self.parser = parser
        self.analyzer = analyzer
        self._geneFamList = list()

    def __len__(self):
        return self.get_number_of_fams(singletons=True)

    def __getitem__(self, key):
        return self.geneFamDict[key]

    def __iter__(self):
        return iter(self.geneFamList)

    @property
    def geneFamList(self):
        return self._geneFamList

    @geneFamList.setter
    def geneFamList(self, gfam_list):
        self._geneFamList = sorted(gfam_list)

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

        self.geneFamDict = {gf.getFamId(): gf for gf in gfamList}
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
        for fam in self:
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
        if singletons:
            return len(self.geneFamList)
        return len([x for x in self if not x.is_singleton()])


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
        self.i1 = iter(fam_history_1.geneFamList)
        self.i2 = iter(fam_history_2.geneFamList)
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
        self.comp.addFamily(FamDupl(self.f1, m))
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
                      else FamNovel)  # this check is probably redundant
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
        self.name = str(fam)
        self.fam = fam

    def __str__(self):
        return "{}: {}\n".format(self.name, self.event)

    def __eq__(self, other):
        return self.name == other.name and self.event == other.event

    def __repr__(self):
        return "{}: {}".format(self.name, self.event)


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
            subfam_names = "; ".join([str(s) for s in subfam])
        else:
            subfam_names = str(subfam)
        self.into = subfam_names
        self.subfams = subfam

    def __str__(self):
        return self.write()

    def __eq__(self, other):
        return super().__eq__(other) and self.into == other.into

    def write(self):
        ''' Construct output for printing. Shows the ids of the duplicated
            GeneFamily and the resulting subfamilies. It also lists the
            members of the original GeneFamily and the new subfamilies.
            Sample output:
            -------
            114 --> 114.2a; 114.2b
            114: 47953; 28082; 11418; 29862; 117; 50097; 13634; 50845; 14300
            114.2a: 50097; 13634
            114.2b: 50845; 14300
            -------
        '''
        output = "-------\n"
        output += ("{} --> {}\n".format(self.fam.getFamId(), self.into))
        members = "; ".join(self.fam.getMemberGenes())
        output += ("{}: {}\n".format(self.fam.getFamId(), members))
        for group in self.subfams:
            members = "; ".join(group.getMemberGenes())
            output += ("{}: {}\n".format(group.getFamId(), members))
        output += ("-------\n")
        return output


class LevelComparisonResult(object):
    """
    Result of comparing two 'FamHistory' objects at different levels
    on the tree.

    """

    groups = {'identical': 0,
              'duplicated': 1,
              'lost': 2,
              'novel': 3,
              'singleton': 4}

    groups_back = {0: 'identical',
                   1: 'duplicated',
                   2: 'lost',
                   3: 'novel',
                   4: 'singleton'}

    @staticmethod
    def sort_key(item):

        if item.name == 'n/a':
            return (MAXINT,)

        return tuple((int(num) if num else alpha) for
                     (num, alpha) in re.findall(r'(\d+)|(\D+)', item.name))

    def group_sort_key(self, item):
        return tuple(itertools.chain((self.group_key(item),),
                                     self.sort_key(item)))

    def group_key(self, item):
        return self.groups[item.event]

    def __init__(self, lev1, lev2):
        self.fams_dict = dict()
        self.lev1 = lev1
        self.lev2 = lev2

    def __str__(self):
        fd = io.StringIO()
        self.write(fd)
        res = fd.getvalue()
        fd.close()
        return res

    def __getitem__(self, key):
        return self.fams_dict[key]

    def __iter__(self):
        return iter(self.fams)

    @property
    def fams(self):
        return sorted(self.fams_dict.values(), key=self.sort_key)

    def addFamily(self, famEvent):
        self.fams_dict[famEvent.name] = famEvent

    def write(self, fd):
        fd.write("\nLevelComparisonResult between taxlevel {} and {}\n".
                 format(self.lev1, self.lev2))
        for fam in self:
            fd.writelines(str(fam))

    def summarise(self):
        summary = {'identical': 0,
                   'novel': 0,
                   'lost': 0,
                   'duplicated': 0,
                   'singleton': 0}

        for family in self.fams_dict.values():
            if family.event in ['identical', 'lost', 'singleton', 'novel']:
                summary[family.event] += 1
            elif family.event == 'duplicated':
                summary['duplicated'] += len(family.into.split('; '))

        self.summary = summary
        return summary

    def filter(self, filters):
        if not isinstance(filters, set):
            filters = {filters}
        if not filters.issubset({'identical', 'lost', 'singleton', 'novel',
                                 'duplicated'}):
            raise Exception('Unexpected filters: {0}'.format(filters))
        return [x for x in self if x.event in filters]

    def group_fams(self):
        return dict([(self.groups_back[num], list(iterator))
                     for (num, iterator)
                     in itertools.groupby(sorted(self.fams_dict.values(),
                                                 key=self.group_sort_key),
                                          self.group_key)])


class GroupAnnotator(object):
    """this class annotates orthologGroup elements with the numbering
    schema presented in the LOFT paper:
    van der Heijden, Snel, van Noort, Huynen
    Orthology prediction at scalable resolution by phylogenetic tree analysis.
    BMC Bioinformatics, 2007, 8, 83

    and adding additional property tags with skipped TaxRange levels."""

    def __init__(self, parser, taxrange_2_taxid=None):
        self.parser = parser
        self.ns = parser.ns
        self.taxrange_2_taxid = {rng: str(val) for rng, val in taxrange_2_taxid.items()} if taxrange_2_taxid is not None else {}

    def _getNextSubId(self, idx):
        """helper method to return the next number at a given depth of
        duplication (idx)"""
        while len(self.dupCnt) < idx:
            self.dupCnt.append(0)
        self.dupCnt[idx - 1] += 1
        return self.dupCnt[idx - 1]

    def _encodeParalogClusterId(self, prefix, nr):
        """encode the paralogGroups at the same level, e.g. 1a, 1b, 1c
        for 3 paralogGroups next to each other. the nr argument
        identifies the individual indices of those 3 paralogGroups."""
        letters = []
        while nr // 26 > 0:
            letters.append(chr(97 + (nr % 26)))
            nr = nr // 26 - 1
        letters.append(chr(97 + (nr % 26)))
        return prefix + ''.join(letters[::-1])  # letters were in reverse order

    def _annotateGroupR(self, node, og, idx=0):
        """create the og attributes at the orthologGroup elements
        according to the naming schema of LOFT. ParalogGroup elements
        do not get own attributes (not possible in the xml schema),
        but propagate their sub-names for the subsequent orthologGroup
        elements."""
        if self.parser.is_ortholog_group(node):
            taxid_node = OrthoXMLQuery.getTaxidNodes(node, recursively=False)
            if len(taxid_node) > 0:
                node.set('id', '{}_{}'.format(og, taxid_node[0].get('value')))
            else:
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

    def _addTaxRangeR(self, node, noUpwardLevels=False):
        """recursive method to add TaxRange property tags."""
        if self.parser.is_ortholog_group(node) or self.parser.is_paralog_group(node) \
                or OrthoXMLQuery.is_geneRef_node(node):
            species_covered, nr_genes = self.parser.get_species_below_node(node, return_gene_total_count=True)
            current_level = self.tax.mrca(species_covered)
            og_tag = '{{{}}}orthologGroup'.format(OrthoXMLQuery.ns['ns0'])

            if self.parser.is_ortholog_group(node):
                comp_score = OrthoXMLQuery.getScoreNodes(node, 'CompletenessScore')
                if len(comp_score) == 0:
                    node.append(self._createCompletnessScoreTag(current_level, species_covered))
                node.append(self._createNrMemberGeneTag(nr_genes))
                taxrange = OrthoXMLQuery.getTaxRangeNodes(node, False)
                taxid = OrthoXMLQuery.getTaxidNodes(node, False)
                if len(taxrange) > 0:
                    # check consistency between current_level and value stored in taxrange
                    if taxrange[0].get('value') != current_level:
                        raise Exception("Inconsistent TaxRange: {} vs current_level {}"
                                        .format(taxrange[0].get('value'), current_level))
                    if len(taxid) > 0:
                        if taxid[0].get('value') != self.taxrange_2_taxid[current_level]:
                            raise Exception("Inconsitency between taxids: {} vs {}"
                                            .format(taxid[0].get('value'), self.taxrange_2_taxid[current_level]))
                    else:
                        try:
                            node.append(self._create_taxid(current_level))
                        except KeyError:
                            pass
                else:
                    node.append(self._createTaxRangeTags(current_level))

            try:  # find the closest ancestral orthogroup if it has a TaxRange property
                parent_orthogroup = next(node.iterancestors(og_tag))
                parent_levels = {z.get('value')
                                 for z in OrthoXMLQuery.getTaxRangeNodes(parent_orthogroup, False)}
            except StopIteration:  # couldn't find a parent with a TaxRange property; no extra annotation possible
                parent_levels = set([])

            if len(parent_levels) > 0:
                most_recent_parent_level = self.tax.mostSpecific(parent_levels)

                # Ortholog Node - append missing tax range(s) as property tags under the current node
                if self.parser.is_ortholog_group(node):
                    self._insertOGs_between(node.getparent(), node, current_level, most_recent_parent_level, nr_genes,
                                            species_covered, include_self=False)

                # Paralog Node - insert ortholog node between self and parent; add missing tax range(s) to new parent
                elif self.parser.is_paralog_group(node):
                    if self.tax.levels_between(most_recent_parent_level, current_level) > 1:
                        self._insertOGs_between(node.getparent(), node, current_level, most_recent_parent_level,
                                                nr_genes, species_covered, include_self=False)

                # GeneRef Node - insert ortholog node between self and parent; add all tax range(s) to new parent
                else:
                    self._insertOGs_between(node.getparent(), node, current_level, most_recent_parent_level,
                                            nr_genes, species_covered, include_self=True)
                    return

            for child in node:
                self._addTaxRangeR(child, noUpwardLevels)

    def _insertOGs_between(self, parent, child, specificLev, beforeLev, nr_genes, covered_species, include_self=True):
        pos = parent.index(child)
        if include_self:
            child = self._insert_one_OG(child, specificLev, covered_species, nr_genes)
        for lev in self.tax.iterParents(specificLev, stopBefore=beforeLev):
            child = self._insert_one_OG(child, lev, covered_species=covered_species, nr_genes=nr_genes)
        parent.insert(pos, child)

    def _insert_one_OG(self, child, level, covered_species, nr_genes):
        el = etree.Element('{{{ns0}}}orthologGroup'.format(**self.parser.ns))
        el.append(self._createCompletnessScoreTag(level, covered_species))
        el.extend(self._createTaxRangeTags(level))
        el.append(self._createNrMemberGeneTag(nr_genes))
        el.append(child)
        return el

    def _createTaxRangeTags(self, lev, **kwargs):
        tags = [self._create_tax_range(lev, **kwargs)]
        try:
            tags.append(self._create_taxid(lev))
        except KeyError:
            pass
        return tags

    def _create_tax_range(self, lev, **kwargs):
        return etree.Element('{{{ns0}}}property'.format(**self.parser.ns),
                             attrib=dict(name='TaxRange', value=lev, **kwargs))

    def _create_taxid(self, lev):
        return etree.Element('{{{ns0}}}property'.format(**self.parser.ns),
                             attrib=dict(name='taxid',
                                         value=self.taxrange_2_taxid[lev]))

    def _createNrMemberGeneTag(self, nr_genes):
        return etree.Element('{{{ns0}}}property'.format(**self.parser.ns),
                             attrib={"name": "NrMemberGenes", "value": str(nr_genes)})

    def completenessScore(self, level, covered_species):
        return "{:.3f}".format(len(covered_species) / len(self.tax.descendents[level]))

    def _createCompletnessScoreTag(self, level, covered_species):
        el = etree.Element('{{{ns0}}}score'.format(**self.parser.ns),
                           attrib={"id": "CompletenessScore",
                                   "value": self.completenessScore(level, covered_species)})
        return el

    def annotateMissingTaxRanges(self, tax, propagate_top=False, verbosity=0):
        """This function adds left-out taxrange property elements to
        the orthologGroup elements in the xml. It will add all the levels
        defined in the 'tax'-Taxonomy between the parents most specific
        level and the current nodes level. If no parent exists, all
        tax-levels above the current one are used."""
        self.tax = tax

        top_level_groups = self.parser.getToplevelGroups()

        if PROGRESSBAR and verbosity > 0:
            pbar = setup_progressbar(
                'Adding missing taxonomy annotation: ',
                len(top_level_groups)
            )
            pbar.start()

        for i, fam in enumerate(top_level_groups, start=1):
            self._addTaxRangeR(fam, noUpwardLevels=not propagate_top)
            if PROGRESSBAR and verbosity > 0:
                pbar.update(i)

        if PROGRESSBAR and verbosity > 0:
            pbar.finish()

        del self.tax

    def annotateDoc(self):
        """apply the LOFT naming schema to all the orthologGroups."""
        for i, fam in enumerate(self.parser.getToplevelGroups()):
            self.dupCnt = list()
            self._annotateGroupR(fam, fam.get('id', str(i)))

    def annotateSingletons(self, verbosity=0):
        """Any input genes that aren't assigned to ortholog groups are
        singletons, which are added to the xml as extra ortholog groups"""
        if PROGRESSBAR and verbosity > 0:
            pbar = setup_progressbar('Adding singletons: ', 1)
            pbar.start()

        highest_group = max(self.parser.getToplevelGroups(),
                            key=lambda x: int(
                                x.get('id')))  # TODO: If top-level orthologNodes have no id field this errors out
        input_genes = set(n.get('id') for n in  # Maybe add code to enumerate OG nodes if ids are missing?
                          OrthoXMLQuery.getInputGenes(self.parser.root))
        grouped_genes = set(n.get('id') for n in
                            OrthoXMLQuery.getGroupedGenes(self.parser.root))
        singletons = input_genes - grouped_genes
        groups_node = OrthoXMLQuery.getSubNodes('groups', self.parser.root)[0]

        fam_num = int(highest_group.get('id')) + 1
        singleton_families = set()

        if PROGRESSBAR and verbosity > 0:
            pbar.maxval = len(singletons)

        for i, gene in enumerate(sorted(singletons), start=1):
            singleton_families.add(str(fam_num))
            species = self.parser.mapGeneToSpecies(gene)
            new_node = etree.Element('{{{ns0}}}orthologGroup'.format(
                **self.parser.ns),
                id=str(fam_num))
            new_node.extend(self._createTaxRangeTags(species))
            new_node.append(etree.Element('{{{ns0}}}geneRef'.format(
                **self.parser.ns),
                id=gene))
            groups_node.append(new_node)
            fam_num += 1
            if PROGRESSBAR and verbosity > 0:
                pbar.update(i)

        self.parser.singletons = singleton_families
        if PROGRESSBAR and verbosity > 0:
            pbar.finish()

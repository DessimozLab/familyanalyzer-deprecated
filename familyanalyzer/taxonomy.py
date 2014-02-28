import io
import itertools
import os
import re
from .orthoxmlquery import OrthoXMLQuery
from .newick import NewickLexer, Streamer
from .tools import PROGRESSBAR, setup_progressbar


class TaxonomyInconsistencyError(Exception):
    pass


class ParseError(Exception):
    pass


class Taxonomy(object):
    def __init__(self):
        raise NotImplementedError("abstract class")

    def __iter__(self):
        return self.hierarchy[self.root].iterDescendents()

    def iterParents(self, node, stopBefore=None):
        """iterates over all the taxonomy nodes towards the root
        which are above 'node' and below 'stopBefore'."""

        if node == stopBefore:
            return
        tn = self.hierarchy[node]
        while tn.up is not None and tn.up.name != stopBefore:
            tn = tn.up
            yield tn.name

    def is_ancestor_of(self, anc, desc):
        """ Returns True if `anc' is an ancestor of `desc'"""
        return anc in self.iterParents(desc)

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

    def get_node(self, node_name):
        return self.hierarchy[node_name]

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


class MiniTaxonomy(Taxonomy):
    """Builds a taxonomy from three FamHistory objects, which are in
    a taxonomic hierarchy - history1 is ancestral to history2, which is
    ancestral to history3. This can be checked against an external Taxonomy
    using the `check_against_taxonomy` method. Comparison1 compares history1
    to history2, and comparison2 compares history2 to history3. These can be
    supplied or generated."""


    def __init__(self, history1, history2, history3, comparison1=None,
                 comparison2=None):

        if comparison1 is None:
            comparison1 = history1.compare(history2)
        if comparison2 is None:
            comparison2 = history2.compare(history3)

        self.check_comparison_consistency(comparison1, comparison2)

        bottom = TaxNode(history1.analyzedLevel)
        middle = TaxNode(history2.analyzedLevel)
        top    = TaxNode(history3.analyzedLevel)

        bottom.attachFamHistory(history1)
        middle.attachFamHistory(history2)
        top.attachFamHistory(history3)

        middle.attachLevelComparisonResult(comparison1)
        top.attachLevelComparisonResult(comparison2)

        bottom.addChild(middle)
        middle.addChild(top)

        middle.addParent(bottom)
        top.addParent(middle)

        self.hierarchy = dict()
        self.root = bottom.name
        self.hierarchy = dict((node.name, node)
                              for node in bottom.iterDescendents())

    def check_against_taxonomy(self, tax):
        label1, label2, label3 = (n.history.analyzedLevel for n in self)
        if not {label1, label2}.issubset(tax.iterParents(label3)):
            raise TaxonomyInconsistencyError('Taxonomic hierarchy is violated')

    def check_comparison_consistency(self, comparison1, comparison2):
        if not comparison1.lev2 == comparison2.lev1:
            raise TaxonomyInconsistencyError('Comparisons don\'t overlap')


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
        if self.parser.is_ortholog_group(grp):
            levels = [l.get('value') for l in grp.findall(
                './{{{ns0}}}property[@name="TaxRange"]'
                .format(**self.parser.ns))]
        directChildNodes = list(grp)
        children = [child for child in directChildNodes
                    if self.parser.is_evolutionary_node(child)]
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

        NHX = self._get_node_NHX() + self._get_edge_NHX()
        if self.isLeaf():
            return '{0}{1}'.format(label,
                            ('[&&NHX{0}]'.format(NHX) if NHX > '' else ''))
        subtree = ', '.join(str(ch) for ch in self.down)
        return '({0}){1}{2}'.format(subtree, label,
                            ('[&&NHX{0}]'.format(NHX) if NHX > '' else ''))

    def _get_node_NHX(self):
        return (':Genes={}'.format(len(self.history))
                if self.history else '')

    def _get_edge_NHX(self):
        if self.comparison:
            if not hasattr(self.comparison, 'summary'):
                self.comparison.summarise()

            n_ident     = self.comparison.summary['identical']
            n_dupl      = self.comparison.summary['duplicated']
            n_lost      = self.comparison.summary['lost']
            n_novel     = self.comparison.summary['novel']
            n_singleton = self.comparison.summary['singleton']

            NHX = (':Identical={0}'
                   ':Duplicated={1}'
                   ':Lost={2}'.format(n_ident,
                                      n_dupl,
                                      n_lost))
            if self.isLeaf():
                NHX += ':Singleton={}'.format(n_singleton)
            else:
                NHX += ':Novel={}'.format(n_novel)

            return NHX
        return ''

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


class NewickTaxonomy(Taxonomy):

    """ Create a taxonomy from a file in newick format. The file should contain
    one tree (further trees are ignored). Only the tree topology is used -
    branch lengths and bootstrap support values are thown away.
    The leaf labels should match those in the orthoXML. Inner labels
    should match too, but for OMA XML will be automatically generated if
    auto_annotate == True """

    def __init__(self, filename):
        if not os.path.exists(filename):
            raise Exception('File not found: {0}'.format(filename))
        self.lexer = NewickLexer(Streamer(open(filename)))
        self.nodes = set()
        self.hierarchy = {}
        self.stack = []
        self.parse()

    def _get_label(self, tokens):
        """ Get the node data attributes 'label' and 'length'. Assumes these
        will be the next tokens in the stream. Throws ParseError if they are
        not. """
        label = next(self.lexer)
        if label.typ not in (tokens.LABEL, tokens.SUPPORT):
            raise ParseError(
                'Expected a label or a support value, found {0}'.format(
                    label))

        length = next(self.lexer)
        if length.typ != tokens.LENGTH:
            raise ParseError('Expected a length, found {0}'.format(
                length))

        return (label.val if label.typ == tokens.LABEL else None)

    def annotate_from_orthoxml(self, xmlparser):
        """ Transfers internal node names from OMA orthoxml """
        def _annotate(self, node, levels_dict):
            if not node.down:
                self.nodes.add(node.name)
                self.hierarchy[node.name] = node
                return
            key = frozenset((n.name for n in node.iterLeaves()))
            node.name = levels_dict[key]
            self.nodes.add(node.name)
            self.hierarchy[node.name] = node
            for child in node.down:
                _annotate(self, child, levels_dict)

        levels = xmlparser.getLevels()
        levels_dict = dict((frozenset(x.split('/')), x) for x in levels)
        root_node = self.hierarchy[self.root]
        root_node.name = 'LUCA'
        self.hierarchy = {'LUCA': root_node}
        self.nodes = set('LUCA')
        self.root = 'LUCA'

        for child in root_node.down:
            _annotate(self, child, levels_dict)

    def populate(self, root):
        if len(root.down) == 0:
            self.hierarchy[root.name] = root
            self.nodes.add(root.name)
            return

        if not root.up:
            self.root = (root.name or 'LUCA')

        self.hierarchy[root.name] = root
        self.nodes.add(root.name)

        for child in root.down:
            self.populate(child)

    def parse(self):
        tmp_name = 1
        tokens = self.lexer.tokens

        for token in self.lexer:
            if token.typ == tokens.EOF:
                return

            elif token.typ == tokens.TREE:
                # invent a name for the node
                n = TaxNode(str(tmp_name))
                tmp_name += 1

                # push it onto the stack
                self.stack.append(n)

                # set as root
                self.root = n

            elif token.typ == tokens.SUBTREE:
                # invent a name and make a node
                n = TaxNode(str(tmp_name))
                tmp_name += 1

                # get parent
                p = self.stack[-1]
                p.addChild(n)
                n.addParent(p)

                # push subtree onto the stack
                self.stack.append(n)

            elif token.typ == tokens.LEAF:
                label = self._get_label(tokens)
                self.nodes.add(label)
                l = TaxNode(label)

                # get parent from stack
                p = self.stack[-1]
                p.addChild(l)
                l.addParent(p)

            elif token.typ == tokens.ENDSUB:
                label = self._get_label(tokens)

                # retrieve node from stack
                subtree = self.stack.pop()

                # update node name
                if isinstance(label, str):
                    subtree.name = label

            elif token.typ == tokens.ENDTREE:  # trigger for tree-finalising functions
                self.populate(self.root)
                return

            elif token.typ in (tokens.LABEL, tokens.LENGTH, tokens.SUPPORT):
                raise ParseError('Unexpected token in stream: {0}'.format(token))

            else:
                raise ParseError('Not sure what happened')

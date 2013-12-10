#!/usr/bin/env python

from familyanalyzer.newick import NewickTaxonomy
from familyanalyzer.familyanalyzer import *


"""
Work in progress - will merge into FamilyAnalyzer
"""

def get_histories(op, tax):

    histories = {}
    for level in tax.hierarchy:
        node = tax.hierarchy[level]
        if len(node.down) > 0:
            history = op.getFamHistory()
            history.analyzeLevel(level)
            histories[level] = history

    return histories

def get_comparisons(op, tax, histories=None):

    if histories is None:
        histories = get_histories(op, tax)

    comparisons = {}
    root_node = tax.hierarchy[tax.root]
    for node in root_node.iterInnerNodes():
        children = node.down
        for child in children:
            # print ('Curr = {}, Desc = {} [{}]'.format(node.name, child.name, ('Leaf' if child.isLeaf() else 'Inner')))
            if child.isLeaf():
                comp = histories[node.name].compareLeaf(child.name)
            else:
                comp = histories[node.name].compareFast(histories[child.name])
            comparisons[(node.name, child.name)] = comp
    return comparisons


class CompEdge(object):

    """
    Connects gene family objects to each other
    Stores  tail_id, tail_level,
            head_id, head_level,
            type
    """
    pass

class CompNode(object):
    def __init__(self, level, gene_family):
        self.level = level
        self.gene_family = gene_family
        self.parental_nodes = dict()
        self.child_nodes = dict()

    def add_parental_node(self, level, gene_family):
        parent_node = CompNode(level, gene_family)
        self.parental_nodes[level] = parent_node

    def add_child_node(self, level, gene_family):
        child_node = CompNode(level, gene_family)
        self.child_nodes.setdefault(level, list()).append(child_node)




class GeneFamilyDummy(object):

    def __init__(self, famId):
        self.famId = famId

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

    def getFamId(self):
        return self.famId

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



op = OrthoXMLParser('/homes/kgori/projects/ferret/family_analyzer/XML/augmented.xml')
tax = NewickTaxonomy('/homes/kgori/projects/ferret/family_analyzer/taxonomy.nwk')
tax.annotate_from_orthoxml(op)

placental = op.getFamHistory()
boreotheria = op.getFamHistory()

placental.analyzeLevel('GORGO/NOMLE/OTOGA/PANTR/MUSPF/RABIT/TARSY/FELCA/RATNO/PIGXX/CANFA/MOUSE/CAVPO/AILME/MACMU/ODORO/HORSE/ERIEU/CALJA/LOXAF/PONAB/HUMAN/MICMU/BOVIN')
boreotheria.analyzeLevel('GORGO/NOMLE/OTOGA/PANTR/MUSPF/RABIT/TARSY/FELCA/RATNO/PIGXX/CANFA/MOUSE/CAVPO/AILME/MACMU/ODORO/HORSE/ERIEU/CALJA/PONAB/HUMAN/MICMU/BOVIN')

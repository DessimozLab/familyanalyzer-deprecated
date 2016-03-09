from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future.builtins import range
from future.builtins import dict
from future.builtins import int
from future import standard_library
standard_library.install_hooks()

import unittest
import familyanalyzer as fa
fa.PROGRESSBAR = False
fa.familyanalyzer.PROGRESSBAR = False


class SetupHelper(object):
    @staticmethod
    def createOrthoXMLParserFromSimpleEx():
        filename = "test/simpleEx.orthoxml"
        return fa.OrthoXMLParser(filename)


class OrthoXMLParserTest(unittest.TestCase):
    def setUp(self):
        self._op = SetupHelper.createOrthoXMLParserFromSimpleEx()

    def test_nrOfToplevelFamilies(self):
        self.assertEqual(len(self._op.getToplevelGroups()), 3)

    def test_returnedSpecies(self):
        expectedSpecies = {'HUMAN', 'PANTR', 'MOUSE', 'RATNO',
                           'CANFA', 'XENTR'}
        self.assertSetEqual(self._op.getSpeciesSet(),
                            expectedSpecies)

    def test_numberOfGenesPerSpecies(self):
        expectedCnts = dict(HUMAN=4, PANTR=4, MOUSE=4, RATNO=2,
                            CANFA=3, XENTR=2)
        allGenes = self._op.getGeneIds()
        for species in expectedCnts.keys():
            self.assertEqual(
                len([gid for gid in allGenes
                    if self._op.mapGeneToSpecies(gid) == species]),
                expectedCnts[species],
                "number of genes not correct for "+species)

    def test_numberOfGenesPerSpecies_specFilter(self):
        expectedCnts = dict(HUMAN=4, PANTR=4, MOUSE=4, RATNO=2,
                            CANFA=3, XENTR=2)
        param_list = [{'HUMAN'}, {'XENTR'}, {'PANTR', 'MOUSE'}, {}]
        for param in param_list:
            expected = sum([expectedCnts[z] for z in param])
            returned = len(self._op.getGeneIds(speciesFilter=param))
            self.assertEqual(expected, returned, 'failed with {}'.format(param))

    def test_xrefMapping(self):
        xreftags = dict(protId='', geneId='g')
        allGenes = self._op.getGeneIds()
        for gid in allGenes:
            species = self._op.mapGeneToSpecies(gid)
            for xreftag, prefix in xreftags.items():
                expectedId = "{}{}{}".format(species, prefix, int(gid) % 10)
                xref = self._op.mapGeneToXRef(gid, xreftag)
                self.assertEqual(xref, expectedId,
                                 "xrefmapping failed for {}: {} vs {}"
                                 .format(gid, xref, expectedId))


class TaxNodeTest(unittest.TestCase):
    def setUp(self):
        self.root = fa.TaxNode('root')
        self.child = fa.TaxNode('child')

        # For the iteration tests
        self.iter_test_root = fa.TaxNode('root')
        self.left = fa.TaxNode('left')
        self.right = fa.TaxNode('right')
        self.leftleft = fa.TaxNode('left-left')
        self.leftright = fa.TaxNode('left-right')
        self.left.add_child(self.leftleft)
        self.left.add_child(self.leftright)
        self.iter_test_root.add_child(self.left)
        self.iter_test_root.add_child(self.right)

    def tearDown(self):
        self.root = None
        self.child = None
        self.iter_test_root = None
        self.left = None
        self.right = None
        self.leftleft = None
        self.leftright = None

    def test_add_child(self):
        self.root.add_child(self.child)
        self.assertIn(self.child, self.root.down)

    def test_nameProperlySet(self):
        self.assertEqual(self.root.name, 'root')

    def test_setParentOnce(self):
        self.child.add_parent(self.root)
        self.assertEqual(self.child.up, self.root)

    def test_setParentTwiceSame(self):
        self.child.add_parent(self.root)
        self.child.add_parent(self.root)
        self.assertTrue(True)

    def test_failWithTwoDifferentParentsSet(self):
        extraNode = fa.TaxNode('extra')
        self.child.add_parent(self.root)
        self.assertRaises(fa.TaxonomyInconsistencyError,
                          self.child.add_parent,
                          extraNode)

    def test_iter_leaves(self):
        leaves = list(self.iter_test_root.iter_leaves())
        self.assertEqual(leaves, [self.leftleft, self.leftright, self.right])

    def test_iter_preorder(self):
        nodes = list(self.iter_test_root.iter_preorder())
        expected = [self.iter_test_root,
                    self.left,
                    self.leftleft,
                    self.leftright,
                    self.right]
        self.assertEqual(nodes, expected)

    def test_iter_postorder(self):
        nodes = list(self.iter_test_root.iter_postorder())
        expected = [self.leftleft,
                    self.leftright,
                    self.left,
                    self.right,
                    self.iter_test_root]
        self.assertEqual(nodes, expected)

    def test_iter_levelorder(self):
        nodes = list(self.iter_test_root.iter_levelorder())
        expected = [self.iter_test_root,
                    self.left,
                    self.right,
                    self.leftleft,
                    self.leftright]
        self.assertEqual(nodes, expected)


class GeneFamilyTest(unittest.TestCase):

    def getLastExampleFamily(self, parser=None):
        if parser is None:
            parser = SetupHelper.createOrthoXMLParserFromSimpleEx()
        fam = parser.getToplevelGroups()[-1]
        return(fam)

    def test_onlyOrthologGroup(self):
        fam = self.getLastExampleFamily()
        paralogNode = fam.find(".//{{{ns0}}}paralogGroup".format(**fa.OrthoXMLParser.ns))
        self.assertRaises(fa.ElementError, fa.GeneFamily, paralogNode)

    def test_members(self):
        fam = self.getLastExampleFamily()
        expectedMembers = {'3', '13', '23', '33', '53', '14', '34'}
        gf = fa.GeneFamily(fam)
        members = set(gf.getMemberGenes())
        self.assertSetEqual(expectedMembers, members)

    def famMemberGenesAtLevel(self, fam, level):
        gf = fa.GeneFamily(fam)
        levelAnalysis = gf.analyzeLevel(level)
        genes = set()
        for subfam in levelAnalysis:
            genes.update(subfam.getMemberGenes())
        return genes

    def test_famMemberAtSubLevel(self):
        parser = SetupHelper.createOrthoXMLParserFromSimpleEx()
        fam = self.getLastExampleFamily(parser)
        cases = {'Primates': {'3', '13'}, 'Euarchontoglires': {'3', '13', '14', '33', '34'}}
        for lev, expMemb in cases.items():
            genes = self.famMemberGenesAtLevel(fam, lev)
            self.assertSetEqual(genes, expMemb)

    def test_simpleAnalyzeStrategy(self):
        """test the classification of genes for a few testcases"""
        parser = SetupHelper.createOrthoXMLParserFromSimpleEx()
        fam = fa.GeneFamily(self.getLastExampleFamily(parser))
        analyzer = fa.BasicLevelAnalysis(parser)
        summary = analyzer.analyzeGeneFam(fam)
        hum = summary['HUMAN']
        self.assertEqual(hum.typ, "SINGLECOPY")
        self.assertSetEqual(hum.genes, {'3'})
        ptr = summary['PANTR']
        self.assertEqual(ptr.typ, "MULTICOPY")
        self.assertSetEqual(ptr.genes, {'13', '14'})


class TwoLineageComparisons(unittest.TestCase):
    """tests the ability to compare the orthoxml at two
    lineages."""
    maxDiff = None

    def getFamHistory(self, parser, level):
        hist = parser.getFamHistory()
        hist.analyzeLevel(level)
        return hist

    def compareLevels(self, lev1, lev2):
        parser = SetupHelper.createOrthoXMLParserFromSimpleEx()
        tax = fa.TaxonomyFactory.newTaxonomy(parser)
        parser.augmentTaxonomyInfo(tax)
        hist1 = self.getFamHistory(parser, lev1)
        hist2 = self.getFamHistory(parser, lev2)
        return hist1.compare(hist2)

    def compareLevelsSingletonAware(self, lev1, lev2):
        parser = SetupHelper.createOrthoXMLParserFromSimpleEx()
        tax = fa.TaxonomyFactory.newTaxonomy(parser)
        parser.augmentTaxonomyInfo(tax)
        parser.augmentSingletons()
        hist1 = self.getFamHistory(parser, lev1)
        hist2 = self.getFamHistory(parser, lev2)
        return hist1.compare(hist2)

    def test_levels(self):
        levPairs = [('Mammalia', 'Primates'),
                    ('Vertebrata', 'Euarchontoglires'),
                    ('Euarchontoglires', 'Rodents'),
                    ('Vertebrata', 'HUMAN'),
                    ('Primates', 'PANTR'),
                    ('Primates', 'HUMAN')]
        expRes = [[fa.FamIdent('1'), fa.FamIdent('2'), fa.FamDupl('3', ['3.1a', '3.1b'])],
                  [fa.FamIdent('1'), fa.FamNovel('2'), fa.FamDupl('3', ['3.1a', '3.1b'])],
                  [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                  [fa.FamIdent('1'), fa.FamNovel('2'), fa.FamDupl('3', '3.1a')],
                  [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                  [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamLost('3.1b')]]
        for i in range(len(levPairs)):
            lev1, lev2 = levPairs[i]
            comp = self.compareLevels(lev1, lev2)
            comp.fams.sort(key=lambda x: x.name)
            self.assertListEqual(comp.fams, expRes[i], 'failed for {} vs {}'.format(lev1, lev2))

    def test_summarise(self):
        levPairs = [('Mammalia', 'Primates'),
                    ('Vertebrata', 'Euarchontoglires'),
                    ('Euarchontoglires', 'Rodents'),
                    ('Vertebrata', 'HUMAN'),
                    ('Primates', 'PANTR'),
                    ('Primates', 'HUMAN')]
        expRes = [{u'identical': 2, u'novel': 0, u'lost': 0, u'duplicated': 2, u'singleton': 0},
                  {u'identical': 1, u'novel': 1, u'lost': 0, u'duplicated': 2, u'singleton': 0},
                  {u'identical': 4, u'novel': 0, u'lost': 0, u'duplicated': 0, u'singleton': 0},
                  {u'identical': 1, u'novel': 1, u'lost': 0, u'duplicated': 1, u'singleton': 0},
                  {u'identical': 4, u'novel': 0, u'lost': 0, u'duplicated': 0, u'singleton': 0},
                  {u'identical': 3, u'novel': 0, u'lost': 1, u'duplicated': 0, u'singleton': 0}]
        for i in range(len(levPairs)):
            lev1, lev2 = levPairs[i]
            comp = self.compareLevels(lev1, lev2)
            summary = comp.summarise()
            self.assertDictEqual(summary, expRes[i], 'failed for {} vs {}'.format(lev1, lev2))

    def test_levels_singleton_aware(self):
        levPairs = [('Mammalia', 'Primates'),
                    ('Vertebrata', 'Euarchontoglires'),
                    ('Euarchontoglires', 'Rodents'),
                    ('Vertebrata', 'HUMAN'),
                    ('Primates', 'PANTR'),
                    ('Primates', 'HUMAN')]
        expRes = [[fa.FamIdent('1'), fa.FamIdent('2'), fa.FamDupl('3', ['3.1a', '3.1b'])],
                  [fa.FamIdent('1'), fa.FamNovel('2'), fa.FamDupl('3', ['3.1a', '3.1b'])],
                  [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                  [fa.FamIdent('1'), fa.FamNovel('2'), fa.FamDupl('3', '3.1a'), fa.FamSingleton('5')],
                  [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                  [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamLost('3.1b'), fa.FamSingleton('5')]]
        for i in range(len(levPairs)):
            lev1, lev2 = levPairs[i]
            comp = self.compareLevelsSingletonAware(lev1, lev2)
            comp.fams.sort(key=lambda x: x.name)
            self.assertListEqual(comp.fams, expRes[i], 'failed for {} vs {}'.format(lev1, lev2))

    def test_summarise_singleton_aware(self):
        levPairs = [('Mammalia', 'Primates'),
                    ('Vertebrata', 'Euarchontoglires'),
                    ('Euarchontoglires', 'Rodents'),
                    ('Vertebrata', 'HUMAN'),
                    ('Primates', 'PANTR'),
                    ('Primates', 'HUMAN')]
        expRes = [{u'identical': 2, u'singleton': 0, u'novel': 0, u'lost': 0, u'duplicated': 2},
                  {u'identical': 1, u'singleton': 0, u'novel': 1, u'lost': 0, u'duplicated': 2},
                  {u'identical': 4, u'singleton': 0, u'novel': 0, u'lost': 0, u'duplicated': 0},
                  {u'identical': 1, u'singleton': 1, u'novel': 1, u'lost': 0, u'duplicated': 1},
                  {u'identical': 4, u'singleton': 0, u'novel': 0, u'lost': 0, u'duplicated': 0},
                  {u'identical': 3, u'singleton': 1, u'novel': 0, u'lost': 1, u'duplicated': 0}]
        for i in range(len(levPairs)):
            lev1, lev2 = levPairs[i]
            comp = self.compareLevelsSingletonAware(lev1, lev2)
            summary = comp.summarise()
            self.assertDictEqual(summary, expRes[i], 'failed for {} vs {}'.format(lev1, lev2))

    def test_filter(self):
        levPairs = [('Mammalia', 'Primates'),
                    ('Vertebrata', 'Euarchontoglires'),
                    ('Euarchontoglires', 'Rodents'),
                    ('Vertebrata', 'HUMAN'),
                    ('Primates', 'PANTR'),
                    ('Primates', 'HUMAN')]
        expRes = [
                  [
                   [fa.FamIdent('1'), fa.FamIdent('2')],
                   [fa.FamIdent('1')],
                   [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                   [fa.FamIdent('1')],
                   [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                   [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a')]
                  ],
                  [
                   [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamDupl('3', ['3.1a', '3.1b'])],
                   [fa.FamIdent('1'), fa.FamNovel('2'), fa.FamDupl('3', ['3.1a', '3.1b'])],
                   [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                   [fa.FamIdent('1'), fa.FamNovel('2'), fa.FamDupl('3', '3.1a'), fa.FamSingleton('5')],
                   [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamIdent('3.1b')],
                   [fa.FamIdent('1'), fa.FamIdent('2'), fa.FamIdent('3.1a'), fa.FamLost('3.1b'), fa.FamSingleton('5')]
                  ]
                 ]
        filters = ['identical', {'identical', 'lost', 'singleton', 'novel', 'duplicated'}]
        res = [ [], [] ]
        for i in range(len(levPairs)):
            lev1, lev2 = levPairs[i]
            comp = self.compareLevelsSingletonAware(lev1, lev2)
            for j in range(len(filters)):
                fams = comp.filter(filters[j])
                res[j].append(fams)
                self.assertListEqual(res[j][i], expRes[j][i], "failed for {} vs {}".format(lev1, lev2))


class TaxonomyFactoryTest(unittest.TestCase):

    def test_xmlTaxonomyNotImpl(self):
        self.assertRaises(NotImplementedError,
                          fa.TaxonomyFactory.newTaxonomy,
                          "test/taxEx.xml")

    def test_taxFromOrthoXMLParser(self):
        p = SetupHelper.createOrthoXMLParserFromSimpleEx()
        tax = fa.TaxonomyFactory.newTaxonomy(p)
        expectedLevels = p.getLevels().union(p.getSpeciesSet())
        self.assertSetEqual(set(tax.hierarchy.keys()), expectedLevels)


class GeneQueryTest(unittest.TestCase):
    """ tests the OrthoXMLQuery methods getInputGenes and getGroupedGenes,
    with and without species filters
    """

    def setUp(self):
        self._op = SetupHelper.createOrthoXMLParserFromSimpleEx()
        self._tax = fa.TaxonomyFactory.newTaxonomy(self._op)
        self._op.augmentTaxonomyInfo(self._tax)
        self._input_query = fa.OrthoXMLQuery.getInputGenes
        self._grouped_query = fa.OrthoXMLQuery.getGroupedGenes

    def test_all_input_genes(self):
        expected_ids = {'1',  '2',  '3',  '5',  '11', '12', '13', '14',
                        '21', '22', '23', '31', '32', '33', '34', '41',
                        '43', '51', '53'}

        result_nodes = self._input_query(self._op.root)
        result_ids = {n.get('id') for n in result_nodes}
        self.assertEqual(result_ids, expected_ids)

    def test_filtered_input_genes(self):
        filters = ['HUMAN', 'PANTR', 'CANFA', 'MOUSE', 'RATNO', 'XENTR']
        expected_ids = dict(HUMAN={'1',  '2',  '3',  '5'},
                            PANTR={'11', '12', '13', '14'},
                            CANFA={'21', '22', '23'},
                            MOUSE={'31', '32', '33', '34'},
                            RATNO={'41', '43'},
                            XENTR={'51', '53'})

        for filter_ in filters:
            result_nodes = self._input_query(self._op.root, filter_)
            result_ids = {n.get('id') for n in result_nodes}
            self.assertEqual(result_ids, expected_ids[filter_],
                             'failed with {}'.format(filter_))

    def test_all_grouped_genes(self):
        expected_ids = {'1',  '2',  '3',  '11', '12', '13',
                        '14', '21', '22', '23', '31', '32',
                        '33', '34', '41', '51', '53'}

        result_nodes = self._grouped_query(self._op.root)
        result_ids = {n.get('id') for n in result_nodes}
        self.assertEqual(result_ids, expected_ids)

    def test_filtered_grouped_genes(self):
        filters = ['HUMAN', 'PANTR', 'CANFA', 'MOUSE', 'RATNO', 'XENTR']
        expected_ids = dict(HUMAN={'1',  '2',  '3'},
                            PANTR={'11', '12', '13', '14'},
                            CANFA={'21', '22', '23'},
                            MOUSE={'31', '32', '33', '34'},
                            RATNO={'41'},
                            XENTR={'51', '53'})

        for filter_ in filters:
            result_nodes = self._grouped_query(self._op.root, filter_)
            result_ids = {n.get('id') for n in result_nodes}
            self.assertEqual(result_ids, expected_ids[filter_],
                             'failed with {}'.format(filter_))


class TaxonomyNewickTest(unittest.TestCase):
    """Tests that the Taxonomy can write itself as newick """

    def setUp(self):
        self._op = SetupHelper.createOrthoXMLParserFromSimpleEx()
        self._tax = fa.TaxonomyFactory.newTaxonomy(self._op)

    def test_plain_newick(self):
        equivalent_trees = {
            ('((((RATNO, MOUSE)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),
            ('((((MOUSE, RATNO)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),
            ('((((RATNO, MOUSE)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),
            ('((((MOUSE, RATNO)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),

            ('((CANFA, ((RATNO, MOUSE)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),
            ('((CANFA, ((MOUSE, RATNO)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),
            ('((CANFA, ((RATNO, MOUSE)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),
            ('((CANFA, ((MOUSE, RATNO)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),

            ('(XENTR, (((RATNO, MOUSE)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),
            ('(XENTR, (((MOUSE, RATNO)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),
            ('(XENTR, (((RATNO, MOUSE)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),
            ('(XENTR, (((MOUSE, RATNO)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),

            ('(XENTR, (CANFA, ((RATNO, MOUSE)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires)Mammalia)Vertebrata;'),
            ('(XENTR, (CANFA, ((MOUSE, RATNO)Rodents, (PANTR, HUMAN)Primates)Euarchontoglires)Mammalia)Vertebrata;'),
            ('(XENTR, (CANFA, ((RATNO, MOUSE)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires)Mammalia)Vertebrata;'),
            ('(XENTR, (CANFA, ((MOUSE, RATNO)Rodents, (HUMAN, PANTR)Primates)Euarchontoglires)Mammalia)Vertebrata;'),

            ('((((PANTR, HUMAN)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),
            ('((((PANTR, HUMAN)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),
            ('((((HUMAN, PANTR)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),
            ('((((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, CANFA)Mammalia, XENTR)Vertebrata;'),

            ('((CANFA, ((PANTR, HUMAN)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),
            ('((CANFA, ((PANTR, HUMAN)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),
            ('((CANFA, ((HUMAN, PANTR)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),
            ('((CANFA, ((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires)Mammalia, XENTR)Vertebrata;'),

            ('(XENTR, (((PANTR, HUMAN)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),
            ('(XENTR, (((PANTR, HUMAN)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),
            ('(XENTR, (((HUMAN, PANTR)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),
            ('(XENTR, (((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires, CANFA)Mammalia)Vertebrata;'),

            ('(XENTR, (CANFA, ((PANTR, HUMAN)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires)Mammalia)Vertebrata;'),
            ('(XENTR, (CANFA, ((PANTR, HUMAN)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires)Mammalia)Vertebrata;'),
            ('(XENTR, (CANFA, ((HUMAN, PANTR)Primates, (RATNO, MOUSE)Rodents)Euarchontoglires)Mammalia)Vertebrata;'),
            ('(XENTR, (CANFA, ((HUMAN, PANTR)Primates, (MOUSE, RATNO)Rodents)Euarchontoglires)Mammalia)Vertebrata;'),
        }

        result = self._tax.newick()
        self.assertIn(result, ' '.join(equivalent_trees))

    def test_decorated_newick(self):
        pass

class GeneTreeNodeTest(unittest.TestCase):

    def setUp(self):
        n = fa.GeneTreeNode('root', 'speciation')
        n.add_child(fa.GeneTreeNode('child1', 'speciation'))
        n.add_child(fa.GeneTreeNode('child2', 'speciation'))
        n.children[0].add_child(fa.GeneTreeNode('grandchild1', 'leaf'))
        n.children[0].add_child(fa.GeneTreeNode('grandchild2', 'leaf'))
        dup = fa.GeneTreeNode('D', 'duplication')
        n.children[1].add_child(dup)
        dup.add_child(fa.GeneTreeNode('grandchild3_A', 'leaf'))
        dup.add_child(fa.GeneTreeNode('grandchild3_B', 'leaf'))
        n.children[1].add_child(fa.GeneTreeNode('grandchild4', 'leaf'))
        self.root = n

    def test_iter(self):
        expected = ('root child1 grandchild1 grandchild2 child2 D'
                    ' grandchild3_A grandchild3_B grandchild4')
        result = ' '.join(x.name for x in self.root)
        self.assertEqual(expected, result)

    def test_delete(self):
        d = self.root.children[1]
        d.delete()
        expected_n = ('root child1 grandchild1 grandchild2 D grandchild3_A'
                      ' grandchild3_B grandchild4')
        expected_d = ('child2')
        result_n = ' '.join(x.name for x in self.root)
        result_d = ' '.join(x.name for x in d)
        self.assertEqual((expected_n, expected_d), (result_n, result_d))

    def test_delete_recursive(self):
        self.setUp()
        d = self.root.children[1]
        d.delete(True)
        expected_n = ('root child1 grandchild1 grandchild2')
        expected_d = ('child2')
        result_n = ' '.join(x.name for x in self.root)
        result_d = ' '.join(x.name for x in d)
        self.assertEqual((expected_n, expected_d), (result_n, result_d))

    def test_reroot(self):
        self.setUp()
        d = self.root.children[1]
        d.reroot()
        expected = ('child2 D grandchild3_A grandchild3_B grandchild4 root'
                    ' child1 grandchild1 grandchild2')
        result = ' '.join(x.name for x in d)
        self.assertEqual(expected, result)

import unittest
import familyanalyzer as fa


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

    def tearDown(self):
        self.root = None
        self.child = None

    def test_addChild(self):
        self.root.addChild(self.child)
        self.assertIn(self.child, self.root.down)

    def test_nameProperlySet(self):
        self.assertEqual(self.root.name, 'root')

    def test_setParentOnce(self):
        self.child.addParent(self.root)
        self.assertEqual(self.child.up, self.root)

    def test_setParentTwiceSame(self):
        self.child.addParent(self.root)
        self.child.addParent(self.root)
        self.assertTrue(True)

    def test_failWithTwoDifferentParentsSet(self):
        extraNode = fa.TaxNode('extra')
        self.child.addParent(self.root)
        self.assertRaises(fa.TaxonomyInconsistencyError,
                          self.child.addParent,
                          extraNode)


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

    def membSubSetsAtLevel(self, fam, level, parser):
        gf = fa.GeneFamily(fam)
        levelAnalysis = gf.analyzeLevel(level, parser)
        return levelAnalysis.geneClasses()

#    def test_humanMemberSubSetsAtPrimates(self):
#        parser = SetupHelper.createOrthoXMLParserFromSimpleEx()
#        fam = self.getLastExampleFamily(parser)
#        geneClasses = self.membSubSetsAtLevel(fam, 'Primates', parser)
#        self.assertFalse(True,"test not yet finished")


class TaxonomyFactoryTest(unittest.TestCase):

    def test_xmlTaxonomyNotImpl(self):
        self.assertRaises(NotImplementedError,
                          fa.TaxonomyFactory.newTaxonomy,
                          "test/taxEx.xml")

    def test_taxFromOrthoXMLParser(self):
        p = SetupHelper.createOrthoXMLParserFromSimpleEx()
        tax = fa.TaxonomyFactory.newTaxonomy(p)
        self.assertSetEqual(set(tax.hierarchy.keys()), p.getLevels())


class SimpleTaxonomyTest(unittest.TestCase):
    def setUp(self):
        pass

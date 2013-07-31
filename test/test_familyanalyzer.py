import unittest
import familyanalyzer as fa

#class ToyExampleTestCase(unittest.TestCase):
#    _fname=None
#
#    @classmethod
#    def setUpClass(cls):
#        cls._fname = tempfile.mktemp(suffix=".h5")
#        setupH5File(cls._fname)
#
#    @classmethod
#    def tearDownClass(cls):
#        destroyH5File(cls._fname)
#
#    def setUp(self):
#        self._reader = siblings.Reader(self._fname)
#
#    def tearDown(self):
#        self._reader.close()


class OrthoXMLParserTest(unittest.TestCase):
    def setUp(self):
        testfile = "test/simpleEx.orthoxml"
        self._op = fa.OrthoXMLParser(testfile)

    def test_nrOfToplevelFamilies(self):
        self.assertEqual(len(self._op.getToplevelGroups()), 3)

    def test_returnedSpecies(self):
        expectedSpecies = {'HUMAN','PANTR','MOUSE','RATNO',
            'CANFA','XENTR'}
        self.assertSetEqual( self._op.getSpeciesSet(), 
            expectedSpecies)

    def test_numberOfGenesPerSpecies(self):
        expectedCnts = dict(HUMAN=4,PANTR=4,MOUSE=4,RATNO=2,CANFA=3,XENTR=2)
        allGenes = self._op.getGeneIds()
        for species in expectedCnts.keys():
            self.assertEqual( 
                len([gid for gid in allGenes if self._op.mapGeneToSpecies(gid)==species]),
                expectedCnts[species],
                "number of genes not correct for "+species)

    def test_xrefMapping(self):
        xreftags = dict(protId='',geneId='g')
        specToOffs = dict(HUMAN=0,PANTR=10,MOUSE=30,RATNO=40,CANFA=20,XENTR=50)
        allGenes = self._op.getGeneIds()
        for gid in allGenes:
            species = self._op.mapGeneToSpecies(gid)
            for xreftag, prefix in xreftags.items():
                expectedId = "{}{}{}".format(species, prefix, int(gid)%10)
                xref = self._op.mapGeneToXRef(gid,xreftag) 
                self.assertEqual( xref, expectedId,
                    "xrefmapping failed for {}: {} vs {}".format(gid, xref, expectedId) )

class TaxNodeTest(unittest.TestCase):
    def setUp(self):
        self.root = fa.TaxNode('root')
        self.child = fa.TaxNode('child')

    def tearDown(self):
        self.root=None
        self.child=None

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
        extraNode=fa.TaxNode('extra')
        self.child.addParent(self.root)
        self.assertRaises( fa.TaxonomyInconsistencyError, 
            self.child.addParent, 
            extraNode)

class SimpleTaxonomyTest(unittest.TestCase):
    def setUp(self):
       pass 

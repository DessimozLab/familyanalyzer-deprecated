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
        self._op = fa.OrthoXMLParser(testfile)

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

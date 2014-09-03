from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_hooks()

import unittest
import familyanalyzer.taxonomy as tax
import io

class NewickTaxonomyTest(unittest.TestCase):
    def test_binary_tree(self):
        newick_tree = "((Leaf1,Leaf2)Internal1_2,Leaf3)Internal1_3;"
        self.check_root_and_return_rootnode(newick_tree)

    def test_multifurcating_tree(self):
        newick_tree = "((Leaf1,Leaf2)Internal1_2,Leaf3,Leaf4)Internal1_4;"
        root = self.check_root_and_return_rootnode(newick_tree) 
        self.assertEqual(sorted([x.name for x in root.down]), 
                         sorted(['Internal1_2','Leaf3','Leaf4']),
                         'children of root node not as expected')

    def check_root_and_return_rootnode(self, tree):
        fp = io.StringIO(tree)
        taxonomy = tax.NewickTaxonomy(fp)
        root_name = tree[tree.rfind(')')+1:tree.rfind(';')]
        root = taxonomy[root_name]
        self.assertIsNone(root.up, 'not root element')
        self.assertListEqual(list(taxonomy.iterParents(root_name)),[],'not root element')
        return root

import collections
import io
import os
import unittest
from unittest.mock import MagicMock
import familyanalyzer.familyanalyzer as fa
import json
import itertools
import tempfile
import shutil
from string import Template
import re

__author__ = 'adriaal'


class OGLevelMapper(object):
    def __init__(self, ogs):
        self.levels = collections.defaultdict(list)
        self.id2pos = {}
        for og in ogs:
            for lev in fa.OrthoXMLQuery.getLevels(og):
                self.levels[lev].append(og)
                pos = len(self.levels[lev]) - 1
                try:
                    assert og.get('og') is not None
                    pos_before = self.id2pos[og.get('og')]
                    if pos != pos_before:
                        raise HOGError(
                            'HOG {:!r} with several levels has inconsistent positions. This should not happen')
                except KeyError:
                    self.id2pos[og.get('og')] = pos

    def index(self, og):
        return self.id2pos[og.get('og')]

    def nr_subhogs_on_level(self, level):
        return len(self.levels[level])


class HOGVisExtractor(object):
    def __init__(self, fname):
        self.parser = fa.OrthoXMLParser(fname)
        self.re_non_char = re.compile(r'\W')  # matches non-standard chars (not A-Za-z0-9 or _)

    def setup_extractor(self):
        self.tax = fa.TaxRangeOrthoXMLTaxonomy(self.parser)
        self.parser.augmentTaxonomyInfo(self.tax)
        anno = fa.GroupAnnotator(self.parser)
        anno.annotateDoc()

    def get_per_species_structure(self, fam):
        genes_per_species = self.genenodes_per_species(fam)
        ogs_mapper = self.get_groups_to_level_mapper(fam)
        per_species = {}
        for org, genelist in genes_per_species.items():
            per_species[org] = cur_map = {}
            cur_map[org] = [[] for _ in range(len(genelist))]
            for pos, gene in enumerate(genelist):
                cur_map[org][pos].append(int(gene.get('id')))
                parent = gene.getparent()
                while self.parser.is_evolutionary_node(parent):
                    if self.parser.is_ortholog_group(parent):
                        for lev in fa.OrthoXMLQuery.getLevels(parent):
                            lev_repr = lev
                            if self.re_non_char.search(lev):
                                lev_repr = '"{}"'.format(lev)
                            if lev_repr not in cur_map:
                                cur_map[lev_repr] = [[] for _ in range(ogs_mapper.nr_subhogs_on_level(lev))]
                            if lev != org:
                                cur_map[lev_repr][ogs_mapper.index(parent)].append(int(gene.get('id')))
                    parent = parent.getparent()
        return per_species

    def get_groups_to_level_mapper(self, fam):
        ogs = fa.OrthoXMLQuery.getSubNodes('orthologGroup', fam)
        ogs.append(fam)
        return OGLevelMapper(ogs)

    def genenodes_per_species(self, fam):
        genes = collections.defaultdict(list)
        geneRefs = fa.OrthoXMLQuery.getSubNodes("geneRef", fam)
        for gref in geneRefs:
            gid = gref.get('id')
            sp = self.parser.mapGeneToSpecies(gid)
            genes[sp].append(gref)
        return genes

    def iter_families(self):
        for fam in self.parser.getToplevelGroups():
            topology = self.get_topology(fam)
            per_species = self.get_per_species_structure(fam)
            xrefs = self.get_xrefs(fam)
            yield topology, per_species, xrefs, int(fam.get('id'))

    def get_topology(self, fam):
        lev_of_famroot = fa.OrthoXMLQuery.getLevels(fam)[0]
        return str(self.tax.hierarchy[lev_of_famroot])+';'

    def get_xrefs(self, fam):
        geneRefs = fa.OrthoXMLQuery.getSubNodes("geneRef", fam)
        xrefs = {}
        for gene in geneRefs:
            gid = gene.get('id')
            xrefs[int(gid)] = dict(fa.OrthoXMLQuery.getGeneFromId(gid, self.parser.root).attrib)
        return xrefs


class HOGError(Exception):
    pass


class Writer(object):
    html_template = Template("""<!DOCTYPE html>
<html>
<head>
 <title>$name</title>
 <!-- D3 -->
 <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>

 <!-- Underscore -->
 <script src="http://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.6.0/underscore-min.js"></script>

 <!-- TnT -->
 <link rel="stylesheet" href="http://cbrg-oma.ethz.ch/static/css/tnt.css" type="text/css" />
 <script src="http://cbrg-oma.ethz.ch/static/js/tnt.js"></script>
 <script src="http://cbrg-oma.ethz.ch/static/js/hog.js"></script>

  <h1>HOG Viewer for $name</h1>
  <div id="hog_tree"></div>

  <script>

    (function () {
      var tree = tnt.tree.parse_newick('$species_tree');
        var query;
        var per_species = $per_species ;
        var tooltip_data = $xrefs ;
        var options = {'show_internal_labels': "true",
                       'oma_info_url_template': ""};

        var viz = tnt.tree_annot();
        var theme = tnt_theme_tree_hog();
        theme (viz, document.getElementById("hog_tree"), query, per_species, tree, tooltip_data, options);
    }) ();

  </script>
</body>
</html>""")

    def __init__(self, folder, extractor):
        if os.path.exists(folder):
            if not os.path.isdir(folder):
                raise FileExistsError(folder)
        else:
            os.mkdir(folder)
        self.outdir = folder
        self.extractor = extractor

    def dump_hogs(self):
        for famdata in self.extractor.iter_families():
            newick, per_species, xrefs, fam_id = famdata
            self.write_hog(fam_id, newick, per_species, xrefs)

    def write_hog(self, fam_id, species_tree, per_species, xrefs):
        with open(os.path.join(self.outdir, "hog{:06d}.html".format(fam_id)), 'w') as fh:
            fh.write(self.html_template.safe_substitute({'name': "HOG{:06d}".format(fam_id),
                                                         'species_tree': species_tree,
                                                         'xrefs': json.dumps(xrefs),
                                                         'per_species': json.dumps(per_species)}))


class HOGVisExtractorTest(unittest.TestCase):
    def setUp(self):
        file = io.BytesIO("""<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="Family Analyzer Testcase" originVersion="0.2">
 <notes>
Example Notes without meaning
 </notes>
 <species name="HUMAN" NCBITaxId="9601">
  <database name="HUMANfake" version="0.1">
   <genes>
    <gene id="3" protId="HUMAN3" geneId="HUMANg3" />
    <gene id="5" protId="HUMAN5" geneId="HUMANg5" />
   </genes>
  </database>
 </species>
 <species name="PANTR" NCBITaxId="9483">
  <database name="PANTRfake" version="0.1">
   <genes>
    <gene id="13" protId="PANTR3" geneId="PANTRg3" />
    <gene id="14" protId="PANTR4" geneId="PANTRg4" />
   </genes>
  </database>
 </species>
 <species name="CANFA" NCBITaxId="9615">
  <database name="CANFAfake" version="0.1">
   <genes>
    <gene id="23" protId="CANFA3" geneId="CANFAg3" />
   </genes>
  </database>
 </species>
 <species name="MOUSE" NCBITaxId="10090">
  <database name="MOUSEfake" version="0.1">
   <genes>
    <gene id="33" protId="MOUSE3" geneId="MOUSEg3" />
    <gene id="34" protId="MOUSE4" geneId="MOUSEg4" />
   </genes>
  </database>
 </species>
 <species name="RATNO" NCBITaxId="10116">
  <database name="RATNOfake" version="0.1">
   <genes>
    <gene id="41" protId="RATNO1" geneId="RATNOg1" />
    <gene id="43" protId="RATNO3" geneId="RATNOg3" />
   </genes>
  </database>
 </species>
 <species name="XENTR" NCBITaxId="1">
  <database name="XENTRfake" version="0.1">
   <genes>
    <gene id="51" protId="XENTR1" geneId="XENTRg1" />
    <gene id="53" protId="XENTR3" geneId="XENTRg3" />
   </genes>
  </database>
 </species>
 <groups>
  <orthologGroup id="3">
   <property name="TaxRange" value="Vertebrata"/>
   <geneRef id="53"/>
   <orthologGroup >
    <property name="TaxRange" value="Mammalia"/>
    <geneRef id="23"/>
    <paralogGroup>
     <orthologGroup >
      <property name="TaxRange" value="Euarchontoglires" />
      <geneRef id="33"/>
      <orthologGroup>
       <property name="TaxRange" value="Primates"/>
       <paralogGroup>
        <geneRef id="3"/>
        <geneRef id="5"/>
       </paralogGroup>
       <geneRef id="13"/>
      </orthologGroup>
     </orthologGroup>
     <orthologGroup >
      <property name="TaxRange" value="Euarchontoglires" />
      <geneRef id="34" />
      <geneRef id="14" />
     </orthologGroup>
    </paralogGroup>
   </orthologGroup>
  </orthologGroup>
 </groups>
</orthoXML>""".encode('utf-8'))
        self.hog_extractor = HOGVisExtractor(file)
        self.hog_extractor.setup_extractor()
        self.fam = self.hog_extractor.parser.getToplevelGroups()[0]
        self.map = self.hog_extractor.get_groups_to_level_mapper(self.fam)

    def test_levelMap(self):
        self.assertEqual(self.map.nr_subhogs_on_level('Euarchontoglires'), 2)
        self.assertEqual(self.map.nr_subhogs_on_level('Primates'), 2)
        self.assertEqual(self.map.nr_subhogs_on_level('Mammalia'), 1)
        self.assertEqual(self.map.nr_subhogs_on_level('LUCA'), 0)

    def test_og_id_equality(self):
        primate_node = self.hog_extractor.parser.getSubFamilies('Primates')[0]
        primate_from_map = self.map.levels['Primates'][0]
        glires = self.map.levels['Euarchontoglires']
        self.assertEqual(primate_node.get('og'), primate_from_map.get('og'))
        self.assertNotEqual(glires[0].get('og'), glires[1].get('og'))

    def test_per_species_mapping(self):
        mapping = self.hog_extractor.genenodes_per_species(self.fam)
        mouse_geneids = [int(z.get('id')) for z in mapping['MOUSE']]
        self.assertEqual([33, 34], mouse_geneids)

    def test_mapping(self):
        expected_HUMAN = {'HUMAN': [[3], [5]], 'Primates': [[3, 5],[]],
                          'Euarchontoglires': [[3, 5], []],
                          'Mammalia': [[3, 5]], 'Vertebrata': [[3, 5]]}
        expected_MOUSE = {'MOUSE': [[33], [34]], 'Euarchontoglires': [[33], [34]],
                          'Mammalia': [[33, 34]], 'Vertebrata': [[33, 34]]}
        per_species = self.hog_extractor.get_per_species_structure(self.fam)
        self.assertDictEqual(expected_HUMAN, per_species['HUMAN'])
        self.assertDictEqual(expected_MOUSE, per_species['MOUSE'])

    @unittest.skip  # newick string is not invariant against ordering...
    def test_topology(self):
        expected = "(XENTR, ((MOUSE, (HUMAN, PANTR)Primates)Euarchontoglires, CANFA)Mammalia)Vertebrata;"
        topology = self.hog_extractor.get_topology(self.fam)
        self.assertEqual(expected, topology)

    def test_xrefs(self):
        exprected_of_human5 = {'id': "5", 'protId': "HUMAN5", 'geneId': "HUMANg5"}
        xrefs = self.hog_extractor.get_xrefs(self.fam)
        self.assertEqual(exprected_of_human5, xrefs[5])


class WriterTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.extractor_mock = HOGVisExtractor(io.BytesIO(b"""<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="Family Analyzer Testcase" originVersion="0.2"></orthoXML>"""))
        self.extractor_mock.iter_families = MagicMock(return_value=(('(HUMAN,MOUSE)CANFA;',
                                                                     {'HUMAN': [2, 4]}, {2: {'id': 2}}, 1),))

    def test_write(self):
        writer = Writer(self.dir, self.extractor_mock)
        writer.dump_hogs()
        expected_filename = os.path.join(self.dir, 'hog{:06d}.html'.format(1))
        self.assertTrue(os.path.exists(expected_filename))
        with open(expected_filename) as fh:
            data = fh.read()
        self.assertIn('<title>HOG000001</title>', data)
        self.assertIn('(HUMAN,MOUSE)CANFA;', data)

    def tearDown(self):
        shutil.rmtree(self.dir)


def handle_args():
    import argparse
    parser = argparse.ArgumentParser(prog='StandaloneHogVis',
                                     description="Tool to prepare HogVis html pages to analyze the HOGs graphically")
    parser.add_argument('-o', '--outdir', default='hogvis', help="directory where to store the html files."
                                                                 "(defaults: 'vis/')")
    parser.add_argument('orthoxml', help="path to the orthoxml file that should be converted")

    conf = parser.parse_args()

    ex = HOGVisExtractor(conf.orthoxml)
    ex.setup_extractor()
    wr = Writer(conf.outdir, ex)
    wr.dump_hogs()

if __name__ == "__main__":
    handle_args()
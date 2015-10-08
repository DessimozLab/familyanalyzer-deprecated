import collections
import io
import unittest
import familyanalyzer.familyanalyzer as fa
import json

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
                    pos_before = self.id2pos[og.get('og')]
                    if pos != pos_before:
                        raise HOGError(
                            'HOG {!r} with several levels has inconsistent positions. This should not happen')
                except KeyError:
                    self.id2pos[og.get('og')] = pos

    def index(self, og):
        return self.id2pos[og.get('og')]

    def nr_subhogs_on_level(self, level):
        return len(self.levels[level])


class HOGVisExtractor(object):
    def __init__(self, fname):
        self.parser = fa.OrthoXMLParser(fname)

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
                            if lev not in cur_map:
                                cur_map[lev] = [[] for _ in range(ogs_mapper.nr_subhogs_on_level(lev))]
                            cur_map[lev][ogs_mapper.index(parent)].append(int(gene.get('id')))
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


class HOGError(Exception):
    pass


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
        anno = fa.GroupAnnotator(self.hog_extractor.parser)
        anno.annotateDoc()

        self.fam = self.hog_extractor.parser.getToplevelGroups()[0]
        self.map = self.hog_extractor.get_groups_to_level_mapper(self.fam)

    def test_levelMap(self):
        self.assertEqual(self.map.nr_subhogs_on_level('Euarchontoglires'), 2)
        self.assertEqual(self.map.nr_subhogs_on_level('Primates'), 1)
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
        expected_HUMAN = {'HUMAN': [[3], [5]], 'Primates': [[3, 5]],
                          'Euarchontoglires': [[3, 5], []],
                          'Mammalia': [[3, 5]], 'Vertebrata': [[3, 5]]}
        expected_MOUSE = {'MOUSE': [[33], [34]], 'Euarchontoglires': [[33], [34]],
                          'Mammalia': [[33, 34]], 'Vertebrata': [[33, 34]]}
        per_species = self.hog_extractor.get_per_species_structure(self.fam)
        self.assertDictEqual(expected_HUMAN, per_species['HUMAN'])
        self.assertDictEqual(expected_MOUSE, per_species['MOUSE'])

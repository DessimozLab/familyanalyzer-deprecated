from .messaging import PROGRESSBAR, setup_progressbar
import re

class GeneTree(object):
    def __init__(self, root):
        self.root = root

    def __str__(self):
        return str(self.root) + ';'

    def __iter__(self):
        return iter(self.root)

    def write(self, NHX=True):
        return self.root.write(NHX) + ';'


class GeneTreeNode(object):

    events = {'speciation': ('S', 'F'),
              'duplication': ('D', 'T'),
              'loss': ('L', 'F'),
              'leaf': ('', 'F')}

    reg = re.compile(r'\W') # matches anything that's NOT a-z, A-Z, 0-9 or _

    def __init__(self, name, node_type, level=''):
        self.name            = name
        self.node_type       = node_type
        self.taxonomic_level = level
        self.children        = list()
        self.parent          = None
        self.genes           = list()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def node_type(self):
        return self._node_type

    @node_type.setter
    def node_type(self, node_type):
        if not node_type in ['speciation', 'duplication', 'loss', 'leaf']:
            raise Exception('Unexpected node_type: {}'.format(node_type))
        self._node_type = node_type

    def add_child(self, child):
        self.children.append(child)
        child.parent = self

    def is_leaf(self):
        return len(self.children) == 0

    def __iter__(self):
        yield self
        for ch in self.children:
            for element in ch:
                yield element

    def __str__(self):
        return self.write(NHX=False)

    def build_NHX_string(self):
        NHX_string = '&&NHX'
        event, duplication = self.events[self.node_type]

        if self.node_type is not 'leaf':
            NHX_string += ':Ev={0}'.format(event)
            if self.node_type is not 'loss':
                NHX_string += ':D={0}'.format(duplication)
        if self.taxonomic_level > '':
            NHX_string += ':S={0}'.format(self.taxonomic_level)
        return '[' + NHX_string + ']'

    def write(self, NHX=True):
        name = self.name
        if self.reg.search(name):
            name = '"{0}"'.format(name)
        NHX_string = self.build_NHX_string()
        if self.is_leaf():
            return '{0}{1}'.format(name, (NHX_string if NHX else ''))
        else:
            subtree = ', '.join(child.write(NHX) for child in self.children)
            return '({0}){1}{2}'.format(subtree, name, (NHX_string if NHX else ''))


class GeneTreeTracer(object):
    def __init__(self, parser, taxonomy):
        self.parser = parser
        self.tax = taxonomy
        self.trees = list()
        self._dup_counter = 0
        self._loss_counter = 0

    def _reset_counters(self):
        self._dup_counter = self._loss_counter = 0

    def _get_job_list(self, add_genelists=False):
        job_list = list()
        for n in self.tax:
            if n.up is None:
                for family in n.history.geneFamList:
                    job_list.append((n, family, None, add_genelists))
            else:
                for family in n.comparison.filter('novel'):
                    job_list.append((n, n.history[family.fam], None,
                        add_genelists))

        return sorted(job_list, key=lambda x: x[1])

    def trace_gene_families(self, add_genelists=False):
        job_list = self._get_job_list(add_genelists)

        if PROGRESSBAR:
            pbar = setup_progressbar('Extracting {0} genetrees: '
                                     .format(len(job_list)), len(job_list))
            pbar.start()

        for i, job in enumerate(job_list, start=1):
            self.trace_gene_family(*job)
            if PROGRESSBAR:
                pbar.update(i)

        if PROGRESSBAR:
            pbar.finish()

    def _create_parent(self, name, node_type, level):
        parent = GeneTreeNode(name, node_type, level)
        self.trees.append(GeneTree(parent))
        return parent

    def trace_gene_family(self,
                          taxnode,
                          family,
                          parent=None,
                          add_genelists=False,
                          ):

        fam_id = family.getFamId()

        if taxnode.isLeaf():
            return

        for child in taxnode.down:
            comparison = child.comparison[fam_id]

            if comparison.event == 'lost':
                self._loss_counter += 1
                node = GeneTreeNode('LOSS_{0}'.format(self._loss_counter),
                            node_type='loss',
                            level=child.name)
                parent.add_child(node)

            elif comparison.event == 'identical':
                if parent is None:
                    parent = self._create_parent(taxnode.name, 'speciation',
                        taxnode.name)
                    if add_genelists:
                        parent.genes = fam.getMemberGenes()
                next_fam = child.history[fam_id]
                if child.isLeaf():
                    gene = next_fam.getMemberGenes()[0]
                    label = self.parser.mapGeneToXRef(gene)
                    node = GeneTreeNode(label, node_type='leaf',
                        level=child.name)
                else:
                    node = GeneTreeNode(child.name, node_type='speciation',
                        level=child.name)
                if add_genelists:
                    node.genes = next_fam.getMemberGenes()
                parent.add_child(node)
                self.trace_gene_family(child, next_fam, node)

            elif comparison.event == 'duplicated':
                if parent is None:
                    parent = self._create_parent(taxnode.name, 'duplication',
                        taxnode.name)
                self._dup_counter += 1
                dup = GeneTreeNode('DUP_{0}'.format(self._dup_counter),
                    node_type='duplication', level='')
                parent.add_child(dup)
                for dup_id in comparison.into.split('; '):
                    next_fam = child.history[dup_id]
                    if child.isLeaf():
                        gene = next_fam.getMemberGenes()[0]
                        label = self.parser.mapGeneToXRef(gene)
                        node = GeneTreeNode(label, node_type='leaf', level=child.name)
                    else:
                        node = GeneTreeNode(child.name, node_type='speciation',
                            level=child.name)
                    if add_genelists:
                        node.genes = next_fam.getMemberGenes()
                    dup.add_child(node)
                    self.trace_gene_family(child, next_fam, node)

        return parent

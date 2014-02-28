try:
    from progressbar import ProgressBar, Percentage, Timer, ETA, Bar
    PROGRESSBAR = True
except ImportError:
    PROGRESSBAR = False

def setup_progressbar(msg, size):
    if not msg.endswith(': '):
        msg += ': '

    widgets = [msg,
               Percentage(), ' ',
               Bar(), ' ',
               Timer(), ' ',
               ETA()]

    pbar = ProgressBar(widgets=widgets, maxval=size)
    return pbar

def enum(*sequential, **named):
    """creates an Enum type with given values"""
    enums = dict(zip(sequential, range(len(sequential))), **named)
    enums['reverse'] = dict((value, key) for key, value in enums.items())
    return type('Enum', (object, ), enums)

def handle_args():
    import argparse, pkg_resources

    description = 'Analyze Hierarchical OrthoXML families.'

    help_messages = {
        'xreftag': ("xref tag of genes to report. OrthoXML allows to "
                    "store multiple ids and xref annotations per gene "
                    "as attributes in the species section. If not set, "
                    "the internal (purely numerical) ids are reported."),

        'show_levels': ('print the levels and species found in the orthoXML'
                        ' file and quit'),

        'taxonomy': ("Taxonomy used to reconstruct intermediate levels. "
                     "Has to be either 'implicit' (default) or a path to "
                     "a file in Newick format. The taxonomy might be "
                     "multifurcating. If set to 'implicit', the "
                     "taxonomy is extracted from the input OrthoXML file. "
                     "The orthoXML level do not have to cover all the "
                     "levels for all families. In order to infer gene losses "
                     "Family-Analyzer needs to infer these skipped levels "
                     "and reconcile each family with the complete taxonomy."),

        'propagate_top': ("propagate taxonomy levels up to the toplevel. As "
                          "an illustration, consider a gene family in a "
                          "eukaryotic analysis that has only mammalian genes. "
                          "Its topmost taxonomic level will therefor be "
                          "'Mammalia' and an ancestral gene was gained at that"
                          " level. However, if '--propagete-top' is set, the "
                          "family is assumed to have already be present in the"
                          " topmost taxonomic level, i.e. Eukaryota in this "
                          "example, and non-mammalian species have all lost"
                          " this gene."),

        'show_taxonomy': 'write the taxonomy used to standard output.',

        'store_augmented_xml': ("filename to which the input orthoxml file with"
                              " augmented annotations is written. The augmented"
                              " annotations include for example the additional "
                              "taxonomic levels of orthologGroup and unique HOG"
                              " IDs."),

        'add_singletons': ("Take singletons - genes from the input set that "
                              "weren't assigned to orthologous groups and "
                              "appear as novel gains in the taxonomy leaves - "
                              "and add them to the xml as single-member "
                              "orthologGroups"),

        'compare_second_level': ("Compare secondary level with primary one, "
                                 "i.e. report what happend between the "
                                 "secondary and primary level to the individual"
                                 " histories. Note that the Second level needs"
                                 " to be younger than the primary."),

        'orthoxml': 'path to orthoxml file to be analyzed',

        'level': 'taxonomic level at which analysis should be done',

        'species': ("(list of) species to be analyzed. "
                    "Note that only genes of the selected species are "
                    "reported. In order for the output to make sense, "
                    "the selected species all must be part of the lineages "
                    "specified in 'level' (and --compare_second_level).")

    }

    parser = argparse.ArgumentParser(prog='FamilyAnalyzer',
                                     description=description)
    parser.add_argument('--xreftag', default=None,
                        help=help_messages['xreftag'])
    parser.add_argument('--show_levels', action='store_true',
                        help=help_messages['show_levels'])
    parser.add_argument('--taxonomy', default='implicit',
                        help=help_messages['taxonomy'])
    parser.add_argument('--propagate_top', action='store_true',
                        help=help_messages['propagate_top'])
    parser.add_argument('--show_taxonomy', action='store_true',
                        help=help_messages['show_taxonomy'])
    parser.add_argument('--store_augmented_xml', default=None,
                        help=help_messages['store_augmented_xml'])
    parser.add_argument('--add_singletons', action='store_true',
                        help=help_messages['add_singletons'])
    parser.add_argument('--compare_second_level', default=None,
                        help=help_messages['compare_second_level'])
    parser.add_argument('orthoxml', help=help_messages['orthoxml'])
    parser.add_argument('level', help=help_messages['level'])
    parser.add_argument('species', nargs="+", help=help_messages['species'])

    return parser.parse_args()

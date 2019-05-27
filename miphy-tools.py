import os
from optparse import OptionParser
from miphy_resources import phylo
from miphy import __author__, __version__


def group_names(gene_tree_file, prefix_size):
    tree = phylo.load_newick(gene_tree_file)
    grouped_names = {}
    for node in tree.leaves:
        name = node.name
        prefix = name[:prefix_size].lower()
        grouped_names.setdefault(prefix, []).append(name)
    return grouped_names, len(tree.leaves)

def save_info_file(output_file, grouped_names):
    species = sorted(list(grouped_names.keys()))
    buff = ['[species tree]', '({});\n'.format(','.join(species)), '[species assignments]']
    for spc in species:
        genes = ','.join(sorted(grouped_names[spc]))
        line = '%s = %s' % (spc, genes)
        buff.append(line)
    with open(output_file, 'w') as f:
        f.write('\n'.join(buff))
    return len(species)

def clean_tree(gene_tree_file, output_file):
    tree = phylo.load_newick(gene_tree_file)
    if output_file == None:
        output_file = gene_tree_file
    tree.save_newick(output_file, support_values=False, comments=False, internal_names=False, support_as_comment=False)

def midpoint_root(gene_tree_file, output_file):
    tree = phylo.load_newick(gene_tree_file)
    tree.root_midpoint()
    if output_file == None:
        output_file = gene_tree_file
    tree.save_newick(output_file, support_values=False, comments=False, internal_names=False, support_as_comment=False)


def setup_parser():
    usage_str = "python %prog COMMAND ARGUMENTS\n" + \
        "Valid commands:\n" + \
        "python %prog info GENE_TREE.nwk PREFIX_SIZE OUTPUT_FILE.txt\n\t-- Parse a phylogenetic tree in newick format, extracting the gene names, grouping those sharing a prefix of the given size, and generating an info file.\n" + \
        "python %prog clean GENE_TREE.nwk [OUTPUT_FILE.nwk]\n\t-- Validates a phylogenetic tree in newick format, ensuring it is a binary tree, and sanitizing any forbidden characters from the gene names. If no OUTPUT_FILE is given, the original tree will be overwritten.\n" + \
        "python %prog midpoint GENE_TREE.nwk [OUTPUT_FILE.nwk]\n\t-- Re-root the phylogenetic tree using the midpoint algorithm. If no OUTPUT_FILE is given, the original tree will be overwritten."
    version_str = "%%prog %s" % __version__
    parser = OptionParser(usage=usage_str, version=version_str)
    return parser
def validate_args(parser, args):
    commands = ['info', 'clean', 'midpoint']
    if len(args) == 0:
        parser.error('incorrect arguments.')
    command = args[0].lower()
    if command not in commands:
        parser.error('invalid COMMAND "%s"' % command)
    return command
def validate_info_args(parser, args):
    if not len(args) == 4:
        parser.error('incorrect number of arguments for "info". You must supply a tree in newick format, the size of the prefix to group the names with, and an output name to save the file.')
    gene_tree_file, prefix_size, output_file = args[1:]
    gene_tree_file = os.path.realpath(gene_tree_file)
    if not os.path.isfile(gene_tree_file):
        parser.error('could not locate the gene tree at %s' % gene_tree_file)
    try:
        prefix_size = int(prefix_size)
    except:
        parser.error('could not convert PREFIX_SIZE "%s" to an integer.'%(prefix_size))
    if not output_file.endswith('.txt'):
        output_file += '.txt'
    return gene_tree_file, prefix_size, output_file
def validate_clean_args(parser, args):
    if not 2 <= len(args) <= 3:
        parser.error('incorrect number of arguments for "clean". You must supply a tree in newick format, and optionally an output file name.')
    gene_tree_file = os.path.realpath(args[1])
    if not os.path.isfile(gene_tree_file):
        parser.error('could not locate the gene tree at %s' % gene_tree_file)
    if len(args) == 3:
        output_file = args[2]
        if not output_file.endswith('.nwk'):
            output_file += '.nwk'
    else:
        output_file = None
    return gene_tree_file, output_file
def validate_midpoint_args(parser, args):
    if not 2 <= len(args) <= 3:
        parser.error('incorrect number of arguments for "midpoint". You must supply a tree in newick format, and optionally an output file name.')
    gene_tree_file = os.path.realpath(args[1])
    if not os.path.isfile(gene_tree_file):
        parser.error('could not locate the gene tree at %s' % gene_tree_file)
    if len(args) == 3:
        output_file = args[2]
        if not output_file.endswith('.nwk'):
            output_file += '.nwk'
    else:
        output_file = None
    return gene_tree_file, output_file


if __name__ == '__main__':
    parser = setup_parser()
    (options, args) = parser.parse_args()
    command = validate_args(parser, args)
    if command == 'info':
        gene_tree_file, prefix_size, output_file = validate_info_args(parser, args)
        grouped_names, num_genes = group_names(gene_tree_file, prefix_size)
        num_groups = save_info_file(output_file, grouped_names)
        print('Wrote info file for %i genes from %i groups to "%s"' % (num_genes, num_groups, output_file))
    elif command == 'clean':
        gene_tree_file, output_file = validate_clean_args(parser, args)
        clean_tree(gene_tree_file, output_file)
    elif command == 'midpoint':
        gene_tree_file, output_file = validate_midpoint_args(parser, args)
        midpoint_root(gene_tree_file, output_file)

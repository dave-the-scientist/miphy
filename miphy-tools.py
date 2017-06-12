import os
from optparse import OptionParser
from scripts import newick_to_coords
from miphy import __author__, __version__
# for 'info' option, should also indicate how many species there are. It can then bin all seqs by first x letters, then combine all of the small groups into one 'misc' category.


def group_names(gene_tree_file, prefix_size):
    with open(gene_tree_file) as f:
        _, leaves, _, _, _ = newick_to_coords.parse_nodes_edges(f.read().strip())
    grouped_names = {}
    for name in leaves:
        prefix = name.replace('"','').replace("'",'').strip()
        prefix = prefix[:prefix_size].lower()
        grouped_names.setdefault(prefix, []).append(name)
    return grouped_names, len(leaves)
def save_info_file(output_file, grouped_names):
    species = sorted(list(grouped_names.keys()))
    buff = ['[species tree]','(%s);\n'%(','.join(species)),'[species assignments]']
    for spc in species:
        genes = ','.join(sorted(grouped_names[spc]))
        line = '%s = %s' % (spc, genes)
        buff.append(line)
    with open(output_file, 'w') as f:
        f.write('\n'.join(buff))
    return len(species)

def clean_tree(gene_tree_file, output_file):
    with open(gene_tree_file) as f:
        newick_str = f.read()
    _, leaves, _, _, _ = newick_to_coords.parse_nodes_edges(newick_str)
    print len(leaves), 'leaves'
    buff = []
    j, not_found = -1, float('inf')
    while True:
        j = min(x for x in [not_found]+[newick_str.find(c, j+1) for c in ',)'] if x>=0)
        if j == not_found:
            buff.append(newick_str)
            break
        i = max(x for x in [newick_str[:j].rfind(c) for c in '(,']) + 1
        if i == -1:
            continue
        name, _, branch = newick_str[i:j].partition(':')
        name = clean_name(name)
        buff.append('%s%s:%s%s' % (newick_str[:i], name, branch, newick_str[j]))
        newick_str = newick_str[j+1:]
        j = -1
    if output_file == None:
        output_file = gene_tree_file
    with open(output_file, 'w') as f:
        f.write(''.join(buff))
def clean_name(name):
    name = name.replace("'",'').replace('"','')
    name = name.replace(' ','')
    return name

def clean_tree2(gene_tree_file, output_file):
    with open(gene_tree_file) as f:
        newick_str = f.read()
    root_node = newick_str[newick_str.rindex(')')+1:-1]
    if not root_node: root_node = 'root'
    nodes_start = newick_str.find('(')
    node_gen_int = 0
    nodes, parents = set(), set()
    edges, parent = {}, {}
    while True:
        i, j = find_parentheses(newick_str)
        #k = newick_str.find(':',j)
        #[k for k in [newick_str.find(c,j+1) for c in ':,)'] if k>j]
        if i == nodes_start: pnode = root_node
        else:
            k = min(x for x in [newick_str.find(c,j+1) for c in ':,)'] if x>j)
            pnode = newick_str[j+1 : k]
        if not pnode or pnode in parents or pnode in parent:
            pnode += '_%i' % node_gen_int
            node_gen_int += 1
        nodes.add(pnode); parents.add(pnode)
        for datum in newick_str[i+1:j].split(','):
            node, _, weight = datum.partition(':')
            if not node or node in parent:
                node += '_%i' % node_gen_int
                node_gen_int += 1
            #node = node.replace("'",'')
            if weight:
                edges[pnode, node] = float(weight)
            parent[node] = pnode
            nodes.add(node)
        if pnode == root_node: break
        newick_str = '%s%s%s' % (newick_str[:i], pnode, newick_str[k:]) ###
        # Somewhere here add to a buffer, to change the node names.
    leaves = sorted([node.replace(',','') for node in nodes if node not in parents])
    return root_node, leaves, nodes, edges, parent
def find_parentheses(newick_str):
    j = newick_str.find(')')
    i = newick_str[:j].rfind('(')
    return i, j


def setup_parser():
    usage_str = "python %prog COMMAND ARGUMENTS\n" + \
        "Valid commands:\n" + \
        "python %prog info GENE_TREE.nwk PREFIX_SIZE OUTPUT_FILE.txt\n\t-- Parse a phylogenetic tree in newick format, extracting the gene names, grouping those sharing a prefix of the given size, and generating an info file.\n" + \
        "python %prog clean GENE_TREE.nwk [OUTPUT_FILE.nwk]\n\t-- Validates a phylogenetic tree in newick format, ensuring it is a binary tree, and sanitizing any forbidden characters from the gene names. If no OUTPUT_FILE is given, the original tree will be overwritten."
    #usage_str = "python %prog GENE_TREE.nwk PREFIX_SIZE OUTPUT_FILE.txt\n\nParse a phylogenetic tree in newick format, extracting the gene names, and grouping those sharing a prefix of the given size."
    version_str = "%%prog %s" % __version__
    parser = OptionParser(usage=usage_str, version=version_str)
    return parser
def validate_args(parser, args):
    commands = ['info', 'clean']
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
    else: output_file = None
    return gene_tree_file, output_file


if __name__ == '__main__':
    parser = setup_parser()
    (options, args) = parser.parse_args()
    command = validate_args(parser, args)
    if command == 'info':
        gene_tree_file, prefix_size, output_file = validate_info_args(parser, args)
        grouped_names, num_genes = group_names(gene_tree_file, prefix_size)
        num_species = save_info_file(output_file, grouped_names)
        print('Wrote %i genes from %i species to "%s"' % (num_genes, num_species, output_file))
    elif command == 'clean':
        gene_tree_file, output_file = validate_clean_args(parser, args)
        clean_tree(gene_tree_file, output_file)

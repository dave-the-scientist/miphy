import itertools, time
import numpy as np

"""
I found the algorithm described at
http://math.stackexchange.com/questions/156161/finding-the-coordinates-of-points-from-distance-matrix/423898#423898

This module parses the pairwise-distances between sequences from the given
phylogenetic tree, and then produces coordinates that satisfy those
distances constraints using the minimum number of dimensions. If the distance
matrix is spatially incompatible (ex. violating the triangle inequality), it
will also return coordinates that minimize the error. The distance matrix must
be a true distance matrix, and symmetrical. If it's not, results are impacted
by which side of the diagonal the inconsistency is found. Ensure that the
average distance between points * number of points > 0.05. Preferebly much
greater; better number might be 1.0 at a minumum. If it is less, floating point
obscures the signal and the coordinates will be noisily random.

Defines one callable function:
    extract(newick_str, [max_dimensions]) => leaves, coords
-- Newick_str is a Newick-formated phylogenetic tree, and max_dimensions (if
given) allows you to place an upper bound on the number of dimensions returned.
I don't think this is an optimized dimensionality reduction; you would likely
lose less information using a proper method. But does a pretty good job. The
returned values are 'leaves', which is an ordered list of the leaf nodes from
the tree, and 'coords' which is a n x d numpy matrix, where n is the number of
leaf nodes and d is the number of dimensions used. So each row contains the
coordinate points for the corresponding leaf node from 'leaves'.

Defines a few global variables:
    root_node -- A string identifying the root given by the newick file.
    nodes -- A set of all nodes used in the tree, the leaves plus all of the
internal nodes.
    sqrd_dist -- A n x n numpy matrix of dtype='float', where n is the number
of leaves, and D[i,j] is the squared distance between leaf i and leaf j. The
order of i through n is the same as that in the list 'leaves'.
"""

# # # # #  Globally accessible variables:
root_node, nodes, sqrd_dist = None, set(), np.zeros((1,1), dtype='float')
default_weight = 0.01 # For trees that have no branch lengths.

# # # # #  Callable function:
def extract(newick_str, max_dimensions=0):
    global root_node, nodes, sqrd_dist
    newick_str = newick_str.strip()
    root_node, leaves, nodes, edges, parent = parse_nodes_edges(newick_str)
    print('Parsed tree with %i terminal nodes.\nGenerating distance matrix...' % len(leaves))
    leaf_paths = get_leaf_paths(root_node, leaves, parent)
    sqrd_dist = generate_dist_matrix(leaves, leaf_paths, edges)
    print('Distance matrix complete.\nCalculating coordinates...')
    m_mat = calculate_M(sqrd_dist)
    coords = calculate_coords(m_mat, max_dimensions)
    print('Coordinates found using %i dimensions.' % coords.shape[1])
    return leaves, coords

def extract_coords(tree_file, max_dimensions=0):
    with open(tree_file) as f:
        newick_str = f.read().strip()
    root_node, leaves, nodes, edges, parent, leaf_paths = parse_newick_data(newick_str)
    print('Parsed tree with %i terminal nodes.\nCalculating coordinates...' % len(leaves))
    coords = calculate_coordinates(leaves, leaf_paths, edges, max_dimensions)
    print('Coordinates found using %i dimensions.' % coords.shape[1])
    return leaves, coords

def parse_newick_data(newick_str):
    root_node, leaves, nodes, edges, parent = parse_nodes_edges(newick_str)
    leaf_paths = get_leaf_paths(root_node, nodes, parent)
    return root_node, leaves, nodes, edges, parent, leaf_paths
def calculate_coordinates(leaves, leaf_paths, edges, max_dimensions, verbose=False):
    t0 = time.time()
    if verbose: print('-- Generating distance matrix...')
    sqrd_dist = generate_dist_matrix(leaves, leaf_paths, edges)
    if verbose: print('-- Complete. Preparing the matrix for MDS...')
    m_mat = calculate_M(sqrd_dist)
    if verbose: print('-- Complete. Performing multi-dimensional scaling...')
    coords = calculate_coords(m_mat, max_dimensions)
    if verbose: print('-- Complete. Calculated coordinate points for %i sequences using %i dimensions in %.2f seconds.' % (coords.shape[0], coords.shape[1], time.time()-t0))
    return coords

# # # # #  Main functions:
def parse_nodes_edges(newick_str):
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
        # pnode in parents or pnode in parent doesn't seem to happen on normal trees; what am I checking for?
        if not pnode or pnode in parents or pnode in parent:
            pnode += '_%i' % node_gen_int
            node_gen_int += 1
        nodes.add(pnode); parents.add(pnode)
        for datum in newick_str[i+1:j].split(','):
            node, _, weight = datum.partition(':')
            if not node or node in parent:
                node += '_%i' % node_gen_int
                node_gen_int += 1
            if not weight:
                weight = default_weight
            else:
                weight = float(weight)
            edges[pnode, node] = weight
            parent[node] = pnode
            nodes.add(node)
        if pnode == root_node: break
        newick_str = '%s%s%s' % (newick_str[:i], pnode, newick_str[k:]) ###
    leaves = sorted([node.replace(',','') for node in nodes if node not in parents])
    return root_node, leaves, nodes, edges, parent

def find_children(parent):
    children = {}
    for cnode, pnode in parent.items():
        children.setdefault(pnode, []).append(cnode)
    return(children)

def get_leaf_paths(root_node, leaves, parent):
    leaf_paths = {}
    for leaf in leaves:
        pnode, path = leaf, [leaf]
        while pnode != root_node:
            pnode = parent[pnode]
            path.append(pnode)
        path.reverse()
        leaf_paths[leaf] = path
    return leaf_paths

def generate_dist_matrix(leaves, leaf_paths, edges):
    l = len(leaves)
    sqrd_dist = np.zeros((l, l) ,dtype='float')
    for (node1, node2), (i, j) in zip(
        itertools.combinations(leaves, 2), itertools.combinations(range(len(leaves)), 2)):
        d = calc_leaf_dist(leaf_paths[node1], leaf_paths[node2], edges)
        sqrd_dist[i,j] = d*d
        sqrd_dist[j,i] = d*d
    return sqrd_dist

def calculate_M(sqrd_dist):
    l = sqrd_dist.shape[0]
    M = np.zeros((l, l) ,dtype='float')
    for i in range(l):
        di1 = sqrd_dist[0,i]
        for j in range(i, l):
            m = di1 + sqrd_dist[0,j] - sqrd_dist[i,j]
            M[i,j] = m
            M[j,i] = m
    return M/2.0

def calculate_coords(M, max_dimensions=0):
    """M is a positive semi-definite square matrix, so its eigen values are
    greater or equal to 0. A square semi-definite matrix is symmetric, and a
    symmetric square has no complex eigen values. So the ones I was finding
    are a result of rounding error. This was fixed by using linalg.eigh instead
    of .eig. The first columns are the least significant dimensions."""
    zero = 1e-5 # floating point 0 check.
    vals, vecs = np.linalg.eigh(M)
    #total_val = sum(v for v in vals if v > zero)
    if max_dimensions: tokeep = max(len(vals) - max_dimensions, 0)
    else: tokeep = 0
    vals, vecs = vals[tokeep:], vecs[:,range(tokeep, len(vals))]
    return np.column_stack(vecs[:,i]*np.sqrt(val) for i, val in enumerate(vals) if val > zero)


# # # # #  Utility functions:
def find_parentheses(newick_str):
    j = newick_str.find(')')
    i = newick_str[:j].rfind('(')
    return i, j

def calc_path_dist(path, edges):
    dist = 0
    node1 = path[0]
    for node2 in path[1:]:
        dist += edges[node1, node2]
        node1 = node2
    return dist
def calc_leaf_dist(path1, path2, edges):
    i = -1
    for n1, n2 in zip(path1, path2):
        if n1 == n2: i += 1
        else: break
    return calc_path_dist(path1[i:], edges) + calc_path_dist(path2[i:], edges)

"""MIPhy class object that performs the clustering and stability analysis.
"""

import os
import numpy as np
from scripts import newick_to_coords
import time


class Clusterer(object):
    def __init__(self, gene_tree_data, species_tree_data, species_map, use_coords, coords_file, verbose=False):
        self.species_map = species_map
        self.use_coords = use_coords
        self.verbose = verbose
        # #  Options for spread calculation
        self.spread_calc = self._cluster_std_dev # Swap method if desired.
        self.singleton_spread = 0.0
        self.relative_avg = 'median' # One of: 'mean' or 'median'. Specifies how the average spread is calculated.
        # Could add an option to iterate refinement.
        # #
        if verbose: print('Setting up the clusterer...')
        self.parse_gene_tree(gene_tree_data, coords_file)
        self.parse_species_tree(species_tree_data)
        if verbose: print('Finished setting up the clusterer.')
        #self.validate_data()

    # # # # #  Main methods:
    def cluster(self, i_coefficient=0.5, d_coefficient=1.0, l_coefficient=1.0, spread_coefficient=1.0):
        self.nodes = {}
        params = {'i_coef':float(i_coefficient), 'd_coef':float(d_coefficient),
            'l_coef':float(l_coefficient), 'm_coef':float(l_coefficient),
            'spread_coef':float(spread_coefficient)}
        self._calc_clusters(self.gene_root, params)
        if self.use_coords:
            self._refine_clusters(params)
        self.nodes[self.gene_root].cluster_roots.sort(key=lambda n: self.nodes[n].total_score)
        self.nodes[self.gene_root].clusters = [sorted(self.nodes[node].clusters[0]) \
            for node in self.nodes[self.gene_root].cluster_roots]
        clusters = [clstr[:] for clstr in self.nodes[self.gene_root].clusters]
        scores = {}
        for clstr_root in self.nodes[self.gene_root].cluster_roots:
            for name in self.nodes[clstr_root].clusters[0]:
                n = self.nodes[clstr_root]
                scores[name] = (self.nodes[clstr_root].total_score, [n.i, n.d, n.l+n.m, n.spread])

        clusters.sort(key=lambda c: -scores[c[0]][0])
        return clusters, scores


    # # # # #  Option parsing and checking methods:
    def parse_gene_tree(self, gene_tree_data, coords_file):
        nwk_str = gene_tree_data.strip()
        gene_root, gene_leaves, _, edges, gene_parent, gene_paths = newick_to_coords.parse_newick_data(nwk_str)
        gene_children = newick_to_coords.find_children(gene_parent)
        if self.use_coords:
            if not coords_file or not os.path.isfile(coords_file):
                self.coords = newick_to_coords.calculate_coordinates(gene_leaves, gene_paths, edges, 0, self.verbose)
            else:
                if self.verbose: print('-- Loading coordinate points from %s...' % coords_file)
                try:
                    self.coords = self._load_coords(coords_file)
                except Exception as err:
                    self.report_error('could not load the distance matrix at %s; the error: %s' % (matrix_file, str(err)))
                if self.coords.shape[0] != len(gene_leaves):
                    self.report_error('the loaded coords do not match the given gene tree')
                if self.verbose: print('-- Complete. Loaded coordinates for %i sequences using %i dimensions' % (self.coords.shape[0], self.coords.shape[1]))
            if coords_file and not os.path.isfile(coords_file):
                self._save_coords(self.coords, coords_file)
                print('Saved coordinate points to %s' % coords_file)
        else:
            self.coords = None
        self.gene_root = gene_root
        self.gene_leaves = gene_leaves
        self.gene_paths = gene_leaves
        self.gene_children = gene_children
    def parse_species_tree(self, species_tree_data):
        _, species, _,_,_, species_paths = newick_to_coords.parse_newick_data(species_tree_data)
        self.species = species
        self.species_paths = species_paths
        self.num_species = len(species)
    def validate_data(self):
        # ensure gene tree is binary.
        for gene in self.gene_leaves:
            if self.gene_leaves.count(gene) > 1:
                self.report_error('%i sequences named "%s" were found in the given gene tree' % (gene_leaves.count(gene), gene))
            elif gene not in self.species_map:
                self.report_error('A sequence named "%s" was not mapped to a species in the information file' % (gene))
            elif self.species_map[gene] not in self.species:
                self.report_error('Species "%s" was not present in the given species tree from the information file' % (spc))


    # # # # #  Misc methods:
    def report_error(self, msg):
        print('\nError: %s.\n' % msg)
        exit()
    def recent_common_ancestor(self, species):
        if len(set(species)) == 1:
            return species[0]
        for nodes in zip(*map(self.species_paths.get, species)):
            if len(set(nodes)) != 1: break
            ancestor = nodes[0]
        return ancestor

    # # # # #  matrix saving / loading methods:
    def _save_coords(self, coords, coords_file):
        with open(coords_file, 'wb') as f:
            np.save(f, coords)
    def _load_coords(self, coords_file):
        coords = np.load(coords_file)
        return coords
    # # # # #  Private methods:
    def _calc_clusters(self, node, params):
        if node in self.gene_leaves:
            species = self.species_map[node]
            n = _PhyloST_node(node, species)
            n.m = self._calc_loss_events((species,))
            n.total_score = n.event_score = n.m * params['m_coef']
            n._separate_event_score = n._combined_event_score = n.event_score
            self.nodes[node] = n
            return n
        child1_name, child2_name = self.gene_children[node]
        child1 = self._calc_clusters(child1_name, params)
        child2 = self._calc_clusters(child2_name, params)
        rca = self.recent_common_ancestor((child1.rca, child2.rca))
        comp = child1.species_comp | child2.species_comp
        d, i, l = self._count_events(child1, child2)
        m = self._calc_loss_events(comp)
        n = _PhyloST_node(node, rca)
        n.species_comp = comp
        n.leaves = child1.leaves | child2.leaves
        n.d, n.i, n.l, n.m = d, i, l, m
        comb_event_score = d*params['d_coef'] + i*params['i_coef'] + l*params['l_coef'] + m*params['m_coef']
        n._combined_event_score = comb_event_score
        n._separate_event_score = child1.event_score + child2.event_score
        ### TEST < vs <= on the below if statement.
        if child1.event_score + child2.event_score < comb_event_score: # Keep as seperate clusters.
            n.total_score = n.event_score = child1.event_score + child2.event_score
            n.clusters = child1.clusters + child2.clusters
            n.cluster_roots = child1.cluster_roots + child2.cluster_roots
        else:  # Merge into one large cluster.
            n.total_score = n.event_score = comb_event_score
            n.clusters = [list(n.leaves)]
            n.cluster_roots = [node]
        self.nodes[node] = n
        return n
    def _count_events(self, child1, child2):
        d, i, l = child1.d+child2.d, child1.i+child2.i, child1.l+child2.l  # The counts of the children.
        if not child1.species_comp & child2.species_comp:  # Some sort of speciation.
            if child1.rca in self.species_paths[child2.rca] \
                or child2.rca in self.species_paths[child1.rca]:
                i += 1  # incongruence event.
            #else: a standard speciation event. Not counted.
        else:  # a duplication event.
            d += 1
            c1_missing = child2.species_comp - child1.species_comp
            c2_missing = child1.species_comp - child2.species_comp
            l += self._calc_loss_events(child1.species_comp, c1_missing)
            l += self._calc_loss_events(child2.species_comp, c2_missing)
        return d, i, l
    def _calc_loss_events(self, present, missing=None):
        """Given a set of ancestors, calculates the mininum number of loss events
        required to explain that set. If 'missing'=None, it is assumed to be all
        species not in 'present'."""
        if missing == None: missing = self.species
        if len(missing) <= 1: return len(missing)
        return max( \
            len(set( self.recent_common_ancestor((g, m)) for m in missing if m not in present)) \
            for g in present)
    def _refine_clusters(self, params):
        spreads = [self.spread_calc(clstr) for clstr in self.nodes[self.gene_root].clusters if len(clstr) != 1]
        if self.relative_avg == 'median':
            spreads.sort()
            med_i = (len(spreads)+1)//2-1
            avg_spread = (spreads[med_i] + spreads[-med_i-1])/2.0
        elif self.relative_avg == 'mean':
            avg_spread = sum(spreads) / float(len(spreads))
        else:
            self.report_error('Incorrect value (%s) given for the "relative_avg" attribute; must be either "mean" or "median".' % self.relative_avg)
        self._recalc_clusters(self.gene_root, params, avg_spread)
    def _recalc_clusters(self, node, params, avg_spread):
        """Overwrites the .clusters and .cluster_roots atts of the nodes, and
        modifies total_score."""
        n = self.nodes[node]
        if node in self.gene_leaves:
            n.clusters = [[node]]
            n.cluster_roots = [node]
            spread_score = self.spread_calc([node], avg_spread)
            n.total_score = n.event_score + spread_score*params['spread_coef']
            return n
        child1, child2 = self.gene_children[node]
        child1 = self._recalc_clusters(child1, params, avg_spread)
        child2 = self._recalc_clusters(child2, params, avg_spread)
        comb_spread_score = self.spread_calc(n.leaves, avg_spread)
        comb_total_score = n._combined_event_score + comb_spread_score*params['spread_coef']
        n._combined_total_score = comb_total_score
        n._separate_total_score = child1.total_score + child2.total_score
        if child1.total_score + child2.total_score < comb_total_score:
            n.total_score = child1.total_score + child2.total_score
            n.clusters = child1.clusters + child2.clusters
            n.cluster_roots = child1.cluster_roots + child2.cluster_roots
            n.spread = None
        else:
            n.total_score = comb_total_score
            n.clusters = [list(n.leaves)]
            n.cluster_roots = [node]
            n.spread = comb_spread_score
        return n

    def _cluster_std_dev(self, cluster, avg_std_dev=None):
        if len(cluster) == 1:
            return self.singleton_spread
        else:
            inds = [self.gene_leaves.index(c) for c in cluster]
            pts = self.coords[inds,:]
            mean = np.average(pts, axis=0)
            std_dev = np.sqrt(np.sum( np.sum((pts-mean)**2, axis=1) ) / float(pts.shape[0]))
        if avg_std_dev == None: # Return the spread score itself.
            return std_dev
        else: # Return the normalized spread score.
            return std_dev / avg_std_dev - 1.0

class _PhyloST_node(object):
    def __init__(self, node, species):
        self.total_score = 0.0 # Score reported by MIPhy.
        self.event_score = 0.0 # Score before refinement step.
        self.node = node
        self.clusters = [[node]]
        self.cluster_roots = [node]
        self.d, self.i, self.l, self.m = 0, 0, 0, 0
        self.spread = None
        self.rca = species
        self.species_comp = set((species,))
        self.leaves = set((node,))
        self._separate_event_score = 0.0 # Not currently used, potentially useful in the future.
        self._combined_event_score = 0.0
        self._separate_total_score = 0.0
        self._combined_total_score = 0.0

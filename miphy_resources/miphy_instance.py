import time
import xml.etree.ElementTree as ET
from miphy_resources.clusterer import Clusterer
from miphy_resources.miphy_common import MiphyValidationError, MiphyRuntimeError
from miphy_resources import phylo


class MiphyInstance(object):
    def __init__(self, gene_tree_data, info_data, gene_tree_format, allowed_wait, use_coords, coords_file, verbose, refine_limit=None):
        self.clusters, self.scores, self.cluster_list, self.init_weights = {}, {}, {}, []
        self.use_coords = use_coords
        self.verbose = verbose
        self.species_tree_data, self.species_mapping = self.parse_species_tree_mapping(info_data)
        if self.verbose: print('Finished parsing the info file.')
        self.species = sorted(list(set(self.species_mapping.values())))
        if gene_tree_format == 'auto':
            gene_tree = phylo.load_tree_string(gene_tree_data)
        elif gene_tree_format == 'newick':
            gene_tree = phylo.load_newick_string(gene_tree_data)
        elif gene_tree_format == 'nexus':
            gene_tree = phylo.load_nexus_string(gene_tree_data)
        elif gene_tree_format == 'phyloxml':
            gene_tree = phylo.load_phyloxml_string(gene_tree_data)
        elif gene_tree_format == 'nexml':
            gene_tree = phylo.load_nexml_string(gene_tree_data)
        self.num_sequences = len(gene_tree.leaves)
        self.tree_data = gene_tree.phyloxml_string(support_values=False, comments=False, internal_names=False)
        self.clusterer = Clusterer(gene_tree, self.species_tree_data, self.species_mapping, use_coords, coords_file, self.verbose)
        self.sequence_names = self.clusterer.gene_leaves
        # # Code to clean dead instances:
        self.been_processed, self.html_loaded = False, False
        self._allowed_wait = allowed_wait # Only used by collect_garbage() in miphy_daemon.py
        self.last_maintained = time.time()

    # # #  Timeout methods:
    def processed(self, params):
        # params should be a tuple of floats in this order: (ils, dup, loss, spread).
        if self.been_processed or self.html_loaded:
            raise MiphyRuntimeError('miphy instance cannot be processed twice, and must be processed before being loaded by the results page.')
        self.init_weights = list(params)
        self.cluster(params)
        self.been_processed = True
        self.last_maintained = time.time()
    def page_loaded(self):
        if not self.been_processed:
            raise MiphyRuntimeError('miphy instance cannot be loaded by the results page before being processed.')
        self.html_loaded = True
        self.last_maintained = time.time()
    def maintain(self):
        if not self.been_processed or not self.html_loaded:
            raise MiphyRuntimeError('miphy instance should not be maintained before being processed and loaded by the results page.')
        self.last_maintained = time.time()
    def cluster(self, params):
        # params should be a tuple of floats in this order: (ils, dup, loss, spread).
        if params not in self.clusters:
            t0 = time.time()
            self.clusters[params], self.scores[params] = self.clusterer.cluster(*params)
            self.cluster_list[params] = []
            for clstr in self.clusters[params]: # Already sorted by descending instability
                clstr_score, clstr_events = self.scores[params][clstr[0]]
                self.cluster_list[params].append([clstr_score, clstr_events, len(clstr), clstr])
            if self.verbose: print('Clustering took %.2f seconds' % (time.time()-t0))
        elif self.verbose:
            print('Clustering pattern retrieved from cache')
    def still_alive(self):
        age = time.time() - self.last_maintained
        if not self.been_processed:
            if age >= self._allowed_wait['after_instance']:
                return False
        elif not self.html_loaded:
            if age >= self._allowed_wait['page_load']:
                return False
        else:
            if age >= self._allowed_wait['between_checks']:
                return False
        return True

    # # #  Data parsing:
    def parse_species_tree_mapping(self, info_data):
        info = self.parse_info_data(info_data)
        species_tree_data = ''.join(info['species tree'])
        species_mapping = {}
        for line in info['species assignments']:
            line = line.strip()
            if not line:
                continue
            elif '=' in line:
                spc, _, genes = line.partition('=')
                spc = spc.strip()
                if spc != ''.join(spc.split()):
                    raise MiphyValidationError("detected a blank in the species name '{}' in the information file. Newick trees cannot contain blanks.".format(spc))
            else:
                genes = line
            for gene in genes.split(','):
                gene = gene.strip()
                if not gene: continue
                if gene in species_mapping: # don't quit, just return error code. miphy.py can quit; server daemon never does.
                    print('In info file %s, gene "%s" was assigned to more than 1 species.' % (info_data, gene))
                    exit()
                species_mapping[gene] = spc
        return species_tree_data, species_mapping
    def parse_info_data(self, info_data):
        info = {}
        group, buff = '', []
        for line in info_data.splitlines():
            line = line.strip()
            if line.startswith('['):
                if group:
                    info[group] = buff
                    buff = []
                group = line[1:-1]
            elif line:
                buff.append(line)
        info[group] = buff
        return info

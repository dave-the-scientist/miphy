import time
import xml.etree.ElementTree as ET
from scripts.convert import NewickToPhyloxml
from scripts.clusterer import Clusterer

# If there is a species in the info file, but the name is different between tree and mapping, it's not caught until way too late. check here.

class MiphyInstance(object):
    def __init__(self, gene_tree_data, info_data, allowed_wait, use_coords, coords_file, verbose, refine_limit=None):
        self.clusters, self.scores, self.cluster_list, self.init_weights = {}, {}, {}, []
        self.use_coords = use_coords
        self.verbose = verbose
        self.species_tree_data, self.species_mapping = self.parse_species_tree_mapping(info_data)
        if self.verbose: print('Finished parsing the info file.')
        self.species = sorted(list(set(self.species_mapping.values())))
        self.tree = Tree(gene_tree_data, self.species_mapping)
        self.num_sequences = self.tree.size
        if refine_limit and self.num_sequences > refine_limit:
            self.use_coords = False
        self.tree_data = self.tree.data
        self.clusterer = Clusterer(gene_tree_data, self.species_tree_data, self.species_mapping, use_coords, coords_file, self.verbose)
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
        species_tree_data = info['species tree'][0]
        species_mapping = {}
        for line in info['species assignments']:
            line = line.strip()
            if not line:
                continue
            elif '=' in line:
                spc, _, genes = line.partition('=')
                spc = spc.strip()
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


class Tree(object):
    def __init__(self, gene_tree_data, species_mapping):
        converter = NewickToPhyloxml(gene_tree_data)
        self.validateTree(converter, species_mapping)
        self.configureTree(converter.etree, species_mapping)
        self.data = converter.tostring()
        self.size = len(converter.leaves)

    def configureTree(self, tree_xml, species_mapping):
        for child in tree_xml:
            tag = child.tag.lower()
            if 'phylogeny' in tag:
                phylogeny = child
                ns = child.tag[ : tag.find('phylogeny')] # xml namespace
                break
        else:
            print('Error parsing phyloxml file of tree')
            exit()
        # Clean out render, charts, styles, if any data was present.
        render = self.subElement(phylogeny, 'render')
        charts = self.subElement(render, 'charts')
        styles = self.subElement(render, 'styles')
        for clade in tree_xml.findall(".//%sname/.." % ns):
            seqID = clade.find("%sname" % ns).text
            spcs = species_mapping[seqID]
            for child in clade:
                if child.tag in ('annotation', "%sannotation" % ns):
                    clade.remove(child)
            chrt = self.subElement(clade, 'chart')
            species = self.subElement(chrt, 'species')
            species.text = spcs+'_style'
    def validateTree(self, converter, species_mapping):
        missing_genes = []
        for gene in converter.leaves:
            count = converter.leaves.count(gene)
            if count != 1:
                raise MiphyValidationError('a sequence named "%s" was not unique in the tree data; it was found %i times' % (gene, count))
            if gene not in species_mapping:
                if gene[0] == gene[-1] == "'" and gene[1:-1] in species_mapping:
                    species_mapping[gene] = species_mapping[gene[1:-1]]
                    del species_mapping[gene[1:-1]]
                else:
                    missing_genes.append(gene)
        if missing_genes:
            missing_genes.sort()
            raise MiphyValidationError('%i sequences from the tree were not found in the species mapping information file:\n%s' % (len(missing_genes), ','.join(missing_genes)))
        if converter.degree[0] != 2 or converter.degree[1] != 2:
            raise MiphyValidationError('the gene tree must be binary; it currently ranges from degree %i to %i. Unrooted trees may appear to be non-binary, so rooting your tree may solve this issue.' % (converter.degree[0],converter.degree[1]))
        if converter.named_internal_nodes == True:
            raise MiphyValidationError('the gene tree included non-terminal nodes with names; these names should be removed.')
    # # #  Private methods:
    def subElement(self, parent, tag):
        for child in parent:
            if child.tag == tag:
                return child
        return ET.SubElement(parent, tag)
    def newChildWithAttribs(self, parent, tag, attribs):
        child = self.subElement(parent, tag)
        for key, val in attribs:
            child.set(key, val)
        return child


class MiphyValidationError(ValueError):
    def __init__(self, *args, **kwargs):
        ValueError.__init__(self, *args, **kwargs)
class MiphyRuntimeError(RuntimeError):
    def __init__(self, *args, **kwargs):
        RuntimeError.__init__(self, *args, **kwargs)

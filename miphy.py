#!/usr/bin/env python
"""
MIPhy - Minimizing Instability in Phylogenetics.

This program reconciles a multi-species phylogenetic tree of some gene family
with the underlying species tree, and clusters the gene tree into phylogenetic
groups. These groups each receive an instability score, used as a proxy of
adaptive evolution, indicating those that are most likely to be interacting with
the environment.
"""
import sys, os, subprocess, threading, webbrowser, socket
from optparse import OptionParser, OptionGroup
from scripts import miphy_daemon
from scripts.miphy_instance import MiphyInstance

import time

__author__ = 'David Curran'
__version__ = '0.8.0'

def setup_parser():
    usage_str = "python %prog GENE_TREE.nwk INFO_FILE.txt [OPTIONS]\n\nPerform clustering on a gene tree, and score each cluster in terms of phylogenetic stability."
    version_str = "%%prog %s" % __version__
    parser = OptionParser(usage=usage_str, version=version_str)
    parser.set_defaults(dup_weight=1.0, ils_weight=0.5, loss_weight=1.0, spread_weight=1.0,
        use_coords=True, coords_file='', results_file='', only_species='',
        manual_browser=False, server_port=0, test=False, verbose=False)
    parser.add_option('-i', '--ILS_weight', dest='ils_weight', type='float',
        help='Cost of an incomplete lineage sorting event [default: %default]')
    parser.add_option('-d', '--duplication_weight', dest='dup_weight', type='float',
        help='Cost of a duplication event [default: %default]')
    parser.add_option('-l', '--loss_weight', dest='loss_weight', type='float',
        help='Cost of a gene loss event [default: %default]')
    parser.add_option('-s', '--spread_weight', dest='spread_weight', type='float',
        help='Weight given to the spread of a cluster [default: %default]')
    parser.add_option('-n', '--no_coords', dest='use_coords', action='store_false',
        help="Don't calculate the full pairwise distance matrix to generate coordinate points. This will cause any SPREAD_WEIGHTs to be ignored")
    parser.add_option('-f', '--coords_file', dest='coords_file', type='string',
        help="Load this file instead of calculating the full pairwise distance matrix and coordinate points. If COORDS_FILE doesn't exist, the coords will be calculated and stored here. IMPORTANT: this must be recalculated if any sequences are added or removed to the gene tree")
    parser.add_option('-r', '--results_file', dest='results_file', type='string',
        help="Save the clustering patterns and instability scores to this file, instead of visualizing the results in a web browser. Use the --only_species option to filter the results")
    parser.add_option('-o', '--only_species', dest='only_species', type='string',
        help='A comma-separated list of species names, exactly as they appear in INFO_FILE. The clustering information from only these species will be saved to RESULTS_FILE; if not supplied all species will be saved. IMPORTANT: --results_file must be supplied with this option.')
    parser.add_option('-m', '--manual_browser', dest='manual_browser', action='store_true',
        help="Starts the MIPhy daemon, but doesn't automatically open the results with your default web browser. A URL will be printed allowing you to access the results using your web browser of choice")
    parser.add_option('-p', '--port', dest='server_port', type='int',
        help='Port used to communicate between MIPhy and the visualization page; setting to 0 will cause your OS to pick one at random [default: %default]')
    parser.add_option('--test', dest='test', action='store_true',
        help="Run MIPhy on included test data and exit")
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true',
        help='Print MIPhy checkpoints [default: %default]')
    return parser
def validate_args(parser, args):
    if not len(args) == 2:
        parser.error('incorrect number of arguments. You must supply a tree in newick format and an information file.')
    gene_tree_file, info_file = os.path.realpath(args[0]), os.path.realpath(args[1])
    if not os.path.isfile(gene_tree_file):
        parser.error('could not locate the gene tree at %s.' % gene_tree_file)
    elif not os.path.isfile(info_file):
        parser.error('could not locate the info file at %s.' % info_file)
    return gene_tree_file, info_file
def validate_options(parser, opts):
    # Validate weights.
    d_weight, i_weight = opts.dup_weight, opts.ils_weight
    l_weight, spread_weight = opts.loss_weight, opts.spread_weight
    if d_weight<0 or i_weight<0 or l_weight<0 or spread_weight<0:
        parser.error('all weight values must be non-negative.')
    # Validate server port.
    server_port = opts.server_port
    if server_port == 0:
        server_port = new_random_port()
    elif not 1024 <= server_port <= 65535:
        parser.error('invalid port; it must fall between 1024 and 65535.')
    elif not port_available(server_port):
        parser.error('unable to start the MIPhy server, as port %i appears to be in use. You may be able to use the command "lsof -i :%i" to identify what is using that port, or specify a different port using the -p flag.' % (server_port, server_port))
    # Validate coords options.
    coords_file = opts.coords_file.strip()
    use_coords = opts.use_coords
    if coords_file:
        coords_file = os.path.abspath(coords_file)
        if not use_coords:
            parser.error('if you specify -n, you should not supply a file with -f')
    if not use_coords:
        spread_weight = 0.0
    # Validate results file options.
    results_file = opts.results_file.strip()
    if results_file:
        results_file = os.path.abspath(results_file)
    only_species = opts.only_species.strip()
    if only_species:
        if not results_file:
            parser.error('if --only_species is given, a results file must be specified with --results_file')
        only_species = only_species.split(',')
    only_species = set(only_species)
    return {'params':(i_weight,d_weight,l_weight,spread_weight), 'server_port':server_port,
        'coords_file':coords_file, 'use_coords':use_coords,
        'results_file':results_file, 'only_species':only_species}

def generate_csv_OLD(only_species, species, species_mapping, scores):
    # scores = {'gene_name': (inst_score, [num_ILS, num_dup, num_loss, rel_spread]), ...}
    if only_species - set(species):
        print('\nError: these species specified from --only_species were not found in the information file: %s' % (', '.join(only_species - set(species))))
        exit()
    to_keep = []
    for name, scrs in scores.items():
        spc = species_mapping[name]
        if only_species and spc not in only_species:
            continue
        to_keep.append((name, spc, scrs[0]))
    to_keep.sort(key=lambda d: d[0])
    to_keep.sort(key=lambda d: d[2])
    buff = ['%s,%s,%.2f'%d for d in to_keep]
    return '\n'.join(buff)
def generate_csv(only_species, mi, params):
    unknown_species = only_species - set(mi.species)
    if unknown_species:
        print('\nError: these species specified by --only_species were not found in the information file: %s' % (', '.join(unknown_species)))
        exit()
    csv_data = []
    for i, clstr_info in enumerate(mi.cluster_list[params]):
        clustID = 'cluster_%i' % i
        clust_score, _, _, seqIDs = clstr_info
        for seqID in sorted(seqIDs, reverse=True):
            spc = mi.species_mapping[seqID]
            if not only_species or spc in only_species:
                csv_data.append('%s,%s,%s,%.2f' % (seqID, spc, clustID, round(clust_score, 2)))
    return '\n'.join(csv_data[::-1])

def port_available(server_port):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('127.0.0.1', server_port))
    if result == 0:
        return False
    else:
        return True
def new_random_port():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('127.0.0.1', 0))
    ip, prt = sock.getsockname()
    if result == 0:
        print('\nError trying to open a random port. Returned code: %i, ip: %s, port: %i' % (result, ip, prt))
        return False
    else:
        return prt

def test_miphy():
    # repeated gene name - test_info.txt & test_unique_names.nwk -> Cre-ugt-18
    # gene missing from info - test_missing_name.txt & test_tree.nwk -> Cbn-ugt-23
    # non-binary tree -
    print('in test. will be implemented soon.')

# Apparently the branch lengths can sometimes be negative in trees.
# - Add a check for this on the upload page.

# Flag to save cluster pattern to file.
# Rearrange help output into option groups.
# Running on the 4cp_ugt trees online, every cluster seems to have a negative relative spread.
    # How can the average be so high when none have high spread?
    # I think I modified some of the spread code. The local version is totally fucked right now.
# When run, check number of sequences, print message warning about long run time if it's large.
    # Also that the web browser is unlikely to ever load. If that is true.
# Create a test module. Begun with a flag, tests setup and imports, tests calculations, tests server, tests visualization.
# Add function to miphy-tools that sanitizes sequence names in nwk file; removes ' marks that FigTree puts there if there is a / in the seq name.
# Make damned sure tree-parsing software is correct and robust.
    # Crashes badly on non-binary trees.
# Make sure this is robust. If it can't handle non-binary trees, ensure the user's trees are binary.
# Local server using IPv4, not IPv6. This a problem? 47044
# Perhaps use some graphical front-end to allow a user to hand-draw their own species tree.
# Clean up newick_to_coords file. Rename, move code into classes, ensure it is useful standalone.
# FOR V2:
# - Add option to specify singleton spread?


if __name__ == '__main__':
    parser = setup_parser()
    (opts, args) = parser.parse_args()
    if opts.test == True:
        test_miphy()
        exit()
    gene_tree_file, info_file = validate_args(parser, args)
    gene_tree_data = open(gene_tree_file).read().strip()
    info_data = open(info_file).read()
    options = validate_options(parser, opts)

    if options['results_file']: # Don't need to start the MIPhy server.
        mi = MiphyInstance(gene_tree_data, info_data, allowed_wait={}, use_coords=options['use_coords'], coords_file=options['coords_file'], verbose=opts.verbose)
        mi.processed(options['params'])
        data = generate_csv_OLD(options['only_species'], mi.species, mi.species_mapping, mi.scores[options['params']])

        print 'testing equal?', data == generate_csv(options['only_species'], mi, options['params'])
        print generate_csv(options['only_species'], mi, options['params'])

        with open(options['results_file'], 'wb') as f:
            f.write(data)
        print('\nInstability scores saved to %s' % options['results_file'])
    else: # Start the MIPhy server.
        daemon = miphy_daemon.Daemon(options['server_port'], web_server=False, instance_timeout_inf=opts.manual_browser, verbose=opts.verbose)
        idnum = daemon.new_instance(gene_tree_data, info_data, use_coords=options['use_coords'], coords_file=options['coords_file'])
        daemon.process_instance(idnum, options['params'])
        results_url = 'http://127.0.0.1:%i/results?%s' % (options['server_port'], idnum)
        if opts.manual_browser:
            print('\nMIPhy daemon started. Use the following URL to access the results:\n%s' % results_url)
        else:
            # This code prevents the browser from continually printing warnings and error messages, which are usually confusing and irrelevant. It was taken from http://stackoverflow.com/questions/2323080/how-can-i-disable-the-webbrowser-message-in-python
            # It works by closing stderr only for the brief second it takes to process the "webbrowser.open(results_url)" call, before restoring it.
            old_stderr = os.dup(2)
            os.close(2)
            os.open(os.devnull, os.O_RDWR)
            try:
                webbrowser.open(results_url)
            finally:
                os.dup2(old_stderr, 2)
        daemon.start_server()

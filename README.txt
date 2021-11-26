# MIPhy - Minimizing Instability in Phylogenetics

MIPhy is software to analyze existing phylogenetic trees, specifically trees reconstructed from large gene family homologs (be they orthologs, paralogs, or even xenologs) from multiple species. This gene tree is reconciled against the underlying species tree, in order to calculate the most likely series of genetic events that led to the observed gene tree.

The tree is then clustered such that every sequence is contained in exactly one MIG (minimum instability group), and the summed score of all genetics events in all MIGs is minimized. There are two useful purposes for doing this: 1) genes within a MIG are more likely to share functional characteristics than genes in separate MIGs, and 2) MIGs containing many genetic events (those with high phylogenetic instability scores) are more likely to be under higher levels of adaptive evolution (related to positive selection).

Visit http://miphy.wasmuthlab.org/ to access the online version of MIPhy, which requires no installation. You will also find a tutorial under the "Documentation" tab of that page, as well as further discussion on how to interpret MIPhy results.

## Installation
The online version of MIPhy has some processing limits, so it can also be installed locally.
- Requires at least Python 2.6, and will work with Python 3 as well.
- Requires that you have installed pip for your version of Python.
- Requires flask, numpy and tkinter to be installed for your version of Python.
  - The installation process will automatically install these dependencies for you if they are not found.

To install, download the [latest release](https://github.com/dave-the-scientist/miphy/releases).
- Unpack the tarball and navigate into the new directory.
- Install with the command "sudo pip install ." (don't include the quotes, but ensure you include the period).

This will install MIPhy for your default version of Python, which is Python 3 for many current computers.
- If you want to install MIPhy for some alternate version of Python, simply replace "pip" in the command above with the appropriate version.
- As an example, "pip2" is often the name of the package manager for Python 2, so you would instead use "sudo pip2 install ."

MIPhy should now be usable from anywhere in your system with the command "miphy.py".
- It is now safe to delete the MIPhy-X.Y.Z tarball and directory.

## Running MIPhy
Whether you use the online or local version of MIPhy, running it requires 2 input files.
- The first is the gene tree. This is a rooted tree in Newick, NEXUS, PhlyXML, or NeXML format. It can be generated using any methods you choose.
  - Unrooted trees may be used with the local version if the "-r" option is used, which skips the visual results page.
- The second is the information file. This is a text file with 2 required sections, and an optional 3rd.
  - The underlying species tree.
    - This should be in Newick format, but does not need or use branch lengths.
    - It is specified by the line "[species tree]", followed by the tree on the next line.
  - The species mappings.
    - This associates each sequence name in the gene tree with its species.
    - It is specified by the line "[species assignments]", followed by a line for each species.
    - Each species line has the form "Species1 = gene1,gene2,gene3,..."
    - The name of each species in the species tree must exactly match the species name in the mapping section, and the name of each gene in the gene tree must exactly match the gene name in the mapping section.
  - The species colours are optional. If not given, the software will assign them. They may also be manually changed on the Results page.
    - It is specified by the line "[species colours]" or "[species colors]", followed by a line for each species.
    - Each species line has the form "Species1 = COLOR".
    - COLOR may be a hex string of the form "#FF9800", a shorthand hex string of the form "#F90", or comma-separated RGB integer values of the form "255, 152, 0".

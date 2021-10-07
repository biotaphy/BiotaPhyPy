===========================
Phylogenetic Beta Diversity
===========================

Description
===========
The ancestral distribution tool (`bin/ancestral_distribution.py`) uses a novel
approach developed by S.A. Smith and B. O'Meara. Given a set of histograms for
species, representing occupancy of environmental space in terms of common bins
(i.e., a PNO or predicted niche occupancy profile), this approach reconstructs
ancestral histograms of occupancy of climate space.

This approach is different from those used previously, based either on (1) summary
statistics (mean, median, maximum, 95th percentile, etc.), or (2) sampling
statistically from present day environmental space. Instead of sampling environmental
space, probabilities of climate occupancy per bin are explicitly reconstructed.
Likewise, unlike summary statistic approaches, which result either in a point
estimate (mean/median) or a minimum and maximum constraint on ancestral
reconstructions (min/max coding), a distribution is explicitly reconstructed
here, revealing the potential shape of ancestral climate space. A key advantage
of this approach is the ability to reconstruct multimodal ancestral
distributions, whereas sampling-based approaches tend to result in normally
distributed ancestral reconstructions regardless of extant species
distributions.

Input species data must have common bins or results will be meaningless.


Using
=====

 ::

    usage: phylo_beta_diversity.py [-h] [-n NUMBER_PERMUTATIONS] [-a ALPHA]
                                   in_tree_filename {newick,nexml,nexus}
                                   pam_filename {csv,json,phylip,table}
                                   {sorensen,jaccard} out_foldername

    Computes phylogenetic & ecological beta diversity components for Sorensen and
    Jaccard Indices.

    positional arguments:
      in_tree_filename      Path to the tree file
      {newick,nexml,nexus}  The format of the tree
      pam_filename          Path to file with presence/ absence data (PAM)
      {csv,json,phylip,table}
                            The format of the PAM
      {sorensen,jaccard}    Beta diversity family metric to calculate
      out_foldername        Write the output of beta diversity calculations to
                            this folder

    optional arguments:
      -h, --help            show this help message and exit
      -n NUMBER_PERMUTATIONS, --number_permutations NUMBER_PERMUTATIONS
                            The number of permuatations to calculate
      -a ALPHA, --alpha ALPHA
                            The alpha value to determine significance

Data formats
============

Alignment data can be provided as CSV [pages/format_csv]_, JSON
[pages/format_csv]_, Phylip [pages/format_phylip]_, or an alignment table.
Tree data can be provided as Newick [pages/format_newick]_, NeXML
[pages/format_nexml]_, or Nexus [pages/format_nexus].

CSV
---
For CSV data, the first row can contain headers for the columns in the file.
Each row should have a header for the taxon that it represents.  An example CSV
alignment file looks like ::

    , var_1, var_2, var_3, var_4, var_5, var_6
    A,   0.9, 0.2, 0.2, 0.3, 0.4, 0.4
    B,  0.01, 0.1, 0.2, 0.3, 0.4, 0.4
    C,   0.8, 0.1, 0.2, 0.3, 0.4, 0.4
    D,   0.3, 0.1, 0.2, 0.3, 0.4, 0.4
    E, 0.001, 0.1, 0.2, 0.3, 0.4, 0.4
    F,  0.11, 0.1, 0.2, 0.3, 0.4, 0.4
    G,  0.99, 0.2, 0.2, 0.3, 0.4, 0.4

JSON
----
If you want to provide alignment data in JSON format, the file should have a
key named "headers" that is an array of headers for each column of data.  It
should also include a key named "values" that is an array of objects with keys
for "name" (taxon name) and "values" (an array of data values).  An example
JSON alignment file looks like ::

    {
       "headers" : [
          "var_1",
          "var_2",
          "var_3",
          "var_4",
          "var_5",
          "var_6"
       ],
       "values" : [
          {
             "name" : "A",
             "values" : [0.9, 0.2, 0.2, 0.3, 0.4, 0.4]
          },
          {
             "name" : "B",
             "values" : [0.01, 0.1, 0.2, 0.3, 0.4, 0.4]
          },
          {
             "name" : "C",
             "values" : [0.8, 0.1, 0.2, 0.3, 0.4, 0.4]
          },
          {
             "name" : "D",
             "values" : [0.3, 0.1, 0.2, 0.3, 0.4, 0.4]
          },
          {
             "name" : "E",
             "values" : [0.001, 0.1, 0.2, 0.3, 0.4, 0.4]
          },
          {
             "name" : "F",
             "values" : [0.11, 0.1, 0.2, 0.3, 0.4, 0.4]
          },
          {
             "name" : "G",
             "values" : [0.99, 0.2, 0.2, 0.3, 0.4, 0.4]
          }
       ]
    }

Phylip
------
Phylip data should be formatted as a list of taxa with corresponding values.
An example phylip alignment file looks like ::

    7 6
    A   0.9 0.2 0.2 0.3 0.4 0.4
    B   0.01 0.1 0.2 0.3 0.4 0.4
    C   0.8 0.1 0.2 0.3 0.4 0.4
    D   0.3 0.1 0.2 0.3 0.4 0.4
    E   0.001 0.1 0.2 0.3 0.4 0.4
    F   0.11 0.1 0.2 0.3 0.4 0.4
    G   0.99 0.2 0.2 0.3 0.4 0.4

Table
-----
You can provide your alignment data as a table as well.  This format looks like
Phylip but does not include metadata for the number of taxa or the number of
data values.  It looks like ::

    A   0.9 0.2 0.2 0.3 0.4 0.4
    B   0.01 0.1 0.2 0.3 0.4 0.4
    C   0.8 0.1 0.2 0.3 0.4 0.4
    D   0.3 0.1 0.2 0.3 0.4 0.4
    E   0.001 0.1 0.2 0.3 0.4 0.4
    F   0.11 0.1 0.2 0.3 0.4 0.4
    G   0.99 0.2 0.2 0.3 0.4 0.4

Newick
------
You can provide your tree data as a Newick file.  You can also request that the
resulting tree be formatted as Newick.  An example Newick file looks like ::

    (A:2.9999,((B:0.1,C:0.1):0.1,(G:0.2,(D:0.1,(E:0.1,F:0.1):0.1):0.1):0.1):0.1);

NeXML
-----
You can provide your tree data as a NeXML file.  You can also request that the
resulting tree be formatted as NeXML.  An example NeXML file looks like ::


    <?xml version="1.0" encoding="ISO-8859-1"?>
    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
        xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
    >
        <otus id="d0">
            <otu id="d1" label="A" />
            <otu id="d2" label="B" />
            <otu id="d3" label="C" />
        </otus>
        <trees id="d4" otus="d0">
            <tree id="d5" xsi:type="nex:FloatTree">
                <node id="d6" />
                <node id="d7" otu="d1" />
                <node id="d8" />
                <node id="d9" otu="d2" />
                <node id="d10" otu="d3" />
                <rootedge id="d11" target="d6" />
                <edge id="d12" source="d6" target="d7" />
                <edge id="d13" source="d6" target="d8" />
                <edge id="d14" source="d8" target="d9" />
                <edge id="d15" source="d8" target="d10" />
            </tree>
        </trees>
    </nex:nexml>

Nexus
-----
You can provide your tree data as a Nexus file.  You can also request that the
resulting tree be formatted as Nexus.  An example Nexus file looks like ::

    #NEXUS

    BEGIN TAXA;
        DIMENSIONS NTAX=7;
        TAXLABELS
            A
            B
            C
            G
            D
            E
            F
      ;
    END;

    BEGIN TREES;
        TREE 1 = (A:2.9999,((B:0.1,C:0.1):0.1,(G:0.2,(D:0.1,(E:0.1,F:0.1):0.1):0.1):0.1):0.1);
    END;


Executable
==========
The phylogenetic beta diversity executable can be found at
`bin/phylo_beta_diversity.py`

Output
======
The phylogenetic beta diversity metrics are printed to the console and written
to the specified folder.

import os
import gzip
import csv
import json
import re

import pandas as pd
import numpy as np
from itertools import combinations

import tskit

from Bio import Phylo
from io import StringIO

RUNS = glob_wildcards("trees/{run}.trees")

rule all: 
    input: 
        expand("results/{run}_traits.json", run=RUNS)

rule align: ########### get_alignment.py {input.ts} {output.alignment}
    message: "generating alignment from simulated network"
    input: 
        ts="trees/{run}.trees"
    output:
        alignment="alignment/{run}.fa"
    shell:
        """
        ts = tskit.load({input.ts})
        with open({output.alignment}, "w") as f:
            f.write(f">{id}\n")
            f.write(a)
            f.write("\n")
        """

rule reference: #######
    message: "generating reference .fa from simulated network"
    input: 
        ts="trees/{run}.trees"
    output:
        reference="reference/{run}.fa"
    shell: 
        """
        ts = tskit.load({input.ts})
        seq = tskit.random_nucleotides(length=3e6, seed=123)
        seq_list = list(seq)
        states_list = [site.ancestral_state for site in ts.sites()]
        sites = ts.tables.sites

        sites_df = pd.DataFrame({
            "position": sites.position,
            "ancestral_state": states_list
        })

        sites_df['position'] = sites_df['position'].astype(int)

        for row in sites_df.itertuples():
            seq_list[row.position] = row.ancestral_state

        res = "".join(seq_list)

        with open({output.reference}, "w") as f:
            f.write(">Reference \n")
            f.write(res)
        """

rule gubbins: 
    message: "running gubbins"
    input:
        alignment="alignment/{run}.fa"
    output:
        gtree=touch("status/notified.done")
    conda: "gubbins"
    shell:
        """
        run_gubbins.py --filter-percentage 100.0 --tree-builder iqtree --best-model --prefix out/ {input.aligment}
        """

rule rescale: # rescale_tree.py 3000000 {input.gtree} {output.tree_raw}#
    message: "rescaling gubbins phylogeny"
    input:
        unscaled_tree="out/{run}.node_labelled.final_tree.tre"
    output:
        scaled_tree="rescaled/{run}.rescaled.tre"
    shell:
        """
        tree = Phylo.read({input.unscaled_tree}, 'newick') 
        for clade in tree.find_clades():
            if clade.branch_length is not None:
                clade.branch_length /= float(3e6)
        
        Phylo.write({output.rescaled_tree}, new_filename, 'newick', format_branch_length="%1.10f")
        """

rule metadata: # get_metadata.py {input.ts} {output.metadata}
    message: "getting metadata for simulated network"
    input:
        ts="trees/{run}.trees"
    output:
        metadata="metadata/{run}.csv"
    shell:
        """
        ts = tskit.load({input.ts})
        rows = []
        for sample in ts.samples():
            pop = ts.population(ts.node(sample).population).metadata["name"]
            pop = int(pop[4:])+1
            rows.append({
                "strain": sample,
                "population": pop,
                "date": "2010-01-01"
            })
        
        res =  pd.DataFrame(rows)
        res.to_csv({output.metadata}, index=False, sep = '\t')
        """

rule refine:
    message: "running augur traits"
    input:
        tree_raw="rescaled/{run}.nwk",
        alignment="alignments/{run}.fa",
        vcf_ref = "vcf_stage/{run}.vcf",
        metadata = "metadata/{run}.tsv"
    output:
        node_data="results/{run}_branch_lengths.json",
        tree="results/{run}.nwk",
    shell:
        """
        augur refine \
        --alignment {input.alignment} \
        --tree {input.tree_raw} \
        --output-tree {output.tree} \
        --output-node-data {output.node_data} \
        --vcf-reference {input.vcf_ref} \
        --keep-root \
        --clock-rate 1.0 \
        --timetree \
        --stochastic-resolve \
         --resolve-polytomies \
        --metadata {input.metadata} 
        """

rule traits:
    message: "running augur traits"
    input:
        tree = "results/{run}.nwk",
        metadata = "metadata/{run}.tsv"
    output:
        node_data = "results/{run}_traits.json"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --columns population \
            --confidence \
            --output-node-data {output.node_data}
        """


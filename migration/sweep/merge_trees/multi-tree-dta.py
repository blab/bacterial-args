#!/usr/bin/env python3

"""
Take a set of newick trees and metadata assigning states to their tips and
perform a DTA reconstruction by first joining them into a single multi-tree
with a zero-branch-length polytomy at the root.

Data assumptions: trees must have all nodes labelled, including internal ones,
and labels must be integers; similarly, states (demes) are also integers.
"""

import argparse
import glob
import os
import re
import json
import subprocess
from dataclasses import dataclass
from treetime.wrappers import reconstruct_discrete_traits
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import csv

from Bio import Phylo


@dataclass
class TreeFile:
    filename: str
    index: int
    tree: Phylo.BaseTree.Tree

@dataclass
class Metadata:
    original: dict[str,int]
    alphabetical: dict[str,str]
    multi_tree: dict[str,str]

BASIC_ALPHABET="ABCDEFGHIJKLMONPQRSTUVWXYZ"

def load_metadata(metadata_file: str, trees: list[TreeFile]) -> Metadata:
    """
    Read a simple TSV file which must map node names (string) to integer state (int).
    This is returned as the **original** property
    We then convert this integer state to 0: A, 1: B etc, as this is nicer to work with.
    This is the **alphabetical** property.
    Because the original node names have changed in the multi-tree (so that each
    tip is unique), we create a mapping of multi-tree node-names to alphabetical state.
    This is the **multi_tree** property.
    """
    original: dict[str,int] = {}
    alphabetical: dict[str,str] = {}
    multi_tree: dict[str,str] = {}
    with open(metadata_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            node_name = row[reader.fieldnames[0]]
            state = int(row[reader.fieldnames[1]])
            original[node_name] = state
            alphabetical[node_name] = BASIC_ALPHABET[state]
            for tree in trees:
                multi_node_name = node_name_multi_tree(nice_node_name(node_name), tree.index)
                multi_tree[multi_node_name] = BASIC_ALPHABET[state]
    print(f"Original (tip) states: {f', '.join({str(s) for s in original.values()})}")
    print(f"Alphabetised (tip) states: {f', '.join({s for s in alphabetical.values()})}")

    
    return Metadata(original=original, multi_tree=multi_tree, alphabetical=alphabetical)


def load_trees(data_folder: str) -> list[TreeFile]:
    """
    Read trees. For each tree, assert all nodes are labelled and
    convert to `nodeXXX` naming format to avoid having to deal with
    integer types
    """
    pattern = os.path.join(data_folder, '*.nwk')
    results = []
    for path in sorted(glob.glob(pattern)):
        basename = os.path.basename(path)
        match = re.match(r'^index(\d+)', basename)
        if not match:
            raise Exception("Unexpected tree filename. Must start with 'index' + number")
        index = int(match.group(1))
        tree = Phylo.read(path, 'newick', rooted=True)
        for clade in tree.find_clades():
            if clade.name is None and clade.confidence is not None:
                clade.name = str(int(clade.confidence))
                clade.confidence = None
            assert clade.name is not None, f"Unnamed node found in {basename}"
            clade.name = nice_node_name(clade.name)
        results.append(TreeFile(filename=basename, index=index, tree=tree))
    return results

def nice_node_name(original_name: str|int) -> str:
    """Turn integer node names into nodeXXX"""
    return f"node{int(original_name):>03}"

def node_name_multi_tree(node_name: str, index: int) -> str:
    "Take 'nodeXXX' and turn it to 'indexYY_nodeXXX'"
    return f"tree{index:>02}_{node_name}"

def node_name_ARG(full_node_name:str) -> tuple[int, int]:
    [idx, node_name] = full_node_name.split("_")
    return (int(idx.replace("tree", "")), int(node_name.replace("node", "")))

def modify_node_names(trees: list[TreeFile]) -> dict[str, list[tuple[int, bool]]]:
    """
    Change the name of nodes in the trees to be unique via the `node_name_multi_tree`
    function. This allows us to combine the trees together into a multi-tree and
    have all node-names be unique.
    """
    original_node_names: dict[str, list[tuple[int, bool]]] = {}
    for tf in trees:
        for node in tf.tree.find_clades():
            original_node_names.setdefault(node.name, []).append((tf.index, node.is_terminal()))
            node.name = node_name_multi_tree(node.name, tf.index)

    tree_indices = {tf.index for tf in trees}

    # Ensure terminal nodes are present (and terminal) in every tree
    for name, entries in original_node_names.items():
        if not any(terminal for _, terminal in entries):
            continue
        if not all(terminal for _, terminal in entries):
            raise Exception(f"Node '{name}' is terminal in some trees but internal in others")
        indices_with_node = {idx for idx, _ in entries}
        if indices_with_node != tree_indices:
            missing = tree_indices - indices_with_node
            raise Exception(f"Terminal node '{name}' is missing from tree(s): {sorted(missing)}")

    # ASCII histogram: for each internal node name, how many trees is it present in?
    from collections import Counter
    n_trees = len(trees)
    internal_counts = Counter()
    for name, entries in original_node_names.items():
        if any(terminal for _, terminal in entries):
            continue
        internal_counts[len(entries)] += 1
    max_count = max(internal_counts.values()) if internal_counts else 0
    bar_width = 40
    print("Number of trees each internal node is present in:")
    for n in range(1, n_trees + 1):
        count = internal_counts.get(n, 0)
        bar = '#' * (round(count / max_count * bar_width) if max_count else 0)
        print(f"  {n:>3} tree(s): {bar} {count}")
    print()
    return original_node_names


def check_all_terminal_nodes_have_metadata(tree: Phylo.BaseTree.Tree, metadata: Metadata):
    terminal_names = {node.name for node in tree.get_terminals()}
    metadata_names = set(metadata.multi_tree.keys())
    missing_metadata = terminal_names - metadata_names
    missing_tree = metadata_names - terminal_names    
    if missing_metadata:
        raise Exception(f"Terminal nodes missing from metadata: {sorted(missing_metadata)}")
    if missing_tree:
        raise Exception(f"Metadata entries missing from tree: {sorted(missing_tree)}")


@dataclass
class DTAResults:
    node_states: dict[str, str]
    node_confidence: dict[str, dict[str, float]]
    node_entropy: dict[str, float]
    alphabet: list[str]
    equilibrium_probabilities: list[float]
    transition_matrix: list[list[float]]


TINY = 1e-12

def run_dta(tree: Phylo.BaseTree.Tree, metadata: Metadata) -> DTAResults:
    """
    Run discrete trait reconstruction on the multi-tree using TreeTime.
    """
    print("\nRunning DTA")
    traits = metadata.multi_tree

    # sampling_bias_correction and weights are supported by reconstruct_discrete_traits
    # but not yet wired up here
    tt, letter_to_state, reverse_alphabet = reconstruct_discrete_traits(
        tree, traits,
    )
    if tt is None:
        raise Exception("TreeTime discrete trait reconstruction failed")

    node_states: dict[str, str] = {}
    node_confidence: dict[str, dict[str, float]] = {}
    node_entropy: dict[str, float] = {}

    for node in tt.tree.find_clades():
        node_states[node.name] = letter_to_state[node.cseq[0]]

        pdis = node.marginal_profile[0]
        node_entropy[node.name] = float(-np.sum(pdis * np.log(pdis + TINY)))

        marginal = {
            letter_to_state[tt.gtr.alphabet[i]]: float(pdis[i])
            for i in range(len(tt.gtr.alphabet))
            if pdis[i] > 0.001
        }
        node_confidence[node.name] = marginal

    return DTAResults(
        node_states=node_states,
        node_confidence=node_confidence,
        node_entropy=node_entropy,
        alphabet=[letter_to_state[k] for k in sorted(letter_to_state.keys())],
        equilibrium_probabilities=[float(x) for x in tt.gtr.Pi],
        transition_matrix=[[float(x) for x in row] for row in tt.gtr.W],
    )


def combine_trees(trees: list[TreeFile], fname: str|None) -> Phylo.BaseTree.Tree:
    print(f"Combining {len(trees)} trees together into a multi-tree with a zero-length polytomy as the new root")
    root_clade = Phylo.BaseTree.Clade(name="ROOT")
    for tf in trees:
        tf.tree.root.branch_length = 0.0
        root_clade.clades.append(tf.tree.root)

    multi_tree = Phylo.BaseTree.Tree(root=root_clade, rooted=True)

    if fname:
        print("\tWriting multi-tree to", fname)
        Phylo.write(multi_tree, fname, 'newick')
        
    return multi_tree

def node_data_json(dta: DTAResults, fname: str) -> None:
    print("\nWriting node-data file to", fname)
    nodes = {}
    branches = {}
    for name, state in dta.node_states.items():
        if name=="ROOT":
            tree_idx, original_node_name = ['', '']
        else:
            tree_idx, original_node_name = node_name_ARG(name)
        nodes[name] = {
            "state": state,
            "state_confidence": dta.node_confidence[name],
            "state_entropy": dta.node_entropy[name],
            "tree_idx": f"index_{tree_idx}", # string so auspice uses categorial scale
            "original_node_name": f"node_{original_node_name}", # string so auspice uses categorial scale
        }
        branches[name] = {
            'labels': {
                "multi-tree-name": name,
                "original_node_name": original_node_name,
            }
        }
    data = {
        "model": {
            "alphabet": dta.alphabet,
            "equilibrium_probabilities": dta.equilibrium_probabilities,
            "transition_matrix": dta.transition_matrix,
        },
        "nodes": nodes,
        "branches": branches,
    }
    
    with open(fname, 'w') as fh:
        json.dump(data, fh, indent=2)

def states_tsv(dta: DTAResults, fname: str) -> None:
    print("Writing states TSV to", fname)
    state_to_int = {s: i for i, s in enumerate(BASIC_ALPHABET)}
    with open(fname, 'w') as fh:
        fh.write("TREE_INDEX\tNODE_NAME\tSTATE\tCONFIDENCE\n")
        for name, state in dta.node_states.items():
            if name == "ROOT":
                continue
            tree_idx, node_name = node_name_ARG(name)
            confidence = dta.node_confidence[name].get(state, 0.0)
            fh.write(f"{tree_idx}\t{node_name}\t{state_to_int[state]}\t{confidence}\n")

def basic_stats(dta: DTAResults, multi_tree: Phylo.BaseTree.Tree, fname: str) -> None:
    from collections import defaultdict, Counter
    print("\nBasic stats across trees (internal nodes only)")

    states_by_node: dict[int, list[str]] = defaultdict(list)
    for node in multi_tree.find_clades():
        if node.name == "ROOT" or node.is_terminal():
            continue
        _tree_idx, node_name = node_name_ARG(node.name)
        states_by_node[node_name].append(dta.node_states[node.name])

    agreements = []
    for node_name, states in sorted(states_by_node.items()):
        if len(states) <= 1: # ignore internal nodes which only appear in one tree
            continue
        # most_common(1) returns [(state, count)], so [0][1] extracts the count
        most_common_count = Counter(states).most_common(1)[0][1]
        agreements.append(most_common_count / len(states) * 100)

    print(f"  {len(states_by_node)} unique internal nodes across all trees")
    print(f"  {len(agreements)} present in >1 tree")
    print(f"  Mean agreement: {np.mean(agreements):.1f}%")
    print(f"  Median agreement: {np.median(agreements):.1f}%")

    fig, ax = plt.subplots()
    # 5% wide bins; 105 so that 100% is included in the final bin
    ax.hist(agreements, bins=np.arange(0, 105, 5), edgecolor="black")
    ax.set_xlabel("Agreement (%)")
    ax.set_ylabel("Number of internal nodes")
    ax.set_title("State agreement across trees for internal nodes")
    fig.tight_layout()
    fig.savefig(fname)
    print(f"  Histogram saved to {fname}")
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trees', type=str, required=True, help='Path to folder with newick files (all newick files will be read)')
    parser.add_argument('--output-prefix', type=str, help='Output prefix (multiple output files)')
    parser.add_argument('--metadata', type=str, required=True, help='TSV file with node states')
    args = parser.parse_args()
    output = Path(args.output_prefix)
    tree_file = str(output) + ".multi-tree.nwk"
    node_data_file = str(output) + ".multi-tree-node-data.json"
    auspice_file = str(output) + ".multi-tree-auspice.json"
    states_file = str(output) + ".states.tsv"
    stats_plot_file = str(output) + ".agreement.png"

    trees = load_trees(args.trees)
    metadata = load_metadata(args.metadata, trees)
    original_node_names = modify_node_names(trees)
    multi_tree = combine_trees(trees, tree_file)
    check_all_terminal_nodes_have_metadata(multi_tree, metadata)
    dta = run_dta(multi_tree, metadata)
    node_data_json(dta, node_data_file)
    states_tsv(dta, states_file)
    basic_stats(dta, multi_tree, stats_plot_file)

    print("\nNow using `augur export` to create an auspice dataset to visualise the multi-tree")
    subprocess.run([
        "augur", "export", "v2",
        "--tree", tree_file,
        "--node-data", node_data_file,
        "--output", auspice_file,
    ], check=True)

    
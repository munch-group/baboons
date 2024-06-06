import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('retina', 'png')

from ete3 import Tree#, NodeStyle, TreeStyle, TextFace


sns.set()
sns.set_style("ticks")

## UPGMA ##############################

def find_lowest_cell(table):
    x = 1
    y = 0
    min_val = table[x][y]
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_val:
                min_val = table[i][j]
                x = i
                y = j
    return [x, y]

def link(x, y, wx, wy):
    return (x * wx + y * wy) / (wx + wy)

def update_table(table, a, b, weight_a, weight_b):
    for i in range(0, b):
        table[b][i] = link(table[b][i], table[a][i], weight_b, weight_a)
    for j in range(b+1, a):
        table[j][b] = link(table[j][b], table[a][j], weight_b, weight_a)
    for i in range(a+1, len(table)):
        table[i][b] = link(table[i][b], table[i][a], weight_b, weight_a)
    for i in range(a+1, len(table)):
        del table[i][a]
    del table[a] 

def update_labels(labels, i, j, di, dj):
    labels[j] = "({}:{},{}:{})".format(labels[j], dj, labels[i], di)
    del labels[i]

def upgma(mat, names):

    table = mat[:]
    labels = names[:]
    node_heights = [0 for _ in labels]

    while len(labels) > 1:
        i, j = find_lowest_cell(table)
        
        dist = table[i][j]

        wi = max(1, labels[i].count(':'))
        wj = max(1, labels[j].count(':'))

        new_node_height = dist / 2
        di = new_node_height - node_heights[i]
        dj = new_node_height - node_heights[j]
        
        update_table(table, i, j, wi, wj)
        update_labels(labels, i, j, di, dj)
        node_heights[j] = new_node_height
        del node_heights[i]
        
    return labels[0] + ';'


## Distance matrix ##############################

def jc(seq1, seq2, case_sensitive=True):
    assert len(seq1) == len(seq2), (len(seq1), len(seq2))
    tot, diff = 0, 0
    for x, y in zip(seq1, seq2):
        if not case_sensitive:
            x, y = x.upper(), y.upper()
        if x in 'ATCG' and y in 'ATCG':
            if x != y:
                diff += 1
            tot += 1    
    if not tot:
        return np.nan
    elif diff:
        distance = -3/4 * math.log(1 - 4/3 * diff/tot)
    else:
        distance = 0.0
    return distance


def dist_matrix(seq_list, case_sensitive=True, dist_fun=jc):
    n = len(seq_list)
    mat = np.zeros((n, n))

    upper_trag_idx = list(zip(*np.triu_indices(n, k=1)))
#     seq_pairs = [(seq_list[i], seq_list[j]) for i, j in upper_trag_idx]
    args = [(seq_list[i], seq_list[j], case_sensitive) for i, j in upper_trag_idx]

    jc_distances = []
    for a, b, c in args:
        jc_distances.append(dist_fun(a, b, c))
    
    # with Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as p:
    #     jc_distances = p.starmap(dist_fun, args)

    for (i, j), d in zip(upper_trag_idx, jc_distances):
        mat[i][j] = d
        mat[j][i] = mat[i][j]

    return mat


def prune_nans(mat, name_list):

    # Find rows (and cols) where all off-diagonal entries are nan:
    mask = np.isnan(mat).sum(axis=0) == np.size(mat, axis=0)-1
    delete = [i for (i, delete) in enumerate(mask) if delete]
    # remove those rows and cols:
    mat = np.delete(mat, delete, axis=0)
    mat = np.delete(mat, delete, axis=1)
    # update name list:
    name_list = [name for (i, name) in enumerate(name_list) if i not in delete]
    
    # if there are any nans left we have to remove all row/cols with a nan:
    if np.isnan(mat).any():    
        # in that case we could do this to remove all row/cols with a nan
        mask = np.isnan(mat).sum(axis=0).astype('bool')
        delete = [i for (i, delete) in enumerate(mask) if delete]
        # remove those rows and cols:
        mat = np.delete(mat, delete, axis=0)
        mat = np.delete(mat, delete, axis=1)
        # update name list:
        name_list = [name for (i, name) in enumerate(name_list) if i not in delete]
        
    return mat, name_list

## Construct and manipulate trees ##############################


def tree_newick(name_list, seq_list, case_sensitive=True, dist_fun=jc):

    mat = dist_matrix(seq_list, case_sensitive=case_sensitive, dist_fun=dist_fun)
    mat, name_list = prune_nans(mat, name_list)
    
    lowtri = [lst[:i] for (i, lst) in enumerate(mat.tolist())]
    newick_str = upgma(lowtri, name_list)
    return newick_str

def order_tree(t, key=lambda c: len(c.get_leaves())):
    """
    without key for sorting it makes a comb tree by default
    """
    for node in t.traverse():
        if not node.is_leaf():
            node.children = sorted(node.children, key=key, reverse=False)

def remove_outgroup(name):
    tree.set_outgroup( tree&name )
    all_leaves = tree.get_leaf_names()
    all_leaves.remove(name)
    tree.prune(all_leaves, preserve_branch_length=False)


# # remove chimp outgroup branches
# tree.set_outgroup( tree&"Outgroup" )
# all_leaves = tree.get_leaf_names()
# all_leaves.remove('Outgroup')
# tree.prune(all_leaves, preserve_branch_length=False)


## Plot tree ############################

def plot_tree(t, ax, leaf_colors=None, show_inner_nodes=False, fontsize=10, 
              text_offset=None, margins=(0.5, 1, 0.5, 1)): # top, right, bottom, left

    y_offset = len(t.get_leaves())
    for node in t.traverse("preorder"):
        node.x_offset = node.dist + sum(x.dist for x in node.get_ancestors())
        if node.is_leaf():
            y_offset -= 1
            node.y_offset = y_offset

    for node in t.traverse("postorder"):
        if not node.is_leaf():
            node.y_offset = sum(x.y_offset for x in node.children) / len(node.children)

    horizontal_lines = list()
    vertical_lines = list()
    node_coords = list()
    leaf_coords = list()
    max_x_offset = 0
    for node in t.traverse("postorder"):
        max_x_offset = max(max_x_offset, node.x_offset)
        node_coords.append((node.x_offset, node.y_offset))
        if node.is_leaf():
            leaf_coords.append([node.name, node.x_offset, node.y_offset])
        if not node.is_root():
            y = node.y_offset
            horizontal_lines.append(([node.up.x_offset, node.x_offset], [y, y]))
        if not node.is_leaf():
            c = sorted(node.children, key=lambda x: x.y_offset)
            bottom, top = c[0], c[-1]
            x = node.x_offset
            vertical_lines.append(([x, x],[bottom.y_offset, top.y_offset]))

    
    # shift the tree to put leaves at zero
    for i in range(len(horizontal_lines)):
        horizontal_lines[i][0][0] -= max_x_offset
        horizontal_lines[i][0][1] -= max_x_offset
    for i in range(len(vertical_lines)):
        vertical_lines[i][0][0] -= max_x_offset
        vertical_lines[i][0][1] -= max_x_offset
    for i in range(len(leaf_coords)):
        leaf_coords[i][1] -= max_x_offset
            
    # draw the tree:
    for x in horizontal_lines:
        ax.plot(*x, c='black', linewidth=0.8)
    for x in vertical_lines:
        ax.plot(*x, c='black', linewidth=0.8)

#     for tup in node_coords:
#         ax.plot(*tup, c='black', marker="o")

    if text_offset is None:
        text_offset = max_x_offset / 20
        
    for name, x, y in leaf_coords:
        ax.text(x+text_offset, y, name, fontsize=fontsize,
                verticalalignment='center', horizontalalignment='left')
        if leaf_colors is None:
            color = 'black'
        else:
            color = leaf_colors[name]
        ax.plot(x, y, c=color, marker="o", ms=3)


#     ax.set_xlim(-margins[3], max_x_offset + margins[1])
    ax.set_xlim(-margins[3]-max_x_offset, margins[1])
    ax.set_ylim(-margins[2], len(leaf_coords)-1+margins[0])


    #ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    ax.spines['top'].set_visible(False) 
    ax.spines['left'].set_visible(False) 
    ax.spines['right'].set_visible(False)
    
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    
    return leaf_coords



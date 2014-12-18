
def RF_to_PG(rfgraph, non_dups):
    for cur_node in non_dups:
        DFSExpand(cur_node, rfgraph)


def DFSExpand(cur_node, rfgraph):
    if cur_node.duplicate

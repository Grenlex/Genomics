import msprime
import tskit
import math
from IPython.display import SVG
from scipy.stats import poisson as pois


class Edge:
    def __init__(self, _child, _parent, _child_time, _parent_time, _left, _right):
        self.parent = _parent
        self.child = _child
        self.parent_time = _parent_time
        self.child_time = _child_time
        self.length = _parent_time - _child_time
        self.left = _left
        self.right = _right
        self.mutations = []

    def print(self):
        print(
            "Edge {} -> {}, child time: {:.2f}, parent time: {:.2f}, lenght: {:.2f}, interval: ({:.2f}-{:.2f}), number of mutatuions: {}".format(
                self.child,
                self.parent,
                self.child_time,
                self.parent_time,
                self.length,
                self.left,
                self.right,
                len(self.mutations)
            ))

    def __str__(self):
        return "{}-{}".format(self.child, self.parent)

    def __repr__(self):
        return "{}-{}".format(self.child, self.parent)


def mark_mutations(edges, mutations_and_times):
    edge_indexes_by_child = {}
    for i in range(len(edges)):
        edge = edges[i]
        edge_indexes_by_child.setdefault(edge.child, [])
        edge_indexes_by_child[edge.child].append(i)

    for key in edge_indexes_by_child:
        edge_indexes_by_child[key].sort(key=lambda index: edges[index].left)

    for mutation, time in mutations_and_times:
        key = mutation.node
        indexes = edge_indexes_by_child[key]
        l_ind = 0  # that is true for all time: _edges[indexes[l_ind]].left <= time
        r_ind = len(indexes)  # that is true for all time: _edges[indexes[r_ind]].left > time
        while r_ind - l_ind > 1:
            m_ind = (l_ind + r_ind) // 2
            if edges[indexes[m_ind]].left <= time:
                l_ind = m_ind
            else:
                r_ind = m_ind
        edges[indexes[l_ind]].mutations.append(mutation)


tree_sequence = msprime.simulate(sample_size=10, Ne=1000, length=1e4, recombination_rate=2e-8, mutation_rate=6.83e-8)
print(
    "\n\033[41m                                                            GENOMICS                                                            \033[0m\n")
print("\n\033[43mSimulating trees\033[0m\n")
for tree in tree_sequence.trees():
    print("tree {}: interval = {}".format(tree.index, tree.interval))
    display(SVG(tree.draw(format="svg")))
    print("\033[34mTotal branch length of the tree: {}\033[0m\n".format(tree.total_branch_length))
    print("\n\033[40m" + "-" * 100 + "\033[0m\n")
t = tskit.EdgeTable()

for elem in tree_sequence.edges():
    t.add_row(elem.left, elem.right, elem.parent, elem.child)

t.squash()

edges = []
for edge in t:
    child = edge.child
    parent = edge.parent
    child_time = tree_sequence.node(child).time
    parent_time = tree_sequence.node(parent).time
    e = Edge(
        child, parent,
        child_time, parent_time,
        edge.left, edge.right
    )
    edges.append(e)

print("\n\033[40m\033[37m{}\033[0m".format(t))

mutations_and_times = []
for site in tree_sequence.sites():
    mutations_and_times.extend([(mutation, site.position) for mutation in site.mutations])

mark_mutations(edges, mutations_and_times)

for edge in edges:
    edge.print()
F = 0.0
for edge in edges:
    F += math.log(pois.pmf(len(edge.mutations), [6.83e-8 * edge.length*(edge.right - edge.left)]))
print(F)

for edge in edges:

    edge.print()

for site in tree.sites():

        for mutation in site.mutations:

            print("Mutation @ position {:.2f} over node {}".format(site.position, mutation.node))

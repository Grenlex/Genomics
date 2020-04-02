import msprime

import tskit

from IPython.display import SVG

 

 

class Edge:

    def __init__(self, _child, _parent, _child_time, _parent_time, _left, _right):

        self.parent = _parent

        self.child = _child

        self.parent_time = _parent_time

        self.child_time = _child_time

        self.length = _parent_time - _child_time

        self.left = _left

        self.right = _right

 

    def print(self):

        print(

            "Edge {} -> {}, child time: {:.2f}, parent time: {:.2f}, lenght: {:.2f}, interval: ({:.2f}-{:.2f})".format(

                self.child,

                self.parent,

                self.child_time,

                self.parent_time,

                self.length,

                self.left,

                self.right))

 

    def __str__(self):

        return "{}-{}".format(self.child, self.parent)

 

    def __repr__(self):

        return "{}-{}".format(self.child, self.parent)

 

 

tree_sequence = msprime.simulate(sample_size=10, Ne=1000, length=1e4, recombination_rate=2e-8, mutation_rate=6.83e-8)

 

print("\n\033[41m                                                            GENOMICS                                                            \033[0m\n")

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

    edges.append(Edge(

        child, parent,

        child_time, parent_time,

        edge.left, edge.right

    ))

print("\n\033[40m\033[37m{}\033[0m".format(t))

for edge in edges:

    edge.print()

for site in tree.sites():

        for mutation in site.mutations:

            print("Mutation @ position {:.2f} over node {}".format(site.position, mutation.node))

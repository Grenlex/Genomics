import msprime
import tskit
import math
from IPython.display import SVG
from scipy.stats import poisson as pois

SAMPLE_SIZE = 10
Ne = 1000
LENGTH = 5e4
RECOMBINATION_RATE = 2e-8
MUTATION_RATE = 6.83e-8

try:
    display
except:
    def display(a):
        print(a)


class Node:
    def __init__(self, _id, _real_time):
        self.id = _id
        self.real_time = _real_time
        self.down_edges = []
        self.up_edges = []
        self.time = -1

    def __str__(self):
        return str(self.id)


class Edge:
    def __init__(self, _id, _child, _parent, _left, _right):
        self.id = _id
        self.parent = _parent
        self.child = _child
        self.left = _left
        self.right = _right
        self.mutations = []

    def length(self):
        if self.parent.real_time == -1:
            return 0
        return self.parent.real_time - self.child.real_time

    def im_length(self):
        if self.parent.real_time == -1:
            return 0
        return self.parent.time - self.child.time

    def mutations_number(self):
        return len(self.mutations)

    def print(self):
        print(self)

    def __str__(self):
        return "Edge {} -> {}, child time: {:.2f}, parent time: {:.2f}, lenght: {:.2f}, interval: ({:.2f}-{:.2f}), number of mutations: {}".format(
            self.child,
            self.parent,
            self.child.real_time,
            self.parent.real_time,
            self.length(),
            self.left,
            self.right,
            self.mutations_number()
        )

    def __repr__(self):
        return "{}-{}".format(self.child, self.parent)


def createRoot(nodes, edges):
    rootId = max(nodes) + 1
    freeEdgeId = max(edges) + 1

    root = Node(rootId, -1)
    for node in nodes.values():
        if len(node.up_edges) == 0:
            edge = Edge(freeEdgeId, node, root, 0, LENGTH)
            edges[freeEdgeId] = edge
            freeEdgeId += 1
            node.up_edges.append(edge)
            root.down_edges.append(edge)
    nodes[rootId] = root
    return root


def markTimes(root):
    if root.time != -1:
        return

    for downEdge in root.down_edges:
        child = downEdge.child
        markTimes(child)
        root.time = max(root.time, child.time + 1)

    if root.time == -1:  # leaf
        root.time = 0


def markMutations(edges, nodes, mutationsAndTimes):
    ##if fact for mutation "time" is its place in DNA
    for node in nodes.values():
        node.up_edges.sort(key=lambda edge: edge.left)

    for mutation, time in mutationsAndTimes:
        ##print("Mutation on node {} with time {}".format(mutation.node, time))
        key = mutation.node
        edgesWithCurrentChild = nodes[key].up_edges
        l_ind = 0  # that is true for all time: edgesWithCurrentChild[l_ind].left <= time
        r_ind = len(edgesWithCurrentChild)  # that is true for all time: edgesWithCurrentChild[r_ind].left > time
        while r_ind - l_ind > 1:
            m_ind = (l_ind + r_ind) // 2
            if edgesWithCurrentChild[m_ind].left <= time:
                l_ind = m_ind
            else:
                r_ind = m_ind

        mutatedEdge = edgesWithCurrentChild[l_ind]
        assert mutatedEdge.left <= time < mutatedEdge.right

        mutatedEdge.mutations.append(mutation)
        ##print("Find edge for mutation: {}".format(mutatedEdge))


treeSequence = msprime.simulate(sample_size=SAMPLE_SIZE, Ne=Ne, length=LENGTH, recombination_rate=RECOMBINATION_RATE,
                                mutation_rate=MUTATION_RATE)
print(
    "\n\033[41m                                                            GENOMICS                                                            \033[0m\n")
print("\n\033[43mSimulating trees\033[0m\n")
for tree in treeSequence.trees():
    print("tree {}: interval = {}".format(tree.index, tree.interval))
    display(SVG(tree.draw(format="svg")))
    print("\033[34mTotal branch length of the tree: {}\033[0m\n".format(tree.total_branch_length))
    print("\n\033[40m" + "-" * 100 + "\033[0m\n")

edges = {}
nodes = {}

for node in treeSequence.nodes():
    nodes[node.id] = Node(node.id, node.time)

for edge in treeSequence.edges():
    edgeId = edge.id
    child = nodes[edge.child]
    parent = nodes[edge.parent]
    e = Edge(
        edgeId,
        child, parent,
        edge.left, edge.right
    )
    child.up_edges.append(e)
    parent.down_edges.append(e)
    edges[edgeId] = e

root = createRoot(nodes, edges)

markTimes(root)

mutationsAndTimes = []
for site in treeSequence.sites():
    for mutation in site.mutations:
        ##print("Mutation @ position {:.2f} over node {}".format(site.position, mutation.node))
        mutationsAndTimes.append((mutation, site.position))

markMutations(edges, nodes, mutationsAndTimes)

for edge in edges.values():
    edge.print()

F_real = 0.
for edge in edges.values():
    F_real += pois.logpmf(edge.mutations_number(), MUTATION_RATE * edge.length() * (edge.right - edge.left))
print("Real time function defenition", F_real)

F_im = 0.
for edge in edges.values():
    F_im += pois.logpmf(edge.mutations_number(), MUTATION_RATE * edge.im_length() * (edge.right - edge.left))
print("Imaginary time function defenition", F_im)

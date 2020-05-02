import msprime
from IPython.display import SVG, display
from scipy.stats import poisson as pois

SAMPLE_SIZE = 10
Ne = 1000
LENGTH = 5e4
RECOMBINATION_RATE = 2e-8
MUTATION_RATE = 6.83e-8


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


def f_(node, t0):
    sum1 = 0
    for downEdge in node.down_edges:
        child = downEdge.child
        mvc = downEdge.mutations_number()
        p = 1
        for downEdge_ in node.down_edges:
            child_ = downEdge_.child
            if child_.id != child.id:
                p *= (t0 - child_.time)
        for upEdge_ in node.up_edges:
            parent_ = upEdge_.parent
            p *= (parent_.time - t0)
        sum1 += p * mvc

    sum2 = 0
    for upEdge in node.up_edges:
        parent = upEdge.parent
        mvp = upEdge.mutations_number()
        p = 1
        for upEdge_ in node.up_edges:
            parent_ = upEdge_.parent
            if parent_.id != parent.id:
                p *= (parent_.time - t0)
        for downEdge_ in node.down_edges:
            child_ = downEdge_.child
            p *= (t0 - child_.time)
        sum2 += p * mvp

    P = 1

    for upEdge in node.up_edges:
        parent = upEdge.parent
        P *= (parent.time - t0)

    for downEdge in node.down_edges:
        child = downEdge.child
        P *= (t0 - child.time)

    result = sum1 - sum2 + (len(node.down_edges) - len(node.up_edges)) * P

    return result


def findRoot(node):
    if not node.down_edges or not node.up_edges:
        return -float("inf")

    highestDownEdge = max(node.down_edges, key=lambda elem: elem.child.time)
    lowestUpEdge = min(node.up_edges, key=lambda elem: elem.parent.time)

    l = highestDownEdge.child.time
    r = lowestUpEdge.parent.time

    f_l = f_(node, l)
    f_r = f_(node, r)
    if f_l == 0:
        return l
    if f_r == 0:
        return r

    assert f_l * f_r < 0
    reverse = (f_(node, l) > 0)
    eps = 1e-5
    stepCounter = 0

    while abs(f_(node, l)) > eps and stepCounter < 1e7:
        m = (r + l) / 2
        if (f_(node, m) > 0) ^ reverse:
            r = m
        else:
            l = m
        stepCounter += 1
    print(stepCounter)
    return l


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
print("Real time function definition", F_real)

F_im = 0.
for edge in edges.values():
    F_im += pois.logpmf(edge.mutations_number(), MUTATION_RATE * edge.im_length() * (edge.right - edge.left))
print("Imaginary time function definition", F_im)

for node in nodes.values():
    t_0 = findRoot(node)
    print("======")
    print("node:", node)
    print("node.time:", node.time)
    if t_0 != -float("inf"):
        print("t_0 =", t_0)
        print("f'(t_0) =", f_(node, t_0))
        assert f_(node, t_0) < 1e-4
    else:
        print("Node without parents or child")
    print("\n")

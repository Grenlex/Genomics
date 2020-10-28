import msprime
import random
import math
from IPython.display import SVG, display
from scipy.stats import poisson as pois
import numpy as np

SAMPLE_SIZE = 20
Ne = 1000
LENGTH = 5e2
RECOMBINATION_RATE = 2e-8
MUTATION_RATE = 6.83e-7
Zero_mutation_num=0.01


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
        if len(self.mutations)!=0:
            return len(self.mutations)
        else:
            return Zero_mutation_num

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
        root.time = max(root.time, child.time + 200)

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

def Delete_Root(nodes, edges):
    Rid=len(nodes)-1
    num=len(edges)
    for edge in range(num):
        if edges[edge].parent.id==Rid:
            edges.pop(edge)
    for node in range(len(nodes)):
        if nodes[node].id==Rid:
            nodes.pop(node)
            break
        
def f_(node, t0):
    sum1 = 0
    for downEdge in node.down_edges:
        child = downEdge.child
        mvc = downEdge.mutations_number()
        dna = downEdge.right - downEdge.left
        p = 1
        for downEdge_ in node.down_edges:
            child_ = downEdge_.child
            if child_.id != child.id:
                p *= (t0 - child_.time)
        for upEdge_ in node.up_edges:
            parent_ = upEdge_.parent
            p *= (parent_.time - t0)
        sum1 += p * mvc * dna

    sum2 = 0
    for upEdge in node.up_edges:
        parent = upEdge.parent
        mvp = upEdge.mutations_number()
        dna = upEdge.right - upEdge.left
        p = 1
        for upEdge_ in node.up_edges:
            parent_ = upEdge_.parent
            if parent_.id != parent.id:
                p *= (parent_.time - t0)
        for downEdge_ in node.down_edges:
            child_ = downEdge_.child
            p *= (t0 - child_.time)
        sum2 += p * mvp * dna

    P = 0

    for upEdge in node.up_edges:
        parent = upEdge.parent
        dna = upEdge.right - upEdge.left
        P += (parent.time - t0) * dna

    for downEdge in node.down_edges:
        child = downEdge.child
        dna = downEdge.right - downEdge.left
        P += (t0 - child.time) * dna

    k = 0
    for downEdge in node.down_edges:
        k -= downEdge.right - downEdge.left

    for upEdge in node.up_edges:
        k += upEdge.right - upEdge.left

    result = sum1 - sum2 - MUTATION_RATE * k * P

    return result

def findRoot(node, debug=False):
    from matplotlib import pyplot as plt
    import numpy as np
    if not node.down_edges or not node.up_edges:
        return None

    highestDownEdge = max(node.down_edges, key=lambda elem: elem.child.time)
    lowestUpEdge = min(node.up_edges, key=lambda elem: elem.parent.time)

    if lowestUpEdge.parent.time ==-1:
        r=2*node.time-highestDownEdge.child.time - 1e-5
    else:
        r = lowestUpEdge.parent.time - 1e-5
    l = highestDownEdge.child.time + 1e-5
    
    
     

    if debug:
        print("l, r", l, r)
    f_l = f_(node, l)
    f_r = f_(node, r)

    if debug:
        print("f_l, f_r", f_l, f_r)
        fx, fy = plot_F_IM(node, l, r)

        x = np.array([])
        y = np.array([])
        for x0 in np.linspace(l, r, endpoint=True, num=100):
            x = np.append(x, x0)
            y = np.append(y, f_(node, x0))
        title = str(node) + " " + " parents: " + str(len(node.up_edges)) + " children: " + \
                str(len(node.down_edges)) + "\n IM F: " + str(F_im)

    if f_l * f_r >= 0:
        return None

    reverse = (f_(node, l) > 0)
    eps = 1e-5
    stepCounter = 0
    m = (l + r) / 2
    while abs(f_(node, m)) > eps and stepCounter < 1e4:
        m = (r + l) / 2
        if (f_(node, m) > 0) ^ reverse:
            r = m
        else:
            l = m
        stepCounter += 1
        if stepCounter > 1e5:
            if debug:
                print("StepCounter overflow")
                print("Calculated value is not 0 but", f_(node, l))
    if debug:
        print("root:", m)
        plt.figure()

        plt.subplot(211)
        plt.axvline(x=node.time, color="black")
        plt.plot(fx, fy, color="blue")

        plt.subplot(212)
        plt.title(title)
        plt.axhline(y=0, color="black")
        plt.axvline(x=node.time, color="black")
        plt.plot(x, y, color="green")
        plt.plot(m, f_(node, m), 'ro')


        plt.show()

    return m


def updateTimes(nodes, edges, debug):
    import random
    keys = list(nodes.values())
    random.shuffle(keys)
    c=0

    for node in keys:
        F_imold = 0.
        for edge in edges.values():
            if edge.im_length()!=0:
                F_imold+=math.log((MUTATION_RATE * edge.im_length()* (edge.right - edge.left))**edge.mutations_number()*math.exp((-1)*MUTATION_RATE * edge.im_length()* (edge.right - edge.left)))
        a=node.time
        t_0 = findRoot(node, debug)
        if not (t_0 is None):
            f_t0 = f_(node, t_0)
            if debug:
                print("t_0 =", t_0)
                print("f'(t_0) =", f_(node, t_0))
            node.time = t_0
        else:
            if debug:
                print(node.id)
                print("Time shouldn't be changed")
        if debug:
            print("\n")
        F_imnew = 0.
        for edge in edges.values():
            if edge.im_length()!=0:
                F_imnew+=math.log((MUTATION_RATE * edge.im_length()* (edge.right - edge.left))**edge.mutations_number()*math.exp((-1)*MUTATION_RATE * edge.im_length()* (edge.right - edge.left)))
        if F_imnew<=F_imold:
            node.time=a
            #print("aa")
        else:
            c+=1
    print(c)

    
def F_IM(edges):
    F_im = 0.
    for edge in edges.values():
        if edge.im_length()!=0:
            F_im+=math.log((MUTATION_RATE * edge.im_length()* (edge.right - edge.left))**edge.mutations_number()*math.exp((-1)*MUTATION_RATE * edge.im_length()* (edge.right - edge.left)))
    return F_im

def plot_F_IM(node, l = None, r = None):
    if l is None:
        l = node.time * 0.8
        r = node.time * 1.2
    vals = 50
    x = np.zeros(vals)
    y = np.zeros(vals)
    beg_time = node.time
    for i, val in enumerate(np.linspace(l, r, 50)):
        x[i] = val
        node.time = val
        y[i] = F_IM(edges)

    node.time = beg_time

    return x, y

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

Delete_Root(nodes, edges)

for edge in edges.values():
    edge.print()

#for i in range(len(nodes)):
  #  print(nodes[i].id)
  #  print(nodes[i].time)
  #  print()
nodes[18].time=nodes[18].real_time

for i in range(len(nodes)):
    print(nodes[i].id)
    print(nodes[i].time)
    print()

#print(Zero_mutation_num)
F_re = 0.
for edge in edges.values():
    if  edge.length()!=0:
        F_re+=math.log((MUTATION_RATE * edge.length()* (edge.right - edge.left))**edge.mutations_number()*math.exp(-1*MUTATION_RATE * edge.length()* (edge.right - edge.left)))
print("F_real=",F_re)
F_im=0.
for edge in edges.values():
    if edge.im_length()!=0:
        F_im+=math.log((MUTATION_RATE * edge.im_length()* (edge.right - edge.left))**edge.mutations_number()*math.exp((-1)*MUTATION_RATE * edge.im_length()* (edge.right - edge.left)))
print("F_im=",F_im)
while (1):
    updateTimes(nodes, edges, debug=True)
    print("F_real=",F_re)
    F_im=0.
    for edge in edges.values():
        if edge.im_length()!=0:
            F_im+=math.log((MUTATION_RATE * edge.im_length()* (edge.right - edge.left))**edge.mutations_number()*math.exp((-1)*MUTATION_RATE * edge.im_length()* (edge.right - edge.left)))
    print("F_im=",F_im)
    print()
    

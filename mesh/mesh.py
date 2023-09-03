from typing import List
from mesh.node import *
from mesh.parameters import Parameters


class Mesh:
    def __init__(self, parameters: Parameters):
        self.parameters = parameters
        self.start = MeshStart(self.parameters)
        self.left = MeshLeft(self.parameters)
        self.right = MeshRight(self.parameters)
        self.center = MeshCenter(self.parameters)
        self.finish = MeshFinish(self.parameters)

    def create_mesh(self, revers_time=False):
        self.start.create()

        self.left.create(self.start)
        self.left.connect_start(self.start)

        self.right.create(self.start, self.left)

        self.right.connect_start(self.start)

        self.center.create(self.start, self.left)
        self.center.connect_start(self.start)
        self.center.connect_left(self.left)
        self.center.connect_right(self.right)

        self.left.connect_center(self.center)
        self.right.connect_center(self.center)

        self.finish.create(self.left, self.center, self.right)

        self.finish.connect_center(self.center)
        self.finish.connect_left(self.left)
        self.finish.connect_right(self.right)
        if revers_time:
            self.start.revers_time()
            self.left.revers_time()
            self.right.revers_time()
            self.center.revers_time()
            self.finish.revers_time()

    def get_list_nodes(self):
        list_nodes = []
        left_nodes = self.left.get_array_node()[0]
        right_nodes = self.right.get_array_node()[0]

        for level in self.start.get_array_node():
            for node in level:
                list_nodes.append(node)
        for i, level in enumerate(self.center.get_array_node()):
            list_ = [left_nodes[i]]+level+[right_nodes[i]]
            for node in self.center.sort_by_solved(list_):
                list_nodes.append(node)

        for level in self.finish.get_array_node():
            for node in level:
                list_nodes.append(node)
        return list_nodes

    def get_XY_finish(self, revers_time=False):
        finish = self.finish.get_array_node()

        def approx(t_u, t_d, f_u, f_d):
            if t_d == t_u:
                return f_u
            elif revers_time:
                return (f_d - f_u) * (self.parameters.t0 - t_u) / (t_d - t_u) + f_u
            else:
                return (f_d-f_u)*(self.parameters.t1-t_u)/(t_d-t_u) + f_u

        x_t1 = [approx(u.t, d.t, u.x, d.x) for u, d in zip(finish[0], finish[1])]
        y_t1 = [approx(u.t, d.t, u.y, d.y) for u, d in zip(finish[0], finish[1])]
        s_ = [n.s for n in finish[0]]

        return s_, x_t1, y_t1


class PartMesh(ABC):
    def __init__(self, parameter):
        self.parameters = parameter
        self.nodes = None
        self.revers = False

    @abstractmethod
    def create(self, *args):
        pass

    @abstractmethod
    def connect(self):
        pass

    @abstractmethod
    def get_array_node(self) -> List[List[Node]]:
        pass

    def revers_time(self):
        self.revers = True
        for level in self.get_array_node():
            for node in level:
                node.t = self.parameters.t1-node.t
                node.s = self.parameters.s1-node.s

    def sort_by_solved(self, list_node):  # сначала те, которые можно решить TODO для обратной
        first = [list_node[0]]
        second = []
        for i, node in enumerate(list_node[1:]):
            if not self.revers:
                if node.t > list_node[i].t:
                    second.append(node)
                else:
                    first.append(node)
            else:
                if node.t < list_node[i].t:
                    second.append(node)
                else:
                    first.append(node)
        return first + second


class MeshStart(PartMesh):
    def get_array_node(self) -> List[List[Node]]:
        return [self.nodes[1], self.nodes[2], self.sort_by_solved(self.nodes[0])]

    def create(self):
        nodes = [[], [], []]
        t_ = self.parameters.t0
        s_ = self.parameters.s0

        for i in range(self.parameters.m):
            nodes[0].append(NodeCenter(s_, t_))
            t_ = (t_ + self.parameters.h[0]) % (self.parameters.h[0]+self.parameters.h[1])
            s_ += self.parameters.ds

        for node in nodes[0]:
            s_right = self.parameters.c1 * (node.t - self.parameters.t0) + node.s
            s_left = self.parameters.c2 * (self.parameters.t0 - node.t) + node.s
            nodes[1].append(NodeStart(s_right, self.parameters.t0))
            nodes[2].append(NodeStart(s_left, self.parameters.t0))

        nodes[0].append(NodeRight(s_, t_))
        nodes[1].append(NodeStart(self.parameters.s1, self.parameters.t0))
        s_left = self.parameters.c2 * (self.parameters.t0 - nodes[0][-1].t) + nodes[0][-1].s
        nodes[2].append(NodeStart(s_left, self.parameters.t0))
        self.nodes = nodes
        self.connect()

    def connect(self):
        self.parameters.inds = self.neighbor_is_left()
        for i, ind in enumerate(self.parameters.inds):
            self.nodes[0][i+1].left = self.nodes[0][i] if ind == 1. else self.nodes[2][i+1]
            self.nodes[0][i].right = self.nodes[0][i+1] if ind == 0. else self.nodes[1][i]
        self.nodes[0][0].left = self.nodes[2][0]
        self.nodes[0][-1].right = self.nodes[1][-1]

    def neighbor_is_left(self):
        ind = np.zeros(self.parameters.m)
        t_hist = 0
        t_ = (t_hist + self.parameters.h[0]) % (self.parameters.h[0] + self.parameters.h[1])
        for i in range(self.parameters.m):
            if t_ > t_hist:
                ind[i] = 1.
            t_hist = t_
            t_ = (t_hist + self.parameters.h[0]) % (self.parameters.h[0] + self.parameters.h[1])
        return ind


class MeshCenter(PartMesh):
    def get_array_node(self) -> List[List[Node]]:
        llist_nodes = []
        for level in self.nodes:
            llist_nodes.append(self.sort_by_solved(level))
        return llist_nodes

    def create(self, start_mesh, left_mesh):
        count_level = len(left_mesh.nodes)
        start_nodes = start_mesh.nodes[0][1:-1]
        self.nodes = [[NodeCenter(node.s, node.t + self.parameters.dt * i) for node in start_nodes] for i in range(1, count_level + 1)]
        self.connect()

    def connect(self):
        count_level = len(self.nodes)
        for level in range(1, count_level):
            for i, ind in enumerate(self.parameters.inds[1:-1]):
                self.nodes[level][i + 1].left = self.nodes[level][i] if ind == 1. else self.nodes[level - 1][i]
                self.nodes[level][i].right = self.nodes[level - 1][i + 1] if ind == 1. else self.nodes[level][i + 1]

        for i, ind in enumerate(self.parameters.inds[1:-1]):
            if ind == 1.:
                self.nodes[0][i + 1].left = self.nodes[0][i]
            else:
                self.nodes[0][i].right = self.nodes[0][i + 1]

    def connect_start(self, start_mesh):
        start_nodes = start_mesh.nodes[0]
        for i, ind in enumerate(self.parameters.inds[1:-1]):
            if ind == 1.:
                self.nodes[0][i].right = start_nodes[i + 2]
            else:
                self.nodes[0][i + 1].left = start_nodes[i+1]
        if self.parameters.inds[-1] == 1.:
            self.nodes[0][-1].right = start_nodes[-1]

    def connect_left(self, left_mesh):
        left_nodes = left_mesh.nodes
        for i in range(len(self.nodes)):
            self.nodes[i][0].left = left_nodes[i]

    def connect_right(self, right_mesh):
        right_nodes = right_mesh.nodes
        if self.parameters.inds[-1] == 1.:
            for i in range(1, len(self.nodes)):
                self.nodes[i][-1].right = right_nodes[i-1]
        else:
            for i in range(len(self.nodes)):
                self.nodes[i][-1].right = right_nodes[i]


class MeshLeft(PartMesh):
    def get_array_node(self) -> List[List[Node]]:
        return [self.nodes]

    def create(self, start_mesh):
        start_node = start_mesh.nodes[0][0]
        dt = self.parameters.dt
        t = np.arange(start_node.t+dt, self.parameters.t1-dt, dt)
        self.nodes = [NodeLeft(self.parameters.s0, ti) for ti in t]
        self.connect()

    def connect(self):
        for i, node in enumerate(self.nodes[1:]):
            node.left = self.nodes[i]

    def connect_start(self, start_mesh):
        start_nodes = start_mesh.nodes[0]
        self.nodes[0].left = start_nodes[0]
        self.nodes[0].right = start_nodes[1]

    def connect_center(self, center_mesh):
        center_nodes = center_mesh.nodes
        for i in range(len(self.nodes)-1):
            self.nodes[i+1].right = center_nodes[i][0]


class MeshRight(PartMesh):
    def get_array_node(self) -> List[List[Node]]:
        return [self.nodes]

    def create(self, start_mesh, left_mesh):
        start_node = start_mesh.nodes[0][-1]
        count_nodes = len(left_mesh.nodes)
        dt = self.parameters.dt
        t = [start_node.t+i*dt for i in range(1, count_nodes+1)]
        self.nodes = [NodeRight(self.parameters.s1, ti) for ti in t]
        self.connect()

    def connect(self):
        for i, node in enumerate(self.nodes[1:]):
            node.right = self.nodes[i]

    def connect_start(self, start_mesh):
        start_nodes = start_mesh.nodes[0]
        self.nodes[0].right = start_nodes[-1]
        if self.parameters.inds[-1] == 0:
            self.nodes[0].left = start_nodes[-2]

    def connect_center(self, center_mesh):
        center_nodes = center_mesh.nodes
        if self.parameters.inds[-1] == 1.:
            for i in range(len(self.nodes)):
                self.nodes[i].left = center_nodes[i][-1]
        else:
            for i in range(1, len(self.nodes)):
                self.nodes[i].left = center_nodes[i-1][-1]


class MeshFinish(PartMesh):
    def get_array_node(self) -> List[List[Node]]:
        return [self.sort_by_solved(self.nodes[0]), self.sort_by_solved(self.nodes[1])]

    def create(self, left_mesh, center_mesh, right_mesh):
        left_node = left_mesh.nodes[-1]
        right_node = right_mesh.nodes[-1]
        center_nodes = center_mesh.nodes[-1]
        finish_nodes: List[List[Node]] = [[], []]
        for i in [0, 1]:
            dt = (1+i)*self.parameters.dt
            finish_nodes[i] = ([NodeLeft(left_node.s, left_node.t+dt)] +
                               [NodeCenter(node.s, node.t + dt) for node in center_nodes] +
                               [NodeRight(right_node.s, right_node.t+dt)])

        self.nodes = finish_nodes
        self.connect()

    def connect(self):
        self.nodes[1][1].left = self.nodes[1][0]
        self.nodes[1][0].right = self.nodes[0][1]
        self.nodes[1][0].left = self.nodes[0][0]
        self.nodes[0][1].left = self.nodes[0][0]
        for i, ind in enumerate(self.parameters.inds[1:-1]):
            if ind == 1.:
                self.nodes[0][i + 2].left = self.nodes[0][i + 1]
                self.nodes[1][i + 2].left = self.nodes[1][i + 1]
                self.nodes[1][i + 1].right = self.nodes[0][i + 2]  # соединение с нижним
            else:
                self.nodes[0][i + 1].right = self.nodes[0][i + 2]
                self.nodes[1][i + 1].right = self.nodes[1][i + 2]
                self.nodes[1][i + 2].left = self.nodes[0][i + 1]  # соединение с нижним
        if self.parameters.inds[-1] == 1.:
            self.nodes[0][-1].left = self.nodes[0][-2]
            self.nodes[1][-1].left = self.nodes[1][-2]
            self.nodes[1][-2].right = self.nodes[0][-1]
        else:
            self.nodes[0][-2].right = self.nodes[0][-1]
            self.nodes[1][-2].right = self.nodes[1][-1]
            self.nodes[1][-1].left = self.nodes[0][-2]
        self.nodes[1][-1].right = self.nodes[0][-1]

    def connect_center(self, center_mesh):
        center_nodes = center_mesh.nodes
        self.nodes[0][0].right = center_nodes[-1][0]
        for i, ind in enumerate(self.parameters.inds[1:-1]):
            if ind == 1.:
                self.nodes[0][i + 1].right = center_nodes[-1][i + 1]  # соединение с нижним
            else:
                self.nodes[0][i + 2].left = center_nodes[-1][i]  # соединение с нижним
        if self.parameters.inds[-1] == 0.:
            self.nodes[0][-1].left = center_nodes[-1][-1]

    def connect_left(self, left_mesh):
        left_nodes = left_mesh.nodes
        self.nodes[0][0].left = left_nodes[-1]

    def connect_right(self, right_mesh):
        right_nodes = right_mesh.nodes
        if self.parameters.inds[-1] == 1.:
            self.nodes[0][-2].right = right_nodes[-1]
        self.nodes[0][-1].right = right_nodes[-1]
from mesh.node import Node
import numpy as np


class Solver:
    def __init__(self, parameters):
        self.C1 = parameters.c1
        self.C2 = parameters.c2
        self.S0 = parameters.s0
        self.S1 = parameters.s1
        self.x0 = lambda s: np.exp(s)
        self.y0 = lambda s: 0

        self.B11 = lambda s, t: - self.C1
        self.B12 = lambda s, t: - np.exp(s) / (s + 2)
        self.B21 = lambda s, t: (s + 2) / np.exp(s)
        self.B22 = lambda s, t: self.C2 / (s + 2)

        self.G11 = lambda t: np.exp(self.S1) * np.sin(t) / ((self.S1 + 2) * np.sin(t) - np.exp(self.S1) * np.cos(t))
        self.G12 = lambda t: - self.G11(t)
        self.G21 = lambda t: 2 * np.cos(t) / (np.cos(t) - 2 * np.sin(t))
        self.G22 = lambda t: - self.G21(t)

    def solver(self, node: Node):
        type_node = node.get_type()
        if type_node == 1:  # начальные точки
            self.solver_start(node)
        elif type_node == 2:  # поверхность
            self.solver_center(node)
        elif type_node == 3:  # левая граница в прямой задачи и правая в обратной
            self.solver_border1(node)
        elif type_node == 4:  # правая граница в прямой задачи и левая в обратной
            self.solver_border2(node)
        elif type_node == 5:
            self.solver_center(node)

    def solver_start(self, node: Node):
        node.x = self.x0(node.s)
        node.y =  self.y0(node.s)
        node.is_resolved = True

    def solver_center(self, node: Node):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = node.get_old_point()

        A = [[1 - hr / 2 * self.B11(node.s, node.t), -hr / 2 * self.B12(node.s, node.t)],
             [-hl / 2 * self.B21(node.s, node.t), 1 - hl / 2 * self.B22(node.s, node.t)]]
        b = [xr + hr / 2 * (self.B11(sr, tr) * xr + self.B12(sr, tr) * yr),
             yl + hl / 2 * (self.B21(sl, tl) * xl + self.B22(sl, tl) * yl)]

        node.x, node.y = np.linalg.solve(A, b)
        node.is_resolved = True

    def solver_border1(self, node: Node):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = node.get_old_point()

        A = [[1 - hr / 2 * self.B11(node.s, node.t), -hr / 2 * self.B12(node.s, node.t)],
             [-hl / 2 * self.G21(node.t), 1 - hl / 2 * self.G22(node.t)]]
        b = [xr + hr / 2 * (self.B11(sr, tr) * xr + self.B12(sr, tr) * yr),
             yl + hl / 2 * (self.G21(tl) * xl + self.G22(tl) * yl)]
        node.x, node.y = np.linalg.solve(A, b)

        node.is_resolved = True

    def solver_border2(self, node: Node):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = node.get_old_point()

        A = [[1 - hr / 2 * self.G11(node.t), -hr / 2 * self.G12(node.t)],
             [-hl / 2 * self.B21(node.s, node.t), 1 - hl / 2 * self.B22(node.s, node.t)]]
        b = [xr + hr / 2 * (self.G11(tr) * xr + self.G12(tr) * yr),
             yl + hl / 2 * (self.B21(sl, tl) * xl + self.B22(sl, tl) * yl)]
        node.x, node.y = np.linalg.solve(A, b)

        node.is_resolved = True

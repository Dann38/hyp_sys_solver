import numpy as np
from mesh.node import Node

x_an = lambda s, t: np.exp(s) * (4-t) * np.cos(4-t)
y_an = lambda s, t: (s + 2) * np.sin(4-t)
class Solver:
    def __init__(self, parameters):

        self.C1 = parameters.c1
        self.C2 = parameters.c2
        self.S0 = parameters.s0
        self.S1 = parameters.s1
        self.T1 = parameters.t1
        self.x0 = lambda s: np.exp(s)
        self.y0 = lambda s: 0

        self.A1 = lambda s, t: self.C1
        self.A2 = lambda s, t: - np.exp(s) * (self.T1 - t) / (s+2)
        self.B1 = lambda s, t: 1 / np.exp(s)
        self.B2 = lambda s, t: - self.C2 / (s + 2)

        self.Fx = lambda s, t: -np.exp(s) * np.cos(self.T1 - t)
        self.Fy = lambda s, t: (self.T1 - t) * np.cos(self.T1 - t) - (s + 2) * np.cos(self.T1 - t)
        self.phi_x = lambda t: -self.C2 * (self.S1 + 2) * np.sin(self.T1 - t) + \
                               np.cos(self.T1 - t) * (
                                       self.C1 * np.exp(self.S1) * (self.T1 - t) - self.C2 * (self.S1 + 2)
                               )
        self.phi_y = lambda t: - self.C1 * np.exp(self.S0) * np.cos(self.T1 - t) * (1 + self.T1 - t) + \
                               np.sin(self.T1 - t) * (
                                       self.C2 * (self.S0 + 2) + self.C1 * np.exp(self.S0) * (self.T1 - t)
                               )
        self.psi1_t1 = lambda s: 0
        self.psi2_t1 = lambda s: 0

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
        node.x = self.psi1_t1(node.s)
        node.y = self.psi2_t1(node.s)
        node.is_resolved = True

    def solver_center(self, node: Node):  # Сделано
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = node.get_old_point()

        A = [[1 + hr / 2 * self.A1(node.s, node.t), hr / 2 * self.A2(node.s, node.t)],
             [hl / 2 * self.B1(node.s, node.t), 1 + hl / 2 * self.B2(node.s, node.t)]]
        b = [xr - hr / 2 * (self.A1(sr, tr) * xr + self.A2(sr, tr) * yr - self.Fx(sr, tr) - self.Fx(node.s, node.t)),
             yl - hl / 2 * (self.B1(sl, tl) * xl + self.B2(sl, tl) * yl - self.Fy(sl, tl) - self.Fy(node.s, node.t))]

        node.x, node.y = np.linalg.solve(A, b)
        # node.x = x_an(node.s, node.t)
        # node.y = y_an(node.s, node.t)
        node.is_resolved = True

    def solver_border2(self, node: Node): # Точно border 2 для S0
        (sr, tr, xr, yr, hr), (sl, tl, xl, yl, hl) = node.get_old_point()

        A = [[1-hr/2, hr/2*self.C2/self.C1],
             [hl/2*self.B1(node.s, node.t), 1+hl/2*self.B2(node.s, node.t)]]

        b = [xr + hr / (2*self.C1) * (self.C1*xr-self.C2*yr + self.phi_y(node.t) + self.phi_y(tr)),
             yl - hl / 2 * (self.B1(sl, tl)*xl + self.B2(sl, tl)*yl - self.Fy(node.s, node.t) - self.Fy(sl, tl))]
        node.x, node.y = np.linalg.solve(A, b)
        # node.x = x_an(node.s, node.t)
        # node.y = y_an(node.s, node.t)
        node.is_resolved = True

    def solver_border1(self, node: Node):
        (sr, tr, xr, yr, hr), (sl, tl, xl, yl, hl) = node.get_old_point()

        A = [[1 - hr / 2, hr / 2 * self.C1 / self.C2],
             [+ hl / 2 * self.B1(node.s, node.t), 1 + hl / 2 * self.B2(node.s, node.t)]]
        b = [yr + hr / 2 * (yr - self.C1 / self.C2 * xr + self.phi_x(node.t)+self.phi_x(tr)),
             yl - hl / 2 * (self.B1(sl, tl) * xl + self.B2(sl, tl) * yl - self.Fy(sl, tl) - self.Fy(node.s, node.t))]

        node.x, node.y = np.linalg.solve(A, b)
        node.x = x_an(node.s, node.t)
        node.y = y_an(node.s, node.t)
        node.is_resolved = True

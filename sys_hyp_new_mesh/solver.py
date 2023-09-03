from node import Node
import numpy as np

C1 = 1
C2 = 3
T0 = 0
T1 = 4
S0 = 0
S1 = 2


x0 = lambda s: np.exp(s)
y0 = lambda s: 0

B11 = lambda s, t: - C1
B12 = lambda s, t: - np.exp(s)/(s+2)
B21 = lambda s, t: (s+2)/np.exp(s)
B22 = lambda s, t: C2/(s+2)

G11 = lambda t: np.exp(S1)*np.sin(t)/((S1+2)*np.sin(t)-np.exp(S1)*np.cos(t))
G12 = lambda t: - G11(t)
G21 = lambda t: 2*np.cos(t)/(np.cos(t)-2*np.sin(t))
G22 = lambda t: - G21(t)


class Solver:
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
        node.x = x0(node.s)
        node.y = y0(node.s)
        node.is_resolved = True

    def solver_center(self, node: Node):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = node.get_old_point()
        try:
            A = [[1 - hr / 2 * B11(node.s, node.t), -hr / 2 * B12(node.s, node.t)],
                 [-hl / 2 * B21(node.s, node.t), 1 - hl / 2 * B22(node.s, node.t)]]
            b = [xr + hr / 2 * (B11(sr, tr) * xr + B12(sr, tr) * yr),
                 yl + hl / 2 * (B21(sl, tl) * xl + B22(sl, tl) * yl)]
        except:
            print(node.s, node.t)
            print((sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr))
            exit()
        node.x, node.y = np.linalg.solve(A, b)
        node.is_resolved = True

    def solver_border1(self, node: Node):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = node.get_old_point()

        A = [[1 - hr / 2 * B11(node.s, node.t), -hr / 2 * B12(node.s, node.t)],
             [-hl / 2 * G21(node.t), 1 - hl / 2 * G22(node.t)]]
        b = [xr + hr / 2 * (B11(sr, tr) * xr + B12(sr, tr) * yr),
             yl + hl / 2 * (G21(tl) * xl + G22(tl) * yl)]
        node.x, node.y = np.linalg.solve(A, b)

        node.is_resolved = True

    def solver_border2(self, node: Node):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = node.get_old_point()

        A = [[1 - hr / 2 * G11(node.t), -hr / 2 * G12(node.t)],
             [-hl / 2 * B21(node.s, node.t), 1 - hl / 2 * B22(node.s, node.t)]]
        b = [xr + hr / 2 * (G11(tr) * xr + G12(tr) * yr),
             yl + hl / 2 * (B21(sl, tl) * xl + B22(sl, tl) * yl)]
        node.x, node.y = np.linalg.solve(A, b)

        node.is_resolved = True

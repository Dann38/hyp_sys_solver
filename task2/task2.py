from mesh.mesh import Mesh
from mesh.parameters import Parameters
from solver import Solver
import numpy as np

from mesh.drawer import Drawer

x_an = lambda s, t: np.exp(s) * (4-t) * np.cos(4-t)
y_an = lambda s, t: (s + 2) * np.sin(4-t)

if __name__ == '__main__':
    parameters = Parameters(c1=1, c2=3, s0=0, s1=2, t0=0, t1=4, count_node=20)
    mesh = Mesh(parameters)
    solver = Solver(parameters)
    drawer = Drawer()

    mesh.create_mesh(revers_time=True)
    nodes = mesh.get_list_nodes()
    for node in nodes:
        solver.solver(node)
        x = x_an(node.s, node.t)
        y = y_an(node.s, node.t)
        print(node.get_type(), f"({node.s:3.2f},{node.t:3.2f})\tреш ({node.x:3.4f}, {node.y:3.4f})\t aн ({x:3.4f}, {y:3.4f})")

    drawer.plot_mesh(nodes)
    s, x, y = mesh.get_XY_finish(revers_time=True)
    drawer.plot_final(s, lambda si: x_an(si, mesh.parameters.t1), lambda si: y_an(si, mesh.parameters.t1), x, y)

    drawer.plot_border(mesh.left.get_array_node()[0], mesh.right.get_array_node()[0],
                       lambda t: x_an(parameters.s1, t), lambda t: y_an(parameters.s1, t),
                       lambda t: x_an(parameters.s0, t), lambda t: y_an(parameters.s0, t))
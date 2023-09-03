from mesh.mesh import Mesh
from mesh.parameters import Parameters
from task1.solver import Solver
import numpy as np

from mesh.drawer import Drawer

x_an = lambda s, t: np.exp(s) * (4-t) * np.cos(4-t)
y_an = lambda s, t: (s + 2) * np.sin(4-t)

if __name__ == '__main__':
    parameters = Parameters(c1=1, c2=3, s0=0, s1=2, t0=0, t1=4, count_node=30)
    mesh = Mesh(parameters)
    solver = Solver(parameters)
    drawer = Drawer()

    mesh.create_mesh(revers_time=True)
    nodes = mesh.get_list_nodes()
    for node in nodes:
        solver.solver(node)

    s, x, y = mesh.get_XY_t1()
    drawer.plot_final(s, lambda si: x_an(si, mesh.parameters.t1), lambda si: y_an(si, mesh.parameters.t1), x, y)
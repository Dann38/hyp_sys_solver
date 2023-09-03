from mesh.mesh import Mesh
from mesh.parameters import Parameters
from solver import Solver
import numpy as np

from mesh.drawer import Drawer

x_an = lambda s, t: np.exp(s) * np.cos(t)
y_an = lambda s, t: (s + 2) * np.sin(t)

if __name__ == '__main__':
    parameters = Parameters(c1=1, c2=3, s0=0, s1=2, t0=0, t1=40, count_node=48)
    mesh = Mesh(parameters)
    solver = Solver(parameters)
    drawer = Drawer()

    mesh.create_mesh(revers_time=False)
    nodes = mesh.get_list_nodes()
    for node in nodes:
        solver.solver(node)

    s, x, y = mesh.get_XY_finish()
    drawer.plot_final(s, lambda si: x_an(si, mesh.parameters.t1), lambda si: y_an(si, mesh.parameters.t1), x, y)
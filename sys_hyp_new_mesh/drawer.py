from typing import List

import matplotlib.pyplot as plt
from node import Node
import numpy as np


class Drawer:
    def plot_mesh(self, list_nodes: List[Node]):

        def plot_arrow_node(node: Node):
            alpha = 0.8
            width_arrow = 0.01
            if node.left is not None:
                s = node.s
                t = node.t
                ds = (node.left.s - node.s) * alpha
                dt = (node.left.t - node.t) * alpha
                plt.arrow(s, t, ds, dt, width=width_arrow, color="b")
            if node.right is not None:
                s = node.s
                t = node.t
                ds = (node.right.s - node.s) * alpha
                dt = (node.right.t - node.t) * alpha
                plt.arrow(s, t, ds, dt, width=width_arrow * 0.7, color="orange")

        def plot_node(node: Node):
            color = ("r" if node.right is None else "b") if node.left is None else (
                "orange" if node.right is None else "g")
            marker = "o" if color == "g" else "o"
            plot_arrow_node(node)
            plt.scatter(node.s, node.t, color=color, marker=marker)

        for node in list_nodes:
            plot_node(node)

        plt.show()

    def plot_border(self, right_nodes, left_nodes, x_an_s1, y_an_s1, x_an_s0, y_an_s0):
        fig, axs = plt.subplots(nrows=2, ncols=2)
        x = [node.x for node in right_nodes]
        y = [node.y for node in right_nodes]
        t = [node.t for node in right_nodes]
        axs[0][0].plot(t, x, "*")
        axs[0][0].plot(np.array(t), x_an_s1(np.array(t)))
        axs[0][0].set_title("Правая граница $x$")
        axs[1][0].plot(t, y, "*")
        axs[1][0].plot(np.array(t), y_an_s1(np.array(t)))
        axs[1][0].set_title("Правая граница $y$")

        x = [node.x for node in left_nodes]
        y = [node.y for node in left_nodes]
        t = [node.t for node in left_nodes]
        axs[0][1].plot(t, x, "*")
        axs[0][1].plot(np.array(t), x_an_s0(np.array(t)))
        axs[0][1].set_title("Левая граница $x$")
        axs[1][1].plot(t, y, "*")
        axs[1][1].plot(np.array(t), y_an_s0(np.array(t)))
        axs[1][1].set_title("Левая граница $y$")
        plt.show()

    def plot_final(self, s, x_an, y_an, x, y):
        s_ = np.sort(s)
        x_an_ = [x_an(si) for si in s_]
        y_an_ = [y_an(si) for si in s_]

        count_node = len(s)
        fig, axs = plt.subplots(nrows=2, ncols=1)
        axs[0].set_title(f"Численное решение (узлов по s: {count_node})")
        axs[0].plot(s_, x_an_, )
        axs[0].plot(s, x, "o")
        axs[0].grid()
        axs[0].legend(["аналитическое решение", "численное решение"])
        axs[0].set_ylabel("$x(s, t_1)$")

        axs[1].plot(s_, y_an_)
        axs[1].plot(s, y, "o")
        axs[1].grid()
        axs[1].legend(["аналитическое решение", "численное решение"])
        axs[1].set_ylabel("$y(s, t_1)$")
        axs[1].set_xlabel("$s$")

        plt.show()

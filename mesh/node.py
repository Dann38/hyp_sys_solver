import numpy as np
from abc import ABC, abstractmethod


class Node(ABC):
    def __init__(self, s, t, left=None, right=None):
        self.s = s
        self.t = t
        self.x = None
        self.y = None
        self.left = left
        self.right = right
        self.is_resolved = False

    def __str__(self):

        # return f"(s: {self.s:.2f}, t: {self.t:.2f})"
        return f"(s: {self.s:.2f}, t: {self.t:.2f}) => x: {self.x: .4f}, y: {self.y: .4f} " \
               # f"\t ан.реш x:{x_an(self.s, self.t): .4f}, y:{y_an(self.s, self.t):.4f}"
               # f"({self.left.s:.2f}, {self.left.t:.2f}), ({self.right.s:.2f}, {self.right.t:.2f})"\
               # f"{type(self)}"
               # f"Node left: ({self.left.s: .2f}, {self.left.t: .2f})]"
               # f" [ Node right: ({self.right.s: .2f},{self.right.s: .2f})" \

    @abstractmethod
    def get_type(self):
        pass

    def get_inf(self):
        return [self.s, self.t, self.x, self.y]

    def get_old_point(self):
        sl, tl, xl, yl = self.left.get_inf()
        hl = abs(self.t - tl)

        sr, tr, xr, yr = self.right.get_inf()
        hr = abs(self.t - tr)
        return (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr)

class NodeStart(Node):
    def get_type(self):
        return 1

class NodeLeft(Node):
    def get_type(self):
        return 3

class NodeRight(Node):
    def get_type(self):
        return 4

class NodeFinish(Node):
    def get_type(self):
        return 5

class NodeCenter(Node):
    def get_type(self):
        return 2
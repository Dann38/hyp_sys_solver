class Parameters:
    def __init__(self, c1, c2, s0, s1, t0, t1, count_node):
        self.c1 = c1
        self.c2 = c2
        self.c = [self.c1, self.c2]
        self.t0 = t0
        self.t1 = t1
        self.t = [self.t0, self.t1]
        self.s0 = s0
        self.s1 = s1
        self.s = [self.s0, self.s1]
        self.m = count_node
        self.ds, self.h_down, self.h_up = self.create_mesh_s()
        self.h = [self.h_down,  self.h_up]
        self.dt = self.h[0]+self.h[1]
        self.inds = None

    def create_mesh_s(self):
        ds = (self.s1 - self.s0) / self.m
        h_down = ds / self.c2
        h_up = ds / self.c1
        return ds, h_down, h_up
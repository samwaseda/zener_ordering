import numpy as np
from pyiron_base.generic.datacontainer import DataContainer
from pyiron_atomistics.atomistics.job.interactivewrapper import InteractiveWrapper
from scipy.spatial import cKDTree


class Metadynamics(InteractiveWrapper):
    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        self.input = DataContainer(table_name='input')
        self.output = DataContainer(table_name='output')
        self.input.update_every_n_steps = 100
        self.input.sigma = 0.01
        self.input.increment = 0.0001
        self.input.n_mesh = 1000
        self.input.n_repeat = None
        self.input.use_gradient = True
        self.output.z_lst = []
        self.gradient = Gradient()
        self._mesh = None
        self._mesh_tree = None
        self._ind_c = None
        self._ind_fe = None

    @property
    def mesh(self):
        if self._mesh is None:
            self._mesh, sp = np.linspace(0, 1, self.input.n_mesh, retstep=True, endpoint=False)
            self._mesh += 0.5 * sp
        return self.mesh

    def get_mesh_index(self, z):
        if self._mesh_tree is None:
            self._mesh_tree = cKDTree(self.mesh)
        return self._mesh_tree.query(z)[1]

    def validate_ready_to_run(self):
        super().validate_ready_to_run()
        if self.input.n_repeat is None:
            raise ValueError('repeat not set')

    def run_static(self):
        self.output.B = np.zeros(self.input.n_mesh)
        self.output.dBds = np.zeros(self.input.n_mesh)
        self.output.ddBdds = np.zeros(self.input.n_mesh)
        self.status.running = True
        self.ref_job_initialize()
        self.ref_job.set_fix_external(self.callback, overload_internal_fix_external=True)
        self.ref_job.run()
        self.status.collect = True
        self.ref_job.interactive_close()
        self.status.finished = True
        self.to_hdf()

    @property
    def cell(self):
        return self.ref_job.interactive_cells_getter() / self.n_repeat

    @property
    def ind_fe(self):
        if self._ind_fe is None:
            self._ind_fe = self.job.structure.select_index('Fe')
        return self._ind_fe

    @property
    def ind_c(self):
        if self._ind_c is None:
            self._ind_c = self.job.structure.select_index('C')
        return self._ind_c

    def callback(self, caller, ntimestep, nlocal, tag, x, fext):
        tags = tag.flatten().argsort()
        index_c = tags[self.ind_c]
        fext.fill(0)
        ext_forces = self.get_forces(x[index_c].copy())
        fext[index_c] += ext_forces
        fext[tags[self.ind_fe]] -= np.sum(ext_forces, axis=0) / len(self.ind_fe)
        if ((ntimestep + 1) % self.input.update_every_n_steps) == 0:
            self.set_s()

    @property
    def order_parameter(self):
        return self.gradient.ord_z

    def get_forces(self, x):
        gradient = self.get_gradient(x)
        index = self.get_mesh_index(self.order_parameter)
        n_dash = gradient * self.output.dBds[index]
        if self.input.use_gradient:
            ds = self.order_parameter - self.mesh[index]
            n_dash += gradient * ds * self.output.ddBdds[index]
        return -n_dash

    def get_gradient(self, x):
        return self.gradient.get_gradient(x, self.cell)

    def set_s(self):
        ds = (self.mesh - self.order_parameter) / self.sigma
        B = self.increment * np.exp(-0.5 * ds**2)
        self.output.z_lst.append(self.gradient.ord_z)
        self.output.B += B
        self.output.dBds -= B * ds / self.sigma
        self.output.ddBdds += B * (ds**2 - 1) / self.sigma**2

    def to_hdf(self, hdf=None, group_name=None):
        super().to_hdf(
            hdf=hdf,
            group_name=group_name
        )
        self.output.to_hdf(hdf=self.project_hdf5, group_name='output')

    def from_hdf(self, hdf=None, group_name=None):
        super().from_hdf(
            hdf=hdf,
            group_name=group_name
        )
        self.output.from_hdf(hdf=self.project_hdf5, group_name='output')

    def write_input(self):
        pass


class Gradient:
    def __init__(self):
        self.x = None
        self._E_mat = None
        self.cell = None
        self._x_to_xi_mat = None
        self._x_rel = None
        self._cos_sq = None
        self._inv_cos_sq_sum = None
        self._ord_z = None
        self._n_vec = None
        self._sine = None

    @property
    def E_mat(self):
        if self.E_mat is None:
            E = np.ones((3, 3)) - np.eye(3)
            self._E_mat = np.append(
                E, np.roll(np.eye(3), 1, axis=0) - np.roll(np.eye(3), -1, axis=0), axis=0
            ).reshape(2, 3, 3)
        return self._E_mat

    @property
    def x_to_xi_mat(self):
        if self._x_to_xi_mat is None:
            self._x_to_xi_mat = np.pi * np.linalg.inv(self.cell)
        return self._x_to_xi_mat

    @property
    def x_rel(self):
        if self._x_rel is None:
            self._x_rel = np.einsum('ij,ni,pkj->npk', self.x_to_xi_mat, self.x, self.E_mat)
        return self._x_rel

    @property
    def cos_sq(self):
        if self._cos_sq is None:
            self._cos_sq = np.cos(self.x_rel)**2
        return self._cos_sq

    @property
    def inv_cos_sq_sum(self):
        if self._inv_cos_sq_sum is None:
            self._inv_cos_sq_sum = 1 / np.einsum('ijk->i', self.cos_sq)
        return self._inv_cos_sq_sum

    def reset(self):
        self._x_to_xi_mat = None
        self._x_rel = None
        self._cos_sq = None
        self._inv_cos_sq_sum = None
        self._ord_z = None
        self._n_vec = None
        self._sine = None

    @property
    def n_vec(self):
        if self._n_vec is None:
            self._n_vec = np.einsum('imj,i->j', self.cos_sq, self.inv_cos_sq_sum)
        return self._n_vec

    @property
    def sine(self):
        if self._sine is None:
            self._sine = np.sin(2 * self.x_rel)
        return self._sine

    @property
    def ord_z(self):
        if self._ord_z is None:
            self._ord_z = np.sqrt(1.5 * np.sum(self.n_vec**2) / np.sum(self.n_vec)**2 - 0.5)
        return self._ord_z

    @property
    def safe_z(self):
        return 1.5 / np.sum(self.n_vec)**2 / (self.ord_z + np.isclose(self.ord_z, 0))

    @property
    def gradient(self):
        n_dash = -np.einsum('pil,api,a->ail', self.E_mat, self.sine, self.inv_cos_sq_sum)
        n_dash += np.einsum(
            'aqi,pjl,apj,a->ail', self.cos_sq, self.E_mat, self.sine, self.inv_cos_sq_sum**2
        )
        return np.einsum('jk,aij,i,->ak', self.x_to_xi_mat, n_dash, self.n_vec, self.safe_z)

    def get_gradient(self, x, cell):
        self.reset()
        self.x = x
        self.cell = cell
        return self.gradient
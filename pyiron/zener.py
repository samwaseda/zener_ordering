import numpy as np
from pyiron_base.storage.datacontainer import DataContainer
from pyiron_atomistics.atomistics.job.interactivewrapper import InteractiveWrapper
from pandas import DataFrame
from pyiron_atomistics import Project as PyironProject
import random
from hashlib import sha1


def get_potential():
    potential = {}
    potential['Config'] = [['pair_style eam/alloy\n', 'pair_coeff * * Fe-C-Bec07.eam Fe C\n']]
    potential['Filename'] = [[
        '/home/jovyan/dev/local_projects/zener_ordering/pyiron/Fe-C-Bec07.eam'
    ]]
    potential['Model'] = ['EAM']
    potential['Name'] = ['Raulot']
    potential['Species'] = [['Fe', 'C']]
    potential = DataFrame(potential)
    return potential


def get_lattice_coeff(temperature):
    coeff = np.array([-3.16238463e-12,  2.06828703e-08,  1.00738345e-05,  2.85520935e+00])
    return np.polyval(coeff, temperature)


def get_bulk(pr, temperature=0, cubic=True):
    return pr.create.structure.bulk('Fe', cubic=cubic, a=get_lattice_coeff(temperature))


def get_job(
    pr,
    structure,
    job_name=None,
    temperature=None,
    pressure=None,
    drag_fix_id=None,
    run=True,
    non_modal=False,
    n_ionic_steps=1e5,
):
    job_name = 'lmp_' + sha1(
        (
            structure.__repr__()
            + str(temperature)
            + str(pressure)
            + str(drag_fix_id)
            + str(job_name)
            + str(n_ionic_steps)
        ).encode()
    ).hexdigest()
    lmp = pr.create.job.Lammps(job_name)
    if lmp.status.finished:
        return lmp
    lmp.potential = get_potential()
    lmp.structure = structure
    if temperature is not None:
        if temperature == 0:
            lmp.calc_minimize(
                pressure=pressure, max_iter=n_ionic_steps, n_print=int(n_ionic_steps / 10)
            )
        else:
            lmp.calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps)
    if drag_fix_id is not None:
        setup_lmp_input(lmp, drag_fix_id)
    if non_modal:
        lmp.server.run_mode.non_modal = True
    if lmp.status.initialized and run:
        lmp.run()
    return lmp


def get_elastic_tensor(pr):
    lmp = get_job(
        pr=pr, structure=get_bulk(pr=pr), job_name="elast", run=False
    )
    lmp.interactive_open()
    elast = lmp.create_job("ElasticTensor", "elast")
    if elast.status.initialized:
        elast.run()
    return elast["output/elastic_tensor"]



def get_value(x, min_dist, structure):
    dist_array = structure.get_distances_array(x, x)
    return np.sum(1 / dist_array[(dist_array > 0) * (dist_array < min_dist)])


class Project(PyironProject):
    @staticmethod
    def get_potential():
        return get_potential()
    
    def get_job(
        self,
        structure,
        job_name=None,
        temperature=None,
        pressure=None,
        drag_fix_id=None,
        run=True,
        n_ionic_steps=1e5,
        non_modal=False
    ):
        return get_job(
            pr=self,
            structure=structure,
            job_name=job_name,
            temperature=temperature,
            n_ionic_steps=n_ionic_steps,
            pressure=pressure,
            drag_fix_id=drag_fix_id,
            run=run,
            non_modal=non_modal
        )

    @staticmethod
    def get_lattice_coeff(temperature):
        return get_lattice_coeff(temperature)

    def get_bulk(self, temperature=0, cubic=True):
        return get_bulk(pr=self, temperature=temperature, cubic=cubic)

    def add_carbon(self, positions, repeat, temperature=0, cubic=True, relative=True):
        bulk = self.get_bulk(temperature=temperature, cubic=cubic)
        if relative:
            x = np.einsum("nj,ji->ni", np.array(positions).reshape(-1, 3), bulk.cell)
        else:
            x = positions
        structure = bulk.repeat(repeat)
        return structure + self.create.structure.atoms(
            elements=len(x) * ['C'], positions=x, cell=structure.cell, pbc=structure.pbc
        )

    def append_carbon(self, structure, c, ordered=True, minimum_dist=1.01, max_iter=1000):
        a_0 = (structure.get_volume(per_atom=True) * 2)**(1 / 3)
        x = structure.positions[np.newaxis, :, :] + 0.5 * np.eye(3)[:, np.newaxis, :] * a_0
        n_C = np.rint(len(structure.positions) * c / (1 - c)).astype(int)
        if ordered:
            x = x.reshape(-1, 3)
            indices = np.isclose(
                np.cos(
                    np.pi * np.einsum('ij,nj->ni', np.ones((3, 3)) - np.eye(3), (x - x[0]) / a_0)
                ).sum(axis=-1),
                3
            )
            x = np.random.permutation(x[indices])[:n_C]
        else:
            min_dist = minimum_dist * a_0
            atoms = np.empty(x.shape[:-1], dtype=bool)
            atoms.fill(False)
            atoms[np.arange(3), :np.round(n_C / 3).astype(int)] = True
            current_value = np.inf
            for iii in range(max_iter):
                i_x, i_n = random.choice(np.stack(np.where(atoms), axis=-1))
                j_n = random.choice(np.where(~atoms[i_x])[0])
                atoms[i_x, i_n] = ~atoms[i_x, i_n]
                atoms[i_x, j_n] = ~atoms[i_x, j_n]
                new_value = get_value(x[atoms], min_dist, structure)
                if current_value < new_value:
                    atoms[i_x, i_n] = ~atoms[i_x, i_n]
                    atoms[i_x, j_n] = ~atoms[i_x, j_n]
                else:
                    current_value = new_value
                if current_value == 0:
                    break
            x = x[atoms]
        return structure + self.create.structure.atoms(
            elements=len(x) * ['C'], positions=x, cell=structure.cell
        )


class Metadynamics(InteractiveWrapper):
    """
    Metadynamics class

    Main inputs:

    - update_every_n_steps (int): After every how many steps to update the
        metadynamics histogram
    - sigma (float): Gaussian smearing width of the metadynamics histogram
    - increment (float): Histogram value change at each update
    - n_mesh (int): Number of bins
    - n_repeat (int): Number of repeats in each direction (optional)
    - use_gradient (bool): Whether to use the gradient values
    - z_lst (list): Initial values of the descriptor values (optional)

    There is in principle no value that has to be changed. Depending on the
    setup, it might be important to change the values of `sigma`, `increment` and `n_mesh`
    """
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
        self.input.z_lst = []
        self.output.z_lst = []
        self.gradient = Gradient()
        self._mesh = None
        self._ind_c = None
        self._ind_fe = None
        self._n_repeat = None

    @property
    def n_repeat(self):
        if self._n_repeat is None:
            if self.input.n_repeat is None:
                self._n_repeat = np.rint((0.5 * len(self.ind_fe))**(1 / 3)).astype(int)
            else:
                self._n_repeat = self.input.n_repeat
        return self._n_repeat

    @property
    def mesh(self):
        if self._mesh is None:
            self._mesh, sp = np.linspace(0, 1, self.input.n_mesh, retstep=True, endpoint=False)
            self._mesh += 0.5 * sp
        return self._mesh

    def get_mesh_index(self, z):
        index = np.maximum(np.rint(z * self.input.n_mesh).astype(int), 0)
        return np.minimum(index, self.input.n_mesh - 1)

    def _initialize_meta(self):
        self.output.B = np.zeros(self.input.n_mesh)
        self.output.dBds = np.zeros(self.input.n_mesh)
        self.output.ddBdds = np.zeros(self.input.n_mesh)
        if len(self.input.z_lst) > 0:
            for z in self.input.z_lst:
                self.set_s(z)

    def run_static(self):
        self._initialize_meta()
        self.status.running = True
        self.ref_job_initialize()
        if 'fixedpoint' not in self.ref_job.input.control['fix___ensemble']:
            self.ref_job.input.control['fix___ensemble'] += ' fixedpoint 0 0 0'
        self.ref_job.set_fix_external(self.callback, overload_internal_fix_external=True)
        self.ref_job.input.control['fix_modify___fix_external'] = ' virial no'
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
            self._ind_fe = self.structure.select_index('Fe')
        return self._ind_fe

    @property
    def ind_c(self):
        if self._ind_c is None:
            self._ind_c = self.structure.select_index('C')
        return self._ind_c

    def callback(self, caller, ntimestep, nlocal, tag, x, fext):
        tags = tag.flatten().argsort()
        index_c = tags[self.ind_c]
        fext.fill(0)
        ext_forces = self.get_forces(x[index_c].copy())
        fext[index_c] += ext_forces
        fext[tags[self.ind_fe]] -= np.sum(ext_forces, axis=0) / len(self.ind_fe)
        if ((ntimestep + 1) % self.input.update_every_n_steps) == 0:
            self.set_s(self.order_parameter)

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

    def set_s(self, s):
        self.output.z_lst.append(s)
        ds = (self.mesh - s) / self.input.sigma
        B = self.input.increment * np.exp(-0.5 * ds**2)
        self.output.B += B
        self.output.dBds -= B * ds / self.input.sigma
        self.output.ddBdds += B * (ds**2 - 1) / self.input.sigma**2

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
        if self._E_mat is None:
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
            self._x_rel = np.einsum(
                'ij,ni,pkj->npk', self.x_to_xi_mat, self.x, self.E_mat, optimize=True
            )
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
            z = 1.5 * np.sum(self.n_vec**2) / np.sum(self.n_vec)**2 - 0.5
            self._ord_z = np.sqrt(np.maximum(z, 0))
        return self._ord_z

    @property
    def safe_z(self):
        return 1.5 / np.sum(self.n_vec)**2 / (self.ord_z + np.isclose(self.ord_z, 0))

    @property
    def gradient(self):
        n_dash = -np.einsum(
            'pil,api,a->ail', self.E_mat, self.sine, self.inv_cos_sq_sum, optimize=True
        )
        n_dash += np.einsum(
            'aqi,pjl,apj,a->ail', self.cos_sq, self.E_mat, self.sine, self.inv_cos_sq_sum**2,
            optimize=True
        )
        return np.einsum(
            'jk,aij,i,->ak', self.x_to_xi_mat, n_dash, self.n_vec, self.safe_z, optimize=True
        )

    def get_gradient(self, x, cell):
        self.reset()
        self.x = x
        self.cell = cell
        return self.gradient

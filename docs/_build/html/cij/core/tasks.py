'''This module analysis the basic units (task) of phonon contribution
calculation to be performed, and the order they are conducted.
'''

import numpy
import itertools
from networkx import nx
from typing import List, Union, NamedTuple, Tuple, Iterable
from collections import UserList, UserDict

from cij.util import c_, C_, ElasticModulusCalculationType
from .phonon_contribution import (
    ShearElasticModulusPhononContribution,
    LongitudinalElasticModulusPhononContribution,
    OffDiagonalElasticModulusPhononContribution
)

import logging
logger = logging.getLogger(__name__)

class PhononContributionTaskParams(NamedTuple):
    '''Parameters for elastic constants phonon contribution calculation task.
    '''
    calc_type: ElasticModulusCalculationType
    params: tuple

    @staticmethod
    def _make_param_by_strain_key(strain: tuple, key: C_):
        # print("!!!!", strain)
        if key.is_shear:
            return strain, key
        else: # i = j, k = l
            i, j, k, l = key.s
            # TODO: sorted
            return (
                strain[:, i-1] / numpy.sum(strain, axis=1),
                strain[:, k-1] / numpy.sum(strain, axis=1)
            )

    @classmethod
    def create(cls, strain: tuple, key: C_) -> 'PhononContributionTaskParams':
        '''Factory method for creating ``PhononContributionTaskParams``.

        :param strain: The strains :math:`e_1/\\Delta`, :math:`e_2/\\Delta`,
            :math:`e_3/\\Delta` in the orthogonal coordinates system under
            hydrostatic pressure.
        :param key: The symbol :math:`c_{ij}` needs to be calculated.
        '''
        # print("!!!! ->", strain)
        return cls(key.calc_type, cls._make_param_by_strain_key(strain, key))
    
    def __eq__(self, other: 'PhononContributionTaskParams') -> bool:
        if self.calc_type != other.calc_type: return False
        if self.calc_type == ElasticModulusCalculationType.SHEAR:
            if self.params[1] != other.params[1]: return False # key comparison
            if not numpy.allclose(self.params[0], other.params[0]): return False # strain comparison
            return True
        else:
            if not numpy.allclose(self.params, other.params): return False
            return True
        
    def __hash__(self):
        if self.calc_type != ElasticModulusCalculationType.SHEAR:
            return hash(self.calc_type) ^ hash(tuple(self.params[0].flatten().tolist())) ^ hash(tuple(self.params[1].flatten().tolist()))
        else:
            return hash(self.calc_type) ^ hash(tuple(self.params[0].flatten().tolist())) ^ hash(self.params[1])


class PhononContributionTaskResults(UserDict):
    '''The results of phonon calculation tasks

        - Key is ``PhononContributionTaskParams``
        - Value is ``numpy.ndarray``
    '''

    def get_results_by_strain_keys(self, strain: tuple, keys: Iterable[C_]) -> dict:
        '''The results for given :math:`c_{ij}` symbol at given set of strain

        :param strain: The strains :math:`e_1/\\Delta`, :math:`e_2/\\Delta`,
            :math:`e_3/\\Delta` in the orthogonal coordinates system under
            hydrostatic pressure.
        :param keys: The symbols :math:`c_{ij}`s are calculated.
        '''
        results = dict()
        for key in keys:
            params = PhononContributionTaskParams.create(strain, key)
            results[key] = self[params]
        return results

    def __setitem__(self, _key, _val):
        if isinstance(_key, PhononContributionTaskParams):
            params = _key
        else:
            strain, key = _key
            params = PhononContributionTaskParams.create(strain, key)
        super().__setitem__(params, _val)

    def __getitem__(self, _key):
        if isinstance(_key, PhononContributionTaskParams):
            params = _key
        else:
            strain, key = _key
            params = PhononContributionTaskParams.create(strain, key)
        return next(v for p, v in self.data.items() if p == params)


class PhononContributionTask:

    def __init__(self, strain: tuple, key: C_, calculator = None):

        self.key = key
        self._task_params = PhononContributionTaskParams.create(strain, key)

        if self.key.is_longitudinal:
            self.calculator = LongitudinalElasticModulusPhononContribution(calculator, self.params)
        elif self.key.is_off_diagonal:
            self.calculator = OffDiagonalElasticModulusPhononContribution(calculator, self.params)
        elif self.key.is_shear:
            self.calculator = ShearElasticModulusPhononContribution(*self.params, calculator)
        
    @property
    def calc_type(self) -> ElasticModulusCalculationType:
        return self.key.calc_type
    
    @property
    def task_params(self) -> PhononContributionTaskParams:
        return self._task_params
    
    @property
    def params(self) -> tuple:
        return self.task_params.params
    
    def get_dependencies(self) -> List[Tuple[Tuple, C_]]:
        if self.calc_type in {
            ElasticModulusCalculationType.LONGITUDINAL,
            ElasticModulusCalculationType.OFF_DIAGONAL
        }: return []
        return list(itertools.product(
                [self.calculator.strain],
                self.calculator.get_modulus_keys()
            )) + list(itertools.product(
                [self.calculator.strain_rotated],
                self.calculator.get_modulus_keys_rotated()
            ))

    def get_modulus_isothermal(self) -> numpy.ndarray:
        '''Get phonon contribution of the isothermal elastic modulus as the 
        result of this task.
        '''
        if self.calc_type == ElasticModulusCalculationType.SHEAR:
            self.calculator.modulus = self.modulus_results
            self.calculator.modulus_rotated = self.modulus_results_rotated
        return self.calculator.value_isothermal

    def get_modulus_adiabatic(self) -> numpy.ndarray:
        '''Get phonon contribution of the adiabatic elastic modulus as the
        result of this task.
        '''
        if self.calc_type == ElasticModulusCalculationType.SHEAR:
            self.calculator.modulus = self.modulus_results
            self.calculator.modulus_rotated = self.modulus_results_rotated
        return self.calculator.value_adiabatic


class PhononContributionTaskList(UserList):

    #type data: List[PhononContributionTask]

    def __init__(self, calculator):

        super().__init__()

        self.calculator = calculator

        self.modulus_isothermal_values = PhononContributionTaskResults()
        self.modulus_adiabatic_values = PhononContributionTaskResults()

    def resolve(self, strain: tuple, keys: Iterable[C_]) -> None:
        '''Create the list of phonon calculation tasks based on the initial
        strain and the :math:`c_{ij}` keys in topological order, i.e. each task
        should have a smaller index than its dependant.

        Example
            :math:`c_{44}` depends on :math:`c_{2'2'}`, :math:`c_{2'3'}` and 
            :math:`c_{3'3'}`, the order of calculation tasks should be:

                [:math:`c_{2'2'}`, :math:`c_{2'3'}`, :math:`c_{3'3'}`,
                :math:`c_{44}`]

        :param strain: The strains :math:`e_1/\\Delta`, :math:`e_2/\\Delta`,
            :math:`e_3/\\Delta` in the orthogonal coordinates system under
            hydrostatic pressure.
        :param keys: The symbols :math:`c_{ij}`s needs to be calculated.
        '''

        self.strain = strain # e1, e2, e3
        self.keys = keys # c11, c22, c33, ...

        q = list(itertools.product([strain], keys, [None])) # strain, key, dependency

        tasks = [] 

        graph = nx.DiGraph()

        while not len(q) == 0:

            strain, key, dep = q.pop()

            # print("202 -> ", strain)
            # print(key, dep, strain[0, :])
            task_params = PhononContributionTaskParams.create(strain, key)
            task = next((t for t in tasks if t.task_params == task_params), None)

            if task is None:
                task = PhononContributionTask(strain, key, self.calculator)
                tasks.append(task)
                curr = len(tasks) - 1
                graph.add_node(curr)
                logger.debug(f"new task -> #{curr}")
            else:
                curr = tasks.index(task)

            if dep is not None:
                graph.add_edge(curr, dep)
                logger.debug(f"new task dependency -> #{dep} depends on #{curr}")

            for _strain, _key in task.get_dependencies():
                q.append((_strain, _key, curr))

        orders = nx.topological_sort(graph)

        self.data = [tasks[idx] for idx in orders]

        for task in self.data:
            logger.debug(str(task.calc_type) + ":" + str(task.params))

        self._graph = graph
        self._tasks = tasks

    def calculate(self) -> None:
        
        for task in self.data:
            # type task: PhononContributionTask
            if task.calc_type == ElasticModulusCalculationType.SHEAR:
                # TODO: do not set value directly
                # print(list(task.calculator.get_modulus_keys()))
                # print(list(task.calculator.get_modulus_keys_rotated()))
                task.modulus_results = self.modulus_isothermal_values.get_results_by_strain_keys(
                    task.calculator.strain,
                    task.calculator.get_modulus_keys()
                )
                task.modulus_results_rotated = self.modulus_isothermal_values.get_results_by_strain_keys(
                    task.calculator.strain_rotated,
                    task.calculator.get_modulus_keys_rotated()
                )

                # print("03 ->", list(task.calculator.get_modulus_keys()))
                # print("03 ->", list(task.calculator.get_modulus_keys_rotated()))

                # print("04 ->", task.modulus_results.keys())
                # print("04 ->", task.modulus_results_rotated.keys())

            self.modulus_isothermal_values[task.task_params] = task.get_modulus_isothermal()
            self.modulus_adiabatic_values[task.task_params] = task.get_modulus_adiabatic()
    
    def get_adiabatic_results(self) -> dict:
        return self.modulus_adiabatic_values.get_results_by_strain_keys(self.strain, self.keys)

    def get_isothermal_results(self) -> dict:
        return self.modulus_isothermal_values.get_results_by_strain_keys(self.strain, self.keys)

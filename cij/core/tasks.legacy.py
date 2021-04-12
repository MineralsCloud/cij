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

from logging import Logger
logger = Logger(__file__)

class PhononContributionTaskParams(NamedTuple):

    calc_type: ElasticModulusCalculationType
    params: dict

    @staticmethod
    def make_param_by_strain_key(strain: tuple, key: C_):
        if key.is_shear:
            return strain, key
        else:
            i, j, k, l = key.s
            return tuple(sorted([
                strain[i-1] / sum(strain),
                strain[j-1] / sum(strain)
            ]))

    @classmethod
    def create(cls, strain: tuple, key: C_) -> cls:
        return cls(key.calc_type, cls.make_param_by_strain_key(strain, key))
    
    def __eq__(self, other: PhononContributionTaskParams) -> bool:
        if self.calc_type != other.calc_type: return False
        if self.calc_type == ElasticModulusCalculationType.SHEAR:
            if self.params[0] != other.params[0]: return False
            if not numpy.allclose(self.params[1], other.params[1]): return False
            return True
        else:
            if not numpy.allclose(self.params, other.params): return False
            return True


def get_param_by_strain_key(strain: tuple, key: C_):
    if key.is_shear:
        return (strain, key)
    else:
        i, j, k, l = key.s
        # print(strain)
        # if not isinstance(strain, tuple):
        #     raise RuntimeError()
        # print(f"strain {i}, {j} -> ", [
        #     strain[i-1] / sum(strain),
        #     strain[j-1] / sum(strain)
        # ])
        return tuple(sorted([
            strain[i-1] / sum(strain),
            strain[j-1] / sum(strain)
        ]))


class PhononContributionResults(UserDict):

    def get_results_by_strain_key(self, strain: tuple, keys: list) -> dict:
        results = dict()
        # print("01 ->", list(keys))
        for key in keys:
            # print("01.5 ->", key)
            results[key] = self[strain, key]
        # print("02 ->", results.keys())
        return results

    def __getitem__(self, _key):
        strain, key = _key
        # print(strain, key)
        idx = index_result_by_strain_key(self.data.keys(), strain, key)
        # print(self.data.keys())
        # print(idx)
        if idx is None:
            raise RuntimeError()
        return list(self.data.values())[idx]


class PhononContributionTask:

    def __init__(self, strain: tuple, key: C_, calculator = None):

        self.key = key

        self.params = get_param_by_strain_key(strain, key)

        if self.key.is_longitudinal:
            self.calculator = LongitudinalElasticModulusPhononContribution(calculator, self.params)

        elif self.key.is_off_diagonal:
            self.calculator = OffDiagonalElasticModulusPhononContribution(calculator, self.params)

        else:
            self.calculator = ShearElasticModulusPhononContribution(*self.params)
    
    @property
    def calc_type(self) -> ElasticModulusCalculationType:
        return self.key.calc_type
    
    def resolve_dependencies(self):
        if self.calc_type in {
            ElasticModulusCalculationType.LONGITUDINAL,
            ElasticModulusCalculationType.OFF_DIAGONAL
        }: return []
        return list(itertools.product(
                [self.calculator.strain], self.calculator.get_modulus_keys()
            )) + list(itertools.product(
                [self.calculator.strain_rotated], self.calculator.get_modulus_keys_rotated()
            ))
    
    def get_modulus_isothermal(self, results: Union[None, PhononContributionResults] = None):
        if self.calc_type == ElasticModulusCalculationType.SHEAR:
            self.calculator.modulus = self.modulus_results
            self.calculator.modulus_rotated = self.modulus_results_rotated
            # return 
        return self.calculator.value_isothermal

    def get_modulus_adiabatic(self, results: Union[None, PhononContributionResults] = None):
        if self.calc_type == ElasticModulusCalculationType.SHEAR:
            self.calculator.modulus = self.modulus_results
            self.calculator.modulus_rotated = self.modulus_results_rotated
        return self.calculator.value_adiabatic


def index_task_by_strain_key(
    tasks: Iterable[PhononContributionTask], strain: tuple, key: C_
) -> Union[int, None]:
    for i, task in enumerate(tasks):
        if task.calc_type != key.calc_type: continue
        if task.calc_type != ElasticModulusCalculationType.SHEAR:
            params = get_param_by_strain_key(strain, key)
            if numpy.allclose(params, task.params):
                return i
        else:
            params = get_param_by_strain_key(strain, key)
            if params[1] == task.params[1] and numpy.allclose(params[0], task.params[0]):
                    return i
    return None

def index_result_by_strain_key(
    tasks: Iterable[PhononContributionTask], strain: tuple, key: C_
) -> Union[int, None]:
    for i, (calc_type, params) in enumerate(tasks):
        if calc_type != key.calc_type: continue
        if calc_type != ElasticModulusCalculationType.SHEAR:
            _params = get_param_by_strain_key(strain, key)
            if numpy.allclose(params, _params):
                return i
        else:
            _params = get_param_by_strain_key(strain, key)
            if params[1] == _params[1] and numpy.allclose(params[0], _params[0]):
                return i
    return None


class PhononContributionTaskList(UserList):

    #type data: List[PhononContributionTask]

    def __init__(self, calculator):

        super().__init__()

        self.calculator = calculator

        self.modulus_isothermal_values = PhononContributionResults()
        self.modulus_adiabatic_values = PhononContributionResults()

    def resolve(self, strain, keys) -> None:

        self.strain = strain
        self.keys = keys

        q = list(itertools.product([strain], keys, [None]))

        tasks = []

        graph = nx.DiGraph()

        while not len(q) == 0:

            strain, key, dep = q.pop()
            print(key, dep)

            curr = index_task_by_strain_key(tasks, strain, key)

            if curr is None:
                task = PhononContributionTask(strain, key, self.calculator)
                tasks.append(task)
                curr = len(tasks) - 1
                graph.add_node(curr)
                print("add ->", curr)
            else:
                task = tasks[curr]

            if dep is not None:
                graph.add_edge(curr, dep)
                print("edge ->", curr, dep)

            # print(q)
            for _strain, _key in task.resolve_dependencies():

                print((_strain, _key, curr))
                q.append((_strain, _key, curr))

        orders = nx.topological_sort(graph)

        self.data = [tasks[idx] for idx in orders]

        for task in self.data:
            print(task.calc_type, task.params)

        self._graph = graph
        self._tasks = tasks

    def calculate(self) -> None:
        
        for task in self.data:
            # type task: PhononContributionTask
            if task.calc_type == ElasticModulusCalculationType.SHEAR:
                # print(list(task.calculator.get_modulus_keys()))
                # print(list(task.calculator.get_modulus_keys_rotated()))
                task.modulus_results = self.modulus_isothermal_values.get_results_by_strain_key(
                    task.calculator.strain,
                    task.calculator.get_modulus_keys()
                )
                task.modulus_results_rotated = self.modulus_isothermal_values.get_results_by_strain_key(
                    task.calculator.strain_rotated,
                    task.calculator.get_modulus_keys_rotated()
                )

                print("03 ->", list(task.calculator.get_modulus_keys()))
                print("03 ->", list(task.calculator.get_modulus_keys_rotated()))

                print("04 ->", task.modulus_results.keys())
                print("04 ->", task.modulus_results_rotated.keys())

            self.modulus_isothermal_values[task.calc_type, task.params] = task.get_modulus_isothermal()
            self.modulus_adiabatic_values[task.calc_type, task.params] = task.get_modulus_adiabatic()
    
    def get_adiabatic_results(self) -> dict:
        return self.modulus_adiabatic_values.get_results_by_strain_key(self.strain, self.keys)

    def get_isothermal_results(self) -> dict:
        return self.modulus_isothermal_values.get_results_by_strain_key(self.strain, self.keys)

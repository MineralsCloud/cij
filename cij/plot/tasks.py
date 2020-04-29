import networkx as nx
from matplotlib.axes import Axes
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from typing import Union

from cij.core.calculator import Calculator
from cij.core.tasks import PhononContributionTaskList, PhononContributionTask
from cij.util import C_, ElasticModulusCalculationType


_default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def _make_task_color(task: PhononContributionTask, colors: list = _default_colors) -> str:
    calc_type = task.calc_type
    if calc_type == ElasticModulusCalculationType.LONGITUDINAL:
        return colors[0]
    elif calc_type == ElasticModulusCalculationType.OFF_DIAGONAL:
        return colors[1]
    else:
        key = task.key
        if key.i == key.j:
            return colors[2]
        else:
            return colors[3]


def _make_task_label(task: PhononContributionTask) -> str:
    calc_type = task.calc_type
    format_strain = lambda c: "(" + ",".join([("%.2f" % x).lstrip("0") for x in c]) + ")"
    if calc_type == ElasticModulusCalculationType.SHEAR:
        return "$c_{%d%d, %s}$" % (*task.params[1].v, format_strain(task.params[0]))
    else:
        return format_strain(task.params)


def plot_phonon_contribution_dependencies(calculator: Calculator, **kwargs):
    return plot_tasklist_dependencies(
        calculator._full_modulus._phonon_contribution_task_list,
        **kwargs
    )


def plot_tasklist_dependencies(task_list: PhononContributionTaskList, **kwargs):

    if "node_color" not in kwargs.keys():
        kwargs.update({"node_color": [
            _make_task_color(task)
            for task in task_list._tasks
        ]})

    if "labels" not in kwargs.keys():
        kwargs.update({"labels": dict(
            (idx, _make_task_label(task))
            for idx, task in enumerate(task_list._tasks)
        )})

    if "font_size" not in kwargs.keys():
        kwargs.update({"font_size": 10})

    return nx.draw(task_list._graph, **kwargs)


def make_legend(ax: Union[Axes, None] = None, node_color=_default_colors):
    if ax is None: ax = plt.gca()
    return ax.legend([
        Line2D([], [], color="none", marker='o', markerfacecolor=c, markeredgecolor="none")
        for c in node_color[:4]
    ], [
        r"$c_{iiii}$ (longitud)",
        r"$c_{iijj}$ (off-diag)",
        r"$c_{ijij}$ (shear)",
        r"$c_{ijkl}$ (shear)"
    ])
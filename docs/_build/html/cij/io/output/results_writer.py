from typing import Union, Callable, TypeVar, Iterable, NamedTuple, Optional, List, Dict
from enum import Enum, auto
import itertools
import yaml
from logging import Logger

from cij.data import get_data_fname
from cij.util import convert_unit, units, C_
# from cij.core.calculator import CijPressureBaseInterface, CijVolumeBaseInterface

from ..traditional.qha_output import save_x_pt


_T = TypeVar("_T")
logger = Logger(__name__)

def _load_writer_rules_file(fname: Optional[str] = None):
    if fname is None:
        fname = get_data_fname("output/writer_rules.yml")
    with open(fname) as fp:
        return yaml.safe_load(fp)

DEFAULT_WRITER_RULES = _load_writer_rules_file()

class ResultsWriterRule(NamedTuple):

    keywords: List[str]
    fname_pattern: str
    prop: str
    unit: Union[units.Unit, str]
    unit_internal: Union[units.Unit, str]
    var_type: str

    @classmethod
    def create(cls, rule: dict) -> 'ResultsWriterRule':
        return cls(
            keywords=rule["keywords"],
            fname_pattern=rule["fname_pattern"],
            prop=rule["prop"],
            unit=rule["unit"],
            unit_internal=rule["unit_internal"],
            var_type=rule["var_type"]
        )
    
    @staticmethod
    def _format_ij(key: C_) -> str:
        return "%d%d" % key.v
    
    def write_variable(self, base, config: Optional[dict] = None):

        _config = self._asdict()
        if config is not None:
            _config.update(config)

        convert = convert_unit(_config["unit_internal"], _config["unit"])

        variable = getattr(base, self.prop)

        if "fname" in _config:
            fname = config["fname"]
        else:
            fname = self.fname_pattern.format(base=base._base_name)

        logger.info(f"Writing output <{fname}>.")
        base.write_table(fname, convert(variable))

    def write_ij_variable(self, base, config: Optional[dict] = None):

        _config = self._asdict()
        if config is not None:
            _config.update(config)

        convert = convert_unit(_config["unit_internal"], _config["unit"])

        variable = getattr(base, self.prop)

        for k, v in variable.items():
            if "fname" in _config:
                fname = config["fname"]
            else:
                fname = self.fname_pattern.format(
                    base=base._base_name,
                    ij=self._format_ij(k)
                )

            logger.info(f"Writing output <{fname}>.")
            base.write_table(fname, convert(v))

    def write(self, base, config: Optional[dict] = None):
        if self.var_type == 'value':
            self.write_variable(base, config)
        elif self.var_type == 'ij_value':
            self.write_ij_variable(base, config)
        else:
            raise NotImplementedError(f"Unknown variable type {self.var_type}")

class ResultsWriter:

    def __init__(self, base, rules: Optional[dict] = None):
        self.registry = {} # type: Dict[str, ResultsWriterRule]
        if rules is None:
            rules = DEFAULT_WRITER_RULES
        self._init_rules(rules)
        self.base = base
    
    def _init_rules(self, rules):
        for rule in rules:
            rule = ResultsWriterRule.create(rule)
            for keyword in rule.keywords:
                self.registry[keyword] = rule
    
    def write(self, config: Union[str, dict]):
        if isinstance(config, str):
            config = { "keyword": config }
        self.registry[config["keyword"]].write(self.base, config)
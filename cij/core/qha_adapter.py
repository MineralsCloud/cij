import logging
import json
import numpy
import copy

import qha.calculator
import qha.basic_io
from qha.settings import DEFAULT_SETTINGS

import cij.io.traditional.models

logger = logging.getLogger(__name__)

class QHACalculator(qha.calculator.Calculator):

    def __init__(self, settings: dict):

        super().__init__(settings)

    def read_input(self, qha_input: cij.io.traditional.models.QHAInputData) -> None:

        self._formula_unit_number: int = qha_input.nm

        self._volumes = numpy.array([
            volume.volume for volume in qha_input.volumes
        ])
        self._static_energies = numpy.array([
            volume.energy for volume in qha_input.volumes
        ]) 
        self._frequencies = numpy.array([
            [modes for q_coord, modes in volume.q_points]
            for volume in qha_input.volumes
        ])
        self._q_weights = numpy.array([
            weight for q_coord, weight in qha_input.weights
        ])

    def desired_pressure_status(self) -> None:

        logger.info(
           "The pressure range can be dealt with: "
           "[{0:6.2f} to {1:6.2f}] GPa".format(
               self.p_tv_gpa[:, 0].max(), self.p_tv_gpa[:, -1].min()
        ))

        if self.p_tv_gpa[:, -1].min() < self.desired_pressures_gpa.max():
            ntv_max = int(
                (self.p_tv_gpa[:, -1].min() - self.desired_pressures_gpa.min()) / self.settings['DELTA_P'])

            logger.error(
                "!!!ATTENTION!!!\n"
                "DESIRED PRESSURE is too high (NTV is too large)!\n"
                "QHA results might not be right!\n"
                "Please reduce the NTV accordingly, for example,"
                "try to set NTV < {:4d}.\n".format(ntv_max)
            )

            raise ValueError(
                "DESIRED PRESSURE is too high (NTV is too large), qha results might not be right!")

class QHACalculatorAdapter():

    def __init__(self, settings, qha_input: cij.io.traditional.models.QHAInputData):

        self.calculator = self._load_qha_calculator(settings, qha_input)
        self.volume_base_results = QHAVolumeBaseInterface(self.calculator)
        self.pressure_base_results = QHAPressureBaseInterface(self.calculator)

    @staticmethod
    def _load_qha_calculator(settings: dict, qha_input: cij.io.traditional.models.QHAInputData):
        
        user_settings = copy.copy(DEFAULT_SETTINGS)
        user_settings.update(settings)

        calculator = QHACalculator(user_settings)

        logger.debug("QHA info %s" % json.dumps(calculator.settings))
        logger.debug(
            qha.basic_io.out.make_tp_info(
                calculator.temperature_array[0],
                calculator.temperature_array[-1 - 4],
                calculator.desired_pressures_gpa[0],
                calculator.desired_pressures_gpa[-1]
            )
        )

        calculator.read_input(qha_input)

        logger.info(
            "Caution: If negative frequencies found, "
            "they are currently treated as 0!")

        tmp = calculator.where_negative_frequencies
        # Don't delete this parenthesis!
        if tmp is not None and not (tmp.T[-1].max() <= 2):
            for indices in tmp:
                logging.info(
                    "Found negative frequency in {0}th volume "
                    "{1}th q-point {2}th band".format(*tuple(indices + 1)))

        calculator.refine_grid()

        logger.info(
            "The volume range used in this calculation expanded"
            "x {0:6.4f}".format(calculator.v_ratio))

        calculator.desired_pressure_status()

        temperature_array = calculator.temperature_array
        desired_pressures_gpa = calculator.desired_pressures_gpa
        temperature_sample = calculator.temperature_sample_array
        p_sample_gpa = calculator.pressure_sample_array

        return calculator

    @property
    def v_array(self):
        return self.calculator.finer_volumes_bohr3

    @property
    def t_array(self):
        return self.calculator.temperature_array

    @property
    def t_sample_array(self):
        return self.calculator.temperature_sample_array

    @property
    def ntv(self):
        return len(self.v_array)
    
    def v2p(self):
        pass

    @property
    def volume_base(self):
        return self.volume_base_results

    @property
    def pressure_base(self):
        return self.pressure_base_results

class QHAVolumeBaseInterface:

    def __init__(self, calculator: QHACalculator):
        self.calculator = calculator

    @property
    def v_array(self):
        return self.calculator.finer_volumes_bohr3

    @property
    def t_array(self):
        return self.calculator.temperature_array
    @property
    def t_sample_array(self):
        return self.calculator.temperature_sample_array

    @property
    def gibbs_free_energies(self):
        return self.calculator.g_tv_ry
    @property
    def helmholtz_free_energies(self):
        return self.calculator.f_tv_ry
    @property
    def enthalpies(self):
        return self.calculator.h_tv_ry

    @property
    def entropies(self):
        raise NotImplementedError()
    @property
    def pressures(self):
        return self.calculator.p_tv_au

    @property
    def thermal_expansivities(self):
        return self.calculator.alpha_tv
    @property
    def bulk_modulus(self):
        return self.calculator.bt_tv_au
    @property
    def bulk_modulus_isothermal(self):
        return self.calculator.bs_tv_au
    @property
    def heat_capacity(self):
        return self.calculator.cv_tv_au

class QHAPressureBaseInterface:

    def __init__(self, calculator):
        self.calculator = calculator

    @property
    def p_array(self):
        return self.calculator.desired_pressures
    @property
    def t_array(self):
        return self.calculator.temperature_array
    @property
    def t_sample_array(self):
        return self.calculator.temperature_sample_array

    @property
    def gibbs_free_energies(self):
        return self.calculator.g_tp_ry
    @property
    def helmholtz_free_energies(self):
        return self.calculator.f_tp_ry
    @property
    def enthalpies(self):
        return self.calculator.h_tp_ry

    @property
    def entropies(self):
        raise NotImplementedError()
    @property
    def volumes(self):
        return self.calculator.v_tp_bohr3

    @property
    def thermal_expansivities(self):
        return self.calculator.alpha_tp
    @property
    def bulk_modulus(self):
        return self.calculator.bt_tp_au
    @property
    def bulk_modulus_isothermal(self):
        return self.calculator.bs_tp_au
    @property
    def heat_capacity(self):
        return self.calculator.cv_tp_au
import click

from logging import getLogger

logger = getLogger(__name__)


@click.command("run-static", help="Fit equation of state based on E, V given in INPUT01, and also elastic moduli and acoustic velocities if static Cij table INPUT02 (elast.dat) is given.")
@click.argument("input01", type=click.Path(exists=True))
@click.argument("input02", required=False, type=click.Path(exists=True))
@click.option("-I", "--interp", type=click.Choice(["none", "pressure", "volume"]), default="none", show_default=True, help="Interpolate result vs. volume or pressure or don't interpolate (none).")
@click.option("-n", "--ntv", type=click.INT, default=201, show_default=True, help="Number of volumes (or equivivalently, pressures) on the grid when interpolating results.")
@click.option("--p-min", type=click.FLOAT, default=0, show_default=True, help="Lower bound of pressure range when interpolate result vs. P, in GPa unit.")
@click.option("--delta-p", type=click.FLOAT, default=1.0, show_default=True, help="The interval between two nearest pressures on the grid when interpolate result vs. P, in GPa unit.")
@click.option("--delta-p-sample", type=click.FLOAT, default=None, help="The interval between two nearest pressures on the grid when interpolate result vs. P, in GPa unit, same as DELTA_P by default.")
@click.option("--cellmass", type=click.FLOAT, help="Mass of the unitcell, value in amu unit.")
@click.option("--v-ratio", type=click.FLOAT, default=1.2, show_default=True, help="Ratio to expand the volume range when extrapolate equation of states. (from (vmin, vmax) to (vmin / volume_ratio, vmax * volume_ratio)).")
@click.option("-s", "--system", default=None, help="Name of the crystal system whose symmetry is applied to fill the missing elastic tensor components. Should be Should be one of: triclinic, monoclinic, hexagonal, trigonal6, trigonal7, orthorhombic, tetragonal6, tetragonal7, cubic.")
def main(input01: str, input02: str, interp: str, ntv: int, cellmass: float, v_ratio: float, p_min: float, delta_p: float, delta_p_sample: float, system: str = None):

    import pandas
    import numpy
    import sys
    import itertools
    from pathlib import Path

    from scipy.interpolate import InterpolatedUnivariateSpline

    from qha.fitting import polynomial_least_square_fitting
    from qha.grid_interpolation import calculate_eulerian_strain
    from qha.v2p import v2p

    from cij.io.traditional import read_elast_data, read_energy
    from cij.util import c_
    from cij.util.units import convert_unit, _from_gpa, _to_gpa, _from_ang3, _to_ang3, _to_gcm3, _to_ev, _to_kms
    from cij.util.fill import fill_cij
    from cij.data import get_data_fname

    def fit_modulus(volumes: numpy.ndarray, v_array: numpy.ndarray, moduli: numpy.ndarray, order: int = 2) -> numpy.ndarray:
        '''Interpolate static elastic constants :math:`c^\\text{st}_{ij}(V)` as
        a function of volume with polynomial least square fitting

        :param moduli: The static elastic moduli :math:`c^\\text{st}_{ij}(V)`
            to be fitted
        :param order: The order for least square fitting

        :returns: The interpolated modulus :math:`c^\\text{st}_{ij}(V)`
        '''

        strains = calculate_eulerian_strain(volumes[0], volumes)
        strain_array = calculate_eulerian_strain(volumes[0], v_array)
        _, modulus_array = polynomial_least_square_fitting(
            strains, moduli, strain_array,
            order=order
        )
        return modulus_array

    def v2p1d(x_old: numpy.array, p_old: numpy.array, p_new: numpy.array):
        from numpy import newaxis as nax
        return v2p(x_old[nax, ::-1], p_old[nax, ::-1], p_new)[0]


    input01 = read_energy(input01)
    if input02:
        input02 = read_elast_data(input02)

    df = pandas.DataFrame(index=range(input01.nv))

    for i in df.index:
        df.loc[i, "V"] = input01.volumes[i].volume
        df.loc[i, "F"] = input01.volumes[i].energy

    volumes = df.loc[:, "V"].to_numpy()

    # fit, pressure

    v_array = numpy.linspace(numpy.min(volumes) / v_ratio, numpy.max(volumes) * v_ratio, ntv)

    energies = df.loc[:, "F"].to_numpy()
    f_array = fit_modulus(volumes, v_array, energies)
    p_array = - numpy.gradient(f_array) / numpy.gradient(v_array)

    if interp == "none":
        df.loc[:, "P"] = InterpolatedUnivariateSpline(v_array, p_array)(volumes)

    elif interp == "volume":
        df = pandas.DataFrame(index=range(ntv))
        df.loc[:, "V"] = v_array
        df.loc[:, "F"] = f_array
        df.loc[:, "P"] = p_array

    elif interp == "pressure":
        _p_array = numpy.linspace(
            _from_gpa(p_min), _from_gpa(p_min + delta_p * (ntv - 1)), ntv)
        _v_array = v2p1d(v_array, p_array, _p_array)
        _f_array = v2p1d(v_array, p_array, _p_array)

        df = pandas.DataFrame(index=range(ntv))
        df.loc[:, "V"] = _v_array
        df.loc[:, "F"] = _f_array
        df.loc[:, "P"] = _p_array

    v_array = df.loc[:, "V"].to_numpy()

    # cij

    if input02:

        df.loc[:, "density"] = input02.cellmass / df.loc[:, "V"]

        for key in input02.volumes[0].static_elastic_modulus.keys():
            volumes, moduli = numpy.array([
                (volume.volume, volume.static_elastic_modulus[key])
                for volume in input02.volumes
            ]).T
            strkey = "c%d%d" % key.voigt
            df.loc[:, strkey] = fit_modulus(volumes, v_array, moduli)

    # fill all non-zero terms based on symmetry constraints
    
    if system == None:
        logger.warning(f"Symmetry constraints check not performed! Make sure to fill in all non-zero terms for correct VRH averages!")
    else:
        df = fill_cij(df, system)
    
    if cellmass:
        df.loc[:, "density"] = cellmass / df.loc[:, "V"]
    
    # VRH

    if input02:

        cij = numpy.zeros((df.shape[0], 6, 6))
        for i, j in itertools.product(range(6), range(6)):
            key = "c%d%d" % tuple(sorted((i+1, j+1)))
            if key in df.columns:
                cij[:, i, j] = df.loc[:, key]
        
        sij = numpy.linalg.inv(cij)

        c = numpy.zeros((cij.shape[0], 7, 7))
        s = numpy.zeros((sij.shape[0], 7, 7))

        c[:, 1:, 1:] = cij[:, :, :]
        s[:, 1:, 1:] = sij[:, :, :]

        df.loc[:, "bm_V"] = (c[:,1,1] + c[:,2,2] + c[:,3,3] \
                    + 2 * (c[:,1,2] + c[:,2,3] + c[:,1,3])) / 9
        df.loc[:, "bm_R"] = 1 / (s[:,1,1] + s[:,2,2] + s[:,3,3] \
                    + 2 * (s[:,1,2] + s[:,2,3] + s[:,1,3]))
        df.loc[:, "bm_VRH"] = (df.loc[:, "bm_V"] + df.loc[:, "bm_R"]) / 2

        df.loc[:, "G_V"] = ( \
                    + (c[:,1,1] + c[:,2,2] + c[:,3,3]) \
                    - (c[:,1,2] + c[:,2,3] + c[:,1,3])
                    + 3 * (c[:,4,4] + c[:,5,5] + c[:,6,6])) / 15
        df.loc[:, "G_R"] = 15 / ( \
                + 4 * (s[:,1,1] + s[:,2,2] + s[:,3,3]) \
                - 4 * (s[:,1,2] + s[:,2,3] + s[:,1,3]) \
                + 3 * (s[:,4,4] + s[:,5,5] + s[:,6,6]))
        df.loc[:, "G_VRH"] = (df.loc[:, "G_V"] + df.loc[:, "G_R"]) / 2


    # Unit conversion

    df["V"] = _to_ang3(df["V"].to_numpy())
    df["F"] = _to_ev(df["F"].to_numpy())
    df["P"] = _to_gpa(df["P"].to_numpy())

    if "density" in df.columns:
        df["density"] = _to_gcm3(df["density"].to_numpy())

    # velocity

    if input02:
        df.loc[:, "v_p"] = numpy.sqrt((df.loc[:, "bm_VRH"] + 4 / 3 * df.loc[:, "G_VRH"]) /  df.loc[:, "density"])
        df.loc[:, "v_s"] = numpy.sqrt(df.loc[:, "G_VRH"] / df.loc[:, "density"])
        df.loc[:, "v_phi"] = numpy.sqrt(df.loc[:, "bm_VRH"] / df.loc[:, "density"])

        df["v_p"] = _to_kms(df["v_p"].to_numpy())
        df["v_s"] = _to_kms(df["v_s"].to_numpy())
        df["v_phi"] = _to_kms(df["v_phi"].to_numpy())

    # sample
    if interp == "pressure" and delta_p_sample:
        step = round(delta_p_sample / delta_p)
        df = df.iloc[::step, :]

    sys.stdout.write(df.to_string())

if __name__ == "__main__":
    main()

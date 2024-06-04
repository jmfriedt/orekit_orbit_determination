# https://github.com/GorgiAstro/tle-fitting/blob/master/01-tle-fitting-example.ipynb

sc_mass = 400.0 # kg
sc_cross_section = 0.3 # m2
cd_drag_coeff = 2.0 
cr_radiation_pressure = 1.0

from datetime import datetime
date_start = datetime(2023, 5, 20)
fitting_duration_d = 1 # days

import numpy as np
a = 7000.0e3  # meters
e = 0.001
i = float(np.deg2rad(98.0))  # Conversion to Python float is required for Orekit
pa = float(np.deg2rad(42.0))
raan = float(np.deg2rad(42.0))
ma = float(np.deg2rad(42.0))  # Mean anomaly

satellite_number = 99999
classification = 'X'
launch_year = 2018
launch_number = 42
launch_piece = 'F'
ephemeris_type = 0
element_number = 999
revolution_number = 100

dt = 60.0  # s, period at which the spacecraft states are saved to fit the TLE

prop_min_step = 0.001 # s
prop_max_step = 300.0 # s
prop_position_error = 10.0 # m

# Estimator parameters
estimator_position_scale = 1.0 # m
estimator_convergence_thres = 1e-3
estimator_max_iterations = 25
estimator_max_evaluations = 35

import orekit
orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir
orekit_data_dir = 'orekit-data-master' # unzipped from orekit-data.zip
setup_orekit_curdir(orekit_data_dir)

from org.orekit.frames import FramesFactory, ITRFVersion
from org.orekit.utils import IERSConventions
gcrf = FramesFactory.getGCRF()
teme = FramesFactory.getTEME()
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)

from org.orekit.models.earth import ReferenceEllipsoid
wgs84_ellipsoid = ReferenceEllipsoid.getWgs84(itrf)

from org.orekit.bodies import CelestialBodyFactory
moon = CelestialBodyFactory.getMoon()
sun = CelestialBodyFactory.getSun()

from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.utils import Constants as orekit_constants
from orekit.pyhelpers import datetime_to_absolutedate
date_start_orekit = datetime_to_absolutedate(date_start)
keplerian_orbit = KeplerianOrbit(a, e, i, pa, raan, ma, PositionAngleType.MEAN, 
                                 gcrf, date_start_orekit, orekit_constants.EIGEN5C_EARTH_MU)

from org.orekit.propagation.analytical.tle import TLE
mean_motion = float(np.sqrt(orekit_constants.EIGEN5C_EARTH_MU / np.power(a, 3)))
mean_motion_first_derivative = 0.0
mean_motion_second_derivative = 0.0
b_star_first_guess = 1e-5  # Does not play any role, because it is a free parameter when fitting the TLE

tle_first_guess = TLE(satellite_number, 
                        classification,
                        launch_year,
                        launch_number,
                        launch_piece,
                        ephemeris_type,
                        element_number,
                        date_start_orekit,
                        mean_motion,
                        mean_motion_first_derivative, 
                        mean_motion_second_derivative,
                        e,
                        i,
                        pa,
                        raan,
                        ma,
                        revolution_number,
                        b_star_first_guess)

print(tle_first_guess)

from org.orekit.attitudes import NadirPointing
nadir_pointing = NadirPointing(gcrf, wgs84_ellipsoid)

from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
integrator_builder = DormandPrince853IntegratorBuilder(prop_min_step, prop_max_step, prop_position_error)

from org.orekit.propagation.conversion import NumericalPropagatorBuilder
propagator_builder = NumericalPropagatorBuilder(keplerian_orbit,
                                               integrator_builder, PositionAngleType.MEAN, estimator_position_scale)
propagator_builder.setMass(sc_mass)
propagator_builder.setAttitudeProvider(nadir_pointing)

# Earth gravity field with degree 64 and order 64
from org.orekit.forces.gravity.potential import GravityFieldFactory
gravity_provider = GravityFieldFactory.getConstantNormalizedProvider(64, 64, date_start_orekit) # update JMF
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
gravity_attraction_model = HolmesFeatherstoneAttractionModel(itrf, gravity_provider)
propagator_builder.addForceModel(gravity_attraction_model)

# Moon and Sun perturbations
from org.orekit.forces.gravity import ThirdBodyAttraction
moon_3dbodyattraction = ThirdBodyAttraction(moon)
propagator_builder.addForceModel(moon_3dbodyattraction)
sun_3dbodyattraction = ThirdBodyAttraction(sun)
propagator_builder.addForceModel(sun_3dbodyattraction)

# Solar radiation pressure
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
isotropic_radiation_single_coeff = IsotropicRadiationSingleCoefficient(sc_cross_section, cr_radiation_pressure);
from org.orekit.forces.radiation import SolarRadiationPressure
solar_radiation_pressure = SolarRadiationPressure(sun, wgs84_ellipsoid, isotropic_radiation_single_coeff)
propagator_builder.addForceModel(solar_radiation_pressure)

# Atmospheric drag
# from org.orekit.forces.drag.atmosphere.data import MarshallSolarActivityFutureEstimation JMF error
from org.orekit.models.earth.atmosphere.data import MarshallSolarActivityFutureEstimation 
# msafe = MarshallSolarActivityFutureEstimation('(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\p{Digit}\p{Digit}\p{Digit}\p{Digit}F10\.(?:txt|TXT)', MarshallSolarActivityFutureEstimation.StrengthLevel.AVERAGE) JMF error
msafe = MarshallSolarActivityFutureEstimation(MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES, MarshallSolarActivityFutureEstimation.StrengthLevel.AVERAGE)
from org.orekit.data import DataProvidersManager
#DM = DataProvidersManager.getProviders() # DataProvidersManager.getInstance() XXX JMF
#DM.feed(msafe.getSupportedNames(), msafe) # Feeding the F10.7 bulletins to Orekit's data manager JMF

from org.orekit.models.earth.atmosphere import NRLMSISE00 # change JMF
atmosphere = NRLMSISE00(msafe, sun, wgs84_ellipsoid)
from org.orekit.forces.drag import IsotropicDrag
isotropic_drag = IsotropicDrag(sc_cross_section, cd_drag_coeff)
from org.orekit.forces.drag import DragForce
drag_force = DragForce(atmosphere, isotropic_drag)
propagator_builder.addForceModel(drag_force)

propagator = propagator_builder.buildPropagator([a, e, i, pa, raan, ma])

from org.orekit.propagation import SpacecraftState
initial_state = SpacecraftState(keplerian_orbit, sc_mass)
propagator.resetInitialState(initial_state)
# propagator.setEphemerisMode() # JMF obsolete
date_end_orekit = date_start_orekit.shiftedBy(fitting_duration_d * 86400.0)
# state_end = propagator.propagate(date_end_orekit) JMF obsolete
generator=propagator.getEphemerisGenerator();
propagator.propagate(date_start_orekit, date_end_orekit);
state_end=generator.getGeneratedEphemeris();  # JMF problem

from java.util import ArrayList
states_list = ArrayList()
bounded_propagator = propagator.getGeneratedEphemeris()

date_current = date_start_orekit
while date_current.compareTo(date_end_orekit) <= 0:
    spacecraft_state = bounded_propagator.propagate(date_current)
    states_list.add(spacecraft_state)
    date_current = date_current.shiftedBy(dt)

from org.orekit.propagation.conversion import TLEPropagatorBuilder, FiniteDifferencePropagatorConverter
from org.orekit.propagation.analytical.tle import TLEPropagator
threshold = 1.0  # "absolute threshold for optimization algorithm", but no idea about its impact
tle_builder = TLEPropagatorBuilder(tle_first_guess, PositionAngleType.MEAN, 1.0)
fitter = FiniteDifferencePropagatorConverter(tle_builder, threshold, 1000)
fitter.convert(states_list, False, 'BSTAR')  # Setting BSTAR as free parameter
tle_propagator = TLEPropagator.cast_(fitter.getAdaptedPropagator())
tle_fitted = tle_propagator.getTLE()

print(tle_first_guess)
print('')
print(tle_fitted)

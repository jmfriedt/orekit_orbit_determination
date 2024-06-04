# JMF simulation ends at date of TLE => use *historical* TLE, *not* latest TLE

sat_list = {
    'envisat': {
        'norad_id': 27386,  # For Space-Track TLE queries
        'cospar_id': '0200901',  # For laser ranging data queries
        'sic_id': '6179',  # For writing in CPF files
        'mass': 8000.0, # kg; TODO: compute proper value
        'cross_section': 100.0, # m2; TODO: compute proper value
        'cd': 2.0, # TODO: compute proper value
        'cr': 1.0  # TODO: compute proper value
    },
    'lageos2': {
        'norad_id': 22195,
        'cospar_id': '9207002',
        'sic_id': '5986',
        'mass': 405.0, # kg
        'cross_section': 0.2827, # m2
        'cd': 2.0, # TODO: compute proper value
        'cr': 1.0  # TODO: compute proper value
    },
    'technosat': {
        'norad_id': 42829,
        'cospar_id': '1704205',
        'sic_id': '6203',
        'mass': 20.0, # kg
        'cross_section': 0.10, # m2,
        'cd': 2.0, # TODO: compute proper value
        'cr': 1.0  # TODO: compute proper value
    },
    'snet1': {
        'norad_id': 43189,
        'cospar_id': '1801410',
        'sic_id': '6204',
        'mass': 8.0, # kg
        'cross_section': 0.07,
        'cd': 2.0, # TODO: compute proper value
        'cr': 1.0  # TODO: compute proper value
    }
}

sc_name = 'lageos2'  # Change the name to select a different satellite in the dict

"""
NPT: Normal point data. Recommended option. The data is pre-filtered by the laser data providers
FRD: Full-rate data. Warning, these are a lot of data points (potentially tens of thousands per day),
    the execution time could be greatly increased
"""
laser_data_type = 'NPT'

range_weight = 1.0 # Will be normalized later (i.e divided by the number of observations)
range_sigma = 1.0 # Estimated covariance of the range measurements, in meters

import numpy as np
az_weight = 0.1  # Do not weigh the Az/El measurements too much because they are much less accurate than ranges
el_weight = 0.1
az_sigma = float(np.deg2rad(0.01))
el_sigma = float(np.deg2rad(0.01))

from datetime import datetime
odDate = datetime(2023, 6, 8) # Beginning of the orbit determination -> JMF: be consistent with downloaded TLE
collectionDuration = 2 # day
from datetime import timedelta
startCollectionDate = odDate + timedelta(days=-collectionDuration)

# Orbit propagator parameters
prop_min_step = 0.001 # s
prop_max_step = 300.0 # s
prop_position_error = 10.0 # m

# Estimator parameters
estimator_position_scale = 1.0 # m
estimator_convergence_thres = 1e-2
estimator_max_iterations = 25
estimator_max_evaluations = 35

import orekit
orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir
orekit_data_dir = 'orekit-data-master'
setup_orekit_curdir(orekit_data_dir)

from org.orekit.utils import Constants as orekit_constants

from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions
tod = FramesFactory.getTOD(IERSConventions.IERS_2010, False) # Taking tidal effects into account when interpolating EOP parameters
gcrf = FramesFactory.getGCRF()
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)
# Selecting frames to use for OD
eci = gcrf
ecef = itrf

from org.orekit.models.earth import ReferenceEllipsoid
wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(ecef)
from org.orekit.bodies import CelestialBodyFactory
moon = CelestialBodyFactory.getMoon()
sun = CelestialBodyFactory.getSun()

from org.orekit.time import AbsoluteDate, TimeScalesFactory
utc = TimeScalesFactory.getUTC()
mjd_utc_epoch = AbsoluteDate(1858, 11, 17, 0, 0, 0.0, utc)

stationFile = 'SLRF2020_POS+VEL_2023.10.02.snx'
stationEccFile = 'ecc_xyz.snx'
from org.orekit.files.sinex import SinexLoader, Station
from org.orekit.data import DataSource
stations_map = SinexLoader(DataSource(stationFile)).getStations()
ecc_map = SinexLoader(DataSource(stationEccFile)).getStations()

from org.orekit.estimation.measurements import GroundStation
from org.orekit.frames import TopocentricFrame

import pandas as pd
from orekit.pyhelpers import datetime_to_absolutedate
import numpy as np

station_keys = stations_map.keySet()

n_errors = 0
station_df = pd.DataFrame(columns=['lat_deg', 'lon_deg', 'alt_m', 'GroundStation'])
for key in station_keys:
    station_data = stations_map.get(key)
    ecc_data = ecc_map.get(key)
    if ecc_data.getEccRefSystem() != Station.ReferenceSystem.XYZ:
        print('Error, eccentricity coordinate system not XYZ')

    epoch_velocity = station_data.getEpoch()
    durationSinceEpoch = datetime_to_absolutedate(odDate).durationFrom(epoch_velocity)  # seconds

    # Computing current station position using velocity data
    station_pos_at_epoch = station_data.getPosition()
    vel = station_data.getVelocity()  # m/s
    station_pos_current = station_pos_at_epoch.add(vel.scalarMultiply(durationSinceEpoch))

    # Adding eccentricity
    try:
        station_pos_current = station_pos_current.add(ecc_data.getEccentricities(datetime_to_absolutedate(odDate)))
        # Converting to ground station object
        geodeticPoint = wgs84Ellipsoid.transform(station_pos_current, itrf, datetime_to_absolutedate(odDate))
        lon_deg = np.rad2deg(geodeticPoint.getLongitude())
        lat_deg = np.rad2deg(geodeticPoint.getLatitude())
        alt_m = geodeticPoint.getAltitude()
        topocentricFrame = TopocentricFrame(wgs84Ellipsoid, geodeticPoint, key)
        groundStation = GroundStation(topocentricFrame)
        station_df.loc[key] = [lat_deg, lon_deg, alt_m, groundStation]
    except:
        # And exception is thrown when the odDate is not in the date range of the eccentricity entry for this station
        # This is simply for stations which do not exist anymore at odDate
        n_errors += 1

station_df = station_df.sort_index()
print(station_df)

# rawTle = st.tle(norad_cat_id=sat_list[sc_name]['norad_id'], epoch='<{}'.format(odDate), orderby='epoch desc', limit=1, format='tle')
# tleLine1 = rawTle.split('\n')[0]
# tleLine2 = rawTle.split('\n')[1]
#tleLine1 ="1 22195U 92070B   23153.01560229 -.00000009  00000+0  00000+0 0  9995"  # JMF dummy date 2023 to stop the simulation early
#tleLine2 ="2 22195  52.6797  18.2088 0137585 220.6419 290.8914  6.47293385747319"  #     -> orbit model from odDate to tle date
tleLine1 ="1 22195U 92070B   23166.70909887 -.00000009  00000-0  00000-0 0  9993"
tleLine2 ="2 22195  52.6449 240.2410 0138828  67.4256 295.6508  6.47294014724379"
#JMF from https://www.satcat.com/sats/22195 (15 Jun 2023)

print(tleLine1)
print(tleLine2)

from org.orekit.propagation.analytical.tle import TLE
orekitTle = TLE(tleLine1, tleLine2)

from org.orekit.attitudes import NadirPointing
nadirPointing = NadirPointing(eci, wgs84Ellipsoid)

from org.orekit.propagation.analytical.tle import SGP4
sgp4Propagator = SGP4(orekitTle, nadirPointing, sat_list[sc_name]['mass'])

tleInitialState = sgp4Propagator.getInitialState()
tleEpoch = tleInitialState.getDate()
tleOrbit_TEME = tleInitialState.getOrbit()
tlePV_ECI = tleOrbit_TEME.getPVCoordinates(eci)

from org.orekit.orbits import CartesianOrbit
tleOrbit_ECI = CartesianOrbit(tlePV_ECI, eci, wgs84Ellipsoid.getGM())

from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
integratorBuilder = DormandPrince853IntegratorBuilder(prop_min_step, prop_max_step, prop_position_error)

from org.orekit.propagation.conversion import NumericalPropagatorBuilder
from org.orekit.orbits import PositionAngleType
propagatorBuilder = NumericalPropagatorBuilder(tleOrbit_ECI,
                                               integratorBuilder, PositionAngleType.MEAN, estimator_position_scale)
propagatorBuilder.setMass(sat_list[sc_name]['mass'])
propagatorBuilder.setAttitudeProvider(nadirPointing)

# Earth gravity field with degree 64 and order 64
from org.orekit.forces.gravity.potential import GravityFieldFactory
gravityProvider = GravityFieldFactory.getNormalizedProvider(64, 64)
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
gravityAttractionModel = HolmesFeatherstoneAttractionModel(ecef, gravityProvider)
propagatorBuilder.addForceModel(gravityAttractionModel)

# Moon and Sun perturbations
from org.orekit.forces.gravity import ThirdBodyAttraction
moon_3dbodyattraction = ThirdBodyAttraction(moon)
propagatorBuilder.addForceModel(moon_3dbodyattraction)
sun_3dbodyattraction = ThirdBodyAttraction(sun)
propagatorBuilder.addForceModel(sun_3dbodyattraction)

# Solar radiation pressure
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(sat_list[sc_name]['cross_section'], sat_list[sc_name]['cr']);
from org.orekit.forces.radiation import SolarRadiationPressure
solarRadiationPressure = SolarRadiationPressure(sun, wgs84Ellipsoid,
                                                isotropicRadiationSingleCoeff)
propagatorBuilder.addForceModel(solarRadiationPressure)

# Relativity
from org.orekit.forces.gravity import Relativity
relativity = Relativity(orekit_constants.EIGEN5C_EARTH_MU)
propagatorBuilder.addForceModel(relativity)

# Atmospheric drag
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
cswl = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

from org.orekit.models.earth.atmosphere import NRLMSISE00
atmosphere = NRLMSISE00(cswl, sun, wgs84Ellipsoid)
#from org.orekit.forces.drag.atmosphere import DTM2000
#atmosphere = DTM2000(msafe, sun, wgs84Ellipsoid)
from org.orekit.forces.drag import IsotropicDrag
isotropicDrag = IsotropicDrag(sat_list[sc_name]['cross_section'], sat_list[sc_name]['cd'])
from org.orekit.forces.drag import DragForce
dragForce = DragForce(atmosphere, isotropicDrag)
propagatorBuilder.addForceModel(dragForce)

from org.hipparchus.linear import QRDecomposer
matrixDecomposer = QRDecomposer(1e-11)
from org.hipparchus.optim.nonlinear.vector.leastsquares import GaussNewtonOptimizer
optimizer = GaussNewtonOptimizer(matrixDecomposer, False)

from org.orekit.estimation.leastsquares import BatchLSEstimator
estimator = BatchLSEstimator(optimizer, propagatorBuilder)
estimator.setParametersConvergenceThreshold(estimator_convergence_thres)
estimator.setMaxIterations(estimator_max_iterations)
estimator.setMaxEvaluations(estimator_max_evaluations)

#############

## Looking for laser ranging data prior to the OD date.
from slrDataUtils import SlrDlManager
slr_dl_manager = SlrDlManager(username_edc="jmfriedt", password_edc="jmfriedt")
url = 'https://edc.dgfi.tum.de/api/v1/'

 
laserDatasetList = slr_dl_manager.querySlrData(laser_data_type,
                                              sat_list[sc_name]['cospar_id'],
                                              startCollectionDate, odDate)
print(laserDatasetList)
slrDataFrame = slr_dl_manager.dlAndParseSlrData(laser_data_type, laserDatasetList, 'slr-data')
print(slrDataFrame)

# JMF https://edc.dgfi.tum.de/pub/slr/data/npt_crd/lageos2/2024
# 
from orekit.pyhelpers import datetime_to_absolutedate
from org.orekit.estimation.measurements import Range, AngularAzEl, ObservableSatellite
from org.orekit.models.earth.troposphere import MendesPavlisModel
from org.orekit.estimation.measurements.modifiers import RangeTroposphericDelayModifier

observableSatellite = ObservableSatellite(0) # Propagator index = 0

for receiveTime, slrData in slrDataFrame.iterrows():
    if slrData['station-id'] in station_df.index: # Checking if station exists in the stations list, because it might not be up-to-date
        if not np.isnan(slrData['range']):  # If this data point contains a valid range measurement
            orekitRange = Range(station_df.loc[slrData['station-id'], 'GroundStation'],
                                True, # Two-way measurement
                                receiveTime,
                                slrData['range'],
                                range_sigma,
                                range_weight,
                                observableSatellite
                               ) # Uses date of signal reception; https://www.orekit.org/static/apidocs/org/orekit/estimation/measurements/Range.html

            range_tropo_delay_modifier = RangeTroposphericDelayModifier(
                MendesPavlisModel(slrData['temperature_K'],
                                  slrData['pressure_mbar'],
                                  slrData['humidity'],
                                  slrData['wavelength_microm'])
            )
            orekitRange.addModifier(range_tropo_delay_modifier)

            estimator.addMeasurement(orekitRange)
        #if not np.isnan(slrData['az']):  # If this data point contains a valid angles measurement
        #    orekitAzEl = AngularAzEl(station_df.loc[slrData['station-id'], 'GroundStation'],
        #                            receiveTime,
        #                            JArray('double')([slrData['az'], slrData['el']]),
        #                            JArray('double')([az_sigma, el_sigma]),
        #                            JArray('double')([az_weight, el_weight]),
        #                            observableSatellite)
        #    estimator.addMeasurement(orekitAzEl)

estimatedPropagatorArray = estimator.estimate()

dt = 300.0
date_start = datetime_to_absolutedate(startCollectionDate)
date_start = date_start.shiftedBy(-86400.0)
# date_end = datetime_to_absolutedate(odDate) JMF 
# date_end = date_end.shiftedBy(86400.0) # Stopping 1 day after OD date JMF
date_end = date_start.shiftedBy(4*86400.0) # Stopping 1 day after OD date JMF : 4 days

# First propagating in ephemeris mode
estimatedPropagator = estimatedPropagatorArray[0]
estimatedInitialState = estimatedPropagator.getInitialState()
actualOdDate = estimatedInitialState.getDate()
estimatedPropagator.resetInitialState(estimatedInitialState)
eph_generator = estimatedPropagator.getEphemerisGenerator()

# Propagating from 1 day before data collection
# To 1 week after orbit determination (for CPF generation)
estimatedPropagator.propagate(date_start, datetime_to_absolutedate(odDate).shiftedBy(7 * 86400.0))
bounded_propagator = eph_generator.getGeneratedEphemeris()

# Creating the LVLH frame
# It must be associated to the bounded propagator, not the original numerical propagator
from org.orekit.frames import LocalOrbitalFrame
from org.orekit.frames import LOFType
lvlh = LocalOrbitalFrame(eci, LOFType.LVLH, bounded_propagator, 'LVLH')

# Getting covariance matrix in ECI frame
covMat_eci_java = estimator.getPhysicalCovariances(1.0e-10)

# JMF ERROR
#    covMat_eci_java = estimator.getPhysicalCovariances(1.0e-10)
#                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#orekit.JavaError: <super: <class 'JavaError'>, <JavaError object>>
#    Java stacktrace:
#org.orekit.errors.OrekitException: matrix is singular
# JMF : missing observables

# Converting matrix to LVLH frame
# Getting an inertial frame aligned with the LVLH frame at this instant
# The LVLH is normally not inertial, but this should not affect results too much
# Reference: David Vallado, Covariance Transformations for Satellite Flight Dynamics Operations, 2003
eci2lvlh_frozen = eci.getTransformTo(lvlh, actualOdDate).freeze()

# Computing Jacobian
from org.orekit.utils import CartesianDerivativesFilter
from orekit.pyhelpers import JArray_double2D
jacobianDoubleArray = JArray_double2D(6, 6)
eci2lvlh_frozen.getJacobian(CartesianDerivativesFilter.USE_PV, jacobianDoubleArray)
from org.hipparchus.linear import Array2DRowRealMatrix
jacobian = Array2DRowRealMatrix(jacobianDoubleArray)
# Applying Jacobian to convert matrix to lvlh
covMat_lvlh_java = jacobian.multiply(
    covMat_eci_java.multiply(jacobian.transpose()))

# Converting the Java matrices to numpy
import numpy as np
covarianceMat_eci = np.matrix([covMat_eci_java.getRow(iRow)
                              for iRow in range(0, covMat_eci_java.getRowDimension())])
covarianceMat_lvlh = np.matrix([covMat_lvlh_java.getRow(iRow)
                              for iRow in range(0, covMat_lvlh_java.getRowDimension())])

pos_std_crossTrack = np.sqrt(covarianceMat_lvlh[0,0])
pos_std_alongTrack = np.sqrt(covarianceMat_lvlh[1,1])
pos_std_outOfPlane = np.sqrt(covarianceMat_lvlh[2,2])
print(f'Position std: cross-track {pos_std_crossTrack:.3e} m, along-track {pos_std_alongTrack:.3e} m, out-of-plane {pos_std_outOfPlane:.3e} m')

vel_std_crossTrack = np.sqrt(covarianceMat_lvlh[3,3])
vel_std_alongTrack = np.sqrt(covarianceMat_lvlh[4,4])
vel_std_outOfPlane = np.sqrt(covarianceMat_lvlh[5,5])
print(f'Velocity std: cross-track {vel_std_crossTrack:.3e} m/s, along-track {vel_std_alongTrack:.3e} m/s, out-of-plane {vel_std_outOfPlane:.3e} m/s')

sat_properties = {
     'mass': sat_list[sc_name]['mass'],
     'solar_rad_area': sat_list[sc_name]['cross_section'],
     'solar_rad_coeff': sat_list[sc_name]['cd'],
     'drag_area': sat_list[sc_name]['cross_section'],
     'drag_coeff': sat_list[sc_name]['cr']
}

from ccsdsUtils import Ccsds
ccsds_writer = Ccsds(originator='GOR', object_name=sc_name, object_id=sat_list[sc_name]['norad_id'], sat_properties=sat_properties)

pv_eci_init = estimatedInitialState.getPVCoordinates()
pos_eci_init = np.array(pv_eci_init.getPosition().toArray())
vel_eci_init = np.array(pv_eci_init.getVelocity().toArray())

from orekit.pyhelpers import absolutedate_to_datetime

ccsds_writer.write_opm('OPM.txt', absolutedate_to_datetime(actualOdDate), pos_eci_init, vel_eci_init, covarianceMat_eci, 'EARTH', 'GCRF')

propagatorParameters   = estimator.getPropagatorParametersDrivers(True)
measurementsParameters = estimator.getMeasurementsParametersDrivers(True)

lastEstimations = estimator.getLastEstimations()
valueSet = lastEstimations.values()
estimatedMeasurements = valueSet.toArray()
keySet = lastEstimations.keySet()
realMeasurements = keySet.toArray()

from org.orekit.estimation.measurements import EstimatedMeasurement

import pandas as pd
range_residuals = pd.DataFrame(columns=['range'])
azel_residuals = pd.DataFrame(columns=['az', 'el'])

for estMeas, realMeas in zip(estimatedMeasurements, realMeasurements):
    estMeas = EstimatedMeasurement.cast_(estMeas)
    estimatedValue = estMeas.getEstimatedValue()
    pyDateTime = absolutedate_to_datetime(estMeas.date)

    if Range.instance_(realMeas):
        observedValue = Range.cast_(realMeas).getObservedValue()
        range_residuals.loc[pyDateTime] = np.array(observedValue) - np.array(estimatedValue)
    elif AngularAzEl.instance_(realMeas):
        observedValue = AngularAzEl.cast_(realMeas).getObservedValue()
        azel_residuals.loc[pyDateTime] = np.array(observedValue) - np.array(estimatedValue)

print(range_residuals)
print(azel_residuals)

import plotly.graph_objs as go

trace = go.Scattergl(
    x=range_residuals.index, y=range_residuals['range'],
    mode='markers',
    name='Range'
)

data = [trace]

layout = go.Layout(
    title = 'Range residuals',
    xaxis = dict(
        title = 'Datetime UTC'
    ),
    yaxis = dict(
        title = 'Range residual (m)'
    )
)

fig = dict(data=data, layout=layout)

pio.show(fig)

import plotly.graph_objs as go

trace_az = go.Scattergl(
    x=azel_residuals.index, y=np.rad2deg(azel_residuals['az']),
    mode='markers',
    name='Azimuth'
)

trace_el = go.Scattergl(
    x=azel_residuals.index, y=np.rad2deg(azel_residuals['el']),
    mode='markers',
    name='Elevation'
)

data = [trace_az, trace_el]

layout = go.Layout(
    title = 'Angle residuals',
    xaxis = dict(
        title = 'Datetime UTC'
    ),
    yaxis = dict(
        title = 'Angle residual (deg)'
    )
)

fig = dict(data=data, layout=layout)

pio.show(fig)

# Propagating the bounded propagator to retrieve the intermediate states

deltaPV_tle_lvlh_dict = {}
deltaPV_cpf_lvlh_dict = {}

from org.hipparchus.geometry.euclidean.threed import Vector3D

date_current = date_start
while date_current.compareTo(date_end) <= 0:
    datetime_current = absolutedate_to_datetime(date_current)
    spacecraftState = bounded_propagator.propagate(date_current)

    '''
    When getting PV coordinates using the SGP4 propagator in LVLH frame,
    it is actually a "delta" from the PV coordinates resulting from the orbit determination
    because this LVLH frame is centered on the satellite's current position based on the orbit determination
    '''
    deltaPV_lvlh = sgp4Propagator.getPVCoordinates(date_current, lvlh)
    deltaPV_tle_lvlh_dict[datetime_current] = [deltaPV_lvlh.getPosition().getX(),
                                                 deltaPV_lvlh.getPosition().getY(),
                                                 deltaPV_lvlh.getPosition().getZ(),
                                                 deltaPV_lvlh.getPosition().getNorm(),
                                                 deltaPV_lvlh.getVelocity().getX(),
                                                 deltaPV_lvlh.getVelocity().getY(),
                                                 deltaPV_lvlh.getVelocity().getZ(),
                                                 deltaPV_lvlh.getVelocity().getNorm()]

    pos_cpf_ecef = cpfDataFrame.loc[datetime_current]
    ecef2lvlh = ecef.getStaticTransformTo(lvlh, date_current)
    delta_pos_cpf_lvlh_vector = ecef2lvlh.transformPosition(Vector3D(float(pos_cpf_ecef['x']),
                                                                     float(pos_cpf_ecef['y']),
                                                                     float(pos_cpf_ecef['z'])))
    deltaPV_cpf_lvlh_dict[datetime_current] = [delta_pos_cpf_lvlh_vector.getX(),
                                                 delta_pos_cpf_lvlh_vector.getY(),
                                                 delta_pos_cpf_lvlh_vector.getZ(),
                                                 delta_pos_cpf_lvlh_vector.getNorm()]

    date_current = date_current.shiftedBy(dt)

deltaPV_tle_lvlh_df = pd.DataFrame.from_dict(deltaPV_tle_lvlh_dict,
                                             columns=['x', 'y', 'z', 'pos_norm', 'vx', 'vy', 'vz', 'vel_norm'],
                                             orient='index')

deltaPV_cpf_lvlh_df = pd.DataFrame.from_dict(deltaPV_cpf_lvlh_dict,
                                             columns=['x', 'y', 'z', 'norm'],
                                             orient='index')


#output:
#
#$ python3 02-orbit-determination-example.py 
#        lat_deg     lon_deg        alt_m                                      GroundStation
#1824  50.363134   30.495889   211.735985  org.orekit.estimation.measurements.GroundStati...
#1831  49.917575   23.954423   359.744840  org.orekit.estimation.measurements.GroundStati...
#1863  38.685740   66.943098  2713.972551  org.orekit.estimation.measurements.GroundStati...
#1864  38.684896   66.943095  2713.612400  org.orekit.estimation.measurements.GroundStati...
#1868  50.694601  136.743837   269.556377  org.orekit.estimation.measurements.GroundStati...
#...         ...         ...          ...                                                ...
#7841  52.383015   13.061439   127.338569  org.orekit.estimation.measurements.GroundStati...
#7845  43.754636    6.921579  1323.338301  org.orekit.estimation.measurements.GroundStati...
#7865  38.499195  -77.371110    32.021791  org.orekit.estimation.measurements.GroundStati...
#7941  40.648675   16.704617   536.982724  org.orekit.estimation.measurements.GroundStati...
#8834  49.144421   12.878015   665.382135  org.orekit.estimation.measurements.GroundStati...
#
#[67 rows x 4 columns]
#1 22195U 92070B   23166.70909887 -.00000009  00000-0  00000-0 0  9993
#2 22195  52.6449 240.2410 0138828  67.4256 295.6508  6.47294014724379
#               end_data_date errors        incoming_date                           incoming_filename  ...   station status version wavelength
#2845295  2023-06-06 04:17:29         2023-11-07 10:46:54         7840_lageos2_crd_20230606_04_01.npt  ...  78403501  valid      01    532.080
#2845296  2023-06-06 08:20:07         2023-11-07 10:46:56         7840_lageos2_crd_20230606_07_01.npt  ...  78403501  valid      01    532.080
#2828898  2023-06-06 14:44:12         2023-06-06 15:06:33         7821_lageos2_crd_20230606_14_00.npt  ...  78212801  valid      00    532.000
#2828974  2023-06-06 14:47:15         2023-06-06 20:43:51         7396_lageos2_crd_20230606_14_00.npt  ...  73964701  valid      00    532.000
#2828942  2023-06-06 18:48:32         2023-06-06 19:22:58         7821_lageos2_crd_20230606_18_00.npt  ...  78212801  valid      00    532.000
#2828975  2023-06-06 18:39:06         2023-06-06 20:45:33         7396_lageos2_crd_20230606_18_00.npt  ...  73964701  valid      00    532.000
#2829159  2023-06-06 18:31:32         2023-06-07 20:39:10       1890_lageos2_crd_20230606_1819_00.npt  ...  18900901  valid      00    532.000
#2829034  2023-06-07 01:53:17         2023-06-07 03:34:49  8834_lageos2_crd_20230607_0148_00.det7.npt  ...  88341001  valid      00   1064.000
#2829035  2023-06-07 02:04:59         2023-06-07 03:34:48  8834_lageos2_crd_20230607_0202_00.det7.npt  ...  88341001  valid      00   1064.000
#2829036  2023-06-07 02:25:46         2023-06-07 03:34:46  8834_lageos2_crd_20230607_0223_00.det7.npt  ...  88341001  valid      00   1064.000
#2829074  2023-06-07 05:56:01         2023-06-07 09:33:55  8834_lageos2_crd_20230607_0551_00.det7.npt  ...  88341001  valid      00   1064.000
#2829075  2023-06-07 06:24:17         2023-06-07 09:33:51  8834_lageos2_crd_20230607_0621_00.det7.npt  ...  88341001  valid      00   1064.000
#2845297  2023-06-07 10:24:41         2023-11-07 10:46:56         7840_lageos2_crd_20230607_09_01.npt  ...  78403501  valid      01    532.080
#2829231  2023-06-07 20:27:26         2023-06-08 05:38:04       1889_lageos2_crd_20230607_2009_00.npt  ...  18899901  valid      00    532.000
#
#[14 rows x 11 columns]
#                                    station-id         range  wavelength_microm  pressure_mbar  temperature_K  humidity
#2023-06-06T04:15:55.84851870000057Z       7840  7.313215e+06            0.53208    1018.901004     281.600000  0.960000
#2023-06-06T04:16:44.16199469999992Z       7840  7.429026e+06            0.53208    1018.893016     281.600000  0.960000
#2023-06-06T07:45:41.67338469999959Z       7840  6.632704e+06            0.53208    1018.680940     287.461881  0.733812
#2023-06-06T07:47:01.26051970000117Z       7840  6.484162e+06            0.53208    1018.691906     287.483813  0.731619
#2023-06-06T07:49:01.54969469999924Z       7840  6.280238e+06            0.53208    1018.683039     287.500000  0.726608
#...                                        ...           ...                ...            ...            ...       ...
#2023-06-07T10:24:26.35970310000266Z       7840  7.629123e+06            0.53208    1017.422583     289.758058  0.730000
#2023-06-07T20:10:14.10909815999912Z       1889  7.303253e+06            0.53200     892.980000     283.850000  0.970000
#2023-06-07T20:19:28.31912308000028Z       1889  6.924117e+06            0.53200     892.980000     283.850000  0.970000
#2023-06-07T20:20:57.68387671999517Z       1889  6.921368e+06            0.53200     892.980000     283.850000  0.970000
#2023-06-07T20:27:09.33502076999866Z       1889  7.090397e+06            0.53200     892.980000     283.850000  0.970000
#
#[76 rows x 6 columns]

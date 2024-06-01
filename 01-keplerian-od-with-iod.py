## Parameters
sigma_position = 100e3  # Noise (in terms of standard deviation of gaussian distribution) of input position data in meters
sigma_velocity = 100.0  # Noise of input velocity data in meters per second

# Estimator parameters
estimator_position_scale = 1.0 # m
estimator_convergence_thres = 1e-2
estimator_max_iterations = 25
estimator_max_evaluations = 35

# Orbit propagator parameters
prop_min_step = 0.001 # s
prop_max_step = 300.0 # s
prop_position_error = 10.0 # m

## Importing generated position/velocity data
import pandas as pd
points_df = pd.read_csv('pos_vel_data_gcrf_with_noise.csv', index_col=0, parse_dates=True)
points_true_df = pd.read_csv('pos_vel_data_gcrf_without_noise.csv', index_col=0, parse_dates=True)
print(points_df)

## Firing up a JVM for Orekit
import orekit
orekit.initVM()

## Downloading and importing the Orekit data ZIP
from orekit.pyhelpers import download_orekit_data_curdir, setup_orekit_curdir
# download_orekit_data_curdir()
setup_orekit_curdir()

from org.orekit.frames import FramesFactory
gcrf = FramesFactory.getGCRF()

from org.orekit.time import TimeScalesFactory
utc = TimeScalesFactory.getUTC()

## Picking 3 vectors for orbit determination. These vectors are selected evenly-spaced in the dataframe
import math
i = math.ceil(len(points_df.index) / 3)
points_for_iod = points_df.iloc[::i, :]
print(points_for_iod)

from org.hipparchus.geometry.euclidean.threed import Vector3D
from orekit.pyhelpers import datetime_to_absolutedate
pos_1 = points_for_iod.iloc[0]
vector_1 = Vector3D(pos_1[['x', 'y', 'z']].to_list())
date_1 = datetime_to_absolutedate(points_for_iod.index[0])

pos_2 = points_for_iod.iloc[1]
vector_2 = Vector3D(pos_2[['x', 'y', 'z']].to_list())
date_2 = datetime_to_absolutedate(points_for_iod.index[1])

pos_3 = points_for_iod.iloc[2]
vector_3 = Vector3D(pos_3[['x', 'y', 'z']].to_list())
date_3 = datetime_to_absolutedate(points_for_iod.index[2])

## Performing the Initial Orbit Determination using Gibb's method. It assumes that at least 3 data points are available
from org.orekit.estimation.iod import IodGibbs
from org.orekit.utils import Constants as orekit_constants
iod_gibbs = IodGibbs(orekit_constants.EIGEN5C_EARTH_MU)
orbit_first_guess = iod_gibbs.estimate(gcrf,
                                      vector_1, date_1,
                                      vector_2, date_2,
                                      vector_3, date_3)
from org.orekit.propagation.analytical import KeplerianPropagator
kepler_propagator_iod = KeplerianPropagator(orbit_first_guess)

print(orbit_first_guess)

## Setting up a numerical propagator. It is not possible in Orekit to perform orbit determination with a Keplerian propagator
from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
integratorBuilder = DormandPrince853IntegratorBuilder(prop_min_step, prop_max_step, prop_position_error)

from org.orekit.propagation.conversion import NumericalPropagatorBuilder
from org.orekit.orbits import PositionAngleType
propagatorBuilder = NumericalPropagatorBuilder(orbit_first_guess,
                                               integratorBuilder, PositionAngleType.TRUE, estimator_position_scale)
from org.hipparchus.linear import QRDecomposer
matrix_decomposer = QRDecomposer(1e-11)
from org.hipparchus.optim.nonlinear.vector.leastsquares import GaussNewtonOptimizer
optimizer = GaussNewtonOptimizer(matrix_decomposer, False)

from org.orekit.estimation.leastsquares import BatchLSEstimator
estimator = BatchLSEstimator(optimizer, propagatorBuilder)
estimator.setParametersConvergenceThreshold(estimator_convergence_thres)
estimator.setMaxIterations(estimator_max_iterations)
estimator.setMaxEvaluations(estimator_max_evaluations)

## Feeding position measurements to the estimator
from orekit.pyhelpers import datetime_to_absolutedate
from org.orekit.estimation.measurements import Position, ObservableSatellite

observableSatellite = ObservableSatellite(0) # Propagator index = 0

for timestamp, pv_gcrf in points_df.iterrows():
    orekit_position = Position(
        datetime_to_absolutedate(timestamp),
        Vector3D(pv_gcrf[['x', 'y', 'z']].to_list()),
        sigma_position,
        1.0,  # Base weight
        observableSatellite
    )
    estimator.addMeasurement(orekit_position)

## Performing the orbit determination
estimatedPropagatorArray = estimator.estimate()

dt = 60.0
date_start = datetime_to_absolutedate(points_df.index[0])
date_end = datetime_to_absolutedate(points_df.index[-1])

# First propagating in ephemeris mode
estimatedPropagator = estimatedPropagatorArray[0]
estimatedInitialState = estimatedPropagator.getInitialState()
print(estimatedInitialState.getOrbit())
estimatedPropagator.resetInitialState(estimatedInitialState)

# estimatedPropagator.setEphemerisMode() # JMF obsolete see https://gitlab.orekit.org/orekit-labs/python-wrapper/-/issues/441
generator = estimatedPropagator.getEphemerisGenerator();
estimatedPropagator.propagate(date_start, date_end);
bounded_propagator= generator.getGeneratedEphemeris();
print(bounded_propagator)

## Propagating the bounded propagator to retrieve the intermediate states
# BLS = batch least squares
import numpy as np
from orekit.pyhelpers import absolutedate_to_datetime
pv_bls_df = pd.DataFrame(columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])
pv_iod_df = pd.DataFrame(columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])
    
date_current = date_start
while date_current.compareTo(date_end) <= 0:
    datetime_current = absolutedate_to_datetime(date_current)    
    spacecraftState = bounded_propagator.propagate(date_current)
    
    pv_bls = spacecraftState.getPVCoordinates(gcrf)
    pos_bls = np.array(pv_bls.getPosition().toArray())
    pv_bls_df.loc[datetime_current] = np.concatenate(
        (pos_bls,
         np.array(pv_bls.getVelocity().toArray()))
    )

    pv_iod = kepler_propagator_iod.getPVCoordinates(date_current, gcrf)
    pos_iod_gcrf = np.array(pv_iod.getPosition().toArray())
    pv_iod_df.loc[datetime_current] = np.concatenate(
        (pos_iod_gcrf,
         np.array(pv_iod.getVelocity().toArray()))
    )
    date_current = date_current.shiftedBy(dt) 

## Plotting orbit in 3D
import plotly.graph_objects as go
fig_data = data=[go.Scatter3d(x=pv_iod_df['x'], y=pv_iod_df['y'], z=pv_iod_df['z'],
                              mode='lines',
                              name='IOD solution'),
                 go.Scatter3d(x=pv_bls_df['x'], y=pv_bls_df['y'], z=pv_bls_df['z'],
                              mode='lines',
                              name='Batch least squares solution'),
                 go.Scatter3d(x=points_df['x'], y=points_df['y'], z=points_df['z'],
                             mode='markers',
                             name='Measurements'),
                 go.Scatter3d(x=points_for_iod['x'], y=points_for_iod['y'], z=points_for_iod['z'],
                             mode='markers',
                             name='Measurements used for IOD')]
scene=dict(aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
          )
layout = dict(
    scene=scene
)
fig = go.Figure(data=fig_data,
               layout=layout)
fig.show()

## Computing residuals
residuals_df = pd.DataFrame(columns=['bls_minus_measurements_norm', 'iod_minus_measurements_norm', 'bls_minus_truth_norm', 'iod_minus_truth_norm'])

for timestamp, pv_gcrf in points_df.iterrows():    
    date_current = datetime_to_absolutedate(timestamp)
    pv_bls = bounded_propagator.getPVCoordinates(date_current, gcrf)
    pos_bls = np.array(pv_bls.getPosition().toArray())

    pv_iod = kepler_propagator_iod.getPVCoordinates(date_current, gcrf)
    pos_iod = np.array(pv_iod.getPosition().toArray())
    
    pv_measurements = points_df.loc[timestamp]
    pos_measurements = pv_measurements[['x', 'y', 'z']]
    
    pv_true = points_true_df.loc[timestamp]
    pos_true = pv_true[['x', 'y', 'z']]
    
    bls_minus_measurements = np.linalg.norm(pos_bls - pos_measurements)
    iod_minus_measurements = np.linalg.norm(pos_iod - pos_measurements)
    bls_minus_truth = np.linalg.norm(pos_bls - pos_true)
    iod_minus_truth = np.linalg.norm(pos_iod - pos_true)
    
    residuals_df.loc[timestamp] = [
        np.linalg.norm(pos_bls - pos_measurements),
        np.linalg.norm(pos_iod - pos_measurements),
        np.linalg.norm(pos_bls - pos_true),
        np.linalg.norm(pos_iod - pos_true),
    ]
    
print(residuals_df)

## Showing the residuals, i.e. the distance between the measurement points and the estimated positions (by the IOD and the BLS respectively).
## This tells how well the estimation fits the measurements, but not how well it fits the "true" orbit, because the measurements contain significant noise here

fig = go.Figure(data=[
    go.Scatter(x=residuals_df.index, y=residuals_df['bls_minus_measurements_norm'], 
               name='BLS - measurements', 
               mode='markers+lines'),
    go.Scatter(x=residuals_df.index, y=residuals_df['iod_minus_measurements_norm'], 
               name='IOD - measurements', 
               mode='markers+lines')    
])

fig.show()

## Showing the "estimation error", i.e. the difference between the truth (used to generate the measurement data) and the estimation.
## The batch least squares estimation gives better results than the measurement noise: it is able to "filter out" the noise a little bit as it fits the ellipse

fig = go.Figure(data=[
    go.Scatter(x=residuals_df.index, y=residuals_df['bls_minus_truth_norm'], name='BLS - truth', mode='markers+lines'),
    go.Scatter(x=residuals_df.index, y=residuals_df['iod_minus_truth_norm'], name='IOD - truth', mode='markers+lines'),
])

fig.show()

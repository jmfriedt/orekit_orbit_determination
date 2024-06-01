## Parameters for data generation
import numpy as np
import datetime
orbit_epoch = datetime.datetime(2020, 1, 1)
sma = 10000e3  # Semi-major axis in meters
ecc = 0.2  # Eccentricity
inc = float(np.deg2rad(60.0))  # Inclination in radians
raan = 0.0  # Right-ascension of ascending node in radians
aop = 0.0  # Argument of perigee in radians
ta = 0.0  # Initial true anomaly in radians

sigma_position = 100e3  # Noise (in terms of standard deviation of gaussian distribution) of output position data in meters
sigma_velocity = 100.0  # Noise of output velocity data in meters per second

n_orbits = 1.0  # Number of orbits to generate
T = 500.0  # Sample time of output data in seconds

orbit_period = 2*np.pi*np.sqrt(sma**3 / 3.98600e14)
duration = float(n_orbits*orbit_period)
print(duration)

## Firing up a JVM for Orekit
import orekit
orekit.initVM()

## Downloading and importing the Orekit data ZIP
from orekit.pyhelpers import download_orekit_data_curdir, setup_orekit_curdir
# download_orekit_data_curdir()
setup_orekit_curdir()

## Setting up models (frames, timescales)
from org.orekit.frames import FramesFactory
gcrf = FramesFactory.getGCRF()



from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from orekit.pyhelpers import datetime_to_absolutedate
from org.orekit.utils import Constants as orekit_constants
orbit = KeplerianOrbit(sma, ecc, inc, aop, raan, ta,
                       PositionAngleType.TRUE, 
                       gcrf,
                       datetime_to_absolutedate(orbit_epoch),
                       orekit_constants.EIGEN5C_EARTH_MU)

from org.orekit.propagation.analytical import KeplerianPropagator
propagator = KeplerianPropagator(orbit)

## Orekit can generate measurements for all types of measurements that are used for OD. However for position/velocity, it is easy to retrieve the measurements by hand.
## Propagating, adding noise and saving data to a pandas dataframe
import pandas as pd
position_df = pd.DataFrame(columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])
position_without_noise_df = pd.DataFrame(columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])

from orekit.pyhelpers import absolutedate_to_datetime
date_start = datetime_to_absolutedate(orbit_epoch)
date_end = date_start.shiftedBy(duration)
date_current = date_start

while date_current.compareTo(date_end) < 0:
    pv_gcrf = propagator.getPVCoordinates(date_current, gcrf)
    
    # Adding noise to position and velocity
    pos_without_noise = np.array(pv_gcrf.getPosition().toArray())
    pos_with_noise = pos_without_noise + np.random.normal(0, sigma_position, len(pos_without_noise))
    vel_without_noise = np.array(pv_gcrf.getVelocity().toArray())
    vel_with_noise = vel_without_noise + np.random.normal(0, sigma_velocity, len(vel_without_noise))
    
    position_df.loc[absolutedate_to_datetime(date_current)] = np.concatenate(
        (pos_with_noise, vel_with_noise)
    )
    position_without_noise_df.loc[absolutedate_to_datetime(date_current)] = np.concatenate(
        (pos_without_noise, vel_without_noise)
    )
    date_current = date_current.shiftedBy(T)
    
print(position_df)

## Finally, saving position data to CSV file
position_df.to_csv('pos_vel_data_gcrf_with_noise.csv')
position_without_noise_df.to_csv('pos_vel_data_gcrf_without_noise.csv')

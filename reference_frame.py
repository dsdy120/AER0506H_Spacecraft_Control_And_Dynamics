'''
Handle reference-frame changes between different coordinate systems.
'''
import numpy as np

OBLATENESS_EARTH = 0.0033528
EARTH_RADIUS = 6378.0


def observation_site_ECI_vector(latitude_deg,datum_elevation_km,loc_sidereal_deg):
    '''
    Returns position vector of observation site in ECI frame.
    '''
    if latitude_deg < -90.0 or latitude_deg > 90.0:
        raise ValueError("Latitude out of range")
    
    loc_sidereal_deg = loc_sidereal_deg % 360.0
    
    latitude_rad = np.radians(latitude_deg)
    loc_sidereal_rad = np.radians(loc_sidereal_deg)

    mag_R_phi = EARTH_RADIUS / np.sqrt(1 - (2*OBLATENESS_EARTH - OBLATENESS_EARTH**2) * np.sin(latitude_rad)**2)

    I_component = (mag_R_phi + datum_elevation_km) \
        * np.cos(latitude_rad) \
            * (np.cos(loc_sidereal_rad))
    J_component = (mag_R_phi + datum_elevation_km) \
        * np.cos(latitude_rad) \
            * (np.sin(loc_sidereal_rad))
    K_component = (mag_R_phi * (1 - OBLATENESS_EARTH)**2 + datum_elevation_km) \
        * np.sin(latitude_rad)

    return np.array(
        [I_component, J_component, K_component]
    )

def target_track_ENZ_vector(azimuth_deg,elevation_deg,range_km):
    '''
    Returns position vector of target track in TEN frame.
    '''
    if elevation_deg < -90.0 or elevation_deg > 90.0:
        raise ValueError("Elevation out of range")
    azimuth_deg = azimuth_deg % 360.0

    azimuth_rad = np.radians(azimuth_deg)
    elevation_rad = np.radians(elevation_deg)

    I_component = range_km * np.cos(elevation_rad) * np.sin(azimuth_rad)
    J_component = range_km * np.cos(elevation_rad) * np.cos(azimuth_rad)
    K_component = range_km * np.sin(elevation_rad)

    return np.array(
        [I_component, J_component, K_component]
    )

def ECI_to_ENZ_matrix(latitude_deg,loc_sidereal_deg):
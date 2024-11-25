import reference_frame
import orbit_determination

LATITUDE_DEG = -20.0
DATUM_ELEVATION_KM = 0.5

r1 = reference_frame.observation_to_ECI_position_vector(
    LATITUDE_DEG,DATUM_ELEVATION_KM,
    60.0,
    165.931,
    9.53549,
    1214.89
)

r2 = reference_frame.observation_to_ECI_position_vector(
    LATITUDE_DEG,DATUM_ELEVATION_KM,
    60.5014,
    145.967,
    45.7711,
    421.441
)

r3 = reference_frame.observation_to_ECI_position_vector(
    LATITUDE_DEG,DATUM_ELEVATION_KM,
    61.0027,
    2.40962,
    21.8825,
    732.079
)

orbital_elements = orbit_determination.gibbs_solver(r1,r2,r3)

print(orbital_elements)
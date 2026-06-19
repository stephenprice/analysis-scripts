from polaris_rio import polar_rio

# case 1: 
# mean thickness = 0.5 m
# concentration = 0.5
# actual thickness in areas where there are ice = 0.5 m / 0.5 = 1.0 m
# thickness used to calculate RIV = 1.0 m

# case 2: 
# mean thickness = 0.5 m
# concentration = 0.5
# thickness used to calculate RIV = 0.5 m

# Which one is correct in the context of sea ice daily stats from mpas-si?
# They give quite diff. RIO values (For NIS: case 1 RIO=-10; case 2 RIO=0)

rio, ice_type, level = polar_rio(
    vessel_ice_class="NIS",
    ice_thickness_m=0.5,
    concentration_tenths=5
)

print(f"RIO = {rio}")
print(f"Ice type used = {ice_type}")
print(f"Operational guidance = {level}")

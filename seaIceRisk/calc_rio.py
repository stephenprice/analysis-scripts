from polaris_rio import polar_rio

rio, ice_type, level = polar_rio(
    vessel_ice_class="PC3",
    ice_thickness_m=0.8,
    concentration_tenths=10
)

print(f"RIO = {rio}")
print(f"Ice type used = {ice_type}")
print(f"Operational guidance = {level}")

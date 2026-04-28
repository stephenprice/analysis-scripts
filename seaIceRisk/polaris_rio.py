# polaris_rio.py
from typing import Tuple, Optional

ICE_TYPE_KEYS = [
    "ice_free", "new_ice", "grey_ice", "grey_white_ice",
    "thin_first_year_1st_stage", "thin_first_year_2nd_stage", "thin_first_year",
    "medium_first_year_lt_1m", "medium_first_year", "thick_first_year",
    "second_year_ice", "light_multi_year_lt2p5", "heavy_multi_year"
]

# ## orginal values - ChatGPT-dervied here:
# RIV_TABLE = {
#     "PC1":  [3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1],
#     "PC2":  [3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 0],
#     "PC3":  [3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 0, -1, -1],
#     "PC4":  [3, 3, 3, 3, 2, 2, 2, 2, 1, 0, -1, -2, -2],
#     "PC5":  [3, 3, 3, 3, 2, 2, 1, 1, 0, -1, -2, -2, -2],
#     "PC6":  [3, 2, 2, 2, 2, 1, 1, 0, -1, -2, -3, -3, -3],
#     "PC7":  [3, 2, 2, 2, 1, 1, 0, -1, -2, -3, -3, -3, -3],
#     "IAsuper": [3, 2, 2, 2, 2, 1, 0, -1, -2, -3, -4, -4, -4],
#     "IA":   [3, 2, 2, 2, 1, 0, -1, -2, -3, -4, -5, -5, -5],
#     "IB":   [3, 2, 2, 1, 0, -1, -2, -3, -4, -5, -6, -6, -6],
#     "IC":   [3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -8],
#     "NIS": [3, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -8, -8],  # NIS = 'Not Ice Strengthened'
}

# Modified table here - relative to above, drops rightmost entry
# PC* values correspond to entries in Appendix D of our 2025 report.
RIV_TABLE = {
    "PC1":  [3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1 ],
    "PC2":  [3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 0 ],
    "PC3":  [3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 0, -1 ],
    "PC4":  [3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 0, -1, -2 ],
    "PC5":  [3, 3, 3, 3, 3, 2, 2, 1, 1, 0, -1, -2, -2 ],
    "PC6":  [3, 3, 2, 2, 2, 2, 1, 1, 0, -1, -2, -3, -3 ],
    "PC7":  [3, 3, 2, 2, 2, 1, 1, 0, -1, -2, -3, -3, -3 ],
    "IAsuper": [3, 3, 2, 2, 2, 2, 1, 0, -1, -2, -3, -4, -4 ],
    "IA":   [3, 3, 2, 2, 2, 1, 0, -1, -2, -3, -4, -5, -5 ],
    "IB":   [3, 3, 2, 2, 1, 0, -1, -2, -3, -4, -5, -6, -6 ],
    "IC":   [3, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8 ],
    "NIS": [3, 3, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -8 ],  # NIS = 'Not Ice Strengthened'
}

def _approx_thickness_to_ice_type(thickness_m: float) -> str:

# ## orginal values - ChatGPT-dervied based on values in tables here:
# ## https://www.nautinst.org/static/uploaded/2f01665c-04f7-4488-802552e5b5db62d9.pdf
#     if thickness_m <= 0.0:
#         return "ice_free"
#     if thickness_m <= 0.05:
#         return "new_ice"
#     if thickness_m <= 0.30:
#         return "grey_ice"
#     if thickness_m <= 0.50:
#         return "grey_white_ice"
#     if thickness_m <= 0.70:
#         return "thin_first_year_1st_stage"
#     if thickness_m <= 1.00:
#         return "thin_first_year_2nd_stage"
#     if thickness_m <= 1.25:
#         return "thin_first_year"
#     if thickness_m <= 1.75:
#         return "medium_first_year_lt_1m"
#     if thickness_m <= 2.00:
#         return "medium_first_year"
#     if thickness_m <= 2.25:
#         return "thick_first_year"
#     if thickness_m <= 2.75:
#         return "second_year_ice"
#     if thickness_m <= 2.5:
#         return "light_multi_year_lt2p5"
#     return "heavy_multi_year"

# Altered - based on ranges used in prev. projects AND accounting for
# what looks like a mistake in the above ("<=" should be "<" ??).
# See Table in Appendix D of our 2025 report.
    if thickness_m <= 0.0:
        return "ice_free"
    if thickness_m > 0 and thickness_m < 0.05:
        return "new_ice"
    if thickness_m >= 0.05 and thickness_m < 0.10:
        return "grey_ice"
    if thickness_m >= 0.10 and thickness_m < 0.15:
        return "grey_white_ice"
    if thickness_m >= 0.15 and thickness_m < 0.30:
        return "thin_first_year_1st_stage"
    if thickness_m >= 0.30 and thickness_m < 0.50:
        return "thin_first_year_2nd_stage"
    if thickness_m >= 0.50 and thickness_m < 0.70:
        return "thin_first_year"
    if thickness_m >= 0.70 and thickness_m < 1.00:
        return "medium_first_year_lt_1m"
    if thickness_m >= 1.00 and thickness_m < 1.20:
        return "medium_first_year"
    if thickness_m >= 1.20 and thickness_m < 1.70:
        return "thick_first_year"
    if thickness_m >= 1.70 and thickness_m < 2.0:
        return "second_year_ice"
    if thickness_m >= 2.0 and thickness_m < 2.5:
        return "light_multi_year_lt2p5"
    elif thickness_m >= 2.5:
        return "heavy_multi_year"
    
def polar_rio(
    vessel_ice_class: str,
    ice_thickness_m: Optional[float] = None,
    ice_type: Optional[str] = None,
    concentration_tenths: int = 10
):
    cls = vessel_ice_class.strip()
    if cls not in RIV_TABLE:
        raise ValueError(f"Unknown ice class '{cls}'")

    if ice_type is None:
        if ice_thickness_m is None:
            raise ValueError("Either ice_type or ice_thickness_m must be provided")
        ice_type_used = _approx_thickness_to_ice_type(float(ice_thickness_m))
    else:
        ice_type_used = ice_type
        
    idx = ICE_TYPE_KEYS.index(ice_type_used)
    riv = RIV_TABLE[cls][idx]

    if not 0 <= concentration_tenths <= 10:
        raise ValueError("concentration_tenths must be between 0 and 10")

    # Altered version of originally proposed code to account for RIV associated w/ open water 
    # rio = concentration_tenths * riv  # original code
    #rio = 999 #DEBUG
    rio = concentration_tenths * riv + (10 - concentration_tenths) * 3

    is_pc = cls.startswith("PC") and cls[2:].isdigit() and 1 <= int(cls[2:]) <= 7
    if is_pc:
        if rio >= 0:
            level = "Normal operation"
        elif rio >= -10:
            level = "Elevated operational risk"
        else:
            level = "Operation subject to special consideration"
    else:
        if rio >= 0:
            level = "Normal operation"
        else:
            level = "Operation subject to special consideration"

    return rio, ice_type_used, level


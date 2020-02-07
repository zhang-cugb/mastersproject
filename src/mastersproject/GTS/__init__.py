from GTS.deprecated_files.geological_model_visualization import GeologicalModel
from GTS.deprecated_files.tests import GeologicalModelMethodTests


# Import new model data
from GTS.ISC_data.isc import ISCData, swiss_to_gts

# Transformations
from GTS.ISC_data.isc import borehole_to_global_coords

# Import fracture tools
from GTS.ISC_data.fracture import (
    convex_plane,
    fracture_network
)

# Plane fit tools
from GTS.fit_plane import plane_from_points  # New, improved method
from GTS.fit_plane import fit_normal_to_points  # Temporary for prototype testing

from GTS.fit_plane import convex_hull  # Construct the convex hull

# -------------------------
# --- SETUPS AND MODELS ---
# -------------------------

# --- SETUP AND RUN MODEL---
from GTS.isc_modelling.setup import (
    run_mechanics_model,
    run_biot_model,
)

# --- MODELS ---
# Contact mechanics model
from GTS.isc_modelling.mechanics import (
    ContactMechanicsISC,
    ContactMechanicsIsotropicISC,
)

# Contact Mechanics Biot
from GTS.isc_modelling.contact_mechanics_biot import ContactMechanicsBiotISC

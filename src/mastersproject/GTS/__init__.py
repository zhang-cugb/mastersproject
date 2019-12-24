from GTS.deprecated_files.geological_model_visualization import GeologicalModel
from GTS.deprecated_files.tests import GeologicalModelMethodTests


# Import new model
from GTS.ISC_data.isc import ISCData, swiss_to_gts

# Import fracture tools
from GTS.ISC_data.fracture import convex_plane, fracture_network

# Plane fit tools
from GTS.fit_plane import plane_from_points  # New, improved method

from GTS.fit_plane import convex_hull  # Construct the convex hull


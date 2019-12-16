""" This is a python implementation of the matlab ISC data manipulation script.

Citation:
"Comprehensive geological dataset for a fractured crystalline rock
volume at the Grimsel Test Site" Krietsch et al. 2018
DOI: 10.3929/ethz-b-000243199

link: https://doi.org/10.3929/ethz-b-000243199

"""
from pathlib import Path

import numpy as np
import pandas as pd

class ISCData:

    def __init__(self):
        """ Initialize the class managing data from the ISC project"""
        self.root = Path.cwd() / 'src/mastersproject'  # should be: C:\Users\Haakon\OneDrive\Dokumenter\FORSKNING\mastersproject\src\mastersproject + \src\mastersproject
        self.data_path = self.root / 'GTS/01BasicInputData'  # Path to available data

        # TEMPORARY CHECK FOR LOCAL TESTING:
        assert(Path.cwd() == Path('C:/Users/Haakon/OneDrive/Dokumenter/FORSKNING/mastersproject'))
        assert(self.data_path.is_dir())

        self.gts_coordinates = np.array((667400, 158800, 1700))

        # 1. Step: Define Tunnel =======================================================================================
        # TODO: Consider importing tunnels at a late time.

        # 2. Step: Drill Borehole ======================================================================================
        # Import all boreholes in order to use their coordinates for further calculations.

        # Name of boreholes.
        self.borehole_types = {'FBS': np.array([1, 2, 3]),
                               'SBH': np.array([1, 3, 4]),  # Note the skip of numbering for SBH
                               'INJ': np.array([1, 2]),
                               'PRP': np.array([1, 2, 3]),
                               'GEO': np.array([1, 2, 3, 4])}
        # Name of shearzones
        self.shearzone_types = {'S1': np.array([1, 2, 3]),
                                'S3': np.array([1, 2])}

    # def _import_data(self, folder: str, df_columns, skiprows=0):
    #     """ General method to import ISC data
    #
    #     Parameters:
    #         folder (str): Folder(s) relative to self.data_path (i.e. .../GTS/01BasicInputData/) to fetch data from
    #         df_columns (list): column names
    #         skiprows (int, default 0): Number of rows to skip.
    #
    #     """
    #     file_loc = self.data_path / folder

    def borehole_data(self):
        """ Fetch data with borehole coordinates

        Assumes existence of the following attributes:
            gts_coordinates
            borehole_types

        Returns
        pd.DataFrame: data on location and orientation of each borehole INJ1, FBS1, etc.
            Columns:
        """
        file_loc = self.data_path / "02_Boreholes"
        columns = ['x', 'y', 'z', 'length', 'diameter', 'azimuth', 'upward_gradient']

        data = []
        for parent in self.borehole_types:
            path = file_loc / (parent + '.txt')
            frame = pd.read_csv(path, sep=None, names=columns, engine='python')
            borehole_name = np.array([parent + str(id) for id in self.borehole_types[parent]])
            frame['borehole'] = borehole_name
            data.append(frame)

        df = pd.concat(data, ignore_index=True)
        return df

    def borehole_structure_data(self):
        """ Data on geological structures' intersections to boreholes.

        Returns
        pd.DataFrame: data on depth, orientation, and thickness of intersecting structures to each borehole.

        """
        file_loc = self.data_path / "03_GeologicalMapping" / "02_BoreholeIntersections"
        columns = ['depth', 'azimuth', 'dip', 'aperture', 'type']

        data = []
        for parent in self.borehole_types:  # Loop INJ, SBH, ...
            for p_num in self.borehole_types[parent]:  # Loop INJ1, INJ2, ...
                borehole_name = parent + str(p_num)
                path = file_loc / (borehole_name + '_structures.txt')
                frame = pd.read_csv(path, sep=None, names=columns, skiprows=2, engine='python')
                frame['borehole'] = borehole_name
                data.append(frame)
        df = pd.concat(data, ignore_index=True)
        return df

    def borehole_shearzone_data(self):
        """ Extract shearzones and fracture data from borehole-structure data.

        Use method borehole_structure_data(self), and filter type by shearzones and fractures.

        Returns
        pd.DataFrame: data on depth, orientation, and thickness of intersecting shearzones and
            fractures to each borehole.

        """
        df = self.borehole_structure_data()
        structure_filter = ['Fracture', 'S1 Shear-zone', 'S3 Shear-zone']
        df = df[df.type.isin(structure_filter)]
        return df

    def tunnel_shearzone_data(self):
        """Data on tunnel intersections with shearzones

        """
        file_loc = self.data_path / "03_GeologicalMapping" / "01_TunnelIntersections"
        columns = ['x', 'y', 'z', 'true_dip_direction', 'dip', 'tunnel', 'shearzone']

        path = file_loc / "Tunnel_intersections.txt"
        df = pd.read_csv(path, sep=None, names=columns, engine='python')
        return df

    def shearzone_borehole_data(self):
        """ Import data on shearzone intersections with boreholes

        This data is 'opposite' of the data from borehole_structure_data(self).
        i.e. We extract for every shearzone (S1_1, S1_2, ...) their intersections with boreholes.

        Note: The borehole_shearzone_data contains multiple intersections of shearzones to boreholes,
        whereas the shearzone_borehole_data only contains one (a subset of the former).
        Why this is so, I do not know.
        TODO: Find out why we have single intersections with this data, but not with the other data.
            (see borehole_shearzone_data(self)).

        """
        file_loc = self.data_path / "06_ShearzoneInterpolation"
        columns = ['borehole', 'depth']

        data = []
        for parent in self.shearzone_types:
            for sz_num in self.shearzone_types[parent]:
                sz_name = parent + "_" + str(sz_num)  # e.g. 'S1_1'
                path = file_loc / (sz_name + ".txt")
                frame = pd.read_csv(path, sep=None, names=columns, skiprows=1, engine='python')
                frame['shearzone'] = sz_name
                data.append(frame)
        df = pd.concat(data, ignore_index=True)
        return df


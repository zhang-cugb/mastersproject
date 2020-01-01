""" This is a python implementation of the matlab ISC data manipulation script.

Citation:
"Comprehensive geological dataset for a fractured crystalline rock
volume at the Grimsel Test Site" Krietsch et al. 2018
DOI: 10.3929/ethz-b-000243199

link: https://doi.org/10.3929/ethz-b-000243199

"""
from pathlib import Path
import logging

import numpy as np
import pandas as pd


class ISCData:

    def __init__(self, path=None):
        """ Initialize the class managing data from the ISC project

        Parameters:
            path : str, pathlib.Path
                Path/to/01BasicInputData/
                or 'windows' or 'linux' for default values.
                If None, 'linux' by default.

        """
        if path is None:
            path = 'linux'
        if path == 'linux':
            _root = Path.cwd()  # should be path/to/mastersproject/src/mastersproject
            self.data_path = _root / 'GTS/01BasicInputData'
        elif path == 'windows':
            _root = Path('C:/Users/Haakon/OneDrive/Dokumenter/FORSKNING/mastersproject/src/mastersproject')
            self.data_path = _root / 'GTS/01BasicInputData'
        else:
            self.data_path = Path(path)

        # TEMPORARY CHECK FOR LOCAL TESTING:
        # assert(Path.cwd() == Path('C:/Users/Haakon/OneDrive/Dokumenter/FORSKNING/mastersproject'))
        logging.info(f"Data located at: {self.data_path}.")
        assert (self.data_path.is_dir())

        self.gts_coordinates = np.array((667400, 158800, 1700))

        # TODO: For all assignments below (except 'self.structures').
        #  Find which ones you can remove.

        # Name of boreholes.
        self.borehole_types = {'FBS': np.array([1, 2, 3]),
                               'SBH': np.array([1, 3, 4]),  # Note the skip of numbering for SBH
                               'INJ': np.array([1, 2]),
                               'PRP': np.array([1, 2, 3]),
                               'GEO': np.array([1, 2, 3, 4])}
        # Name of shearzones
        self.shearzone_types = {'S1': np.array([1, 2, 3]),
                                'S3': np.array([1, 2])}
        self.shearzones = [sz_set+'_'+str(sz_num) for sz_set in self.shearzone_types
                           for sz_num in self.shearzone_types[sz_set]]

        # x. Step: Load all available data. ============================================================================
        # Load borehole data
        self.borehole_geometry = self.borehole_data()

        # Load borehole structure data
        self.borehole_structures = self.borehole_structure_data()

        # Filter by shear-zone and fracture structures only
        # self.borehole_frac_sz = self.borehole_shearzone_data()

        # Load tunnel structures (only shear-zones and fractures)
        self.tunnel_structures = self.tunnel_shearzone_data()

        # Load interpolation-ready shear-zone - borehole intersections
        # i.e. 1-1 (-0) mapping between shear-zones and boreholes.
        self.shearzone_borehole_geometry = self.shearzone_borehole_data()

        # y. Step: All characterized structures ========================================================================
        self.structures = self.full_structure_geometry()

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
        df['shearzone'] = df['shearzone'].apply(rename_sz)
        df = df.rename(
            columns={
                'true_dip_direction': 'azimuth_struc',
                'tunnel': 'borehole',
            })
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

    @staticmethod
    def bh_struc_to_global_coords(data: pd.DataFrame, *,
                                  x: str, y: str, z: str, depth: str,
                                  upward_gradient: str, azimuth: str,
                                  ):
        """ Convert coordinates in a borehole to global coordinates

        For all rows in a DataFrame, convert some (x,y,z) coordinates to global
        coordinates, localized to Swiss and/or GTS.

        Parameters:
        data (pd.DataFrame)
        x, y, z, depth, upward_gradient, azimuth (str):
            Column names for the respective quantities
        """

        # Compute angle scalers
        rad = np.pi / 180
        data.loc[:, '_trig_x'] = (data[upward_gradient] * rad).apply(np.cos) * \
                                 (data[azimuth] * rad).apply(np.sin)

        data.loc[:, '_trig_y'] = (data[upward_gradient] * rad).apply(np.cos) * \
                                 (data[azimuth] * rad).apply(np.cos)

        data.loc[:, '_trig_z'] = (data[upward_gradient] * rad).apply(np.sin)

        # Swiss coordinates
        data.loc[:, 'x_swiss'] = data[x] + (data[depth] * data['_trig_x'])
        data.loc[:, 'y_swiss'] = data[y] + (data[depth] * data['_trig_y'])
        data.loc[:, 'z_swiss'] = data[z] + (data[depth] * data['_trig_z'])

        # GTS coordinates
        data[['x_gts', 'y_gts', 'z_gts']] = data[['x_swiss', 'y_swiss', 'z_swiss']] \
            .apply(swiss_to_gts, axis=1, raw=True, result_type='expand')

    def _characterize_shearzones(self):
        """ Classify all structures as specific shear-zones (e.g. S1_1) or none

        Load all borehole structures, and classify each as not a shear-zone,
        or as Sx_y. This method combines structure-borehole data with
        shear-zone -- borehole data.

        Note:
            - S1_2 intersects GEO3 is of type 'Minor ductile Shear-zone'
            - One shearzone is inconsistent (to +- 0.01) across the two datasets.
        These two notes has significantly complicated the code.

        Returns:
        pd.DataFrame: Characterized structures.
        """

        # Import shear-zone -- borehole data.
        sz_bh = self.shearzone_borehole_data()
        sz_bh = sz_bh[sz_bh['depth'].notna()]  # Remove rows with no intersection.

        # Import borehole - structure data.
        strc = self.borehole_structure_data().merge(self.borehole_data(), how='outer', on='borehole',
                                                    suffixes=('_struc', '_bh'), validate='m:1')

        # left merge data structures with simulation-shearzones.
        # Notes:
        # - Absolute tolerance = 0.01
        # - exact match on column = 'borehole'
        _merge = pd.merge_asof(strc.sort_values('depth'),
                               sz_bh.sort_values('depth'),
                               by='borehole',
                               on='depth',
                               tolerance=0.01,
                               direction='nearest')

        # The above merge includes some nearby 'fractures', etc.
        # All shearzones are classified as 'S1 Shear-zone' or 'S3 Shear-zone'
        # except once: S1_2 intersects GEO3 is of type 'Minor ductile Shear-zone'.
        # With this set of 'type' filters, no spurious shearzone intersections are located.
        ### Below is the complete classified shearzone set.
        _mask_nna = _merge.shearzone.notna()
        _mask_szset = _merge.type.isin(['S1 Shear-zone', 'S3 Shear-zone', 'Minor ductile Shear-zone'])
        _shearzones = _merge[_mask_nna & _mask_szset]

        # Step: Get indices indices of non-shearzones:
        _index_other = set(_merge.index) - set(_shearzones.index)

        # Set non-shearzones to nan
        _merge.loc[np.array(list(_index_other)), 'shearzone'] = np.nan

        return _merge

    def full_structure_geometry(self):
        """ Compute geometry of all structures in ISC.

        Geometry of all structures in boreholes and shear-zones in tunnels are located,
        and global coordinates computed.

        Returns:
            pd.DataFrame: Full data set - with all characterized structures in local coordinates.

        """
        # Characterized borehole structures
        borehole_structures = self._characterize_shearzones()

        # Tunnel shearzone data
        tunnel_structures = self.tunnel_shearzone_data()

        structures = pd.concat([
            borehole_structures,
            tunnel_structures],
            ignore_index=True, sort=False)

        # Fill NaN-values in all columns to 0 except in column 'shearzone', for which we do nothing.
        structures = structures.fillna(
            value={
                **{s: 0 for s in borehole_structures},
                **{'shearzone': np.nan}
            })

        mapping = {'x': 'x', 'y': 'y', 'z': 'z', 'depth': 'depth',
                   'upward_gradient': 'upward_gradient', 'azimuth': 'azimuth_bh'}
        self.bh_struc_to_global_coords(structures, **mapping)

        return structures

    def get_shearzone(self, sz: str, coords: str = 'gts'):
        """ Extract shear-zone coordinates for a given shear-zone

        Coordinates extracted will either be 'swiss' or 'gts'.

        Parameters:
        sz (str): Name of shear-zone (S1_1, S1_2, S1_3, S3_1, S3_2)
        coords (str, Default: 'gts'):
            Get coordinates in 'gts' or 'swiss'.

        Returns
        np.ndarray (3, n): Coordinates of shearzone intersections.
        """
        df = self.structures
        assert sz in ['S1_1', 'S1_2', 'S1_3', 'S3_1', 'S3_2'], f"unknown shear-zone {sz}."
        assert coords in ['swiss', 'gts'], f"unknown coordinate system {coords}."
        sz = df.loc[df.shearzone == sz, (f'x_{coords}', f'y_{coords}', f'z_{coords}')]
        return sz.to_numpy().T

    def structures_depth(self, borehole: str, depth: np.ndarray, structure=None, shearzone=None, coords='gts'):
        """ Get structures in a borehole at depth

        For a given borehole, and a given depth (or depth interval),
        get all structures - or a subset of structures, or specific shearzones.

        Parameters:
        borehole (str): name of borehole (INJ1, INJ2, ...)
        depth (np.array): Depth interval in borehole.
        structure (str or list, Optiona): Filter by certain structures
            (Fracture, Minor ductile Shear-zone, S1 Shear-zone, Quartz, ...)
        shearzone (str or list, Optional): Filter by certain shear-zones
            (S1_1, S1_2, ...)
        coords (str, optional): which coordinate system to return

        Returns:
        pd.DataFrame: Filtered dataframe

        """
        assert (depth.shape[0] == 2) and (depth[0] <= depth[1]), "Depth must be given as an interval."
        assert coords in ['swiss', 'gts'], f"unknown coordinate system {coords}."
        df = self.structures

        # Structure mask
        if isinstance(structure, str):
            structure = [structure]
        if structure is not None:
            _mask_struc = df.type.isin(structure)
        else:
            _mask_struc = np.ones(df.shape[0], dtype=bool)

        # Shear-zone mask
        if isinstance(shearzone, str):
            shearzone = [shearzone]
        if shearzone is not None:
            _mask_sz = df.shearzone.isin(shearzone)
        else:
            _mask_sz = np.ones(df.shape[0], dtype=bool)

        # Borehole mask
        _mask_bh = df.borehole == borehole

        # Depth mask
        _mask_depth = (depth[0] <= df.depth) & (df.depth <= depth[1])

        # Full mask #
        _mask = _mask_bh & _mask_depth & _mask_struc & _mask_sz

        # Filter DataFrame #
        _bh = df.loc[_mask, ('depth', 'azimuth_struc', 'dip', 'aperture', 'type',
                             'borehole', 'shearzone',
                             f'x_{coords}', f'y_{coords}', f'z_{coords}')]
        return _bh


def swiss_to_gts(v):
    """ Convert from swiss coordinates to gts coordinates

    GTS coordinates are: (x,y,z) = (667400, 158800, 1700)

    Parameters:
    v (np.array (3,)): Coordinate array

    """
    return v - np.array([667400, 158800, 1700])


def rename_sz(sz):
    """ Rename shearzone on form '12' to 'S1_2'. """
    sz = str(sz)
    sz_set = sz[0]
    num = sz[1]
    return f'S{sz_set}_{num}'

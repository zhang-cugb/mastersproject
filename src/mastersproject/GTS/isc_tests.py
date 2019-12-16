from GTS.ISC_data.isc import ISCData

class test_isc_data_import:

    def __init__(self):
        self.cls = ISCData()

    def test_import_borehole_data(self):
        """ Test the method 'cls.borehole_data(self)'.

        """
        cls = self.cls
        df = cls.borehole_data()
        # assert(df[])
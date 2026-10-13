import hashlib
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file

from lagcheck.models import InputError
from lagcheck.netcdf import inspect_netcdf, load_netcdf_signal


class NetCDFTests(unittest.TestCase):
    def test_equal_length_but_different_dimensions_rejected(self):
        path=self.root/'different_dimensions.nc'
        with netcdf_file(path,'w') as nc:
            nc.createDimension('time_axis',3)
            nc.createDimension('unrelated_axis',3)
            nc.createVariable('clock','d',('time_axis',))[:]=[0,1,2]
            nc.createVariable('Fe','d',('unrelated_axis',))[:]=[10,20,30]
        with self.assertRaises(InputError):
            load_netcdf_signal(path,dict(self.config,signal_variables=['Fe']))
        self.assertEqual(inspect_netcdf(path)['variables']['clock']['dimensions'],['time_axis'])

    def test_conflicting_known_time_units_rejected(self):
        path=self.root/'units.nc'
        with netcdf_file(path,'w') as nc:
            nc.createDimension('sample',3)
            time=nc.createVariable('clock','d',('sample',));time[:]=[0,1,2]
            time.units='seconds since 2026-01-01 00:00:00'
            nc.createVariable('Fe','d',('sample',))[:]=[10,20,30]
        with self.assertRaises(InputError):
            load_netcdf_signal(path,dict(self.config,signal_variables=['Fe']))
        result=load_netcdf_signal(path,dict(self.config,time_unit='s',signal_variables=['Fe']))
        self.assertTrue(result.metadata['time_unit_metadata_checked'])
        self.assertEqual(result.metadata['variable_dimensions']['Fe'],['sample'])

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.config = {"time_variable": "clock", "time_unit": "ms", "signal_variables": ["Fe", "Si"]}

    def tearDown(self):
        self.temp.cleanup()

    def fixture(self, time=(1000, 2000, 3000), fe=(1, 2, 3), si=(4, 5, 6), attrs=None, matrix=False, version=1, time_unit='milliseconds'):
        path = self.root / "test.nc"
        with netcdf_file(path, "w", version=version) as output:
            output.title = "Synthetic fixture; no vendor export"
            output.createDimension("n", len(time))
            clock = output.createVariable("clock", "d", ("n",))
            clock.units = time_unit
            clock[:] = time
            if matrix:
                output.createDimension("mass", 1)
                iron = output.createVariable("Fe", "d", ("n", "mass"))
                iron[:] = np.asarray(fe)[:, None]
            else:
                dim='n' if len(fe)==len(time) else 'nf'
                if dim=='nf': output.createDimension(dim,len(fe))
                iron = output.createVariable("Fe", "d", (dim,))
                iron[:] = fe
            for name, value in (attrs or {}).items():
                setattr(iron, name, value)
            dim='n' if len(si)==len(time) else 'ns'
            if dim=='ns': output.createDimension(dim,len(si))
            silicon = output.createVariable("Si", "d", (dim,))
            silicon[:] = si
        return path

    def test_units_sum_input_hash_and_classic_versions(self):
        for version in (1, 2):
            path = self.fixture(version=version)
            original = path.read_bytes()
            result = load_netcdf_signal(path, self.config)
            np.testing.assert_array_equal(result.time, [0, 1, 2])
            np.testing.assert_array_equal(result.signal, [5, 7, 9])
            self.assertEqual(result.metadata["time_first_seconds"], 1.)
            self.assertEqual(result.metadata["sha256"], hashlib.sha256(original).hexdigest())
            self.assertEqual(path.read_bytes(), original)
        for unit, factor in (("s", 1), ("ms", .001), ("us", .000001), ("min", 60)):
            path=self.fixture(time_unit=unit)
            result = load_netcdf_signal(path, {**self.config, "time_unit": unit})
            np.testing.assert_allclose(result.time, np.array([0, 1000, 2000]) * factor)

    def test_inspection_has_metadata_without_signal_arrays(self):
        path = self.fixture()
        result = inspect_netcdf(path)
        self.assertEqual(result["variables"]["Fe"]["shape"], [3])
        self.assertEqual(result["variables"]["clock"]["units"], "milliseconds")
        self.assertEqual(result["global_attributes"]["title"], "Synthetic fixture; no vendor export")
        self.assertNotIn("data", result["variables"]["Fe"])
        self.assertTrue(result["input_unchanged"])

    def test_explicit_mapping_required(self):
        path = self.fixture()
        for key in self.config:
            config = dict(self.config)
            del config[key]
            with self.assertRaises(InputError):
                load_netcdf_signal(path, config)
        for update in ({"signal_variables": []}, {"signal_variables": ["Fe", "Fe"]},
                       {"signal_variables": ["clock"]}, {"time_variable": "missing"},
                       {"signal_variables": ["missing"]}, {"time_unit": "hours"}):
            with self.assertRaises(InputError):
                load_netcdf_signal(path, {**self.config, **update})

    def test_nan_infinite_and_unordered_time_rejected(self):
        for time in ((0, 0, 1), (0, 2, 1), (0, np.nan, 2), (0, 1, np.inf)):
            with self.assertRaises(InputError):
                load_netcdf_signal(self.fixture(time=time), self.config)
        for fe in ((0, np.nan, 2), (0, 1, np.inf)):
            with self.assertRaises(InputError):
                load_netcdf_signal(self.fixture(fe=fe), self.config)

    def test_matrix_unequal_length_and_short_trace_rejected(self):
        for options in ({"matrix": True}, {"fe": (1, 2)}, {"time": (1, 2), "fe": (1, 2), "si": (1, 2)}):
            with self.assertRaises(InputError):
                load_netcdf_signal(self.fixture(**options), self.config)

    def test_fill_missing_and_packed_variables_rejected(self):
        for attributes in ({"_FillValue": -999.}, {"missing_value": -999.},
                           {"missing_value": np.array([-999., -888.])}):
            with self.assertRaises(InputError):
                load_netcdf_signal(self.fixture(fe=(1, -999, 3), attrs=attributes), self.config)
        with self.assertRaises(InputError):
            load_netcdf_signal(self.fixture(fe=(1, 9.969209968386869e36, 3)), self.config)
        for attrs in ({"scale_factor": .1}, {"add_offset": 10.}, {"_Unsigned": "true"}):
            with self.assertRaises(InputError):
                load_netcdf_signal(self.fixture(attrs=attrs), self.config)

    def test_declared_fill_without_missing_cells_can_be_read(self):
        result = load_netcdf_signal(self.fixture(attrs={"_FillValue": -999.}), self.config)
        np.testing.assert_array_equal(result.signal, [5, 7, 9])

    def test_sum_overflow_and_unsupported_container_rejected(self):
        with self.assertRaises(InputError):
            load_netcdf_signal(self.fixture(fe=(1e308, 1e308, 1e308), si=(1e308, 1e308, 1e308)), self.config)
        path = self.root / "hdf5.nc"
        path.write_bytes(b"\x89HDF\r\n\x1a\n")
        with self.assertRaises(InputError):
            inspect_netcdf(path)
        path.write_bytes(b"CDF\x01")
        with self.assertRaises(InputError):
            inspect_netcdf(path)

    def test_scalar_and_character_variables_rejected(self):
        for dtype, dimensions in (("d", ()), ("c", ("n",))):
            path = self.fixture()
            with netcdf_file(path, "a") as output:
                item = output.createVariable("unsupported", dtype, dimensions)
                item.data[...] = 1 if dtype == "d" else b"a"
            with self.assertRaises(InputError):
                load_netcdf_signal(path, {**self.config, "signal_variables": ["unsupported"]})


if __name__ == "__main__":
    unittest.main()

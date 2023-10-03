import matplotlib
matplotlib.use('agg')
import unittest

#from tempdir.tempfile_ import TemporaryDirectory

class TestPopIIIHburn(unittest.TestCase):

    def test_popIII_abu(self):
        from .test_functions import ppn_abundance_test
        from .test_functions import dataget_ppn
        import os

        dataget_ppn('ppn_PopIIIHburn', ['00001', '00020', '00038'])
        ppn_abundance_test('ppn_PopIIIHburn', ['00001', '00020', '00038'])
        os.system('rm -rf NuPPN/examples/regression_tests/ppn_PopIIIHburn')

class Testc13_pocket_abu(unittest.TestCase):

    def test_c13_abu(self):
        from .test_functions import ppn_abundance_test
        from .test_functions import dataget_ppn
        import os

        #modified cycle from 1,62, 125 to 1,20,38 for now
        dataget_ppn('ppn_C13_pocket', ['00001', '00062', '00125'])
        ppn_abundance_test('ppn_C13_pocket', ['00001', '00062', '00125'])
        os.system('rm -rf NuPPN/examples/regression_tests/ppn_C13_pocket')

class test_mppnp_Hcore_burning(unittest.TestCase):
    '''
    Restart M3.00Z0.020 from cycle 500.
    '''
    def test_abu_profiles(self):
	'''
	Test profile after 20 cycles.
	'''
        from .test_functions import mppnp_H5_test
        from .test_functions import dataget_mppnp
	import os

	files = ['M3.00Z0.020.0000001.out.h5']
        case='mppnp_Hcore_burning'
	test_cycles = [520] 

        dataget_mppnp(case, files=files)
	mppnp_H5_test(test_cycles,case,files)
        #shutil.rmtree('NuPPN/examples/regression_tests/' + case)
        os.system('rm -rf NuPPN/examples/regression_tests/'+case)

class test_mppnp_Hecore_burning(unittest.TestCase):
    '''
    Restart M5.00Z0.020 from cycle 1500.
    '''
    def test_abu_profiles(self):
	'''
	Test profile after 20 cycles.
	'''

        from .test_functions import mppnp_H5_test
        from .test_functions import dataget_mppnp
	import os

	files = ['M5.00Z0.020.0001501.out.h5']
        case='mppnp_Hecore_burning'
	test_cycles = [1520]        

        dataget_mppnp(case, files)
	mppnp_H5_test(test_cycles,case,files)

        os.system('rm -rf NuPPN/examples/regression_tests/'+case)
        #shutil.rmtree('NuPPN/examples/regression_tests/' + case)


class test_mppnp_HBB(unittest.TestCase):
    '''
    Restart M6.00Z.0100 from cycle 4800.
    '''
    def test_abu_profiles(self):
	'''
	Test profile after 20 cycles
	'''

        from .test_functions import mppnp_H5_test
        #import shutil
        from .test_functions import dataget_mppnp
	import os

	files=['M6.00Z.0100.0004801.out.h5']
	case='mppnp_HBB'
	test_cycles = [4820]

	dataget_mppnp(case, files)
        mppnp_H5_test(test_cycles,case,files)

        os.system('rm -rf NuPPN/examples/regression_tests/'+case)
        #shutil.rmtree('NuPPN/examples/regression_tests/' + case)

    def test_surface_abu_evol(self):
	'''
	Test surface evolution of 20 cycles.
	'''

        #import shutil
        from .test_functions import dataget_mppnp
        from .test_functions import mppnp_H5_test
	import os

	files=['M6.00Z.0100.0004801.surf.h5']
	case='mppnp_HBB'
	test_cycles = range(4801,4821,1)

	dataget_mppnp(case, files)
	mppnp_H5_test(test_cycles,case,files)

        os.system('rm -rf NuPPN/examples/regression_tests/'+case)


if __name__ == '__main__':
    unittest.main()

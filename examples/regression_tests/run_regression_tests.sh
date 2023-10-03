#this resembles what is executed in the travis file.
#all output goes into
echo 'For output see file run_regression_tests.out'
cd ../../
docker build -t testsuite . >> examples/regression_tests/run_regression_tests.out 2>&1
#run ppn examples
#PoPIII does not work and therefore commented
#docker run testsuite /bin/sh -c "cp home/NuPPN/frames/ppn/CODE/MAKE_LOCAL/Make.local.docker home/NuPPN/frames/ppn/CODE/Make.local; cd home/NuPPN/examples/ppn_PopIIIHburn; make distclean; make; ./ppn.exe; cd /home; python -m NuPPN.examples.regression_tests.selftest TestPopIIIHburn" 
#docker run testsuite /bin/sh -c "cp home/NuPPN/frames/ppn/CODE/MAKE_LOCAL/Make.local.docker home/NuPPN/frames/ppn/CODE/Make.local; cd home/NuPPN/examples/ppn_C13_pocket; make distclean; make; ./ppn.exe; cd /home; python -m NuPPN.examples.regression_tests.selftest Testc13_pocket_abu" >> examples/regression_tests/run_regression_tests.out 2>&1
#run mppnp examples
docker run testsuite /bin/sh -c "cp home/NuPPN/frames/mppnp/CODE/MAKE_LOCAL/Make.local.docker home/NuPPN/frames/mppnp/CODE/Make.local; cd home/NuPPN/examples/mppnp_Hcore_burning; make distclean; make clean; make; ./setup.sh; ./run_regression_test.sh; cd /home; python -m NuPPN.examples.regression_tests.selftest test_mppnp_Hcore_burning" memory=inf, memory-swap=inf  >> examples/regression_tests/run_regression_tests.out 2>&1
docker run testsuite /bin/sh -c "cp home/NuPPN/frames/mppnp/CODE/MAKE_LOCAL/Make.local.docker home/NuPPN/frames/mppnp/CODE/Make.local; cd home/NuPPN/examples/mppnp_Hecore_burning; make distclean; make clean; make; ./setup.sh; ./run_regression_test.sh; cd /home; python -m NuPPN.examples.regression_tests.selftest test_mppnp_Hecore_burning" memory=inf, memory-swap=inf >> examples/regression_tests/run_regression_tests.out 2>&1
docker run testsuite /bin/sh -c "cp home/NuPPN/frames/mppnp/CODE/MAKE_LOCAL/Make.local.docker home/NuPPN/frames/mppnp/CODE/Make.local; cd home/NuPPN/examples/mppnp_HBB; make distclean; make clean; make; ./setup.sh; ./run_regression_test.sh; cd /home; python -m NuPPN.examples.regression_tests.selftest test_mppnp_HBB" memory=inf, memory-swap=inf >> examples/regression_tests/run_regression_tests.out 2>&1


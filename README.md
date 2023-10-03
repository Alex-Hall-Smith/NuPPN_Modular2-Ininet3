# NuPPN

Changes compared to NuGrid/nuppn-modular2_ininet3 branch:
- edit parameter_physics for new REACLIB size
- reaclib.F90 edited name for new REACLIB files
- NPDATA/REACLIB now contains both REACLIB files

## TODO:

* instantaneous decay at network edges should be optional
* more appropriately-sized networks should be available in the
  `isotopedatabase*.txt` files. Otherwise we do un-necessary work and store
  unnecessary data
* the default output from tppnp should be single precision and the python codes
  should handle this
* network should perform a check each time step to see if the abundances along
  the network edges are growing, and stop the code if the network is too small
  for the problem

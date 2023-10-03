# NuPPN

These are the NuGrid Post-Processing Network (PPN) codes. 

_See the [wiki](https://github.com/NuGrid/NuPPN/wiki) for more details_

This repo starts with the svn rev:6626. svn repo will remain available for archival purpose, no new commits to svn.

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

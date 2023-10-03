# TPPNP Python modules

__[TL;DR]__ - there is a crash course at the bottom of this page

The core of the python framework is the `particle_set` method in `tppnp.py`. It
is a class providing a data structure and methods for a set of tracer
particles. 

The particle set stores time-evolution and meta-data ('trajectories') for
tracer particles along wth initial or final abundances for those particles.

A particle set can be either an 'input' or 'output' particle set, pertaining to
the TPPNP code. An 'input' particle set would contain time-evolution data
(time, temperature, density), 'meta-data' (particle mass, initial position etc)
and initial abundances. An 'output' particle set would contain meta-data and
final abundances.

Once trajectories and initial abundances have been loaded for a particle set,
we can write a binary file that is readable by the post-processing code. One
can also load the results from a post-processing calculation as a particle set.

Abundances, either initial or final, are stored for each particle as an
instance of the abundances class (provided by the abundances module in the file
`abundances.py`). The abundances class has three properties: `a`, `z` and
`xiso`, which are the mass number, charge number and mass fractions of the
isotopes in the object. There are methods to do things such as calculate the
electron fraction or return elemental masses, and so on.

## Loading an input particle set and write it to the TPPNP input format

Define a new particle set by:

```python
import tppnp as t
p = t.particle_set()
```

Tracer particle information can be loaded by:

```python
p.load_trajectories(*args,**kwargs)
```

Currently we can load particles from Chris Fryer's 1d CCSN code or particles
from Dean Townsley's 2d/3d FLASH simulations.

Initial abundances can be loaded a few different ways. Look at the
`load_initial_abundances` docstring for more information. You can currently:
1. provide the composition of a progenitor star as 1d profiles and have the method
interpolate it onto the tracer particles for you
2. provide the initial composition of each particle
3. provide one initial composition to be used for all particles

e.g. to load initial composition for the particles from a progenitor model with
mass, abundances, a and z saved in a `.npz` (numpy zipped binary) file:

```python
d = np.load("M15.npz")
p.load_initial_abundances(mode=1, mass=d['m'], yps=d['yps'], a=d['a'], z=d['z'])
```

Once trajectories and initial abundances have been loaded, write the input file
for tppnp like so:

```python
p.write_tppnp_input("tp")
```

where "tp" is the filename you want to write to (will actually be appended with
the suffix appropriate for the endianness of your system, so tp.lab or tp.msb).

## Loading results from a tppnp post-processing calculation

Define a new particle set by:

```python
p = t.particle_set()
```

If your tppnp run directory was called "./shazbat", you can load the result of
a TPPNP calculation by:

```python
p.load_final_abundances("./shazbat/output")
```

## Data storage and accessing data once loaded into a particle set object

The bulk of the data are stored in a dictionary `p.data`. You can see the
objects in here by getting the dictionary keys by the usual `p.data.keys()`
call.

Directly accessing the data array may not be the easiest/most efficient way to
get the data you want. Additionally, there is no guarantee that a particle set
is complete, because there is no requirement that the particle numbers be a
continuous set of integers. This is accounted for in the various get methods
listed below.

There are multiple get methods of the particle set class that try to get you
the right data. They are (where n is the number of the particle you want the
data for; input/output denotes for what type of particle set this method should
work):

```python
p.xini(n)    # initial abundances (input)
p.xfin(n)    # final abundances (output)
p.temp(n)    # temperature evolution (input)
p.dens(n)    # density evolution (input)
p.time(n)    # time (input)
p.dm(n)      # particle mass (input/output)
p.posini(n)  # initial position (input/output)
p.posfin(n)  # final position (input/output)
p.tpeak(n)   # peak temperature (input/output)
p.rhopeak(n) # peak density (input/output)
p.velfin(n)  # final velocity (input/output)
p.pname(n)   # particle name (input/output)
p.nstp(n)    # number of timesteps (input)
p.nisoini(n) # number of isotopes in initial abundances (input)
p.nisofin(n) # number of isotopes in final abundances (output)
p.afin(n)    # mass number for final abundances (output)
p.aini(n)    # mass number for initial abundances (input)
p.zini(n)    # charge number for initial abundances (input)
p.zfin(n)    # charge number for final abundances (output)
p.yeini(n)   # initial electron fraction (input)
p.yefin(n)   # final electron fraction (output)
p.xisoini(n) # initial mass fractions of isotopes (input)
p.xisofin(n) # final mass fractions of isotopes(output)
p.isosini(n) # isotope names (input)
p.isosfin(n) # isotope names (output)
```

There's also some methods that return data for all particles at once (the names
should make sense by referencing the list above):

```python
p.tpeak_all():
p.posfin_all(dim=3):
p.posini_all():
p.dm_all():
p.xfin_all():
p.xini_all():
p.rhopeak_all():
p.pname_all():
p.pnum_all():
p.m()        # mass coordinate array - only makes sense for 1d lagrangian data (input/output)
```

Some other useful methods with their docstrings:

```python
p.yields(mode = 0, chop_mass = -1, el = False):
"""
mode 0: total mass of isotopes in final particle set
mode 1: total mass of isotopes in supernova explosion (including
        progenitor abundances).
el: True if you want elemental yields, False for isotopic
chop_mass: delimiting mass between processed and unprocessed ejecta. By default
it is set to be the edge of the particle set domain, but can be user-specified"""


p.isoyield(iso="PROT ", chop_mass = -1):
"""return yield (in solar masses) of given isotope"""

p.elyield(Z = 1, chop_mass = -1):
"""return yield (in solar masses) of given element"""

p.abundance(iso="PROT "):
""" return final abundance of given iso for every particle"""

```

## Loading a progenitor star - the `progenitor` class

For computing CCSN yields, we're going to need a progenitor star so that we can
add its composition to the shock-heated particles to make a complete set of
ejected isotopic masses.

There is a `progenitor` class provided by the file `progenitor.py` that stores
this data in a simple manner that can be used in conjunction with a particle
set to calculate yields.

Initialising a progenitor object requires the mass coordinate and composition
of the 1d progenitor star. This can be done either by passing the data or by
passing the path to a npz file containing the necessary data, e.g.:

```python
from progenitor import progenitor

# load from a npz file
star = progenitor(file="./M15.npz")

# or pass the data yourself:
d = np.load("M15.npz")
star = progenitor(m=d['m'], yps=d['yps'], a=d['a'], z=d['z'])
```

You can then attach a progenitor object to a particle set (i.e. "this is the
progenitor star of this explosion") by:

```python
p.load_progenitor(star)
```

The data stored for a progenitor are:

```python
star.a        # mass number
star.z        # charge number
star.m        # mass coordinate
star.nzones   # number of cells in 1d structure
star.nspecies # number of isotopes
star.x        # mass fractions (nzones,nspecies)
```

There are also methods to get other data, some of which are trivially computed
but would take up memory unnecessarily:

```python
star.dm(): # return mass of each zone
star.zuq(): # return list of charge number for the elements in this model
star.isos(): # return list of isotope names in this model
```

Other methods include:

```python
star.yields(mass_cut = 0., el = False):
"""return the abundances above the mass cut. This could be a physical mass cut or the
dimiting mass between the material processed/not processed by the supernova shock.
Returns integrated yield of each isotope (Y = sum_i(dm_i*x_ij)) in solar masses.
If el == True, elemental yields are returned.
"""

star.abundance(iso = "", mass_cut = 0.):
"""Return mass fractions of isotope as a function of mass coordinate.
If no isotope name is given, all isotopes are returned"""

```

## Crash course for the impatient

To load the results from a TPPNP post-processing simulation:

```python
import tppnp as t
p = t.particle_set()
p.load_final_abundances("<path_to_run_dir>/output")
```

Load the progenitor star and link it to the particle set:

```python
from progenitor import progenitor
star = progenitor("M15.npz")
p.load_progenitor(star)
```

get various types of yields:

```python
# integrated isotopic masses of progenitor star above mass cut:
p.star.yields(mass_cut = p.posini(1)[0])

# integrated isotopic masses of shock-heated particles:
p.yields(mode=0)

# ejected mass from whole supernova (shock-heated particles + pre-SN
# composition):
p.yields(mode=1)

```

Take a look at certain abundances vs mass coordinate:

```python
# mass coordinate and abundance of Titanium-44 in progenitor star:
m, xti44 = p.star.abundance("TI 44")

# mass coordinate and abundance of Titanium-44 in shock-heated particles:
m = p.m()
xti44 = p.abundance("TI 44")
```


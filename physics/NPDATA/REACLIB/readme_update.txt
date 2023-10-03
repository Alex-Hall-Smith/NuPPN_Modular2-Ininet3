This is a short description of the update of the (n,gamma) (Section 4) and 
(gamma,n) reactions (Section 2) performed in November 2006 by I. Dillmann.

Basis was the file reaclib.nosmo available on the nucastro.org website.
This file was updated by inclusion of the most recent stellar (n,gamma)
cross sections from KADoNiS v0.2 (http://nuclear-astrophysics.fzk.de/kadonis).

This replaces the previous "baka" entries (Bao et al. 1987) but also
extends the experimental data up to Bi. While the previous "baka" entries
used only the a0 parameter (a1...a6=0) for an assumed s wave behaviour,
all new "ka02" entries assume the NON-SMOKER (ADNDT 79 (2001) 47)
energy-dependencies.

The NON-SMOKER renormalization factors f are derived from the ratio between
the calculated stellar Maxwellian averaged cross section
(from ADNDT 75 (2000) 1) at kT=30 keV (<sigma*>_30) and the new values
available in KADoNiS v0.2 multiplied with the respective stellar enhancement
factors (which also comes from the calculation):
f= <sigma*>_30(calc) / {<sigma>_30(KADoNiS) * SEF}.

A list with the normalization factors can be found in the file "factors.txt".
It supersedes the factors given in Table 5 of Astrophys. J. 576 (2002) 323
(see also http://download.nucastro.org/astro/stellar).
It should be noted that nuclei with A less than about 40 and nuclei at magic
neutron numbers may have a level density which is too low to allow the
application of the statistical model at 30 keV. However, due to the lack of
other information, the renormalized Hauser-Feshbach result is used nevertheless.
This may explain certain sizeable renormalizations. These cases cannot be
improved by improving the statistical model calculation but have to be
treated by explicitly accounting for resonant and direct contributions.

All updated entries have the new reference "ka02".

Iris Dillmann

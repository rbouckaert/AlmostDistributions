# AlmostDistributions

[![Build Status](https://travis-ci.org/rbouckaert/asc.svg?branch=master)](https://travis-ci.org/rbouckaert/AlmostDistributions)

BEAST 2 package for distributions of infinite support:

AlmostUniform: has very low support for values below 'lower' and above 'upper'. Useful for setting up analyses that are difficult to initialise. The log density of x when x is outside the lower-upper range is calculated as -penalty * (x-centre) * (x-centre) where penalty is determined by the penalty input (defaults to 10000) and centre is the middle of the range between lower and upper.

AlmostNormal: has very low support for values outside 95% HPD of a Normal distribution. For values outside that range, the log density of x is calculated using the same formula as for AlmostUniform.


# Installing the package
AlmostDistributions is a [BEAST2](http://beast2.org) package.
If you have not already done so, you can get BEAST 2 from [here](http://beast2.org).

To install AlmostDistributions, it is easiest to start BEAUti (a program that is part of BEAST), and select the menu File/Manage packages. A package manager dialog pops up, that looks something like this:

![Package Manager](https://github.com/rbouckaert/AlmostDistributions/raw/master/doc/package-manager.png)

If the AlmostDistributions package is listed, just click on it to select it, and hit the Install/Upgrade button.

If the AlmostDistributions package is not listed, you may need to add a package repository by clicking the "Package repositories" button. A window pops up where you can click "Add URL" and add "https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml" (or if you still run BEAST v2.6 "https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra.xml") in the entry. After clicking OK, the dialog should look something like this:

![Package Repositories](https://github.com/rbouckaert/obama/raw/master/doc/package_repos.png)

Click OK and now AlmostDistributions should be listed in the package manager (as in the first dialog above). Select and click Install/Upgrade to install.

# Citing the package

Remco Bouckaert AlmostDistributions package for BEAST 2. 2021 DOI: 10.5281/zenodo.5348449
[![DOI](https://zenodo.org/badge/109171998.svg)](https://zenodo.org/badge/latestdoi/109171998)

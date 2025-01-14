<!-- [![Travis-CI Build Status](https://travis-ci.org/rdiaz02/OncoSimul.svg?branch=master)](https://travis-ci.org/rdiaz02/OncoSimul) -->
<!-- [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rdiaz02/OncoSimul?branch=master&svg=true)](https://ci.appveyor.com/project/rdiaz02/OncoSimul) -->
<!-- [![codecov.io](https://codecov.io/github/rdiaz02/OncoSimul/coverage.svg?branch=master)](https://codecov.io/github/rdiaz02/OncoSimul?branch=master) -->



# OncoSimul

Code for forward population genetic simulation in asexual populations,
with special focus on cancer progression.  Fitness can be an arbitrary
function of genetic interactions between multiple genes or modules of
genes, including epistasis, order restrictions in mutation accumulation,
and order effects. Fitness can also be a function of the relative and
absolute frequencies of other genotypes (i.e., frequency-dependent
fitness).  Mutation rates can differ between genes, and we can include
mutator/antimutator genes (to model mutator phenotypes).  Simulating
multi-species scenarios and therapeutic interventions is also
possible. Simulations use continuous-time models and can include driver
and passenger genes and modules. Also included are functions for:
simulating random DAGs of the type found in Oncogenetic Trees, Conjunctive
Bayesian Networks, and other cancer progression models; plotting and
sampling from single or multiple realizations of the simulations,
including single-cell sampling; plotting the parent-child relationships of
the clones; generating random fitness landscapes (Rough Mount Fuji, House
of Cards, additive NK, Ising, and Eggbox models) and plotting them.

There is also new functionality to allow for frequency-dependent fitness, interventions (including adaptive therapy), user specification of both birth and death rates, and user-defined variables that can be used for interventions.




The /OncoSimulR directory contains the code for the [BioConductor](http://www.bioconductor.org) package
[OncoSimulR](http://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html). The
/miscell-files directory contains additional files so far only related to
the above.


A former version of this code has been used in the paper "Identifying
restrictions in the order of accumulation of mutations during tumor
progression: effects of passengers, evolutionary models, and
sampling",
[BMC Bioinformatics, 2015, 16:41](http://www.biomedcentral.com/1471-2105/16/41).
OncoSimulR has also been used extensively in the simulations reported in
the Bioinformatics 
paper
["Cancer Progression Models And Fitness Landscapes: A Many-To-Many
Relationship"](https://doi.org/10.1093/bioinformatics/btx663) and the PLoS
Computational Biology paper 
["Every which way? On predicting tumor evolution using cancer progression
models"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007246).




You can also find
[OncoSimulR](https://popmodels.cancercontrol.cancer.gov/gsr/packages/oncosimulr/)
on the [Genetic Simulation Resources
catalogue](https://popmodels.cancercontrol.cancer.gov/gsr/).
<!-- <a href="http://popmodels.cancercontrol.cancer.gov/gsr/"><img src="http://popmodels.cancercontrol.cancer.gov/gsr/static/img/gsr_tile.jpg" alt="Catalogued on GSR" width="180" height="60" /></a> -->

# Installation 

The frequency-dependent fitness functionality became available in BioConductor
starting with version 3.13 and is now available in both the stable and devel
versions. To use the most recent code in BioConductor, install the devel version.

```r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("OncoSimulR", version = "devel")

```

If you want to use the stable BioConductor version, remove the `version =
"devel"` part of the invocation. 



<!-- <\!-- that I regard as stable you first need to [use -\-> -->
<!-- <\!-- a development version of Bioconductor](http://www.bioconductor.org/developers/how-to/useDevel/). Most -\-> -->
<!-- <\!-- of the time, from R you only need to do: -\-> -->


<!-- ```r -->
<!--     library(BiocInstaller)  -->
<!--     useDevel() -->

<!-- ``` -->

<!-- Then the next code will install the development version of the package and -->
<!-- its dependencies, if needed: -->


<!-- ```r -->
<!--     source("http://bioconductor.org/biocLite.R") -->
<!--     biocLite("OncoSimulR") -->
<!-- ``` -->

To start using it:

```r
library(OncoSimulR)
```

<!-- The above, however, **will not install the version with frequency -->
<!-- dependent fitness**.  -->
<!-- **To use the frequency-dependent fitness version** read the following. -->


<!-- ## Installing the frequency-dependent fitness branch -->

<!-- The code is available from BioConductor. To use the development version from -->
<!-- BioConductor, do (`BiocManager::install("OncoSimulR", version = "devel")`);; -->

This should work for all operating systems.  If not, or if you want to install
from sources, read on.

### If you use Linux and other Unixes (Macs)

You should install from github as follows: <!-- (and this might be newer
than the BioC code) -->

<!-- ```r -->
<!-- if (!require("devtools")) -->
<!--     install.packages("devtools") ## if you don't have it already -->
<!-- library(devtools) -->
<!-- install_github("rdiaz02/OncoSimul/OncoSimulR") -->
<!-- ```  -->

<!-- If you want to do development with the package, you might want to do, instead -->

```r
if (!require("devtools"))
    install.packages("devtools") ## if you don't have it already
library(devtools)
install_github("rdiaz02/OncoSimul/OncoSimulR", 
               dependencies = TRUE)
``` 
(setting `depèndencies = TRUE` ensures that "Suggests" are also installed).



### If you use Windows

If you use Rtools40, the new toolchain starting R-4.0.0, April 2000 (so a long time ago already)
(https://cran.r-project.org/bin/windows/Rtools/), you can install
OncoSimulR from sources. These are old notes.

<!-- We have not uploaded the changes in the freq-dep-fitness branch to -->
<!-- BioConductor because there are known problems compiling ExprTk with MinGW -->
<!-- (e.g., -->
<!-- https://sourceforge.net/p/mingw-w64/discussion/723797/thread/c6b70624/). -->

<!-- Things work with other toolchains and, eventually, Rtools 40 should become -->
<!-- the default toolchain, and the problem will get solved. In the meantime, -->
<!-- you have two options: -->


<!-- #### Install from a zip file: -->

<!-- (First, install the current OncoSimulR from BioConductor, to resolve all the -->
<!-- dependencies; see above or go to http://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html). -->


<!-- Download the [OncoSimulR_2.17.xyz.zip file we -->
<!--  provide](https://rdiaz02.github.io/OncoSimul/OncoSimulR_2.17.992.zip); for example, from your -->
<!--  web browser, place the mouse on the link, and right click to select "Save -->
<!--  link as", or similar incantation. -->
 
<!-- Now install that zip file from R (e.g., from the menu, go to "Packages", -->
<!-- "Install package(s) from local file(s)", and select the file). This works -->
<!-- with [R-3.6.0](https://cran.r-project.org/bin/windows/base/) (and -->
<!-- [R-3.6.0, -->
<!-- patched](https://cran.r-project.org/bin/windows/base/rtest.html) and -->
<!-- [R-devel, future R-3.7.0](https://cran.r-project.org/bin/windows/base/rdevel.html)) -->

<!-- <\!-- That zip file, for the sake of being small and fast to build and -\-> -->
<!-- <\!-- install, does not contain the vignette. Use the links to the [HTML -\-> -->
<!-- <\!-- vignette](https://rdiaz02.github.io/OncoSimul/OncoSimulR.html) or the [PDF -\-> -->
<!-- <\!-- vignette](https://rdiaz02.github.io/OncoSimul/pdfs/OncoSimulR.pdf). (The -\-> -->
<!-- <\!-- standard help is of course included). -\-> -->


<!-- #### Install using Rtools40 -->

<!-- (This is more work and takes more much more time) -->

<!-- <\!-- **The frequency-dependent fitness version will not install from source in -\-> -->
<!-- <\!-- Windows unless you use Rtools40.** This is a known problem with ExprTk and -\-> -->
<!-- <\!-- MinGW -\-> -->

<!-- We have verified that OncoSimulR (at least as of 2019-05-24) does install -->
<!-- with [Rtools40](https://cran.r-project.org/bin/windows/testing/rtools40.html). -->

How to do it: The standard installation procedure should work, but the
following steps might help.

 1. Install Rtools42 and R as explained in https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html

 <!-- 2. (If you experience issues with igraph, install igraph following these notes: -->
 <!--       https://github.com/r-windows/checks/issues/2 -->
 <!-- 	   ) -->
 2. Now, install OncoSimulR from BioConductor (to resolve all
    dependencies in one go): https://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html
	For a fresh R installation this can take more than one hour.
 3. Clone the git repo and move to that directory.
 4. Go to a MINGW shell console, and install. For example, if you have
    installed R-testing under 'C:\R', you can do
```
/c/R/R-testing/bin/x64/R CMD INSTALL --no-multiarch OncoSimulR
```

   Alternatively, install from a local file, but you need to specify the
   tar.gz (the zip file will not work, of course, since the R-testing that
   ships for/with Rtools will not install from zip files)
  
   Installing from source takes a while (more than 5 minutes). 





<!-- Note that sometimes the latest additions in this repo could be broken (see -->
<!-- [Software status](#software-status)). And you can of course clone this -->
<!-- repo, and install from there. -->


## BioConductor github repository

The github repository for this package is this one:
https://github.com/rdiaz02/OncoSimul . Since mid-2017 BioConductor is
maintained using git, but since this directory contains other files and
directories in addition to the OncoSimulR package itself, I have not used
option ["Sync an existing GitHub repository with
Bioconductor"](https://www.bioconductor.org/developers/how-to/git/sync-existing-repositories). Instead,
I continue using this github repo, but then locally update a
Bioconductor-only repository of just the OncoSimulR code (as explained in
[Maintain a Bioconductor-only repository for an existing
package](https://www.bioconductor.org/developers/how-to/git/maintain-bioc-only/)).

# Documentation


As any R/BioConductor package, OncoSimulR comes with documentation for its
user-visible functions and data sets (using the help is just standard R
usage). From
[OncoSimulR's BioConductor page](https://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html)
you have access to the standard documentation both the manual and overview
---the vignette. The best place to start is the vignette (created from the
`OncoSimulR/vignettes/OncoSimulR.Rnw` file that includes both text and
code).

<!-- you can obtain the  -->
<!-- [PDF reference manual](http://www.bioconductor.org/packages/3.2/bioc/manuals/OncoSimulR/man/OncoSimulR.pdf). A -->
<!-- better place to start, though, is the long vignette, with commented -->
<!-- examples (and created from the `OncoSimulR/vignettes/OncoSimulR.Rnw` file -->
<!-- that includes both text and code). Here is -->
<!-- [the vignette as PDF](http://www.bioconductor.org/packages/devel/bioc/vignettes/OncoSimulR/inst/doc/OncoSimulR.pdf), -->
<!-- from the BioConductor site (the development branch). -->


You can view the vignette from R itself doing


```r
browseVignettes("OncoSimulR")
```

and this gives you access to the HTML, the Rmd file (markdown + R), and the R code.

## Documentation: HTML and PDF for this repo's version


From these two links you can also
[browse the HTML vignette](https://rdiaz02.github.io/OncoSimul/OncoSimulR.html)
and [get a PDF version](https://rdiaz02.github.io/OncoSimul/pdfs/OncoSimulR.pdf).

These files correspond to the most recent, github version, of the package
(i.e., they might include changes not yet available from the BioConudctor
package). Beware that the might have figures and R code that
do not fit on the page, etc.


## Further documentation

This [paper published in
Bioinformatics](https://doi.org/10.1093/bioinformatics/btx077) gives a
quick overview of OncoSimulR (a former version is available as a [bioRxiv
preprint](http://biorxiv.org/content/early/2016/08/14/069500)). You can
also take a look at this [poster presented at ECCB
2016](http://dx.doi.org/10.7490/f1000research.1112860.1). A chapter titled
"Simulating evolution in asexual populations with epistasis" will appear
in K.-C. Wong (ed.), Epistasis: Methods and Protocols. Methods in
Molecular Biology, later during 2020; here is a [preprint
version](https://repositorio.uam.es/bitstream/handle/10486/690422/simul_evol_epist_preprint_v3_with_figures.pdf).


If you use the package in publications **please cite the Bioinformatics paper**.


## Major contributors and main additions to the code

The original implementation is by Ramon Diaz-Uriarte. 

The frequency-dependent fitness functionality is based on [Sergio Sanchez-Carrilo's Master's thesis](https://repositorio.uam.es/handle/10486/685417) (see also file ['miscell-files/Sergio_Sanchez-Carillo-improvements-post-TFM.pdf'](miscell-files/Sergio_Sanchez-Carillo-improvements-post-TFM.pdf) for additional features that were not described in the original thesis) and additional functionality has been added by Juan Antonio Miguel González (2020). Ramon Diaz-Uriarte further simplified and generalized the specification of fitness under frequency-dependence.

Work by Niklas Endres (late 2019 to early 2020) allows to use T (time) in fitness; fitness could already be made to depend on population sizes, and now, then, also on time. This allows for complex scenarios, as illustrated in the vignette (but adaptive therapy was not really possible as we cannot condition on previous states arbitrarily). 

Work by Javier Muñoz Haro (spring 2021) makes interventions possible on population sizes; these can take place at arbitrary times, be repeated, etc. We can have extremely flexible intervention mechanisms, also as illustrated in the vignette (but adaptive therapy still not possible as we cannot use arbitrary user variables to condition on previous states of the population). 

Work by Alberto González Klein (spring 2021) allows users to specify death rates. Before this work, death rates where not user-modifiable: they were 1 for the exponential model and a fixed function of population size in the McFarland model. Now, users can specify death rates and use frequency-dependent models that affect birth, or death, or both.

Work by Javier López Cano (spring 2022) unifies the additions by Javier Muñoz Haro and Alberto González Klein and also makes user variables possible; now we can emulate adaptive therapy (see vignette). (Note, though, that two minor limitations remain: (a) we cannot change mutation rates; (b) it is unclear that we can modify fitness based on user-variables.)



# Licenses and copyright


The R/BioConductor OncoSimulR package is licensed under the GPLv3
license. The code for the OncoSimulR BioConductor package, except for
functions `plot.stream` and `plot.stacked`, is Copyright 2012-2020 by
Ramon Diaz-Uriarte; the code for frequency dependent fitness is Copyright
2017-2019 Sergio Sanchez-Carrillo and 2019-2020 Juan Antonio Miguel
Gonzalez. `plot.stream` and `plot.stacked` are Copyright 2013-2016 by Marc
Taylor (see also https://github.com/marchtaylor/sinkr and
http://menugget.blogspot.com.es/2013/12/data-mountains-and-streams-stacked-area.html).


The code under `src/FitnessLandascape` is from [MAGELLAN, Maps of
Genetical
Landscapes](http://wwwabi.snv.jussieu.fr/public/magellan/Magellan.help.html).
The authors are S. Brouillet, G. Achaz, S. Matuszewski, H. Annoni, and
L. Ferreti. I downloaded the sources on 2019-06-05 from
http://wwwabi.snv.jussieu.fr/public/magellan/latest.tgz. The code is under
the GPLv3. MAGELLAN is "an integrated tool to visualize and analyze
fitness landscapes of small dimension (up to 7-8 loci)". In OncoSimulR we
use only a very limited subset of the functionality of MAGELLAN (mostly to
generate different types of random fitness landscapes and to compute
statistics of epistasis); the `Makevars` file we use only compiles two of
the executables (`fl_statistics` and `fl_generate`) ---the directory
`src/FitnessLandascape` contains, however, the complete sources. Note also
that the plots of fitness landscapes used in OncoSimulR are actually
blatantly copied in looks from MAGELLAN's plots. 


The code under `OncoSimulR/src/exprtk.h` is from [The C++ Mathematical
Expression Toolkit Library
(ExprTk)](http://www.partow.net/programming/exprtk/index.html). This code
is copyright Arash Partow, and is licensed under "The MIT License (MIT)"
(http://www.opensource.org/licenses/MIT) and is compatible with GPL
(http://directory.fsf.org/wiki/License:Expat). The file was originally
downloaded from http://www.partow.net/programming/exprtk/index.html on
2017-05-15. The most recent version was downloaded again on 2022-09-12,
from the [exprTk repo](https://github.com/ArashPartow/exprtk).<!-- at
commit --> <!--
https://github.com/ArashPartow/exprtk/commit/9fad72832c70348725c073e369a3321781001766). -->
The file was originally named `exprtk.hpp`; to conform to R's
requirements, it was renamed as `exprtk.h`


The code in `OncoSimulR/src/multivariate_hypergeometric.cpp` and `OncoSimulR/src/multivariate_hypergeometric.h` is from [extraDistr](https://github.com/twolodzko/extraDistr). This code is copyright Tymoteusz Wolodzko, and is licensed under the GNU GPL 2.


The code in `miscell-files/randutils.h` is copyright Melissa E. O'Neill,
and is licensed under "The MIT License (MIT)" in the terms explained in
the file itself. This is a license that is
[compatible with the GPL](http://directory.fsf.org/wiki/License:Expat).
The file randutils.hpp was downloaded from
https://gist.github.com/imneme/540829265469e673d045 on 2015-06-20 and is
also referenced from the main article [Ease of Use without Loss of Power]
(http://www.pcg-random.org/posts/ease-of-use-without-loss-of-power.html). I
renamed it to randutils.h to conform to R's requirements (and changed the
`auto exit_func = hash(&_Exit);` line to keep R from complaining about the
Exit function). I had to disable usage of randutils for now, since I could
not get it to compile with gcc-4.6 (since version 3.3 of R, 
the official Rtools for Windows now support C++-11, so I might change 
this in the near future).



The file under gitinfo-hooks is Copyright 2011 Brent Longborough, is
part of the
[gitinfo package](https://www.ctan.org/tex-archive/macros/latex/contrib/gitinfo?lang=en),
and is under the LaTeX Project Public License 1.3, which is
[incompatible with the GPL](http://directory.fsf.org/wiki/License:LPPLv1.3a). Note
this file is not part of the OncoSimulR BioConductor package.


Functions `nem_transitive.reduction` and `nem_transitive.closure` are from the
`nem` package, which was last available for version 3.11 of BioConductor
(https://www.bioconductor.org/packages/3.11/bioc/html/nem.html). The authors of
the package are Holger Froehlich, Florian Markowetz, Achim Tresch, Theresa
Niederberger, Christian Bender, Matthias Maneck, Claudio Lottaz, Tim Beissbarth,
and the package


The files under miscell-files/AParramon_discrete_time are copyright
Alberto Parramon, unless otherwise specified. This is an implementation of
a discrete-time version of OncoSimulR.


# Software status


|             | Bioconductor (multiple platforms)   | Travis CI  (Linux)  | Appveyor (Windows)  |
| ------------- | ------------------- | ------------- | ---------------- |
| R CMD check   | <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/OncoSimulR/"><img border="0" src="http://bioconductor.org/shields/build/release/bioc/OncoSimulR.svg" alt="Build status"></a> (release)</br><a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/OncoSimulR/"><img border="0" src="http://bioconductor.org/shields/build/devel/bioc/OncoSimulR.svg" alt="Build status"></a> (devel) | <a href="https://travis-ci.org/rdiaz02/OncoSimul"><img src="https://travis-ci.org/rdiaz02/OncoSimul.svg?branch=master" alt="Build status"></a> | <a href="https://ci.appveyor.com/project/rdiaz02/OncoSimul"><img src="https://ci.appveyor.com/api/projects/status/github/rdiaz02/OncoSimul?branch=master&svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/github/rdiaz02/OncoSimul?branch=master"><img src="https://codecov.io/github/rdiaz02/OncoSimul/coverage.svg?branch=master" alt="Coverage Status"/></a>   |                  |

(Note: Appveyor and Travis can fail for reasons that have nothing to do with the
package, such as R not being downloaded correctly, etc. Look at the
details of each failure. Similarly, some of the errors in BioConductor,
specially in the development branch, can be caused, specially in Windows,
by some required packages not being yet available, often "car" and
_"igraph". <!--            Finally, some failures in Travis can be caused by timeouts in -->
<!-- the coverage tests ---tests that run in about 5 minutes on my laptop, but -->
<!-- can occasionally apparently take longer than the maximum time for jobs -->
<!-- under Travis. --> Again, look at the details of each failure.)

<!-- Based on https://raw.githubusercontent.com/Bioconductor-mirror/illuminaio/master/README.md -->

<!-- [![Travis-CI Build Status](https://travis-ci.org/rdiaz02/OncoSimul.svg?branch=master)](https://travis-ci.org/rdiaz02/OncoSimul) -->
<!-- [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rdiaz02/OncoSimul?branch=master&svg=true)](https://ci.appveyor.com/project/rdiaz02/OncoSimul) -->
<!-- [![codecov.io](https://codecov.io/github/rdiaz02/OncoSimul/coverage.svg?branch=master)](https://codecov.io/github/rdiaz02/OncoSimul?branch=master) -->

# Funding
Supported by: grant BFU2015-67302-R (MINECO/FEDER, EU) funded by MCIN/AEI/10.13039/501100011033 and by ERDF A way of making Europe to R. Diaz-Uriarte; grant PID2019-111256RB-I00 funded by MCIN/AEI/10.13039/501100011033 to R. Diaz-Uriarte; "Beca de Colaboración" at the Universidad Autónoma de Madrid from Spanish Ministry of Education, 2017-18, to S. Sánchez Carrillo; Comunidad de Madrid's
PEJ16/MED/AI-1709 and PEJ-2019-AI/BMD-13961 to R. Diaz-Uriarte. 

<p align="center">
<img src="logo-micin-aei-uefeder.png" alt="micin-aei-uefeder logo" width="200"/>.
</p>


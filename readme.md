# Overview

`mpqDist` is an R package for computing a family of phylogenetically-informed distances and making a interactive plots based on them.
The MPQ distances are a family of distance with a continuous tuning parameter which controls the extent to which the distance is influenced by basal vs. terminal phylogenetic structure. 
The package allows the user to give a template for a plot which will be filled in with the full spectrum of tuning parameters.
Given a template, the package uses `plotly` to create an interactive plot that allows the user to scan through versions of the plot filled in with all the members of the family.

# Installation

The package can be installed using `devtools::install_github('jfukuyama/mpqDist')`.

# Vignettes

Vignettes showing simulations and real data analysis are available in the `vignettes` folder.
Compiled html documents from those vignettes are available in `docs`.
At this time, the html documents need to be downloaded to your machine (github doesn't display them nicely), but they are self-contained and once they are downloaded you can see all the figures and the interactive plots that are created in the vignettes.

You can also install the package and run `vignette(package = "mpqDist")` to see the vignettes.
For example, one of the vignettes is called `gentry`, which can be accessed using `vignette("gentry", package = "mpqDist")`.

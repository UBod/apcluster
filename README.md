# APCluster - An R Package for Affinity Propagation Clustering
In order to make Affinity Propagation Clustering introduced by Frey and
Dueck (2007;
[DOI:10.1126/science.1136800](http://doi.org/10.1126/science.1136800))
accessible to a wider audience, we ported the Matlab code published by the authors to R.
The algorithms are largely
analogous to the Matlab code published by Frey and Dueck.
The package further provides leveraged affinity propagation and an
algorithm for exemplar-based agglomerative clustering that can also be
used to join clusters obtained from affinity propagation. Various
plotting functions are available for analyzing clustering results.

The package is maintained by Ulrich Bodenhofer. The package itself has grown
over the years in which multiple students have contributed
significant parts: Johannes Palme, Chrats Melkonian, Andreas Kothmeier, and
Nikola Kostic

## Installation

The package can be installed from
[CRAN](https://CRAN.R-project.org/package=apcluster). Therefore, the the simplest way to install the package is to enter
```
install.packages("apcluster")
```
into your R session. If, for what reason ever, you prefer to install the package manually, follow the instructions in the [user manual](https://cran.r-project.org/web/packages/apcluster/vignettes/rococo.pdf).


## User support

If you encounter any issues or if you have any question that might be of interest also for other users, before writing a private message to the package developers/maintainers, please create an issue in this repository and also consider posting to the [R-help Mailing List](https://stat.ethz.ch/mailman/listinfo/r-help) or on [StackOverflow](https://stackoverflow.com/). For other matters regarding the package, please contact the package author.

## Citing this package

If you use this package for research that is published later, you are kindly asked to cite it as follows:

- U. Bodenhofer, A. Kothmeier, and S. Hochreiter (2011). APCluster: an R package for affinity propagation clustering. *Bioinformatics* **27**:2463-2464. DOI: [10.1093/bioinformatics/btr406](http://doi.org/10.1093/bioinformatics/btr406)

Moreover, we insist that, any time you use/cite the package, you also cite the original paper in which affinity propagation has been introduced:

- B. J. Frey and D. Dueck (2007). Clustering by passing messages between data points. *Science* **315**:972-976. DOI: [10.1126/science.1136800](http://doi.org/10.1126/science.1136800)

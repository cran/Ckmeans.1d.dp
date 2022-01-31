---
title: The 'Ckmeans.1d.dp' R package
# bibliography: inst/REFERENCES.bib
---

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Ckmeans.1d.dp)](https://cran.r-project.org/package=Ckmeans.1d.dp)
[![CRAN_latest_release_date](https://www.r-pkg.org/badges/last-release/Ckmeans.1d.dp)](https://cran.r-project.org/package=Ckmeans.1d.dp)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/Ckmeans.1d.dp)](https://cran.r-project.org/package=Ckmeans.1d.dp)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/Ckmeans.1d.dp)](https://cran.r-project.org/package=Ckmeans.1d.dp)

### Overview

The package provides a powerful set of tools for fast, optimal, and reproducible univariate clustering by dynamic programming. It is practical to cluster millions of sample points into a few clusters in seconds using a single core on a typical desktop computer. It solves four types of problem, including univariate $k$-means, $k$-median, $k$-segments, and multi-channel weighted $k$-means. Dynamic programming is used to minimize the (weighted) sum of within-cluster distances using respective metrics. Its advantage over heuristic clustering in efficiency and accuracy is increasingly pronounced as the number of clusters $k$ increases. Weighted $k$-means can also optimally segment time series to perform peak calling. An auxiliary function generates histograms that are adaptive to patterns in data. The package was recently used to map dysregulated zones in genomes of 17 human cancer types (Song and Zhong, 2020).

### The main method

The Ckmeans.1d.dp algorithms cluster (weighted) univariate data given by a numeric vector $x$ into $k$ groups by dynamic programming (Wang and Song, 2011; Zhong, 2019; Song and Zhong, 2020). They guarantee the optimality of clustering---the total of within-cluster sums of squared distances is always minimized given the number of clusters $k$. In contrast, heuristic univariate clustering algorithms may be sub-optimal at large $k$ or inconsistent from run to run. As non-negative weights are supported for each point, the algorithm can also partition a time course using the time points as input and the values at each time point as weight. Utilizing the optimal clusters, a function can generate histograms adaptive to patterns in data.

Excluding the time for sorting $x$, the default weighted univariate clustering algorithm takes a runtime of $O(kn)$ (Song and Zhong, 2020), linear in both sample size $n$ and the number of clusters $k$, using a divide-and-conquer strategy based on a theoretical result on matrix search [(Aggarwal et al., 1987) <doi:10.1007/BF01840359>](https://doi.org/10.1007/BF01840359) implemented via a novel in-place search space reduction method (Song and Zhong, 2020). The space complexity is $O(kn)$. This method is numerically stable.

### When to use the package

As an alternative to popular heuristic clustering methods, this package provides functionality for (weighted) univariate clustering, segmentation, and peak calling with guaranteed optimality and efficiency.

An adaptive histogram based on optimal clusters is also recommended if an equal-bin-width histogram is inadequate to characterize clusters that vary in width.

### To download and install the package
```
install.packages("Ckmeans.1d.dp")
```

### Citing the package

Song M, Zhong H (2020). "Efficient weighted univariate clustering maps outstanding dysregulated genomic zones in human cancers." _Bioinformatics_, 36(20), 5027–5036. 
https://doi.org/10.1093/bioinformatics/btaa613

Wang H, Song M (2011). "Ckmeans.1d.dp: Optimal k-means clustering in one dimension by dynamic programming." _The R Journal_, 3(2), 29–33. https://doi.org/10.32614/RJ-2011-015

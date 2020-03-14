## ---- fig.width=6-------------------------------------------------------------
require("Ckmeans.1d.dp")
x <- c(rnorm(40, mean=-2, sd=0.3),
       rnorm(45, mean=1, sd=0.1),
       rnorm(70, mean=3, sd=0.2))
ahist(x, col="lightblue", sub=paste("n =", length(x)),
      col.stick="darkblue", lwd=2, xlim=c(-4,4),
      main="Example 1. Gaussian mixture model with 3 components\n(one bin per component)\nAdaptive histogram")

## ---- fig.width=6-------------------------------------------------------------
ahist(x, breaks=3, col="lightgreen", sub=paste("n =", length(x)),
      col.stick="forestgreen", lwd=2,
      main="Example 1. Regular histogram")

## ---- fig.width=6-------------------------------------------------------------
ahist(x, k=9, col="lavender", col.stick="navy",
      sub=paste("n =", length(x)), lwd=2,
      main="Example 2. Gaussian mixture model with 3 components\n(on average 3 bins per component)\nAdaptive histogram")

## ---- fig.width=6-------------------------------------------------------------
ahist(x, breaks=9, col="lightgreen", col.stick="forestgreen",
      sub=paste("n =", length(x)), lwd=2,
      main="Example 2. Regular histogram")

## ---- fig.show='hold', fig.width=6--------------------------------------------
data(DNase)
res <- Ckmeans.1d.dp(DNase$density)
kopt <- length(res$size)
ahist(res, data=DNase$density, col=rainbow(kopt), col.stick=rainbow(kopt)[res$cluster],
      sub=paste("n =", length(x)), border="transparent",
      xlab="Optical density of protein DNase",
      main="Example 3. Elisa assay of DNase in rat serum\nAdaptive histogram")

## ---- fig.show='hold', fig.width=6--------------------------------------------
ahist(DNase$density, breaks="Sturges", col="palegreen",
      add.sticks=TRUE, col.stick="darkgreen",
      main="Example 3. Elisa assay of DNase in rat serum\nRegular histogram (equal bin width)",
      xlab="Optical density of protein DNase")

## ---- fig.show='hold', fig.width=6--------------------------------------------
x <- c(1,1,1,1, 3,4,4, 6,6,6)
ahist(x, k=c(2,4), col="gray",
      lwd=2, lwd.stick=6, col.stick="chocolate",
      main="Example 4. Adaptive histogram of repetitive elements")
ahist(x, breaks=3, col="lightgreen",
      lwd=2, lwd.stick=6, col.stick="forestgreen",
      main="Example 4. Regular histogram")



\name{modsymmtest-package}
\alias{modsymmtest-package}
\alias{modsymmtest}
\docType{package}
\title{
  The Wilcoxon and Sign tests for symmetry about an unknown center 
}
\description{
  Provides functionality of testing whether a data set comes from a symmetric 
  distribution when the center of symmetry is unknown, including modified wilcoxon test 
  and modified sign test. In addition, sample size calculated is provided. 
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Jiaojiao Zhou, Xinyu Gao, Albert Vexler

Maintainer: Jiaojiao Zhou <woerge99@gmail.com>
}
\references{
  Vexler, A., Gao, X., & Zhou, J. (2023). How to implement signed-rank wilcox. test () type procedures when a center of symmetry is unknown. Computational Statistics & Data Analysis, 107746. 
  
  Gastwirth, J. L. (1971). On the sign test for symmetry. Journal of the American Statistical Association, 66(336), 821-823.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
    data("plasma.silicon")
    post <- plasma.silicon$postoperative 
    pre <- plasma.silicon$preoperative
    # post <- c(0.21,0.24,0.1,0.12,0.28,0.25,0.22,0.21,0.22,0.23,0.22,0.24,0.45,0.38,
    #           0.23,0.22,0.18,0.15,0.04,0.14,0.24,0.2,0.24,0.18,0.19,0.15,0.26,0.3,0.22,0.24)
    # pre <- c(0.15,0.13,0.39,0.2,0.39,0.42,0.24,0.18,0.26,0.12,0.1,0.11,0.19, 0.15,0.27,
    #          0.28,0.11,0.11,0.18,0.18,0.24,0.48,0.27,0.22,0.18,0.19,0.32,0.31,0.19,0.21)
    mod.symm.test(post, pre)
  }
}

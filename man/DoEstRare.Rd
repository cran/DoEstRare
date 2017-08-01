\name{DoEstRare}
\alias{DoEstRare}
\title{DoEstRare}

\description{
Rare variant association test comparing position density functions and mutation counts between cases and controls.
}
\usage{
DoEstRare(pheno, geno, position, genome.size,
          perm=NULL, alpha=NULL, c=NULL,
          autosomal=TRUE, gender=NULL)
}
\arguments{
      \item{pheno}{a numeric vector of phenotypes. Affected individuals are coded 1 and unaffected individuals are coded 0.}
      \item{geno}{a numeric matrix of genotypes (row: individual, column: variant). Genotypes are coded 0,1 or 2 corresponding to the number of minor alleles. }
      \item{position}{a numeric vector of variant positions. }
      \item{genome.size}{a numeric value corresponding to the length of the analyzed region. }
      \item{perm}{number of permutations. If not NULL,a "standard permutation procedure" is performed to compute the significance. See Details. }
      \item{alpha}{error level. If not NULL,an "adaptive permutation procedure" is performed to compute the significance. See Details.}
      \item{c}{precision of the p-value. If not NULL,an "adaptive permutation procedure" is performed to compute the significance. See Details.}
      \item{autosomal}{boolean. If TRUE, autosomal chromosome; FALSE, X chromosome.}
       \item{gender}{numeric vector. 1=male; 2=female.}
}
\value{
	\item{p.value}{the p-value obtained by the phenotype permutation procedure. }
	\item{stat}{the test statistic. }

}

\details{
Two types of permutations procedures can be defined in the function: the standard permutation procedure and the adaptive permutation procedure.

In the standard permutation procedure, the user specifies, in the argument "perm", the number of permutations to be done. The p-value will be \eqn{(R+1)(B+1)}. With \eqn{R} the number of permutation statistics superior to the observed statistic and \eqn{B} the number of permutations.

In the adaptive permutation procedure, the user specifies, in the argument "alpha", the significance to achieve after multiple testing correction. In the argument "c", the estimation precision of the p-value. In function of these two paremeters, the maximal number of permutations and the maximal number of success to achieve will be computed. If the maximal number of success is reached, the p-value will be \eqn{R/B}. If not, the maximal number of permutations will be used to compute the p-value \eqn{(R+1)(B+1)}.
}

\references{
Persyn E, Karakachoff M, Le Scouarnec S, Le Cl??zio C, Campion D, French Exome Consortium, et al. DoEstRare: A statistical test to identify local enrichments in rare genomic variants associated with disease. Wang K, editor. PLOS ONE. 2017 Jul 24;12(7):e0179364.

Che R, Jack JR, Motsinger-Reif AA, Brown CC. An adaptive permutation approach for genome-wide association study: evaluation and recommendations for use. BioData Min. 2014;7:9.

}

\author{Elodie Persyn, \email{elodie.persyn@univ-nantes.fr}}

\examples{
pheno=rep(c(0,1), 500)
geno=matrix(sample(c(0,1),prob=c(0.7,0.3) ,1000*30, replace=TRUE), ncol=30)
position=sample(1:500,30)
genome.size=500
perm=200

#Autosomal gene
#standard phenotype permutation procedure
DoEstRare(pheno, geno, position, genome.size, perm)
#adaptive phenotype permutation procedure
DoEstRare(pheno, geno, position, genome.size, alpha=0.05, c=0.2)

#X gene
gender=rep(c(1,2), each=500)
#standard phenotype permutation procedure
DoEstRare(pheno, geno, position, genome.size, perm, autosomal=FALSE, gender=gender)
#adaptive phenotype permutation procedure
DoEstRare(pheno, geno, position, genome.size, alpha=0.05, c=0.2, autosomal=FALSE, gender=gender)


}

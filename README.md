# TLegene
# Introduction
TLegene is a R procedure for borrowing the idea of transfer learning to integrate useful genetic information available from external studies for multilocus-based eGene identification method. In TLegene, the identification of eGene consists of two components: the first component represents the indirect influence of the auxiliary study after transfer learning, and the second component represents the direct effect of the target study.

Specifically, let e be a n by 1 vector of gene expression level on n individuals in the target study, X is a n by p matrix for covariates, G is a n by m matrix for genotypes of cis-SNPs for for a given gene in the target study, θ quantifies the association between the gene expression level and the weighted genetic score Gγ which is 
the indirect effect of auxiliary study, b quantifies the association between gene expression level and genotypes G which is the direct effect not completely interpreted 
by auxiliary data and α quantifies a p-vector of fixed effect sizes for clinical covariates. We relate e, Gγ, X and G by a linear mixed model:
<p align="center">
e= Xα + (Gγ) × θ + Gb, , b ~ N(0, τ)
</p>
Above, τ is the genetic variance which is the direct effect not completely interpreted by auxiliary data.

TLegene examines the association of G and Gγ with e (while controlling for X) by testing for:
<p align="center">
H0: θ = 0 and b = 0 <==> H0: θ = 0 and τ = 0
</p>
This is a joint test which requires simultaneously assessing the significance of both fixed effects and random effects: the first part of H0 evaluates the indirect influence of auxiliary samples, whereas the second part assesses the direct impact of target samples.

# R2FixedEffectsGLMM

Code from R2 for Fixed Effects in the GLMM by Jaeger et. al., 2016

These files can be used to generate the results in my paper, 'An R^2 Statistic for Fixed Effects in the GLMM.'

Glimmix_R2_V3 provides a SAS macro that extends R^2_beta to the GLMM using penalized quasi-likelihood (PQL) Estimation. This SAS macro computes model and semi-partial R squared statistics for each fixed effect in the model. Due to a lack of the non-central beta distribution in SAS, confidence limits are not currently supplied. The user may specify denominator degrees of freedom estimation. Kenward-Roger (DDF = KR) is recommended, and is compatible with PQL estimation. 

R2_Dental_Dat_Sim provides code in R to reproduce the results from Jaeger et al., 2016

R2_Sim_Cont_Cnt_Bin provides code in R (but requires results from SAS) to reproduce results from Jaeger et al., 2016

R2_Sim_Cont_Cnt_Bin_Simple provides code in R that doesn't reproduce results from Jaeger et al., but is much easier to run.

Rsq_Simlongdat is the code in SAS that produces the results which must be read into the R2_Sim_Cont_Cnt_Bin simulation.

For more on R^2_\beta, Please see 

Edwards et al., 2008 (http://onlinelibrary.wiley.com/doi/10.1002/sim.3429/abstract) . 

Jaeger et al.,  2016 (http://www.tandfonline.com/doi/full/10.1080/02664763.2016.1193725#.V1lzNpErJhE)


# r2glmm

An R package for computation of model R squared and semi-partial 
R squared (with confidence limits) in linear and generalized linear mixed models

This package computes model and semi partial R squared and 
for the linear and generalized linear mixed model (LMM and GLMM). 
The R squared measure from Edwards et.al (2008) is extended
to the GLMM using penalized quasi-likelihood (PQL) estimation 
(see Jaeger et al. 2016). Three methods of computation are provided:

(1) The Kenward-Roger approach (see note A),

(2) The method introduced by Nakagawa and Schielzeth (2013) (see note B), and

(3) an approach using standardized generalized variance (SGV)
that can be used for both mean model and covariance model selection.

Confidence limits for semi partial R squared and model R squared are
computed for each of the methods listed.

Instructions for installation:

After installing the devtools package, run this code in the console:

devtools::install_github('bcjaeger/r2glmm')

(A.) Due to some inconsistency between the pbkrtest package and the glmmPQL
function, the Kenward-Roger approach in the r2glmm package is limited to
the LMM.


(B.) The r2glmm package only computes marginal R squared for the LMM and
does not generalize the statistic to the GLMM; however, confidence limits
may be computed for this marginal R squared in the LMM using this package


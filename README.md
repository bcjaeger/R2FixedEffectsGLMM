# R2FixedEffectsGLMM

Code from R2 for Fixed Effects in the GLMM by Jaeger et. al., 2016

These files can be used to generate the results in my paper, 'An R^2 Statistic for Fixed Effects in the GLMM.'

Glimmix_R2_V3 provides a SAS macro that extends R^2_beta to the GLMM using penalized quasi-likelihood (PQL) Estimation.
  The macro also 
   (1) computes semi-partial R squared statistics for each fixed effect in the model.
   (2) allows the user to specify denominator degrees of freedom estimation. Kenward-Roger (DDF = KR) is recommended.

R2_Dental_Dat_Sim provides code in R to reproduce the results from Jaeger et al., 2016

R2_Sim_Cont_Cnt_Bin provides code in R (but requires results from SAS) to reproduce results from Jaeger et al., 2016

R2_Sim_Cont_Cnt_Bin_Simple provides code in R that doesn't reproduce results from Jaeger et al., but is much easier to run.

Rsq_Simlongdat is the code in SAS that produces the results which must be read into the R2_Sim_Cont_Cnt_Bin simulation.

For more on R^2_\beta, Please see 

Edwards et al., 2008 (http://onlinelibrary.wiley.com/doi/10.1002/sim.3429/abstract) . 

Jaeger et al.,  2016 (http://www.tandfonline.com/doi/full/10.1080/02664763.2016.1193725#.V1lzNpErJhE)

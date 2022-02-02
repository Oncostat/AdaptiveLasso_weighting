# AdaptiveLasso_weighting


Source code for manuscript "Accounting for grouped predictor variables or pathways in high-dimensional penalized Cox regression models"
Shaima Belhechmi, Riccardo De Bin, Federico Rotolo & Stefan Michiels.

To reproduce the results presented in the article, run the files in the following folder "./functions":
- simdata.R: function implementing the data simulation method.
- runsims.R: function generating the data sets according to the 8 scenarios presented in the article.

- AL.R : function implementing the different weighting strategies proposed for the Adaptive Lasso method.
- analysisAL.R: function applying the AL.R function on the simulated data sets.

- SGL.R, gel.R, cMCP.R, ipflasso.R and ipflasso2.R: functions implementing the different methods SGL, gel, cMCP and IPFlasso with two penalty factors.
- analysisSGL.R, analysisGEL.R, analysiscMCP.R, analysisIPFL.R and analysisIPFL2.R: functions applying the SGL.R, gel.R, cMCP.R, ipflasso.R and ipflasso2.R functions on the simulated datasets.

For questions, comments or remarks about the code please contact S. Michiels (stefan.michiels@gustaveroussy.fr) or S. Belhechmi (belhechmishaima@gmail.com).


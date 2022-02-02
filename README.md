# AdaptiveLasso_weighting


Source code for manuscript "Accounting for grouped predictor variables or pathways in high-dimensional penalized Cox regression models"
Shaima Belhechmi, Riccardo De Bin, Federico Rotolo & Stefan Michiels.

For questions, comments or remarks about the code please contact S. Michiels (stefan.michiels@gustaveroussy.fr) or S. Belhechmi (belhechmishaima@gmail.com).

The code has been written using R. 

To reproduce the results presented in the article, execute the files in the following "./functions" folder:
1- simdata.R: the R function implementing the data simulation method.
2- runsims.R: the R function generating the data sets according to the 8 scenarios presented in the article.
3- AL.R: the R function implementing the different weighting strategies proposed for the Adaptive Lasso method.
4- analysisAL.R: the R function applying the AL.R function on the simulated datasets.
5- SGL.R, gel.R, cMCP.R, ipflasso.R and ipflasso2.R: the R functions implementing the different methods SGL, gel, cMCP and IPFlasso with two penalty factors.
6- analysisSGL.R, analysisGEL.R, analysiscMCP.R, analysisIPFL.R and analysisIPFL2.R: the R functions applying the SGL.R, gel.R, cMCP.R, ipflasso.R and ipflasso2.R functions on the simulated datasets.

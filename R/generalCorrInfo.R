#' generalCorr package description: 
#' 
#' This package provides convenient software tools for causal path determinations
#' using Vinod (2014, 2015, 2018, 2021) and is explained in many package vignettes.
#' \code{causeSummary(mtx)}, \code{causeSummary2(mtx)},\code{causeSum2Blk(mtx)},
#' \code{causeSummBlk} are various versions reporting pair-wise causal 
#' path directions and causal strengths. We fit
#' a kernel regression of X1 on (X2, X3,..Xk) and another flipped regression
#' of X2 on (X1, x3, ..Xk).  We compare the two fits using three sophisticated criteria
#' called Cr1 to Cr3. We rescale the
#' weighted sum of the quantified three criteria to the [-100, 100] range.
#' The sign of the weighted sum gives the direction of the causal path, and 
#' the magnitude of the weighted sum gives the strength of the causal path.
#' A matrix of non-symmetric generalized correlations r*(x|y) is reported by the
#' functions \code{rstar()} and \code{gmcmtx0()}. 
#' \code{sudoCoefParcor()} computes pseudo kernel regression coefficients based on
#' generalized partial correlation coefficients (GPCC)
#' \code{depMeas()} a measure of nonlinear nonparametric dependence between two vectors.
#' \code{parcorVec()} has generalized partial correlation coefficients, Vinod (2021)
#' \code{parcorVecH()} has a hybrid version of the above (using HGPCC).
#' The usual partial correlations r(x,y|z) for regression of y on (x, z) measure
#' the effect of y on x after removing the effect of z, where z can have several variables.
#' Vinod (2021) suggests new generalized partial correlation coefficients (GPCC)
#' using kernel regressions, r*(x,y|z).
#' 
#' The criterion Cr1 uses observable values of standard exogeneity test criterion,
#' namely, (kernel regression residual) times (regressor values)
#' Cr2 computes absolute values kernel regression residuals.
#' The quantification of Cr1 and Cr2 further uses four orders of stochastic 
#' dominance measures.
#' Cr3 compares the R-square of the two fits.
#' The package provides additional tools for matrix algebra, such as 
#' \code{cofactor()}, for outlier detection \code{get0outlier()}, 
#' for numerical integration by the trapezoidal rule, stochastic dominance
#' \code{stochdom2()} and \code{comp_portfo2()}, etc.
#' The package has a function \code{pcause()} for bootstrap-based statistical 
#' inference and another one 
#' for a heuristic t-test called \code{heurist()}.  Pairwise deletion of missing data
#' is done in \code{napair()}, while triplet-wise deletion is in \code{naTriplet()}
#' intended for use when control variable(s) are also present. If one has
#' panel data, functions \code{PanelLag()} and \code{Panel2Lag()} are relevant.
#' \code{pillar3D} provides 3-dimensional plots of data that look
#' more like surfaces, than usual plots with vertical pins.
#' 
#' Recent 2020 additions include \code{canonRho()} for generalized canonical 
#' correlations, and many 
#' functions for Granger causality between lagged time series including
#' \code{GcRsqX12()}, \code{bootGcRsq()} and \code{GcRsqYXc()}.
#' 
#' Recent 2021 additions include several functions for portfolio choice.
#' \code{sudoCoefParcor()} for pseudo regression coefficients for kernel regressions.
#' \code{decileVote()}, \code{momentVote()}, \code{exactSdMtx()} for exact
#' computation of stochastic dominance from ECDF areas. The newer stochastic
#' dominance tools are used in \code{causeSummary2(mtx)},\code{causeSum2Blk(mtx)}
#' \code{dif4mtx()}
#' computes growth, change in growth etc. up-to order 4 differencing of time series.
#'   
#' @note Six vignettes provided with this package at CRAN
#' describe the theory and usage of the package with examples. Read them using
#' the command: 
#' \code{vignette("generalCorr-vignette")} to read the first vignette. 
#' vignettes 2 to 6 can be read by including the vignette number. For 
#' example, 
#' \code{vignette("generalCorr-vignette6")} to read the sixth vignette.
#'
#' 
#' @references Vinod, H. D.'Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics' in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \doi{gffn86} 
#'  
#' @references Vinod, H. D. 'Matrix Algebra Topics in Statistics and Economics
#' Using R', Chapter 4 in 'Handbook of Statistics: Computational Statistics
#' with R', Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' 
#' @references Zheng, S., Shi, N.-Z., and Zhang, Z. (2012). 'Generalized measures 
#'  of correlation for asymmetry, nonlinearity, and beyond,' 
#'  Journal of the American Statistical Association, vol. 107, pp. 1239-1252.
#'  
#' @references Vinod, H. D. (2021) 'Generalized, Partial and Canonical Correlation
#' Coefficients' Computational Economics, 59(1), 1--28.
#'   
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=2982128}
#' 
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2019, pp. 33-64.    
#'  
#' @docType package  
#' @name generalCorrInfo
NULL
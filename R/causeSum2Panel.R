#' Kernel regressions based causal paths in Panel Data.
#' 
#' We assume that panel data have space (space=individual region)
#' and time (e.g., year) dimensions. We use upper case
#' X to denote a common prefix in the panel data.
#' Xs =name of the space variable, e.g., state or individual.
#' The range of values for s is 1 to nspace.
#' Xt =name of the time variable, e.g., year.
#' The range of values for t is 1 to ntime.
#' Xy =the dependent variable(s) value at time t in state s.
#' Since panel data causal analysis can take a long computer time,
#' we allow the user to choose subsets of time and space values
#' called chosenTimes and chosenSpaces, respectively.
#' Various input parameters starting with "nam" specify the names
#' of variables in the panel study.
#' 
#' The algorithm calls some function fn(mtx) where mtx is the data matrix,
#' and fn is causeSummary2NoP(mtx). The causal paths between (y, xj)
#' pairs of variables in mtx are computed following 3 sophisticated
#' criteria involving exact stochastic dominance. Type "?causeSummary2"
#' on the R console to get details (omitted here for brevity).
#' Panel data consist of a time
#' series of cross-sections and are also called longitudinal data.
#' We provide estimates of causal path directions and strengths
#' for both the time-series and cross-sectional
#' views of panel data. Since our regressions are kernel type
#' with no functional forms, fixed effects for time and space
#' are being suppressed when computing the causality.
#' 
#' @param da {panel dat having a named column for space and time}
#' @param fn {an R function causeSummary2NoP(mtx)}
#' @param rowfnout {the number of rows output by fn}
#' @param colfnout {the number of columns output by fn}
#' @param fnoutNames {the column names of output by fn, for example,
#' fnoutNames=c("cause","effect","strength","r","p-val")}
#' @param namXs {title of the column in da having the space variable}
#' @param namXt {title of the column in da having the time variable}
#' @param namXy {title of the column in da having the dependent y variable}
#' @param namXc {title(s) of the column(s) in da 
#' having control variable(s), default=0 means none specified}
#' @param namXjmtx {title(s) of the column(s) in da having regressor(s)}
#' @param chosenTimes {subset of values of time variable chosen for quick
#' results, There are NchosenTimes values chosen in the subset.
#' default=NULL means all time identifiers in the data are included.}
#' @param chosenSpaces {subset of values of space variable chosen for quick
#' results, There are NchosenSpaces values chosen in the subset.
#' default=NULL means all space identifiers are included. The degrees
#' of freedom for Studentized statistic for Granger causality tests 
#' are df=(NchosenSpaces -1).}
#' @param ylag {time lag in Granger causality study of time dimension
#'  the default ylag=0 is not really zero. It means ylag=
#'  min(4, round(NchosenTimes/5,0)),
#'  where NchosenTimes is the length of chosenTimes vector}
#' @param verbo {print detail results along the way, default=FALSE} 
#' 
#' @description
#' The algorithm of this function uses an internal function
#' fminmax=function(x){min(x)==max(x)}.
#' The subsets mtx2 of the original data da for a specific time or space
#' can become degenerate if the columns of mtx2 have no variability.
#' The apply function of R is applied to the columns of mtx2 as follows.
#' "ap1=apply(mtx2,2,fminmax)." Now, "sumap1=sum(ap1)" counts how many
#' columns of the data matrix are degenerate.  We have a degeneracy
#' problem only if sumap1 is >1 or =1. For example, the panel consists
#' of data on 50 United States and 20 years. Now, consumer price index
#' (cpi) data may be common for all states. That is, the min(cpi)
#' equals max(cpi) for all states. Then, the variance of cpi is zero,
#' and we have degeneracy. When this happens, the regressor cpi should
#' not be involved in determining causal paths.
#' We identify degeneracy using "fminmax=function(x){min(x)==max(x)}"
#' 
#' @returns The causeSum2Panel(.) produces many output matrices and vectors. 
#' The first
#' "outt" gives a 3-dimensional array of panel causal path output focused on 
#' time series for each space value using fixed space value.
#' It reports causal path directions, and strengths for (y, xj) pairs.
#' The second output array, called
#' "outs", gives similar 3D panel causal path output focused on 
#' space cross sections using fixed time value.
#' The third output matrix called
#' "outdif" gives causal paths using Granger causality for each
#' pair (y, xj). They are not causal strengths but differences 
#' between Rsquare values of two flipped kernel regressions.
#' The summary of Granger causality answer is an output matrix called
#' grangerAns (first row average of differences in R-squares and
#' second row has its test statistic with degrees of freedom n-1),
#' and grangerStat for related t-statistic for formal inference.
#' based on column means and variances of "outdif". This function
#' also produces a matrix summarizing "outt" and "outs" into two-dimensional
#' matrices reporting averages of signed strengths as  "strentime"
#' and "strenspace", Also, "pearsontime" reports the Pearson
#' correlation coefficients for various time values and their average in the
#' last column. It determines the overall direction of
#' the causal relation between y and xj.
#' For example, a negative average correlation means y and xj are negatively
#' correlated (xj goes up, y goes down). Similarly, "pearsonspace"
#' summarizes "outs" correlations.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY.
#' @seealso See  \code{\link{causeSummary2}}
#' @seealso See  \code{\link{causeSummary}} is subject to trapezoidal approximation. 
#' @references Vinod, H. D. 'Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics' in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \doi{10.1080/03610918.2015.1122048}
#' @references Vinod, H. D. 'New exogeneity tests and causal paths,'
#'  Chapter 2 in 'Handbook of Statistics: Conceptual Econometrics 
#' Using R', Vol.32, co-editors: H. D. Vinod and C.R. Rao. New York:
#' North Holland, Elsevier Science Publishers, 2019, pp. 33-64.
#'  
#' @references Vinod, H. D. Causal Paths and Exogeneity Tests 
#' in {Generalcorr} Package for Air Pollution and Monetary Policy 
#' (June 6, 2017). Available at SSRN: 
#' \url{https://www.ssrn.com/abstract=2982128}    
#' 
#' @references Vinod, Hrishikesh D., R Package GeneralCorr 
#' Functions for Portfolio Choice 
#' (November 11, 2021). Available at SSRN: 
#' https://ssrn.com/abstract=3961683 
#' 
#' @references Vinod, Hrishikesh D., Stochastic Dominance 
#' Without Tears (January 26, 2021). Available at 
#' SSRN: https://ssrn.com/abstract=3773309 
#'  
#' @concept  Granger causality 
#' @concept stochastic dominance orders
#' @concept summary index
#' @note The function prints to the screen some summaries of the three
#' output matrices. It reports how often a variable is a cause in
#' various pairs as time series or as cross sections. It also reports the
#' average strengths of causal paths for "outt" and "outs" matrices.
#' We compute the difference between two R-square values to find which
#' causal direction is more plausible.  This involves kernel regressions
#' of y on its own lags and lags of a regressor. Unlike the usual Granger
#' causality we estimate better-fitting nonlinear kernel regressions.
#' If the averages in "outdif" matrix are negative, the Granger causal
#' paths go from y to xj. This may be unexpected when the model assumes
#' that y depends on x1 to xp, that is, the causal paths go from xj to y.
#' In studying the causal pairs, the function creates mixtures of
#' names y and xj. Character vectors containing the mixed names are
#' are column names or row names depending on the context. For example,
#' the output matrix grangerAns column names help identify the relevant
#' regressor name. The first row of the grangerAns matrix has column averages
#' of outdiff matrix to help get an overall estimate of the Granger-causal paths.
#' The second row of the grangerAns has the Studentized test
#' statistic for formal testing of the significance of Granger causal paths.
#' Collecting the results for the time dimension strengths with
#' suitable sign (negative strength means cause reversal xj->y) is output
#' named strentime.  The corresponding Pearson
#' correlations as an output is named pearsontime.
#' Collecting the results for the space dimension strengths with
#' suitable sign (negative strength means cause reversal xj->y) is output
#' named strenspace.  The corresponding Pearson
#' correlations are named pearsonspace. A grand summary of average
#' strengths and correlations is output matrix named grandsum.
#' It is intended to provide an overall picture of causal paths in Panel
#' data. These paths should not be confused with Granger causal paths which
#' always involve time lags and causes are presumed to precede effects in time.
#' 
#' @examples
#'
#'
#' \dontrun{
#' library(plm);data(Grunfeld)
#' options(np.messages=FALSE)
#' namXs="firm"
#' print("initial values identifying the space variable")
#' head(da[,namXs],3)
#' print(str(da[,namXs]))
#' chosenSpaces=(3:10)                        
#' if(is.numeric(da[,namXs])){
#'   chosenSpaces=as.numeric(chosenSpaces)}
#' if(!is.numeric(da[,namXs])){
#'   chosenSpaces=as.character(chosenSpaces)}
#' 
#' namXt="year"
#' print("initial values identifying the time variable")
#' head(da[,namXt],3)
#' print(str(da[,namXt]))
#' chosenTimes=1940:1949
#' if(is.numeric(da[,namXt])){
#'   chosenTimes=as.numeric(chosenTimes)}
#' if(!is.numeric(da[,namXt])){
#'   chosenTimes=as.character(chosenTimes)}
#' 
#' namXy="inv"
#' namXc=0
#' namXjmtx=c("value","capital")
#' p=length(namXjmtx)
#' fn=causeSummary2NoP
#' fnout=matrix(NA,nrow=p,ncol=5)
#' fnoutNames=c("cause","effect","strength","r","p-val")
#' causeSum2Panel(da, fn=causeSummary2NoP,
#'                rowfnout=p, colfnout=5, 
#'                fnoutNames=c("cause","effect","strength","r","p-val"),
#'                namXs=namXs,
#'                namXt=namXt,
#'                namXy=namXy,
#'                namXc=namXc,
#'                namXjmtx=namXjmtx,
#'                chosenTimes=chosenTimes,
#'                chosenSpaces=chosenSpaces,
#'                verbo=FALSE)
#' }
#' 
#' 
#' @export

causeSum2Panel=function(da, fn=causeSummary2NoP, rowfnout, colfnout, fnoutNames,
                        namXs,
                        namXt,
                        namXy,
                        namXc=0,
                        namXjmtx,
                        chosenTimes=NULL,
                        chosenSpaces=NULL,
                        ylag=0,
                        verbo=FALSE)
{
  Xs=da[,namXs]
  uniXs=unique(Xs)
  nspace=length(uniXs)
  if(length(chosenSpaces)==0) chosenSpaces=uniXs #select All
  Xt=da[,namXt]
  uniXt=unique(Xt)
  if(length(chosenTimes)==0) chosenTimes=uniXt
  ntime=length(uniXt)
  Xy=da[,namXy]
  if(namXc != 0) Xcmtx=da[,namXc]
  Xjmtx=da[,namXjmtx]
  p=length(namXjmtx)
  if (verbo){
    print(c("number of regressors in data=",p))
    print(c("number of (states, individuals) spaces=",nspace))
    print(c("number of times=",ntime))
  }#end if verbo
  
  NchosenSpaces=length(chosenSpaces)
  NchosenTimes=length(chosenTimes)
  if(verbo){
    print("how many spaces and times chosen for analysis")  
    print(c(NchosenSpaces,NchosenTimes))
  }#endif verbo
  
  if (ylag> ntime/3) print("Granger Warning: ylag>ntime/3")
  if(ylag==0) ylag=min(4, round(NchosenTimes/5,0)) 
  print(c("Granger Causality: max time lag used=",ylag))
  
  #now define the space range for analysis
  spaces=1:NchosenSpaces
  
  times=1:NchosenTimes
  
  #create mixed names with y and each regressor
  
  mixnam=rep(NA, p)
  for (j in 1:p){
    mixnam[j]=paste("y-",namXjmtx[j],sep="")
  } #end j loop
  mixnam2s=mixnam #initialize if no degeneracy is found
  mixnam2t=mixnam #initialize if no degeneracy is found
  mixnam0=mixnam  #original mix names
  
  if(verbo){print("Causal paths with fixed-time cross sections")
    print("mixed names of y and regressors")
    print(mixnam)  }
  # places to store output for various space values
  outs=array(NA,c(rowfnout,colfnout,NchosenSpaces))
  dimnames(outs)=list(mixnam2s,fnoutNames,chosenSpaces)
  outt=array(NA,c(rowfnout,colfnout,NchosenTimes))
  dimnames(outt)=list(mixnam2t,fnoutNames,chosenTimes)
  outdif=matrix(NA,nrow=NchosenSpaces,ncol=p) 
  
  
  nammtx=c(namXy, namXjmtx)
  if(verbo)print("Next: Granger type time-lag Causality, fix Xs space value ")
  for (i in 1:NchosenSpaces){
    # i=1 is time series for the first chosen space (state) entity
    mtx=subset(da,da[,namXs]==chosenSpaces[i])
    mtx2=mtx[,nammtx] #select data columns for y x1 x2 .. xp
    if (verbo){
      print(head(mtx,2))
      print(tail(mtx,2))
     print(head(mtx2,2)) #selected data columns for y x1 x2 .. xp
    print(tail(mtx2,2))
    }#endif verbo
    fminmax=function(x){min(x,na.rm=TRUE)==max(x,na.rm=TRUE)} #checking degenerate data
    ap0=apply(mtx2,2,fminmax) # this uses names y x1 x2 .. xp
    ap1=ap0[-1] #oimt the first for y
  #  print(c("ap0",ap0))                ################
  #  print(c("ap1",ap1))                ##################
    sumap1=sum(ap1)
    if(verbo)print(c("how many variables have zero variability min=max",sumap1))
    if(verbo)print(c("which variable has zero variability",nammtx[ap0]))
    if(sumap1 >= 1) {
      mixnam[ap1]="NA"
      mixnam2s=mixnam
 #     dimnames(outs)=list(mixnam2s,fnoutNames,chosenSpaces)
      if(verbo) print(c("variable has no variation",nammtx[ap1]))
    } #endif sumap1>=1, mixnam2s allows for degenerate data
    
    nammtx2=setdiff(nammtx, nammtx[ap0])
  #  nRHS=length(nammtx2)
    nRHS=sum(!is.na(nammtx2))
    mtx3=mtx2[,nammtx2]
    fir=fn(mtx3)
 #   if(i==1) { print(" output of fn for the first chosenSpace")

  #    print(fir)
 #     print(" outs matrix designed to store above output: should match dimensions")
 #     print(outs[c(1:(p-sumap1)),,i])
#    }#endif i==1
    outs[c(1:(p-sumap1)),,i]=fir
 #  print("outs")
 #   print(fir) ################## ##################
    
    if(verbo)print("call Granger Causality function: GcRsqX12")
    for (j in 2:nRHS){
      dif1j=GcRsqX12(mtx3[,1],mtx3[,j],px1=ylag,px2=ylag,pwanted=ylag)[3]
    #  print(dif1j) #third value is difference between Rsquare values
      outdif[i,(j-1)]=as.numeric(dif1j)  
    }#  end j loop for Granger causality
  }# end the for loop for i-th chosen space/indivi/state
  
  colnames(outdif)=mixnam
  if(verbo)print("Next:Cross-section Causality between variables fix time unit")
  
  #begin loop for i over chosen time units for cross sections at each time
  for (i in 1:NchosenTimes){
    mtx=subset(da,da[,namXt]==chosenTimes[i])
    mtx2=mtx[,nammtx]
    if (verbo){
      print(head(mtx,2))
      print(tail(mtx,2))
      print(head(mtx2,2))
      print(tail(mtx2,2))
    }#endif verbo
    fminmax=function(x){min(x)==max(x)}
    ap0=apply(mtx2,2,fminmax)#how many variables have zero variability min=max
    ap1=ap0[-1]
  #  print(c("ap0",ap0))                ################
  #  print(c("ap1",ap1))                ################## 
   sumap1=sum(ap1)
     if(verbo)print(c("how many variables have zero variability min=max",sumap1))
    if(verbo)print(c("which variable has zero variability",nammtx[ap0]))
   
    if(sumap1 >= 1) {
      mixnam[ap1]="NA"
      mixnam2t=mixnam
 #     dimnames(outt)=list(mixnam2t,fnoutNames,chosenTimes)
      print(c("variable(s) have no variation",nammtx[ap0])) } #endif
    nammtx2=setdiff(nammtx, nammtx[ap0])
 #   print(nammtx2) # after omitting the degenerate regressor ############
    mtx3=mtx2[,nammtx2]
 #     print(head(mtx3,2))   ##################
 #     print(tail(mtx3,2))   ##################
    fir=fn(mtx3)
    outt[c(1:(p-sumap1)),,i]=fir
  }# end for loop over i
 # print("outt fir")
 # print(fir) ##################  ##################
  print("How frequent is a variable designated as cause in outs, dim= space")
  mtxScause=outs[,1,]
  print(table(mtxScause))
  
  mtxstren=outs[,3,]
  stravg=rep(NA, p) #average strength of causal relation
  names(stravg)=mixnam0 #mixture of two names at a time
  for ( j in 1:p){
    stre=mean(as.numeric(mtxstren[j,]),na.rm=TRUE)
    stravg[j]=stre
  }
  print("outs:Average causal strengths for mixed names defining variable pairs")
  names(stravg)=mixnam2s
 # print(stravg)
  
 print("How frequent is a variable designated cause in outt?")
  mtxTcause=outt[,1,]
  print(table(mtxTcause))
  print("Cause in Granger causality from following sign(s)")
  print(mixnam0)
  grangerStat=rep(NA,p)
  grangerSD=rep(NA,p)
  grangerAns0=apply(outdif,2,mean,na.rm=TRUE)
  print(grangerAns0)
  print("negative sign means y causes xj by Granger method")
  if (nspace>5){grangerSD=apply(outdif,2,sd,na.rm=TRUE)
  grangerStat=sqrt(nspace)*grangerAns0/grangerSD
  } #endif 
  #  print(c("p=",p)) ###########
#  print(mixnam2t) ######
  mtxstren=outt[,3,]
#  print(head(mtxstren),2)   #############
  stravg=rep(NA, p) #average  Granger strength of causal relation
  names(stravg)=mixnam #mixture of two names at a time
  for ( j in 1:p){
    stre=mean(as.numeric(mtxstren[j,]),na.rm=TRUE)
    stravg[j]=stre
  }
 # print(mixnam2t) ######
  names(stravg)=mixnam2t
print("Granger Average causal strengths for mixed names defining variable pairs")
  print(stravg)
  ######
  ###### begin outt summarization
  ######
  print("Causality for various times next")
  strent=outt[,3,] #third col. has abs(strength) values needs sign info
  pearsont=outt[,4,]# 4-th column has pearson correlations
  dim1=dim(outt)[1]
  dim3=dim(outt)[3]
  dimn1=dimnames(outt)[[1]]
  dimn3=dimnames(outt)[[3]]
#  print("c(dim1,dim3,dimn1,dimn3)")
#  print(c(dim1,dim3,dimn1,dimn3))
  strent2=matrix(NA,nrow=dim1, ncol=dim3) #place holder
  pearsont2=matrix(NA,nrow=dim1, ncol=dim3) #place holder
  causeT=outt[,1,] #If first column cause==dependent variable, minus strength
  for (j in 1:dim3){
  for (i in 1:dim1){
  strent2[i,j]=as.numeric(strent[i,j])
  pearsont2[i,j]=as.numeric(pearsont[i,j])
  }}
  #print("pearson time dimension next")
  #print(pearsont2)
  strenT=strent2 #place to store
  for (j in 1:dim3){
  for (i in 1:dim1){
  strenT[i,j]=as.numeric(strent2[i,j])
##If first column cause==dependent variable namXy then change sign of strength
  if (!is.na(causeT[i,j])){
  if (causeT[i,j]==namXy) strenT[i,j]=-as.numeric(strent2[i,j])
  }}}
rownames(strenT)=dimn1
colnames(strenT)=dimn3
rownames(pearsont2)=dimn1
colnames(pearsont2)=dimn3
 # print(strenT)
  avgstrent=apply(strenT,1,mean,na.rm=TRUE)
  avgPearsont=apply(pearsont2,1,mean,na.rm=TRUE)
  strentime=cbind(strenT,avgstrent)
  pearsontime=cbind(pearsont2,avgPearsont)
colnames(strentime)=c(dimn3,"avg")
colnames(pearsontime)=c(dimn3,"avg")
print("strentime or causal strengths by time dimention")
print(strentime)
print("pearsontime or Pearson correlations by time Last col. has average")
print(pearsontime)
print("Above Causality for various times (averaging over space)")
######
###### begin outs summarization wrt s=space
######
print("Causality for various spaces next (averaging over time)")
strens=outs[,3,]
pearsons=outs[,4,]#
dim1=dim(outs)[1]
dim3=dim(outs)[3]
dimn1=dimnames(outs)[[1]]
dimn3=dimnames(outs)[[3]]
#print("c(dim1,dim3,dimn1,dimn3)")
#print(c(dim1,dim3,dimn1,dimn3))
strens2=matrix(NA,nrow=dim1, ncol=dim3) #place holder
pearsons2=matrix(NA,nrow=dim1, ncol=dim3) #place holder
causeS=outs[,1,]
  for (j in 1:dim3){
  for (i in 1:dim1){
    strens2[i,j]=as.numeric(strens[i,j])
    pearsons2[i,j]=as.numeric(pearsons[i,j])
  }}
#print("pearsons2 next")
#print(pearsons2)
strenS=strens2 #place to store
for (j in 1:dim3){
  for (i in 1:dim1){
    strenS[i,j]=as.numeric(strens2[i,j])
#  print(c(strens[i,j],strens2[i,j],namXy))
  if (!is.na(causeS[i,j])){
    if (causeS[i,j]==namXy) strenS[i,j]=-as.numeric(strens2[i,j])
  }}}
rownames(strenS)=dimn1
colnames(strenS)=dimn3
#print(strenS)
avgstrens=apply(strenS,1,mean,na.rm=TRUE)
avgPearsons=apply(pearsons2,1,mean,na.rm=TRUE)
strenspace=cbind(strenS,avgstrens)
pearsonspace=cbind(pearsons2,avgPearsons)
colnames(strenspace)=c(dimn3,"avg")
colnames(pearsonspace)=c(dimn3,"avg")
rownames(pearsonspace)=dimn1
print("strenspace: causal strengths by space dimention, last col.has average")
print(strenspace)
print("pearsonspace or Pearson correlations by space dimension")
print(pearsonspace)
print("report the grand summary in 4 columns")
strtime=(strentime[,NCOL(strentime)])
cortime=(pearsontime[,NCOL(pearsontime)])

strspace=(strenspace[,NCOL(strenspace)])
corspace=(pearsonspace[,NCOL(pearsonspace)])

grandsum=cbind(strtime,cortime,strspace,corspace)
print(grandsum)

print("Nonlinear Granger Causality average of difference in two R-squares")
print(grangerAns0)
print("Granger Causality test statistic mean of differences in two R-squares")
print(grangerStat)
grangerAns=rbind(grangerAns0,grangerStat)
rownames(grangerAns)=c("mean of differeces in Rsq","test statistic")
print(grangerAns)
list(outt=outt, outs=outs, outdif=outdif, 
       grangerAns=grangerAns,
       strentime=strentime,strenspace=strenspace,
       pearsontime=pearsontime, pearsonspace=pearsonspace,
       grandsum=grandsum)
} #end function


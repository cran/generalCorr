#' Kernel regressions based causal paths in Panel Data.
#' 
#' We assume that panel data have space (space=individual region)
#' and time (e.g., year) dimension, and We use upper case
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
#' and fn is causeSummary2(mtx). The causal paths between (y, xj)
#' pairs of variables in mtx are computed following 3 sophisticated
#' criteria involving exact stochastic dominance. type "?causeSummary2"
#' on the R console to get details (omitted here for brevity).
#' Panel data consist of a time
#' series of cross-sections and are also called longitudinal data.
#' We provide estimates of causal path directions and strengths
#' for both the time-series and cross-sectional
#' views of panel data. Since our regressions are kernel type
#' with not functional forms fixed effect for time and space
#' are being suppressed when computing the causality.
#' 
#' @param da {panel dat having a named column for space and time}
#' @param fn {an R function causeSummary2(mtx)}
#' @param rowfnout {the number of rows output by fn}
#' @param colfnout {the number of columns output by fn}
#' @param fnoutNames {the column names of output by fn, for example,
#' fnoutNames=c("cause","effect","strength","r","p-val")}
#' @param namXs {title of the column in da having the space variable}
#' @param namXt {title of the column in da having the time variable}
#' @param namXy {title of the column in da having the dependent y variable}
#' @param namXc {title(s) of the column(s) in da having control variable(s)}
#' @param namXjmtx {title(s) of the column(s) in da having regressor(s)}
#' @param chosenTimes {subset of values of time variable chosen for quick
#' results, There are NchosenTimes values chosen in the subset.
#' default=NULL means all time identifiers in the data are included.}
#' @param chosenSpaces {subset of values of space variable chosen for quick
#' results, There are NchosenSpaces values chosen in the subset.
#' default=NULL means all space identifiers are included.}
#' @param ylag {time lag in Granger causality study of time dimension
#'  the default ylag=0 means ylag=min(4, round(NchosenTimes/5,0))
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
#' of data on 50 united states and 20 years. Now consumer price index
#' (cpi) data may be common for all states. That is, the min(cpi)
#' equals max(cpi) for all states. Then the variance of cpi is zero,
#' and we have degeneracy. When this happens, the regressor cpi should
#' not be involved in determining causal paths.
#' We identify degeneracy using "fminmax=function(x){min(x)==max(x)}"
#' 
#' @returns The causeSum2Panel() produces three output matrices. The first
#' "outt" gives panel causal path output focused on 
#' time series for each space value using fixed space value.
#' It reports causal path directions, and strengths for (y, xj) pairs.
#' The second output matrix called
#' "outs" gives panel causal path output focused on 
#' space cross sections using fixed time value.
#' The third output matrix called
#' "outdiff" gives causal paths using Granger causality for each
#' pair (y, xj).
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
#' If the averages in "outdiff" matrix are negative, the Granger causal
#' paths go from y to xj. This may be unexpected when the model assumes
#' that y depends on x1 to xp,and One expects the causal paths from xj to y.
#' In studying the causal pairs, the function creates mixtures of
#' names y and xj. Character vectors containing the mixed names used
#' are also a part of the output.  The vectors are named mixnam2t for
#' time series and mixnam2s for cross-section portion of panel data.
#' The two sets of mixed names are needed because degeneracies on the
#' time dimension may be the same as degeneracies on the space dimension.
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
#' fn=causeSummary2
#' fnout=matrix(NA,nrow=p,ncol=5)
#' fnoutNames=c("cause","effect","strength","r","p-val")
#' causeSum2Panel(da, fn=causeSummary2,
#'                rowfnout=p, colfnout=5, 
#'                fnoutNames=c("cause","effect","strength","r","p-val"),
#'                namXs=namXs,
#'                namXt=namXt,
#'                namXy=namXy,
#'                namXc=namXc,
#'                namXjmtx=namXjmtx,
#'                chosenTimes=chosenTimes,
#'                chosenSpaces=chosenSpaces,
#'                verbo=TRUE)
#' }
#' 
#' 
#' @export

causeSum2Panel=function(da, fn, rowfnout, colfnout, fnoutNames,
                        namXs,
                        namXt,
                        namXy,
                        namXc=NULL,
                        namXjmtx,
                        chosenTimes=NULL,
                        chosenSpaces=NULL,
                        ylag=0,
                        verbo=FALSE)
{
  Xs=da[,namXs]
  uniXs=unique(Xs)
  nspace=length(uniXs)
  Xt=da[,namXt]
  uniXt=unique(Xt)
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
  print(c("max time lag in Granger Causality=",ylag))
  
  #now define the space range for analysis
  if (NchosenSpaces==0)  spaces=1:nspace
  if (NchosenSpaces!=0)  spaces=1:NchosenSpaces
  
  #option to have all as chosen times by setiing zero
  if (NchosenTimes==0)  times=1:ntime
  if (NchosenTimes!=0)  times=1:NchosenTimes
  
  #create mixed names with y and each regressor
  
  mixnam=rep(NA, p)
  for (j in 1:p){
    mixnam[j]=paste(namXy,namXjmtx[j],sep="")
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
    fminmax=function(x){min(x)==max(x)} #checking degenerate data
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
    nRHS=length(nammtx2)
    mtx3=mtx2[,nammtx2]
    fir=causeSummary2(mtx3)
    if(i==1) { print(" output of fn for first chosenSpace")
      print(fir)
 #     print(" outs matrix designed to store above output: should match dimensions")
 #     print(outs[c(1:(p-sumap1)),,i])
    }#endif i==1
    outs[c(1:(p-sumap1)),,i]=fir
    
    
    if(verbo)print("call Granger Causality function: GcRsqX12")
    for (j in 2:nRHS){
      dif1j=GcRsqX12(mtx3[,1],mtx3[,j],px1=ylag,px2=ylag,pwanted=ylag)[3]
    #  print(dif1j)
      outdif[i,(j-1)]=as.numeric(dif1j)  
    }#  end j loop for Granger causality
  }# end the for loop for i-th chosen space/indivi/state
  
  
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
    fir=causeSummary2(mtx3)
    outt[c(1:(p-sumap1)),,i]=fir
  }# end for loop over i
  print("How frequent is a variable designated as cause: select cross-section space")
  mtxScause=outs[,1,]
  print(table(mtxScause))
  
  mtxstren=outs[,3,]
  stravg=rep(NA, p) #average strength of causal relation
  names(stravg)=mixnam0 #mixture of two names at a time
  for ( j in 1:p){
    stre=mean(as.numeric(mtxstren[j,]),na.rm=TRUE)
    stravg[j]=stre
  }
  print("outs: Average causal strengths for mixed names defining variable pairs")
  names(stravg)=mixnam2s
  print(stravg)
  
   print("How frequent is a variable designated as cause: selecting time series")
  mtxTcause=outt[,1,]
  print(table(mtxTcause))
  print("Cause in Granger causality from following sign(s)")
  print(mixnam0)
  print(apply(outdif,2,mean,na.rm=TRUE))
  print("negative sign means y causes xj by Granger method")
#  print(c("p=",p)) ###########
#  print(mixnam2t) ######
  mtxstren=outt[,3,]
#  print(head(mtxstren),2)   #############
  stravg=rep(NA, p) #average strength of causal relation
  names(stravg)=mixnam #mixture of two names at a time
  for ( j in 1:p){
    stre=mean(as.numeric(mtxstren[j,]),na.rm=TRUE)
    stravg[j]=stre
  }
 # print(mixnam2t) ######
  names(stravg)=mixnam2t
  print("outt: Average causal strengths for mixed names defining variable pairs")
  print(stravg)

  list(outt=outt, outs=outs, outdif=outdif, 
       mixnam2s=mixnam2s, mixnam2t=mixnam2t)
} #end function


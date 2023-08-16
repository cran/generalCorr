#' Compare out-of-sample portfolio choice algorithms by a leave-percent-out method.
#' 
#' This function randomly leaves out 5 percent (`pctOut'=5 by default) 
#' data and finds portfolio choice by seven different
#' portfolio selection algorithms using the data on the remaining 95 percent (say). 
#' The randomization removes any bias in time series definitions of `out-of-sample' data.
#' For example, the input to \code{outOFsamp(.)} named `mtx' is a matrix with 
#' p columns for p stocks and n returns. Also, let  the maximum number of
#' stocks admitted to belong in the portfolio be four, or `maxChosen=4'.
#' Now \code{outOFsamp} function computes the returns earned by the
#' seven portfolio selection algorithms, called
#' "SD1", "SD2", "SD3", "SD4", "SDAll4", "decile," and "moment," where SDAll4 refers
#' to a weighted sum of SD1 to SD4 algorithms. Each algorithm provides
#' a choice ranking of p stocks with choice values 1,2,3,..,p where stock ranked
#' 1 should get the highest portfolio weight.
#' The \code{outOFsamp} function then calls the
#' function `rank2return,' which uses these rank choice numbers to the selected
#' `maxChosen' stocks.  The allocation is linearly declining. For example, it is
#' 4/10, 3/10, 2/10, and 1/10, with the top choice stock receiving 4/10 of the capital.
#' Each choice of `pctOut' rows of the `mtx' data yields an outOFsamp return for each
#' of the seven portfolio selection algorithms.  These outOFsamp return
#' computations are repeated \code{reps} times. 
#' A new random selection of `pctOut' rows (must be 2 or more) of data is made
#' for each repetition. We set 
#' reps=20 by default. The low default is set
#' to save processing time in early phases, but we recommend reps=100+. 
#'  The final choice of stock-picking algorithm out of seven
#' is suggested by the one yielding largest average out-of-sample 
#' return over the `reps' repetitions.`Its standard deviation
#' measures the variability of performance over the reps repititions.
#' 
#' @param mtx {matrix size n by p of data on n returns from p stocks}
#' @param pctOut {percent of n randomly chosen rows left out as out-of-sample, default=5
#'  percent. One must leave out at least two rows of data}
#' @param reps {number of random repetitions of left-out rows over which we average
#'   the out-of-sample performance of a stock-picking algorithm,
#'   default reps=20}
#' @param seed {seed for random number generation, default =23}
#' @param verbo {logical, TRUE means print details, default=FALSE}
#' @param maxChosen {number of stocks (out of p) with nonzero weights in the portfolio}
#' @importFrom utils head tail
#' @return a matrix called `avgRet' with seven columns for seven stock-picking
#' algorithms "SD1","SD2","SD3","SD4","SDAll4","decile",and "moment," containing
#' out-of-sample average returns for linearly declining allocation in a portfolio.
#' User needs to change rank2return() for alternate portfolio allocations.
#' @note The traditional time-series out-of-sample leaves out the last few
#' time periods, and estimates the stock-picking model using part of the data
#' time periods. The pandemic of 2019 has revealed that the traditional
#' out-of-sample would have a severe bias in favor of pessimistic stock-picking
#' algorithms.  The traditional method is fundamentally flawed since it is
#' sensitive to the trends (ups and downs) in the out-of-sample period. The
#' method proposed here is free from such biases. The stock-picking algorithm
#' recommended by our outOFsamp() is claimed to be robust against such biases.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{\link{rank2return}}
#' @examples
#' \dontrun{
#' x1=c(2,5,6,9,13,18,21,5,11,14,4,7,12,13,6,3,8,1,15,2,10,9)
#' x2=c(3,6,9,12,14,19,27,9,11,2,3,8,1,6,15,10,13,14,5,7,4,12)
#' x3=c(2,6,NA,11,13,25,25,11,9,10,12,6,4,3,2,1,7,8,5,15,14,13)
#' mtx=cbind(x1,x2,x3)
#' mtx=mtx[complete.cases(mtx),]
#' os=outOFsamp(mtx,verbo=FALSE,maxChosen=2, reps=3)
#' apply(os,2,mean)}
#' @export

outOFsamp=function(mtx, pctOut=5, reps=10, seed=23, 
          maxChosen=2, verbo=FALSE){
n=NROW(mtx)
p=NCOL(mtx)
if (pctOut>50) {print("percent left-out exceeds 50, reset as 50")
  pctOut=50}
if (maxChosen>p){ print("outOFsamp function maxChosen>p in portfolio")
print("maxChosen is reset to p")
maxChosen=p}#endif >p
if(verbo){
print(c("outOFsamp mtx: n,p",n,p))
print(c("maxChosen",maxChosen)) }#endif verbo
# select pctOut% of n observations for out-of-sample
n5=max(1,round(n*pctOut/100,0)) #rounded integer to 0 digits for out-of-sample
if(n5==1) n5=2   #leave-out only one fails (a limitation of R syntax).
insamp=n-n5 #number of in-sample oservations
if(verbo){
print(c("nobs in each replicate",insamp))
print(c("leave-out",n5,"observations")) }
#reps number of times we repeat leaving out n5 observations
avgRet=matrix(NA,nrow=reps, ncol=7)#average returns place holder
# first five columns are sd1 to sd4, and their average
# sixth col. is for average of decile ranks
# seventh col. is for 4 moments and Sharpe Ratio based average rank
# get myrank for each criterion for each column of mtx
# Begin replication loop
seed=seed+1
set.seed(seed)

for (irep in (1:reps)){
if(verbo) print(c("replication begins",irep))
#randomize  for out-of-sample
sampN=sample(1:n,replace=FALSE)  #random vector 
myin=sampN[1:insamp] #select the first set for insamp
myout=sampN[(insamp+1):n]
#print(c("myout",myout))
if(length(myout)<1) stop("myout set has too few rows for testing")
esti=mtx[myin,] #chosen subset of rows for estimation
testi=mtx[as.numeric(myout), ]#remaining n5 rows for out-of-sample test
#print("testi")
#print(testi)
if (verbo) {print("dimensions of testi matrix")
print(c(NROW(testi),NCOL(testi)))
print(myout)}
colnames(esti)=colnames(mtx)
if(n5>1) colnames(testi)=colnames(mtx)
  if(irep==1){
if(verbo)  print(c("n5 rows for out-of-samp testiing",head(testi)))} 
  if(irep==1) print(c("head:testi[,1:3]",head(testi[,1:3],2)))
  if(irep==1) print(c("tail:testi[,1:3]",tail(testi[,1:3],2))) 
#if(irep==1) print(esti)
#using chosen rows from bigger mtx find myrank vector
ex1=exactSdMtx(esti)
ex2=summaryRank(ex1$out)
#having estimated all SD1 to SD4 now use testing rows called testi
choice1=ex2[5,]  #fifth row has the choice results for SD1
if (NROW(testi)<1) stop("testing subset rows <1")
if (NCOL(testi)<p) stop("testing subset columns <p")
r1= rank2return(testi,maxChosen=maxChosen, 
    myrank=choice1,verbo=verbo)  
avgRet[irep,1]=r1   #first col of avgRet matrix has SD1 results 

choice2=ex2[6,]  #sixth row has the choice results for SD2
r2= rank2return(testi,maxChosen=maxChosen, 
                myrank=choice2,verbo=verbo)  
avgRet[irep,2]=r2   #second col of avgRet matrix has SD2 results 

choice3=ex2[7,]  #7-th row has the choice results for SD3
r3= rank2return(testi,maxChosen=maxChosen, 
                myrank=choice3,verbo=verbo)  
avgRet[irep,3]=r3   #3-rd col of avgRet matrix has SD3 results 

choice4=ex2[8,]  #f8-th row has the choice results for SD4
r4= rank2return(testi,maxChosen=maxChosen, 
                myrank=choice4,verbo=verbo)  
avgRet[irep,4]=r4   #4-th col of avgRet matrix has SD4 results 

choiceAll4=ex2[10,]  #10-th row has the choice results for sum of all SD1 to SD4
rAll4= rank2return(testi,maxChosen=maxChosen, 
                myrank=choiceAll4,verbo=verbo)  
avgRet[irep,5]=rAll4   #5-th col of avgRet matrix has sum SD1 to SD4 results 

d1=decileVote(esti)
dd1=d1$out
choiceDecile=dd1[NROW(dd1),]
#print(choiceDecile)
rdecile= rank2return(testi,maxChosen=maxChosen, 
                myrank=choiceDecile,verbo=verbo)  
avgRet[irep,6]=rdecile   #6-th col of avgRet matrix has decile results 

m1=momentVote(esti)
mm1=NROW(m1)
choiceMoment=m1[mm1,1:p]
#print("choiceMoment")
#  print(choiceMoment)
rmoment= rank2return(testi,maxChosen=maxChosen, 
                     myrank=choiceMoment,verbo=verbo)  
avgRet[irep,7]=rmoment   #7-th col of avgRet matrix has moment results 


if (irep%%20==1){
print(c("replication No. above",irep)) }
if(irep==1) print(c(r1,r2,r3,r4,rAll4,rdecile,rmoment))
}#end irep loop
#print(c("average return over reps",irep ))
colnames(avgRet)=c("SD1","SD2","SD3","SD4","SDAll4","decile","moment")
if(verbo){
  print(c("comparing average returns over",NROW(avgRet),"out-of-sample rows"))
  am=apply(avgRet,2,summary)
  print(am)} #endif verbo average returns over replications
return(avgRet)
}# end outOFsamp function
#example


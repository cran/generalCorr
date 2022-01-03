#' Function compares Pearson Stats and Sharpe Ratio for a matrix of stock returns
#' 
#' The first step computes mean, std.dev, skewness, kurtosis (kurt),and
#' the Sharpe Ratio (mean/sd) representing risk-adjusted return where sd measures
#' the risk. The input x must be a matrix having p columns (col.names recommended).
#' and n rows as in the data.  If data are missing for some columns, insert NA's.
#' Thus x has p column of data matrix ready for comparison and ranking. 
#' For example, x has a matrix of stock returns.
#' The output matrix produced by this function has p columns for each data
#' column (i.e. for each stock being compared). The output matrix has
#' twelve rows. Top five rows have the magnitudes of 
#' mean, sd, skew, kurt, Sharpe ratios.  Output matrix 
#' rows 6 to 10 have respective ranks of moment stats. 
#' The output 11-th row  reports a weighted sum
#' of ranks with following weights mean=1,sd=-1,skew=0.5,kurt=-0.5,Sharpe Ratio=1.
#' User has the option to change the weights. They measure relative importance.
#' 
#' 
#' Since skewness and kurtosis are measured relatively less
#' reliably (have greater sampling variation due to higher powers) their weight
#' is 0.5. Our ranking gives the smallest number 1 to the most desirable outcome.  
#' The 11-th line of the
#' output matrix has weighted sum of ranks and we suggest higher portfolio weight
#' be given to the column having smallest value (in the bottom line).
#' The 12-th row of output matrix has `choice,' where input weights give
#' the number 1 is for the top choice column of data and all other choice numbers.
#' The (p+1)-th column of the output matrix has the chosen weights.  The argument
#' weight to the `momentVote' function allows one to change these weights.
#' 
#' @param mtx {n by p matrix of data, For example, n stock returns
#' for p stocks. The mtx columns
#' should have some names (ticker symbols)}
#' @param weight {vector of reliability weights. default: mean=1, sd=1,
#'  skew=0.5,kurt=0.5,sharpe=1}
#' @return a matrix with same number of columns as in the input matrix x and
#' eleven rows. Top five rows have moment quantities, next five are their ranks
#' the eleventh row has weighted sum of ranks with the input weights (see default)
#' and the 12-th row has choice numbers (choice=1 is best)
### @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
### @seealso \code{\link{}}
#' @examples
#' x1=c(1,4,7,2,6)
#' x2=c(3,4,8,4,7)
#' momentVote(cbind(x1,x2))
#' 
#' @export

momentVote=function(mtx, weight=c(1,-1,0.5,-0.5,1)){
  n=NROW(mtx)
  p=NCOL(mtx)
  skew=function(x){x=x[!is.na(x)]
  dev=x-mean(x)
  num=sum(dev^3)/length(x)
  var=sum(dev^2)/length(x)
  rat=num/var^(3/2)
  return(rat)
  }
  kurt=function(x){x=x[!is.na(x)]
  dev=x-mean(x)
  num=sum(dev^4)/length(x)
  var=sum(dev^2)/length(x)
  den=var^2
  return(num/den)
  }
  #  if(verbo) print(c("n,p",n,p))
  MN=apply(mtx,2,mean,na.rm=TRUE)
  SD=apply(mtx,2,sd,na.rm=TRUE)

  SK=apply(mtx,2,skew)
#  print(c("SK",SK))
  KU=apply(mtx,2,kurt)  
  SH=rep(NA,p) #sharpe ratio
  for (j in 1:p) if (SD[j]!=0) SH[j]=MN[j]/SD[j]
  #rank function in R gives rank=1 to smallest magnitude
  #most desirable rank=1 is for the largest mean, skewness and Sharpe Ratio
  #but to the lowest sd and lowest kurt
  #so we subtract rank by R from p+1 as in following line
  RMN=(p+1)-rank(MN,na.last = TRUE,  ties.method = "average")
  RSD=rank(SD, na.last = TRUE,  ties.method = "average")
  RSK=(p+1)-rank(SK,na.last = TRUE,  ties.method = "average")
  RKU=rank(KU,na.last = TRUE,  ties.method = "average")
  RSH=(p+1)-rank(SH,na.last = TRUE,  ties.method = "average")
  out=matrix(NA,ncol=p,nrow=10)
  #first 5 rows have magnitudes, next 5 have ranks
  out[1,]=MN
  out[2,]=SD
  out[3,]=SK
  out[4,]=KU 
  #compute Sharpe ratio=(mean/sd) only if denominator is nonzero
  out[5,]=SH
 
  out[6,]=RMN
  out[7,]=RSD
  #compute Sharpe ratio=mean/sd only if denominator is nonzero
  out[8,]=RSK
  out[9,]=RKU 
  out[10,]=RSH
  wtsum=t(t(out[6:10,])%*%weight)
  ord1=rank(wtsum,na.last=TRUE)
  #ord1=sort(wtsum,index.return=TRUE)
 # print(ord1)
  choice=ord1 
 # print(wtsum)
 # print(choice)
  weigh=c(weight,weight,NA,NA)
  ou1=rbind(out,wtsum,choice)
  out2=cbind(ou1,weigh)
rownames(out2)=c("mean", "sd", "skewness","kurtosis","Sharpe Ratio",
"Rank.mean", "Rank.sd", "Rank.skewness","Rnk.kurtosis","Rank.ShRatio",
"wtedSumRanks","choice")
colnames(out2)=c(colnames(mtx),"smpReliWt")
    return(out2)
}#end momentVote function


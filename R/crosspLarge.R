crosspLarg = function(x){
n = NROW(x)
if(NCOL(x)!=n) stop("NCOL does not equal NROW in x")
y=apply(x,2,as.character)
for (i in 1:n){
for(j in 1:(i)){
  if (j==i) next
if(abs(x[i,j])>abs(x[j,i])) y[i,j]=paste(y[i,j],"L",sep="")
  else  y[j,i]=paste(y[j,i],"L",sep="")
}#end j
}# end i
return(y) 
}#end function
#Example
x=matrix(1:16,nrow=4)
x
crosspLarg(x)

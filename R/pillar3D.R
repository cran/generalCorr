#' Create a 3D pillar chart to display (x, y, z) data coordinate surface.
#' 
#' Give data on x, y, z coordinate values of a 3D surface, this function
#' plots them after making pillars near each z value by adding and
#' subtracting small amounts dz. Instead of pins of the height z this
#' creates pillars which better resemble a surface. It uses the
#' wireframe() function of `lattice' package to do the plotting.
#' 
#' For additional plotting features type `pillar3D()' on the R console to
#' get my code and adjust wireframe() function defaults.
#' 
#' @param x {x-coordinate values}
#' @param y {y-coordinate values}
#' @param z {z-coordinate values}
#' @param drape {logical value, default drape=TRUE to give color to heights}
#' @param xlab {default "x" label on the x axis}
#' @param ylab {default "y" label on the y axis}
#' @param zlab {default "z" label on the z axis}
#' @param mymain {default "Pillar Chart" main label on the plot}
#' @importFrom lattice wireframe
#' @return A 3D plot 
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @concept  3D plot
#' @concept wireframe plot
#' @examples
#' 

#' \dontrun{
#' pillar3D())}
#' 
#' @export

pillar3D=function(z=c(657,  936, 1111, 1201),
x=c(280,  542,  722, 1168), 
y=c(162, 214, 186, 246), drape=TRUE,
xlab="y", ylab="x", zlab="z", mymain="Pillar Chart") {
##require(lattice)
nz=length(z)
nx=length(x)
ny=length(y)
if(ny != nz) stop("elements in y and z are not equal")
if(nx != nz) stop("elements in x and z are not equal")
ns=3*nz+2
minz=min(z)
maxz=max(z)
dz=(maxz-minz)/(ns-2)
newz=matrix(0,ns,ns) #initialize
for  (i in 1:nz)#tim=ti minus
{ti=3*i;tim1=ti-1
tip1=ti+1 #tip=ti plus
newz[tim1,tim1]=z[i]-dz
newz[tim1,ti]=z[i]-dz
newz[tim1,tip1]=z[i]-dz
newz[ti,tim1]=z[i]
newz[ti,ti]=z[i]
newz[ti,tip1]=z[i]
newz[tip1,tim1]=z[i]+dz
newz[tip1,ti]=z[i]+dz
newz[tip1,tip1]=z[i]+dz
}
wireframe(newz,screen = list(z = 1, x = -10),drape=drape,
  xlab=xlab, ylab=ylab, zlab=zlab, main=mymain)
}

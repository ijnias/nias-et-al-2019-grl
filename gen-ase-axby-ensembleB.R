# Compute an ensemble with members that have fields like
# C <- C_mb * p[1]*exp(sigma_lnC) [basal friction coefficient]
# D <- D_mb * p[2]*exp(sigma_lnD) [viscosity factor]
# M <- M_mb * p[3]*exp(sigma_lnM) [melt rate]
#
# Create data files used for each of these parameters
# based on a perturbation obtained from Latin hypercube sampling
#
# Create configuration files required to run BISICLES
# based on input.template

require(libamrfile)
require(lhs)

amr.free.all()
#eneumerate components 

#from the *coef* file
CONE <- 0
CTHIRD <- 1
MUCOEF <- 2
MELT <- 3

#from the *balance* file
GLMELT <- 0
AMBMELT <- 1


amr.aXbYplusc <- function(X,Y,a,b=function(x,y,comp){1-a(x,y,comp)},c=function(x,y,comp){0},nghost)
{

  #given two AMR files X and Y , produce an amr hierarchy one whose components are
  #a(x,y,comp) * X + a(b,y,comp) * Y + c(x,y,comp) 
  
  #ASSUMES X,Y HAVE THE SAME DISJOINT BOX LAYOUT!
   xID <- amr.load(X)
   yID <- amr.load(Y)

   maxlev <- amr.query.nlevel(xID) - 1
   maxlevy <- amr.query.nlevel(yID) - 1
   if (maxlev != maxlevy)
     {
       stop ("amr.axby assumes identical meshes but level count differs")
     }

   for (lev in 0:maxlev)
     {
       nfab <- amr.query.nfab(xID,lev) - 1
       nfaby <- amr.query.nfab(yID,lev) - 1
       if (nfab != nfaby)
         {
           stop ("amr.axby assumes identical meshes but fab count differs")
         }

       ncomp <- amr.query.ncomp(xID, lev)
       ncompy <- amr.query.ncomp(yID, lev)
       if (ncomp != ncompy)
         {
           stop ("amr.axby assumes same number of components")
         }

       for (ifab in 0:nfab)
         {
           for (icomp in 0:(ncomp-1))
             {
               #print("lev,ifab,icomp");print(lev); print(ifab); print(icomp)
               
               fab <- amr.read.fab(xID,lev,ifab,icomp,ng=nghost)
               faby <- amr.read.fab(yID,lev,ifab,icomp,ng=nghost)
               #print(names(fab))
               print(dim(fab$v))
               A <-  a(fab$x,fab$y,icomp)
               print("A");print(A)
               data <- A * fab$v +
                 b(fab$x,fab$y,icomp) * faby$v +
                   c(fab$x,fab$y,icomp)
               
               
               amr.write.fab(xID, lev, ifab, data,icomp , ng=nghost)
               
             }
         }
     }
   
   #print(xID)
   #ncomp <- amr.query.ncomp(xID, 0)
   #print(ncomp)
   #print(yID)
   #ncomp <- amr.query.ncomp(xID, 0)
   #print(ncomp)
   #amr.free(yID)
   #print(xID)
   #ncomp <- amr.query.ncomp(xID, 0)
   #print(ncomp)
   xID
 }




apply_pert <- function(p,name)
  {
    system(paste("mkdir -p ensembleB/" , name,sep=""))

    coeffile <- paste("ase-coef-",name,"-1km.2d.hdf5",sep="")
    balfile <- paste("ase-balance-melt-",name,"-1km.2d.hdf5",sep="")

    fullname <- function(name,tagdomain,bflaw,geometry){paste(name,tagdomain,bflaw,geometry,sep="-")}
    
    pathto <- function(file)
      {
        paste("ensembleB/",name,"/",file,sep="")
      }
 
    sedargs <- function(tagdomain,bflaw,geometry)
      {
        #contruct the arguments for sed commands that will
        #replace placeholders (starting with @) in the
        #input template file
        
        
        bfcoef <- bflaw
        bfpower <- ""
        if (bflaw == "Cone") { bfpower <- "1.0"}
        if (bflaw == "Cthird") { bfpower <- "0.333"}
        if (bfpower == ""){stop("unknown basal friction law ")}
       
        geofile <- ""
        if (geometry == "mb") { geofile <- "ase-geometry-mb-1km.2d.hdf5"}
        if (geometry == "zero") { geofile <- "ase-geometry-zero-1km.2d.hdf5"}    
        if (geofile == ""){stop("unknown geometry selected ")}
        
        a <- rbind(c("@NAME",name),
                   c("@FULLNAME",fullname(name,tagdomain,bflaw,geometry)),
                   c("@TAGDOMAIN",tagdomain),
                   c("@BFCOEF",bfcoef),
                   c("@BFPOWER",bfpower),
                   c("@GEOFILE",geofile),
                   c("@GEOMETRY",geometry))
        s <- NULL
        for (i in 1:dim(a)[1])
          {
            s <- paste(s,"-e s:",a[i,1],":",a[i,2],": ",sep="")
          }
        s
      }

    for (tagdomain in c("thg","pig","ase")){
      for (bflaw in c("Cone","Cthird")){
        for (geometry in c("mb","zero")){
        inputfile <- pathto(paste("inputs.",fullname(name,tagdomain,bflaw,geometry),sep=""))
        cmd <- paste("sed", sedargs(tagdomain,bflaw,geometry) ," ensembleB/inputs.template > " ,inputfile)
        system(cmd)
               } 
      }
    }

    print(p)
    
    #p[1,2,3] define C (basal traction coefs Cone * Cthird) & D (muCoef) & melt
    #map 0 to C/2 and 1 to 2*C etc
    acoef <- function(x,y,comp){
      print("acoef");  print(p);print(comp)
      r <- 0.0
      if ( (comp == CONE) | (comp == CTHIRD) )
        {
          r <- 2.0 ^ (2 * p[1]- 1) 
        }
      if (comp == MUCOEF )
        {
          r <- 2.0 ^ (2 * p[2]- 1) 
        }
      if (comp == MELT )
        {
          r <- 2.0 ^ (2 * p[3]- 1) 
        }
      print(r)
      r
    }
    bcoef <- function(x,y,comp){0.0}
    ccoef <- function(x,y,comp){matrix(0.0,length(x),length(y))}

    #UNCOMMENT TO CREATE DATA FILES
    #coef <- amr.aXbYplusc("ase-coef-mb-1km.2d.hdf5","ase-coef-zero-1km.2d.hdf5", acoef,bcoef,ccoef,nghost=1  )  
    #amr.write(coef,pathto(coeffile))
    #amr.free(coef)

    #p[3] also defines the melt rate facrros glmelt and ambmelt
    abal <- function(x,y,comp)
      {  r <- 0.0
         if ( (comp == GLMELT) | (comp == AMBMELT) )
           {
             r <- 2.0 ^ (2 * p[3]- 1) 
           }
         r
       }
    bbal <- function(x,y,comp){0.0}
    cbal <- function(x,y,comp){matrix(0.0,length(x),length(y))}

    #UNCOMMENT TO CREATE DATA FILES
    #bal <- amr.aXbYplusc("ase-balance-melt-1km.2d.hdf5","ase-balance-melt-1km.2d.hdf5", abal,bbal,cbal,nghost=0  ) 
    #amr.write(bal,pathto(balfile))
    #amr.free(bal)
  }


# central and end members
p0 <- c(1/2,1/2,1/2) #unperturbed
p1 <- c(0,1/2,1/2)  # lo traction
p2 <- c(1,1/2,1/2)  # hi traction
p3 <- c(1/2,0,1/2)  # lo viscosity
p4 <- c(1/2,1,1/2)  # hi viscosity
p5 <- c(1/2,1/2,0)  # lo melt rate
p6 <- c(1/2,1/2,1)  # hi melt rate


if (FALSE)
  {
   #create a new LHS for the remaining members
   LHS <- maximinLHS(n=64,k=3)
   #re-order LHS from the outside in, so that if we run the ensemble in
   #order, we get the extreme members first
   LHS <- LHS[rev(order(apply(LHS,FUN=function(x){sum( (x-1/2)^2 )},1))),]
   write.table(LHS,"ensembleB/lhs-B.tab",row.names=FALSE)
} else {
  #read the LHS rather than create a new one
  LHS <- read.table(file="ensembleB/lhs-B.tab",header=TRUE)
}
 

P <- rbind(p0,p1,p2,p3,p4,p5,p6)

# Set end members to have names counting from B0000
for (i in 1:dim(P)[1])
  {
    apply_pert (as.numeric(P[i,]),paste("ase-lhs-B",formatC(i-1,width=4,flag="0"),sep=""))
    amr.free.all()
  }

# Set LHS members to have names counting from B1000
P <- LHS

for (i in 1:dim(P)[1])
  {
    apply_pert (as.numeric(P[i,]),paste("ase-lhs-B",formatC(999+i,width=4,flag="0"),sep=""))
    amr.free.all()
  }

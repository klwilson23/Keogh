library(MARSS)
RShowDoc("Quick_Start",package="MARSS")

stockRec <- read.csv("/Data/Keogh_StockRecruitment.csv",stringsAsFactors = F,header=T)

?MARSS

MLEobj <- MARSS(stockRec, model=list(), ..., fit=TRUE)

# example:

data(lakeWAplankton, package="MARSS")
# lakeWA
fulldat = lakeWAplanktonTrans
years = fulldat[,"Year"]>=1965 & fulldat[,"Year"]<1975
dat = t(fulldat[years,c("Greens", "Bluegreens")])
covariates = t(fulldat[years,c("Temp", "TP")])

# re-standardize the data
the.mean = apply(dat,1,mean,na.rm=TRUE)
the.sigma = sqrt(apply(dat,1,var,na.rm=TRUE))
dat = (dat-the.mean)*(1/the.sigma)

# re-standardize the covariates
the.mean = apply(covariates,1,mean,na.rm=TRUE)
the.sigma = sqrt(apply(covariates,1,var,na.rm=TRUE))
covariates = (covariates-the.mean)*(1/the.sigma)

Q <- U <- x0 <- "zero"; B <- Z <- "identity"
d <- covariates
A <- "zero"
D <- "unconstrained"
y <- dat # to show relationship between dat & the equation
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,D=D,d=d,x0=x0)
kem <- MARSS(y, model=model.list)

# process-error only

R <- A <- U <- "zero"; B <- Z <- "identity"
Q <- "equalvarcov"
C <- "unconstrained"
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=covariates)
kem <- MARSS(dat, model=model.list)

# 
model.list$B <- "diagonal and unequal"
kem <- MARSS(dat, model=model.list)

x0 <- dat[,1,drop=FALSE]
model.list$tinitx <- 1
model.list$x0 <- x0
kem <- MARSS(dat, model=model.list)

# process and observation error
D <- d <- A <- U <- "zero"; Z <- "identity"
B <- "diagonal and unequal"
Q <- "equalvarcov"
C <- "unconstrained"
c <- covariates
R <- diag(0.16,2)
x0 <- "unequal"
tinitx <- 1
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,D=D,d=d,C=C,c=c,x0=x0,tinitx=tinitx)
kem <- MARSS(dat, model=model.list)

# process and observation error, but covariates only affect observation process
C <- c <- A <- U <- "zero"; Z <- "identity"
B <- "diagonal and unequal"
Q <- "equalvarcov"
D <- "unconstrained"
d <- covariates
R <- diag(0.16,2)
x0 <- "unequal"
tinitx <- 1
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,D=D,d=d,C=C,c=c,x0=x0,tinitx=tinitx)
kem <- MARSS(dat, model=model.list)

library(Metrics)

shp <- readOGR("BRMIE235GC_SIR.shp")
adj_mat<-nb2mat(poly2nb(shp_sul,queen=F,row.names=shp_sul$NM_MICRO.x), style="B")
A<-adj_mat
g <- network(A, ignore.eval=FALSE, names.eval="a", directed =FALSE, label=TRUE)
C<-as.matrix(g,matrix.type="incidence")
p=length(C[1,])
n=length(C[,1])
A<-(C)%*%t(C)-diag(diag((C)%*%t(C)))
Ar<-t(C)%*%C-2*diag(rep(1,p))

theta_suave<-  rnorm(n, 0, 0.001)
phi.all<- 0.03 + theta_suave + rnorm(n, 0, 0.001)
mean.phi<-rpois(n, 5)
LP <- as.numeric(phi.all)
prob <- mean.phi*exp(LP)
Y<- rpois(n, prob)

formula <- Y ~ offset(log(mean.phi))
GLM<-glm(formula, family=poisson)

burnin<-35
n.sample<-1000

a_prior_invGamma<-0.001
b_prior_invGamma<-0.001

RENeGe <- poisson.mtv(formula=formula, W=A, Ar=Ar,C=C,  burnin=burnin,n.sample=n.sample)  

RENeGe_gamma <-poisson.mtvSIGMab(formula=formula, W=A, Ar=Ar, burnin=burnin, n.sample=n.sample, a_prior_invGamma=a_prior_invGamma, b_prior_invGamma=b_prior_invGamma) 




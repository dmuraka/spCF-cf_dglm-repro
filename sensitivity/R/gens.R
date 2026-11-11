SIG <- 0.5
gen_panel <- function(n_pop=400,T=40,rho=0.7,range=0.10,seed=1){
  set.seed(seed); s<-cbind(runif(n_pop),runif(n_pop)); W<-exp(-fields::rdist(s)/range); W<-W/sqrt(rowSums(W^2))
  b<-matrix(0,T,n_pop); b[1,]<-as.numeric(W%*%rnorm(n_pop))/sqrt(1-rho^2); for(t in 2:T) b[t,]<-rho*b[t-1,]+as.numeric(W%*%rnorm(n_pop))
  x1<-matrix(0,T,n_pop);x2<-matrix(0,T,n_pop); for(t in 1:T){x1[t,]<-as.numeric(W%*%rnorm(n_pop));x2[t,]<-as.numeric(W%*%rnorm(n_pop))}
  pt<-rep(1:T,each=n_pop);loc<-rep(1:n_pop,T);coords<-s[loc,]; al<-b[cbind(pt,loc)];X1<-x1[cbind(pt,loc)];X2<-x2[cbind(pt,loc)]
  eta<-2*X1-1.5*X2+al; y<-eta+rnorm(n_pop*T,sd=SIG)
  set.seed(seed+10000); te<-sort(sample.int(n_pop,round(n_pop*0.3))); ix<-loc%in%te
  list(y=y,eta=eta,coords=coords,pt=pt,x=cbind(x1=X1,x2=X2),ix=ix)}
gen_moving <- function(n_pop=400,T=40,rho=0.7,range=0.10,seed=1){
  set.seed(seed); anc<-cbind(runif(n_pop),runif(n_pop)); Wa<-exp(-fields::rdist(anc)/range); Wa<-Wa/sqrt(rowSums(Wa^2))
  b<-matrix(0,T,n_pop); b[1,]<-as.numeric(Wa%*%rnorm(n_pop))/sqrt(1-rho^2); for(t in 2:T) b[t,]<-rho*b[t-1,]+as.numeric(Wa%*%rnorm(n_pop))
  N<-n_pop*T; coords<-matrix(0,N,2); al<-numeric(N); X1<-numeric(N); X2<-numeric(N); pt<-rep(1:T,each=n_pop)
  for(t in 1:T){ st<-cbind(runif(n_pop),runif(n_pop)); Wi<-exp(-fields::rdist(st,anc)/range); Wi<-Wi/rowSums(Wi)
    Wx<-exp(-fields::rdist(st)/range); Wx<-Wx/sqrt(rowSums(Wx^2)); idx<-((t-1)*n_pop+1):(t*n_pop); coords[idx,]<-st
    al[idx]<-as.numeric(Wi%*%b[t,]); X1[idx]<-as.numeric(Wx%*%rnorm(n_pop)); X2[idx]<-as.numeric(Wx%*%rnorm(n_pop)) }
  al<-al*(sd(as.vector(b))/sd(al)); eta<-2*X1-1.5*X2+al; y<-eta+rnorm(N,sd=SIG)
  set.seed(seed+10000); ix<-logical(N); ix[sample.int(N,round(N*0.3))]<-TRUE
  list(y=y,eta=eta,coords=coords,pt=pt,x=cbind(x1=X1,x2=X2),ix=ix)}

##Following is the R script associated with the manuscript: "A Bayesian optimal
###escape model reveals bird species differ in their capacity to habituate to humans"
###Note that "alpha," our symbol for priors, will often be refered to as "h" or "H" in this script.


##This script will infer prior distributions and perform goodness of fit tests for comparison of a
###seven waterbird species escape behavior across three levels of human activity. A Fisher Information
###and Wilcoxin rank sum test are also included to test for significant differences in MLE parameters.

##Additionally, the script will calculate the relative differences in inferred prior distributions
###for activity level comparison per species.



##Data set input (by default must be a tab delimited file: see associated
###README for specific formatting instructions). Simply specify path to file in quotes.

FullFIDdataset<-read.delim("~/GitHub/WaterbirdPriors/Sutton et al Waterbird Data.txt")
attach(FullFIDdataset)


##This code is set up to test an optimal escape model based either on distance alone, or on distance and
###group size. To include group size in the model, set FLRF to TRUE. (FLRF is short for flock size risk factor)

FLRF<-FALSE


##Below are species specific energetic parameters. May be adjusted as neccessary depending on species
###being modeled.

################################################################################
###Energetic parameters (species specific)


###Daily energy budgets for each species
BS_E<-5650
DD_E<-1082
DM_E<-531.3
EC_E<-531
PBD_E<-1074
PSH_E<-988
SG_E<-288.5

BMs1<-c(BS_E,DD_E,DM_E,EC_E,PBD_E,PSH_E,SG_E)
Es<-(16.42)*(BMs1^.655)*.239



###Cost to flee for each species
BS_B<-5.650
DD_B<-1.082
DM_B<-.5313
EC_B<-.531
PBD_B<-1.074
PSH_B<-.988
SG_B<-.2885

Bs1<-c(BS_B,DD_B,DM_B,EC_B,PBD_B,PSH_B,SG_B)
Bs<-67.16*(Bs1^(0.825))*0.86


###End of energetic parameters
################################################################################


##Please specify the number of data sets to be generated for the goodness of fit test.
###Smaller values will run faster, but with reduced accuracy.

gof_number<-100

###You may optinally only run the parameter estimation by turning off gof test

#Run gof test? (Yes<-TRUE | No<-FALSE)
GOF<-TRUE


####The script is now ready to run!

################################################################################

system.time({
###Following are the packages used in this script
#nonlinear root solving package
library("rootSolve")
#for handling hypergeometric functions
library("hypergeo")
#for extending the functionality of base beta functions in R
#numerical differentiation
library("numDeriv")
#additional packages
library("gsl")
library("optimx")
library("dplyr")
library("tidyr")


means<-c()
ses<-c()
sds<-c()
vars<-c()
ps<-c()
qs<-c()
hessianps<-c()
hessianqs<-c()
uses<-c()
specs<-c()


FIDFsets<-c()
ADFsets<-c()
FLOCKFsets<-c()
FTFsets<-c()
speciesFsets<-c()
xFsets<-c()

ratiomeans<-c()
ratioses<-c()
FIDmeans<-c()
FIDses<-c()
BDmeans<-c()
BDses<-c()
ADmeans<-c()
ADses<-c()
samples<-c()
FullH<-c()

###Dataset subsetting to remove nonflight cases (FID=0) and FID>=AD cases, extreme cases the current model
##is not designed to handle
FullFIDdataset$FullAD<-FullFIDdataset$AD
FullFIDdataset$FullFID<-FullFIDdataset$Seperation.Distance
FullFIDdataset<-subset(FullFIDdataset,FullAD-FullFIDdataset$Seperation.Distance>-0.01 & FullFID>0.001)

FullFIDdataset$FullDL<-as.numeric(as.list(rep(0, nrow(FullFIDdataset))))

###Distance of closest approach (i.e. closest a pred could get to prey)
FTs<-c()
for(k in 1:nrow(FullFIDdataset)){
  if(FullFIDdataset$On.Shore.In.Water[k]=="In Water" & FullFIDdataset$Seperation.Distance[k] > FullFIDdataset$Distance.from.Edge[k]){
  FT<-FullFIDdataset$Distance.from.Edge[k]
  }else{
  FT<-FullFIDdataset$FullDL[k]
  }
  FTs<-c(FTs,FT)
}

x<-c()
for(k in 1:nrow(FullFIDdataset)){
  if(FullFIDdataset$Towards..Away[k]=="Towards"){
    xs<-0.5
  }else{
    xs<-0
  }
  x<-c(x,xs)
}



FullFIDdataset<-cbind(FullFIDdataset,as.data.frame(FTs),as.data.frame(x))

attach(FullFIDdataset)

###Human activity level binning
usage<-c()
for(u in 1:nrow(FullFIDdataset)){
  if(FullFIDdataset$Walker.Rate[u]<=summary(FullFIDdataset$Walker.Rate)["1st Qu."]){
    Usage<-"low"
  }
  if(FullFIDdataset$Walker.Rate[u]>summary(FullFIDdataset$Walker.Rate)["1st Qu."]&FullFIDdataset$Walker.Rate[u]<=summary(FullFIDdataset$Walker.Rate)["3rd Qu."]){
    Usage<-"mid"
  }
  if(FullFIDdataset$Walker.Rate[u]>summary(FullFIDdataset$Walker.Rate)["3rd Qu."]){
    Usage<-"high"
  }
  usage<-c(usage,Usage)
}
FullFIDdataset$usage<-usage
usages<-levels(as.factor(FullFIDdataset$usage))
usages<-c("low","mid","high")


species<-levels(FullFIDdataset$Species)
FullFIDdataset$FullFID<-FullFIDdataset$Seperation.Distance
FullFIDdataset$FullDW<-FullFIDdataset$Distance.from.Edge

#for(q in 1:nrow(FullFIDdataset)){
#  if(FullFIDdataset$FullFID[q]==FullFIDdataset$FullAD[q]){
#    FullFIDdataset$FullFID[q]<-(FullFIDdataset$FullFID[q]-0.001)
#  }
#}




################################################################################
###Functions used for FID, likelihood , and relative differebce calculations

rBeta <- function(n, shape1=2, shape2 =3, params = list(shape1, shape2),...){
  if(!missing(params)){
    shape1 <- params$shape1
    shape2 <- params$shape2
  }
  rbeta(n, shape1 = shape1, shape2 = shape2)
}


###DECISION MECHANISMS
rfs<-function(FID,AD,x,FLOCK){
  rf1<-(1-(x*(FID/AD)))*(1-(FID/AD))
  if(FLRF==FALSE){
    rf2<-1
  }else{
    rf2<-(0.999/FLOCK)
  }
  rf1*rf2
}

rfsc<-function(FID,AD,x,FLOCK){
  rf1<-(1-(x*(FID/AD)))*(1-(FID/AD))
  if(FLRF==FALSE){
    rf2<-0
  }else{
    rf2<-(0.999/FLOCK)
  }
  (1-rf1)*(1-rf2)
}



###FID function: predicts FID based on risk factors and appraoch path variables
FID1<-function(AD,h,x,FLOCK){
  #FID solution as follows
  flightdist<-function(FID){
    B<-Bs[n]
    E<-Es[n]
    #FID equation
    B-(E*(((rfs(FID,AD,x,FLOCK))*h)/(((rfs(FID,AD,x,FLOCK))*h)+((rfsc(FID,AD,x,FLOCK))*(1-h)))))
  }
  #find roots to get FID (use 'max' as the larger root would occur first in an encounter)
  #root solver checks for roots of 'flightdist' function from radial distances of zero to AD
  ft<-0
  fid<-max(uniroot.all(flightdist,c(ft,AD),n=10000))
  fid
}


###H as a function of FID and approach path variables
H_i<-function(FID,AD,x,FLOCK){
  B<-Bs[n]
  E<-Es[n]
  (B*rfsc(FID,AD,x,FLOCK))/((E*rfs(FID,AD,x,FLOCK))+(B*rfsc(FID,AD,x,FLOCK))-(B*rfs(FID,AD,x,FLOCK)))
}


###dH/dF i.e. the derivative of H wrt F (random sets)
dHFset<-function(FID){
  H_iset<-function(FID){
    AD<-ADFset[i]
    x<-xFset[i]
    FLOCK<-FLOCKFset[i]
    B<-Bs[n]
    E<-Es[n]
    (B*rfsc(FID,AD,x,FLOCK))/((E*rfs(FID,AD,x,FLOCK))+(B*rfsc(FID,AD,x,FLOCK))-(B*rfs(FID,AD,x,FLOCK)))
  }
  dHdF<-grad(H_iset,FID,method="Richardson")
}


###dH/dF (observed data)
dHFobs<-function(FID){
  H_iobs<-function(FID){
    AD<-fdata$FullSD[i]
    x<-fdata$x[i]
    FLOCK<-fdata$Flock.Size[i]
    B<-Bs[n]
    E<-Es[n]
    (B*rfsc(FID,AD,x,FLOCK))/((E*rfs(FID,AD,x,FLOCK))+(B*rfsc(FID,AD,x,FLOCK))-(B*rfs(FID,AD,x,FLOCK)))
  }
  grad(H_iobs,FID,method="Richardson")
}


###Incomplete regularized beta function
ibeta<-function(z,p,q){
  log(beta_inc(p,q,z))
}
ibeta2<-function(z,p,q){
  if(z > 0){
    log(beta_inc(p,q,z))
  }else{0}
}

###derivative of incomplete regularized beta function wrt first shape parameter
dpibeta<-function(p,q,z){
  libeta<-function(p){
    log(hyperg_2F1(p+q,1,p+1,z)) + p*log(z)+q*log(1-z)-log(p) - lbeta(p,q)
  }
  grad(func=libeta,p,method="Richardson")
}

###derivative of incomplete regularized beta function wrt second shape parameter
dqibeta<-function(p,q,z){
  libeta<-function(q){
    log(hyperg_2F1(p+q,1,p+1,z)) + p*log(z)+q*log(1-z)-log(p) - lbeta(p,q)
  }
  grad(func=libeta,q,method="Richardson")
}

###derivative of log beta wrt first shape parameter
dplbeta<-function(p,q){
  lbeta<-function(p){
    log(beta(p,q))
  }
  grad(func=lbeta,p,method="Richardson")
}

###derivative of log beta wrt second shape parameter
dqlbeta<-function(p,q){
  lbeta<-function(q){
    log(beta(p,q))
  }
  grad(func=lbeta,q,method="Richardson")
}


###probability density function for beta distributed variable
rhoalpha<-function(z,p,q,da){
  ((((z^(p-1))*((1-z)^(q-1)))/(beta(p,q)))*da)
}

###dH/dF for use in parameter estimation
dHFobs2<-function(FID,AD,x,FLOCK){
  H_iobs<-function(FID){
    B<-Bs[n]
    E<-Es[n]
    (B*rfsc(FID,AD,x,FLOCK))/((E*rfs(FID,AD,x,FLOCK))+(B*rfsc(FID,AD,x,FLOCK))-(B*rfs(FID,AD,x,FLOCK)))
  }
  grad(H_iobs,FID,method="Richardson")
}


###ll function for beta distributed variable
rhoalphac<-function(FID,AD,x,p,q,FLOCK){
  da<-abs((dHFobs2(FID,AD,x,FLOCK)))
  z<-H_i(FID,AD,x,FLOCK)
  rhoalpha(z,p,q,da)
}

###Beta function
betah<-function(z){
  ((((z)^(p1-1))*((1-z)^(q1-1)))/(exp(lbeta(p1,q1))))
}

###Absolute relative difference calculation for difference between two inferred prior distributions
relativediff<-function(z,p1,q1,p2,q2){
  one<-((((z)^(p1-1))*((1-z)^(q1-1)))/(exp(lbeta(p1,q1))))
  two<-((((z)^(p2-1))*((1-z)^(q2-1)))/(exp(lbeta(p2,q2))))
  abs(one-two)/max(abs(one),abs(two))
}


###End of FID, likelihood, and relative difference functions
################################################################################


################################################################################
###Following is the script for parameter estimation and a goodnes of fit test
################################################################################
percentiles<-c()
LLObss<-c()
MeanLLs<-c()
for(n in 1:length(species)){
  ###Parameters for subsetting of data set
  #variable for subsetting
  ssvar<-FullFIDdataset$Species
  sscond<-species[n]
  print(species[n])
  fullset1<-subset(FullFIDdataset,ssvar == sscond & FullFID>0.001)
  
  for(z in 1:length(usages)){
    ssvar2<-fullset1$usage
    sscond2<-usages[z]
    fullset<-subset(fullset1,ssvar2 == sscond2)
    nfdata<-subset(fullset,FullFID==0)
    fdata<-subset(fullset,FullFID>0)
    
    hcalcs<-c()
    for(g in 1:length(fdata$FullFID)){
      FID<-fdata$FullFID[g]
      AD<-fdata$FullAD[g]
      x<-fdata$x[g]
      FLOCK<-fdata$Flock.Size[g]
    hcalc<-H_i(FID,AD,x,FLOCK)
    hcalcs<-c(hcalcs,hcalc)
    }
    FullH<-c(FullH,hcalcs)
    
    
    BDmean<-mean(fdata$FullAD-fdata$FullFID)
    BDse<-sd(fdata$FullAD-fdata$FullFID)/sqrt(nrow(fdata))
    BDses<-c(BDses,BDse)
    BDmeans<-c(BDmeans,BDmean)
    FIDmean<-mean(fdata$FullFID)
    ratiomean<-mean(fdata$FullFID/fdata$FullAD)
    ratiose<-sd(fdata$FullFID/fdata$FullAD)/sqrt(nrow(fdata))
    ratiomeans<-c(ratiomeans,ratiomean)
    ratioses<-c(ratioses,ratiose)
    FIDse<-sd(fdata$FullFID)/sqrt(nrow(fdata))
    FIDses<-c(FIDses,FIDse)
    FIDmeans<-c(FIDmeans,FIDmean)
    ADmean<-mean(fdata$FullAD)
    ADmeans<-c(ADmeans,ADmean)
    ADse<-sd(fdata$FullAD)/sqrt(nrow(fdata))
    ADses<-c(ADses,ADse)
    sample<-nrow(fdata)
    samples<-c(samples,sample)
    
    
###Following is the script for maximum likelihood parameter estimation
    
    #log likelihood function
    LLfun<-function(x){
      p<-abs(x[1])
      q<-abs(x[2])

      intalphaspar<-c()
      for(i in 1:nrow(fdata)){
        if((fdata$FullFID[i]-min(0.5,fdata$FullFID[i]-(fdata$FTs[i]))) <= 0){
          lower=0.01
        }else{lower=fdata$FullFID[i]-min(0.5,fdata$FullFID[i]-(fdata$FTs[i]))}
        

        upper=fdata$FullFID[i]+min(0.5,fdata$FullAD[i]-fdata$FullFID[i])
        
        intalpha<-integrate(rhoalphac,lower=lower,upper=upper,AD=fdata$FullAD[i],x=fdata$x[i],p=p,q=q,FLOCK=fdata$Flock.Size[i])
        
        intalphaspar<-c(intalphaspar,intalpha$value)
      }
      
      sum(log(intalphaspar))
    }
    shape<-optimx(c(1,1),LLfun,method=c("BFGS"),hessian=TRUE,control=list(maximize=TRUE))
    
    p1<-abs(shape$p1)
    q1<-abs(shape$p2)
    if(FLRF==TRUE){
      cat("Inferred shape parameters for group size risk factor model: ",p1,q1,"\n")
      cat("\n")
    }else{cat("Inferred shape parameters for distance only risk risk factor model",p1,q1,"\n")
      cat("\n")}
    curve(betah,from=0,to=1,xlim=c(0,1),xlab="Alpha",ylab="Probability Density")
    hessianp<-hessian(LLfun,c(p1,q1))[1,1]
    hessianq<-hessian(LLfun,c(p1,q1))[2,2]
    hessianps<-c(hessianps,hessianp)
    hessianqs<-c(hessianqs,hessianq)
    
    
    
    
    
    ###Following is the script for calulating the log likelihood of the observed data based on
    ####parameter estimates
    
    intalphas<-c()
    for(i in 1:nrow(fdata)){
      intalpha<-integrate(rhoalphac,lower=fdata$FullFID[i]-min(0.5,fdata$FullFID[i]-(fdata$FTs[i]+0.01)),upper=fdata$FullFID[i]+min(0.5,fdata$FullAD[i]-fdata$FullFID[i]),AD=fdata$FullAD[i],x=fdata$x[i],p=p1,q=q1,FLOCK=fdata$Flock.Size[i])
      intalphas<-c(intalphas,intalpha$value)
    }
    sumflight<-sum(log(intalphas))
    
    LL_Obs<-sumflight
    LLObss<-c(LLObss,LL_Obs)
    
    
    ###Following is the script for reporting summary stats on the obs H dist
    meanh<-p1/(p1+q1)
    sdh<-sqrt((p1*q1)/(((p1+q1)^2)*(p1+q1+1)))
    seh<-sdh/sqrt(nrow(fullset))
    varh<-(p1*q1)/(((p1+q1)^2)*(p1+q1+1))
    means<-c(means,meanh)
    ses<-c(ses,seh)
    sds<-c(sds,sdh)
    vars<-c(vars,varh)
    ps<-c(ps,p1)
    qs<-c(qs,q1)
    uses<-c(uses,sscond2)
    specs<-c(specs,sscond)
    
    
    if(GOF==TRUE){  

      ###Following is the script for the exact test (goodness of fit)
      ###Main for loop calculates likelihoods for n number of generated datasets
      LL_set<-c()
      for(j in 1:gof_number){
        ###dataset reset after each likelihood calculation
        FIDFset<-c()
        ADFset<-c()
        hFset<-c()
        FTFset<-c()
        speciesFset<-c()
        xFset<-c()
        FLOCKFset<-c()
        fullset<-rbind(fdata,nfdata)
        ###For loop generates data sets of length n
        for(i in 1:(nrow(fullset))){
          ###FIDs calculated based on observed AD and randomnly generated H based on fitted H distribution
          h<-rBeta(1,shape1=p1,shape2=q1)
          randomFID<-FID1(fullset$FullAD[i],h,fullset$x[i],fullset$Flock.Size[i])

            FIDFset<-c(FIDFset,randomFID)
            FIDFdataset<-as.data.frame(FIDFset)
            FTFset<-c(FTFset,fullset$FTs[i])
            ADFset<-c(ADFset,fullset$FullAD[i])
            ADFdataset<-as.data.frame(ADFset)
            FLOCKFset<-c(FLOCKFset,fullset$Flock.Size[i])
            hFset<-c(hFset,h)
            speciesFset<-c(speciesFset,sscond)
            xFset<-c(xFset,fullset$x[i])
          
        }
        
        FIDFsets<-c(FIDFsets,FIDFset)
        ADFsets<-c(ADFsets,ADFset)
        FLOCKFsets<-c(FLOCKFsets,FLOCKFset)
        FTFsets<-c(FTFsets,FTFset)
        speciesFsets<-c(speciesFsets,speciesFset)
        xFsets<-c(xFsets,xFset)

        
        
        ###Fourth sum in likelihood calc
        intalphasset<-c()
        for(i in 1:length(FIDFset)){
          intalpha<-integrate(rhoalphac,lower=FIDFset[i]-min(0.5,FIDFset[i]),upper=FIDFset[i]+min(0.5,ADFset[i]-FIDFset[i]),AD=ADFset[i],x=xFset[i],p=p1,q=q1,FLOCK=FLOCKFset[i])
          intalphasset<-c(intalphasset,intalpha$value)
        }
        sumintalphasset<-sum(log(intalphasset))

        ###Likelihood calculation
        likelihood_set<-sumintalphasset
        LL_set<-c(LL_set,likelihood_set)
        print(cbind(j,species[n],usages[z]))

      }

      MeanLL<-mean(LL_set)
      MeanLLs<-c(MeanLLs,MeanLL)
      SDLL<-(sd(LL_set))
      Percentile<-(ecdf(LL_set)(LL_Obs))
      
      percentiles<-c(percentiles,Percentile)
      
      
    }
  }
}

})

###Reporting results of goodness of fit test
dat<-as.data.frame(cbind(specs,uses,LLObss,MeanLLs,percentiles))
colnames(dat)<-c("Species","Activity Level","Observed log likelihood","Mean generated log likelihood","Percentile")
GOFtest<-dat
View(GOFtest)

cat("Fisher Information Analysis\n")
###Caclulate hessians for FI test
hessianps<-c()
hessianqs<-c()
donors<-c()
for(n in 1:length(species)){
  ###Parameters for subsetting of data set
  #variable for subsetting
  hessp<-tibble(specs,uses,ps)
  hessq<-tibble(specs,uses,qs)
  ssvar<-FullFIDdataset$Species
  sscond<-species[n]
  print(species[n])
  fullset1<-subset(FullFIDdataset,ssvar == sscond & FullFID>0.001)
  
  hessp<-subset(hessp,specs == sscond)
  hessq<-subset(hessq,specs == sscond)

for(z in 1:length(usages)){
    ssvar2<-fullset1$usage
    sscond2<-usages[z]
    fullset<-subset(fullset1,ssvar2 == sscond2)
    nfdata<-subset(fullset,FullFID==0)
    fdata<-subset(fullset,FullFID>0)

    

    
    
    ###Following is the script for maximum likelihood parameter estimation
    
    #log likelihood function
    LLfun<-function(x){
      p<-abs(x[1])
      q<-abs(x[2])

      intalphaspar<-c()
      for(i in 1:nrow(fdata)){
        intalpha<-integrate(rhoalphac,lower=fdata$FullFID[i]-min(0.5,fdata$FullFID[i]-(fdata$FTs[i])),upper=fdata$FullFID[i]+min(0.5,fdata$FullAD[i]-fdata$FullFID[i]),AD=fdata$FullAD[i],x=fdata$x[i],p=p,q=q,FLOCK=fdata$Flock.Size[i])
        intalphaspar<-c(intalphaspar,intalpha$value)
      }
 
      sum(log(intalphaspar))
    }

for(t in 1:3){
  hp<-as.numeric(hessp$ps[t])
  hq<-as.numeric(hessq$qs[t])
    hessianp<-hessian(LLfun,c(hp,hq))[1,1]
    hessianq<-hessian(LLfun,c(hp,hq))[2,2]
    hessianps<-c(hessianps,hessianp)
    hessianqs<-c(hessianqs,hessianq)
    donors<-c(donors,t)
}
}}


###Fisher Information test
{
fishersp<- -hessianps
fishersq<- -hessianqs
FIdat<-tibble(specs,uses,samples,ps,qs)


x <-
  FIdat %>%
  as_tibble %>%
  left_join({
    FIdat %>%
      spread(key=uses,value=samples,fill=0) %>%
      mutate(
        low=as.numeric(low),
        mid=as.numeric(mid),
        high=as.numeric(high)
      ) %>%
      select(c(specs,low,mid,high)) %>%
      left_join(FIdat)
  })


x<-
  x%>%
  mutate(donor=1*sign(low)+2*sign(mid)+3*sign(high))

fx<-
  x%>%
  select(c(specs,uses))

fx<-
  fx%>%
  mutate(fishersp=fishersp,
         fishersq=fishersq,
         donor=donors)

fx<-
  fx%>%
  arrange(donor,fx$donor)

x<-
  x%>%
  arrange(donor,x$donor)

x<-
  x%>%
  mutate(fishersp=fx$fishersp,
         fishersq=fx$fishersq)

x1 <-
  FIdat %>%
  as_tibble %>%
  left_join({
    FIdat %>%
      spread(key=uses,value=ps,fill=0) %>%
      mutate(
        lowp=as.numeric(low),
        midp=as.numeric(mid),
        highp=as.numeric(high)
      ) %>%
      select(c(specs,lowp,midp,highp)) %>%
      left_join(FIdat)
  })
x1<-
  x1%>%
  mutate(donor=1*sign(lowp)+2*sign(midp)+3*sign(highp))

x1<-
  x1%>%
  arrange(donor,x1$donor)

xps<-
  x1 %>%
  as_tibble %>%
  select(c(specs,uses,lowp,midp,highp,donor))
x2 <-
  FIdat %>%
  as_tibble %>%
  left_join({
    FIdat %>%
      spread(key=uses,value=qs,fill=0) %>%
      mutate(
        lowq=as.numeric(low),
        midq=as.numeric(mid),
        highq=as.numeric(high)
      ) %>%
      select(c(specs,lowq,midq,highq)) %>%
      left_join(FIdat)
  })
x2<-
  x2%>%
  mutate(donor=1*sign(lowq)+2*sign(midq)+3*sign(highq))

x2<-
  x2%>%
  arrange(donor,x2$donor)
xqs<-
  x2 %>%
  as_tibble %>%
  select(c(specs,uses,lowq,midq,highq,donor))

x<-
  x %>%
  mutate(
    lowp=xps$lowp,
    lowq=xqs$lowq,
    midp=xps$midp,
    midq=xqs$midq,
    highp=xps$highp,
    highq=xqs$highq)

x<-
  x %>%
  mutate(
    p=(lowp+midp+highp),
    q=(lowq+midq+highq)
  )

x <-
  x %>%
  mutate(
    plower=ps-1.96*sqrt(1/(low+mid+high)*1/fishersp),
    pupper=ps+1.96*sqrt(1/(low+mid+high)*1/fishersp),
    qlower=qs-1.96*sqrt(1/(low+mid+high)*1/fishersq),
    qupper=qs+1.96*sqrt(1/(low+mid+high)*1/fishersq)
  ) %>%
  mutate(donor=1*sign(low)+2*sign(mid)+3*sign(high)) %>%
  mutate(donor=c('low','mid','high')[donor]) %>%
  mutate(test=paste(donor,uses)) %>%
  mutate(resultp = p>=plower & p<=pupper) %>%
  mutate(resultq = q>=qlower & q<=qupper)





FItest<-
  x %>%
  mutate(pCI=paste(round(plower,2),round(pupper,2)),
         qCI=paste(round(qlower,2),round(qupper,2))) %>%
  mutate(p=round(p,2))%>%
  mutate(q=round(q,2))

FItest<-
  FItest%>%
  select(c(specs,test,resultp,resultq,p,pCI,q,qCI))

FItest<-
  FItest%>%
  arrange(specs,FItest$specs)

FItest<-rename(FItest,Species=specs)

View(FItest)}


birdprior<-function(z,p1,q1){
  ((((z)^(p1-1))*((1-z)^(q1-1)))/(exp(lbeta(p1,q1))))
}



###Relative Difference Analysis
reldiffs<-c()
for(u in 1:nrow(x)){
  reldiff<-integrate(relativediff,lower=0,upper=1,p1=(x$lowp[u]+x$midp[u]+x$highp[u]),q1=(x$lowq[u]+x$midq[u]+x$highq[u]),p2=x$ps[u],q2=x$qs[u])$value
  reldiffs<-c(reldiffs,reldiff)
}
diffdat <-
  x<-
  x %>%
  mutate(
    reldiff=reldiffs)
diffdat<-
  diffdat%>%
  select(c(specs,test,reldiff))


diffdat<-diffdat[diffdat$reldiff !=0,]

diffdat<-subset(as.data.frame(diffdat),test=="low mid" | test=="low high" | test=="mid high")
tests<-c("low mid","mid high","low high")
diffdat$test<-factor(diffdat$test,levels=tests)
diffdat$specs<-factor(diffdat$specs,levels=species)

diffdat = diffdat %>% 
  group_by(specs) 




meandiffs<-c()
for(y in 1:length(species)){
  subdiffdat<-subset(diffdat,diffdat$specs==species[y])
  meandiff<-mean(subdiffdat$reldiff)
  meandiffs<-c(meandiffs,meandiff)
}

colnames(diffdat)<-c("Species","Comparison","Relative Difference")

RelDiff<-diffdat

View(RelDiff)

reldiffdat<-cbind(species,meandiffs)
colnames(reldiffdat)<-c("Species","Mean Relative Difference")
MeanDiff<-reldiffdat
View(MeanDiff)

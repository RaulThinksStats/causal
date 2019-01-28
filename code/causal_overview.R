#**************************************************************************
# Causal Inference - causal_overview.R
#
# Author: Raul JTA
#
# Summary: Cursory overview of fundamental ideas of causal inference.
#***************************************************************************


#***************************************************************************
# preliminary settings------------------------------------------------------
#***************************************************************************


library(tidyverse)
library(optmatch)
library(MatchIt)
library(tableone)
library(ggdag)


#***************************************************************************
# setting up data-----------------------------------------------------------
#***************************************************************************


load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))

#keeping only a subset of covariates
ARF<-as.numeric(rhc$cat1=='ARF')
CHF<-as.numeric(rhc$cat1=='CHF')
Cirr<-as.numeric(rhc$cat1=='Cirrhosis')
colcan<-as.numeric(rhc$cat1=='Colon Cancer')
Coma<-as.numeric(rhc$cat1=='Coma')
lungcan<-as.numeric(rhc$cat1=='Lung Cancer')
MOSF<-as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1=='MOSF w/Sepsis')
female<-as.numeric(rhc$sex=='Female')
age<-rhc$age
meanbp1<-rhc$meanbp1
aps1<-rhc$aps1

#treatment indicator
treatment<-as.numeric(rhc$swang1=='RHC')

#outcome
died<-as.numeric(rhc$death=='Yes')

#cleaning things up (colnames, tibble)
cov <- c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis","female","age","meanbp1","aps1")
rhc <- as_tibble(cbind(died,treatment,ARF,CHF,Cirr,colcan,Coma,lungcan,MOSF,sepsis,female,age,meanbp1,aps1))
rm(list=setdiff(ls(), c("rhc","cov")))

write_csv(rhc, "data/rhc.csv")

#***************************************************************************
# exploring open paths  ----------------------------------------------------
#***************************************************************************


#setting up tidy dag object
ex_dag <-  dagify(Y ~ A + X,
                  A ~ X, 
                  labels = c("Y" = "Death",
                             "A" = "Treatment",
                             "X" = "Pretreatment\n Covariates"),
                  exposure = "A",
                  outcome = "Y")

#calls for basic dag, list of open paths, and minimal adjusted set
ggdag(ex_dag, text = FALSE, use_labels = "label")
ggdag_paths(ex_dag, text = FALSE, use_labels = "label")
ggdag_adjustment_set(ex_dag, text = T, use_labels = "label")


#***************************************************************************
# matching methods  --------------------------------------------------------
#***************************************************************************


#Unmatched SMD
table1<- CreateTableOne(vars=cov,strata="treatment", data=rhc, test=FALSE)
unmatched <- print(table1,smd=TRUE)[-1,3]

#nearest neighbor pair matching
m.out <- matchit(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+
                   age+female+meanbp1+aps1, data=rhc, method = "nearest",distance = "mahalanobis")
matched1<- match.data(m.out)
table1<- CreateTableOne(vars=cov,strata="treatment", data=matched1, test=FALSE)
nn <- print(table1,smd=TRUE)[-1,3] #very good

#optimal matching
m.out <- matchit(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+
                   age+female+meanbp1+aps1, data=rhc, method = "optimal",distance = "mahalanobis")
matched2<- match.data(m.out)
table1<- CreateTableOne(vars=cov,strata="treatment", data=matched2, test=FALSE)
optimal <- print(table1,smd=TRUE)[-1,3]

#nearest neighbor pair matching with propensity score caliper
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+female+age+meanbp1+aps1,
             family=binomial(),data=rhc)
pscore<-psmodel$fitted.values

m.out <- matchit(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+
                   age+female+meanbp1+aps1, data=rhc, method = "nearest",distance = "mahalanobis",
                 caliper=.1*sd(pscore))
matched3<- match.data(m.out)
table1<- CreateTableOne(vars=cov,strata="treatment", data=matched3, test=FALSE)
nn_caliper <- print(table1,smd=TRUE)[-1,3] #very good


#***************************************************************************
# identifying and visualizing a good match  --------------------------------
#***************************************************************************


smd <- as.data.frame(cbind(unmatched, nn, optimal, nn_caliper))
row_names <- c(str_sub(rownames(smd), end = -13L), "total_smd")

#small adjustments to be able to add smds
smd <- apply(smd, 2, function(x) as.numeric(ifelse(x == "<0.001", "0", x)))
smd <- rbind(smd, apply(smd, 2, sum))
rownames(smd) <- row_names
smd

#visualizing nn match.
par(mfrow = c(1,2))

#fit a propensity score model for unmatched
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+female+age+meanbp1+aps1,
             family=binomial(),data=rhc)
pscore<-psmodel$fitted.values
hist(pscore[rhc$treatment==1],breaks=16, col=rgb(1,0,0,0.5), main = "Unmatched Balance", xlab="P-Score",ylim = c(0,400))
hist(pscore[rhc$treatment==0], col=rgb(0,0,1,0.5), add=T,breaks=16)
legend("topright", c("Treated", "Untreated","Overlap"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),rgb(.6,0,.6,1)), lwd=4)
box()

#fit a propensity score model for nn pair matched data.
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+female+age+meanbp1+aps1,
             family=binomial(),data=matched1)
pscore<-psmodel$fitted.values
hist(pscore[matched1$treatment==1],breaks=16, col=rgb(1,0,0,0.5), main = "Nearest Neighbor Balance", xlab="P-Score",ylim = c(0,400))
hist(pscore[matched1$treatment==0], col=rgb(0,0,1,0.5), add=T,breaks=16)
legend("topright", c("Treated", "Untreated","Overlap"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),rgb(.6,0,.6,1)), lwd=4)
box()


#***************************************************************************
# estimating causal risk difference with 95% CI  ---------------------------
#***************************************************************************

#outcome analysis
y_trt<-matched1$died[matched1$treatment==1]
y_con<-matched1$died[matched1$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)

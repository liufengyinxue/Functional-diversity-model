
setwd("C:\\Users\\adler\\Dropbox\\coexistence traits papers\\revision\\simulation")

## get functions
source("simulation_functions.r")
source("trait_dispersion_test.r")

## Run one simulation of the tolerance-fecundity trade-off model  ################################### 
# This make take a couple of minutes

# mu.hab is the mean resource availability on the simulated landscape
# sigma.hab is the standard deviation of the resource availability
simout=SimUnlimDispersal(mu.hab=0,sigma.hab=400,save.details=T)    
withinCommData=simout$cell.occupants  # save this data for final figure

#check coexistence criteria (a la Muller-Landau 2010)
N=simout$richness
fs=simout$survivors$fs; hs=simout$survivors$hs
crit1=fs*hs # should show decreasing order
crit2=(fs[2:N]*hs[2:N]-fs[1:(N-1)]*hs[1:(N-1)])/(fs[2:N]-fs[1:(N-1)]) # should show increasing order


## Run factorial simulation experiment ##################################
# This make take a couple of hours

mu=c(-100,0,100)  # community mean stress
sigma=c(0,200,400)  # within community standard deviation in stress
reps=10   # number of replicate runs for each combination of mu and sigma

allResults=matrix(NA,length(sigma)*length(mu)*reps,6)
commList=list(1)
colnames(allResults)=c("mu","sigma","richness","SScwm","rangemin","rangemax")

counter=0
for(imu in mu){
  for(isigma in sigma){
    for(irep in 1:reps){
      counter=counter+1
      simout=SimUnlimDispersal(imu,isigma)
      allResults[counter,]=c(imu,isigma,simout$richness,simout$SS.cwm,simout$range[1],simout$range[2])
      commList[[counter]]=simout$survivors
      print(counter)
    }
  }
}
allResults=as.data.frame(allResults)
meanResults=aggregate(cbind(richness,SScwm)~mu+sigma,data=allResults,mean)
sdResults=aggregate(cbind(richness,SScwm)~mu+sigma,data=allResults,sd)


## test for trait dispersion #############################################
keep=which(allResults$sigma>0)
commList2=commList[keep]  # only use communities that contain spatial heterogeneity
reps=length(commList2)
sppPool=NULL

# use only species that persist (across all communities) as species pool
#for(i in 1:reps) sppPool=c(sppPool,commList2[[i]]$SS)
#sppPool=unique(sppPool)

# use all virtual species between min and max established species
#sppPool=c(min(sppPool):max(sppPool))

# use all virtual species that can tolerate some portion of the landscape
range=c(min(allResults$rangemin),max(allResults$rangemax))
sppPool=seq(range[1],range[2],1)

pool<-data.frame(as.character(sppPool), sppPool)
for(i in 1:reps){
  sp_comm=as.character(commList2[[i]]$SS)
  kraftOut=test_trait_data(sp_comm, pool, log=FALSE, reps=999, abweight=FALSE, verbose=TRUE)
  if(i==1) {nulltests=as.data.frame(kraftOut)}
  else {nulltests=rbind(nulltests,as.data.frame(kraftOut))}
}
write.table(nulltests,"kraft_test_results.csv",row.names=F,sep=",")

## final figure ###########################################################

png("simFigure-new.png",width=3.25,height=9,units="in",res=300)
par(mfrow=c(4,1),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1),oma=c(0.5,0.5,0,0),cex.lab=1.2)
y=matrix(meanResults$SScwm,3,3)

#panel A
matplot(meanResults$mu[1:3],y,type="o",xlab="Mean resource availability",ylab="Mean seed size",
          ylim=c(1950,2800),pch=c(22,16,24),bg=c("white","white","grey"),col="black",lty=1,cex=1.1)
legend("topright",legend=sigma,pch=c(22,16,24),pt.bg=c("white","white","grey"),
        col="black",title=expression(paste("Resource ",sigma)),bty="n",cex=0.9)
mtext("A)",side=3,adj=0,line=0.5)

#panel B
y=matrix(meanResults$richness,3,3)
matplot(meanResults$mu[1:3],y,ylim=c(0,max(y)+1),type="o",xlab="Mean resource availability",ylab="Species richness",
          pch=c(22,16,24),bg=c("white","white","grey"),col="black",lty=1,cex=1.1)
mtext("B)",side=3,adj=0,line=0.5)

#panel C
plot(withinCommData,xlab="Resource availability",ylab="Seed size",xaxt="n")
axis(side=1,at=c(-1000,0,1000))
mtext("C)",side=3,adj=0,line=0.5)

#panel D
plot(obs_SDNDr~test.richness, data=nulltests, xlab="Species richness", ylab="Seed size SDNDr", type='n')
tapply(nulltests$null_SDNDr_025, INDEX=list(nulltests$test.richness), FUN=mean)->lo
tapply(nulltests$null_SDNDr_975, INDEX=list(nulltests$test.richness), FUN=mean)->hi
tapply(nulltests$null_mean_SDNDr, INDEX=list(nulltests$test.richness), FUN=mean)->mean
xp<-as.numeric(names(lo))
polygon(x=c(xp, 18,18, rev(xp), 2,2), y=c(lo,.03, .08, rev(hi), 5, .03), col="lightgray", border=NA)
lines(c(.4, mean, .05)~c(2, xp, 18), lwd=2)
points(obs_SDNDr~test.richness, data=nulltests)
mtext("D)",side=3,adj=0,line=0.5)

dev.off()



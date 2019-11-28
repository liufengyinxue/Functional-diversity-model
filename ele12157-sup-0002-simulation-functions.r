# Peter Adler Nov. 2012

# Unlimited dispersal simulation

# The only reason to simulate is to avoid having to test coexistence criteria for n species

# Stochasticity shows up in initial conditions, choices of sites to disturb, and 
# colonization of those sites (if multiple species have some prob. of colonization)

# setwd("C:\\Users\\adler\\Dropbox\\coexistence traits papers\\revision\\simulation")

# The simulation all wrapped up in one function:
SimUnlimDispersal=function(mu.hab,sigma.hab,save.details=F){
  # mu.hab = mean resource availability
  # sigma.hab = within community standard deviation of resource availability
  # save.details = T: Returns data on cell-level resource and seed size of cell occupants
 
  ## FIXED INPUTS ###################################################
  
  # assign species parameters
  Nspp=4400
  SS=seq(1,Nspp,1)  # assign seed size to n species
  htol=-1*(SS-max(SS)/2)  # assign habitat tolerance, centered on 0
  fs=SS[length(SS):1]   # assign fecundity (in reverse order of seed size)
  traits=data.frame(SS=SS,htol=htol,fs=fs)
  rm(SS,htol,fs)
  m=0.1  # mortality rate
  
  # set up simluation
  totTime=200
  ncells=Nspp*20
  habitat=sort(round(rnorm(ncells,mu.hab,sigma.hab),2))  # normally distributed habitat
  #habitat=round(runif(ncells,mu.hab-sigma.hab,mu.hab+sigma.hab),3) #uniformly distributed habitat
  
  # calculate proportion of habitat each species can occupy
  traits$hs=numeric(Nspp)
  for(i in 1:Nspp){
    traits$hs[i]=sum(habitat>=traits$htol[i])/ncells
  }
  
  # remove species that cannot tolerate any habitat
  traits=subset(traits,hs>0)
  
  ## FUNCTIONS ################################################
  
  # occupancy (p) and seed production function
  stateVars=function(presence,traits){
    counts=table(presence)
    ii=as.numeric(row.names(counts))
    occurrences=rep(0,NROW(traits))
    occurrences[ii]=counts
    seeds=occurrences*traits$fs
    return(list(p=occurrences/length(presence),seeds=seeds))
  }
  
  # competition to colonize an empty cell
  colonize=function(cell.hab,traits,seeds){
    dead=which(traits$htol>cell.hab)
    seeds[dead]=0
    if(sum(seeds)==0) stop("No spp can tolerate cell habitat")
    prob.col=seeds/sum(seeds)
    winner=which(prob.col==max(prob.col))
    if(length(winner)>1){
      winner=sample(winner,1)
    }
    return(winner)
  }
  
  ## MAIN LOOP ###############################################
  
  #initial conditions: start with equally weighted lottery for each site
  presence=numeric(ncells)  # vector to store the ID (row in traits data frame) of each spp in each cell
  for(j in 1:ncells){
    keep=which(traits$htol<=habitat[j]) 
    if(length(keep)==0) stop("No species can tolerate one of the cells")
    presence[j]=sample(keep,1)
  }
  pOut=matrix(NA,NROW(traits),totTime)
  
  tmp=stateVars(presence,traits)
  seeds=tmp$seeds
  pOut[,1]=tmp$p
  
  for(iT in 2:totTime){
    
    # kill cells
    kill=rbinom(ncells,1,m)
    presence[kill==1]=0
    
    # colonize empty cells
    empties=which(presence==0)
    for(j in 1:length(empties)){
      presence[empties[j]]=colonize(habitat[empties[j]],traits,seeds)
    }
    
    # calculate seed production and species abundances
    tmp=stateVars(presence,traits)
    seeds=tmp$seeds
    if(sum(tmp$p>0)==1) iT=totTime  # end simulation when only 1 spp left
    pOut[,iT]=tmp$p
  }
  
  if(save.details==T){
    cell.occupants=cbind(habitat,traits$SS[presence])
  }else{
    cell.occupants=NA
  }
  
  matplot(log10(t(pOut)),xlab="Time",ylab="log P",type="l")
  survivors=which(pOut[,totTime]>0)   
  richness=length(survivors)   # species richness
  SS.cwm=sum(traits$SS[survivors]*pOut[survivors,totTime]) # community weighted mean seed size
  range=c(min(traits$SS),max(traits$SS))
  out=list(richness=richness,SS.cwm=SS.cwm,survivors=traits[survivors,],
           cover=pOut[survivors,totTime],range=range,cell.occupants=cell.occupants)
  return(out)
  

}  # end simulation function


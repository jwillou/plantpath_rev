library(ridge)
library(glmnet)
library(lme4)
library(caret)
setwd("/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/pathogen meta/plantpath_rev/") 
setwd("G:/.shortcut-targets-by-id/1aXvrnwE2PSp9STOlaK560N8LSoVBPM0N/Willoughby lab/projects - active/pathogen meta/plantpath_rev/")

####read in and prep data####
data = read.table("meta-standardized_feb17.csv", header=T, sep=",")
data$aauthor = data$Study_char
table(data$Disease_Variable_Recat)

#calculate Cohen's d and clean data
m1  = data$Mutant_dis_mean_stand
m2  = data$Wildtype_dis_mean_stand
sd1 = data$Mutant_dis_sd_stand
sd2 = data$Wildtype_dis_sd_stand
n1  = data$Mutant_N
n2  = data$Wildtype_N
data$Cohen_d = (m1 - m2)/ (sqrt(((sd1^2 * (n1 - 1)) + (sd2^2 * (n2 - 1))) / (n1 + n2 - 2)))
hist(data$Cohen_d)
filter.data = data[!is.na(data$Cohen_d),]
filter.data$PathGenus = rep(NA, nrow(filter.data))
for(i in 1:nrow(filter.data)){filter.data$PathGenus[i] = strsplit(filter.data$Pathogen[i], " ")[[1]][1]}

#####distribution of effect sizes in categories####
discats = c("lesion size 1D", "lesion size 2D", "colony diameter", "disease incidence","disease index","fungal biomass","percentage survival") #unique(filter.data$Disease_Variable_Recat)
colors7 = c("dodgerblue3", "firebrick3", "darkorange3", "goldenrod2", "chartreuse3", "darkorchid", "navy")
means=NULL
for(i in 1:length(discats)){
  if(i==1){
    hist(filter.data$Cohen_d[filter.data$Disease_Variable_Recat==discats[i]], xlim=c(-100,20), ylim=c(0,200), breaks=seq(-100,20,10), col=alpha(colors7[i], 0.5), xlab="Effect Size", ylab="Frequency", main="")
    
  }else{
    par(new=T)
    hist(filter.data$Cohen_d[filter.data$Disease_Variable_Recat==discats[i]], xlim=c(-100,20), ylim=c(0,200), breaks=seq(-100,20,10), col=alpha(colors7[i], 0.5), xlab="", ylab="", main="", axes=F)
  }
  means = c(means, mean(filter.data$Cohen_d[filter.data$Disease_Variable_Recat==discats[i]]))
}

table(filter.data$Disease_Variable_Recat)
cbind(discats, means)

####single ridge####
diseaseindex = filter.data[!is.na(filter.data$Cohen_d),] #Disease_Variable_Recat + 
trial.all = train(Cohen_d ~ Virulence_factor + PathGenus + Necrotrophic  + Hemibiotrophic + Biotrophic + Saprophytic + Soil + Water +  Air + Leaves.stem + Root.Vascular + Fruits.Flower + Seeds.grains - 1, data=diseaseindex, method="glmnet", tuneGrid=expand.grid(alpha=0,lambda=seq(0,10,by=0.01)),trControl=trainControl(method="repeatedcv",n=10,repeats=100))
trial.all$bestTune$lambda #the best lambda, need to  plug into next lines of code
diseaseindex.out = linearRidge(Cohen_d ~ Virulence_factor + PathGenus + Necrotrophic  + Hemibiotrophic + Biotrophic + Saprophytic + Soil + Water +  Air + Leaves.stem + Root.Vascular + Fruits.Flower + Seeds.grains - 1, data=diseaseindex,lambda=trial.all$bestTune$lambda)
summary(diseaseindex.out)


####histograms####
#geneID by virulence factor
sdata = data[order(data$Virulence_factor, data$Gene_ID), ]
sdata$Gene_ID = as.factor(sdata$Gene_ID)
sdata$Virulence_factor = as.factor(sdata$Virulence_factor)
sdata = data[order(data$Virulence_factor, data$Gene_ID), ]
genes = unique(sdata$Gene_ID)
nbars = length(unique(sdata$Gene_ID))

vircols = data.frame(type=c("Effector", "Enzymes", "TF", "others"), color=c("dodgerblue3", "firebrick3", "goldenrod3", "grey50"))
plot(-100,-100, xlim=c(1,nbars), ylim=c(0,45), xlab="", ylab="", xaxt = "n")
axis(side=1, tick=T, xlim=c(0,nbars), at=1:nbars, labels=genes, las = 2, cex.axis=0.5)
barwidth=0.5
for(r in 1:nbars){
  t = sdata[sdata$Gene_ID==as.character(genes[r]),,drop=F]
  xpts = c((r-barwidth), (r+barwidth), (r+barwidth), (r-barwidth))
  ypts = c(0, 0, nrow(t), nrow(t))
  barcol = vircols$color[vircols$type==t$Virulence_factor[1]]
  polygon(x=xpts, y=ypts, col=barcol)
}
#infection strategy
idata = data[order(data$Infection_Strategy), ]
idata$Infection_Strategy = as.factor(idata$Infection_Strategy)
idata = data[order(data$Infection_Strategy), ]
insvar = c("Biotrophic", "Hemibiotrophic", "Necrotrophic", "Saprophytic") #unique(idata$Infection_Strategy)
vartall = c(132, 212, 497, 4)
nbars = length(insvar)

insvarcols = c("chartreuse3")
plot(-100,-100, xlim=c(0.5,(nbars+0.5)), ylim=c(0,500), xlab="", ylab="", xaxt = "n")
axis(side=1, tick=T, xlim=c(1,(nbars)), at=1:nbars, labels=insvar, las = 2, cex.axis=0.5)
barwidth=0.45
for(r in 1:nbars){
  t = idata[idata$Infection_Strategy==as.character(insvar[r]),,drop=F]
  xpts = c((r-barwidth), (r+barwidth), (r+barwidth), (r-barwidth))
  ypts = c(0, 0, vartall[r], vartall[r])
  barcol = insvarcols[1]
  polygon(x=xpts, y=ypts, col=barcol)
}

#plant host type
hdata = data[order(data$Host.tissue.target), ]
hdata$Host.tissue.target = as.factor(hdata$Host.tissue.target)
hdata = data[order(data$Host.tissue.target), ]
hostvar = c("Fruits and flowers", "Leaves and stems", "Roots and vascular tissue", "Seeds and grains") #unique(hdata$Host.tissue.target)
vartall = c(242, 258, 235, 56)
nbars = length(hostvar)

insvarcols = c("deeppink3")
plot(-100,-100, xlim=c(0.5,(nbars+0.5)), ylim=c(0,300), xlab="", ylab="", xaxt = "n")
axis(side=1, tick=T, xlim=c(1,(nbars)), at=1:nbars, labels=hostvar, las = 2, cex.axis=0.5)
barwidth=0.45
for(r in 1:nbars){
  t = hdata[hdata$Host.tissue.target==as.character(hostvar[r]),,drop=F]
  xpts = c((r-barwidth), (r+barwidth), (r+barwidth), (r-barwidth))
  ypts = c(0, 0, vartall[r],vartall[r])
  barcol = insvarcols[1]
  polygon(x=xpts, y=ypts, col=barcol)
}

#transmission mode
tdata = data[order(data$Mode.of.transmission), ]
tdata$Mode.of.transmission = as.factor(tdata$Mode.of.transmission)
tdata = data[order(data$Mode.of.transmission), ]
transvar = c("Air", "Soil", "Water"   )#unique(tdata$Mode.of.transmission)
transval = c(257, 250, 178)
nbars = length(transvar)

transvarcols = c("darkorange2")
plot(-100,-100, xlim=c(0.5,(nbars+0.5)), ylim=c(0,300), xlab="", ylab="", xaxt = "n")
axis(side=1, tick=T, xlim=c(1,(nbars)), at=1:nbars, labels=transvar, las = 2, cex.axis=0.5)
barwidth=0.45
for(r in 1:nbars){
  t = tdata[tdata$Mode.of.transmission==as.character(transvar[r]),,drop=F]
  xpts = c((r-barwidth), (r+barwidth), (r+barwidth), (r-barwidth))
  ypts = c(0, 0, transval[r], transval[r])
  barcol = transvarcols[1]
  polygon(x=xpts, y=ypts, col=barcol)
}


####Old and not used in final####
#ridge for each disease type/measurements individually
diseasemeasures = unique(data$Disease_Variable_Recat)

#lesionsize 1D group - boot 100
lesionsize1D = filter.data[filter.data$Disease_Variable_Recat=="lesion size 1D",]
lesionsize1D = lesionsize1D[complete.cases(lesionsize1D),]
#remove + soil, no data
trial.lesionsize1D = train(Cohen_d ~ Virulence_factor + PathGenus + Pathogen.lifecycle + Necrotrophic + Saprophytic + Hemibiotrophic + Biotrophic +  Air + Water + Leaves.stem + Fruits.Flower + Root.Vascular + Seeds.grains, data=lesionsize1D, method="glmnet", tuneGrid=expand.grid(alpha=0,lambda=seq(0,10,by=0.01)),trControl=trainControl(method="repeatedcv",n=10,repeats=100))
trial.lesionsize1D$bestTune$lambda #the best lambda, need to  plug into next lines of code
lesionsize1D.out = linearRidge(Cohen_d ~ Virulence_factor + PathGenus + Pathogen.lifecycle + Necrotrophic + Saprophytic + Hemibiotrophic + Biotrophic +  Air + Water + Leaves.stem + Fruits.Flower + Root.Vascular + Seeds.grains -, data=lesionsize1D,lambda=trial.lesionsize1D$bestTune$lambda)
summary(lesionsize1D.out)


#disease index - boot 100 
diseaseindex = filter.data[filter.data$Disease_Variable_Recat=="disease index",]
diseaseindex = diseaseindex[complete.cases(diseaseindex),]
#remove + Saprophytic +  Air + Fruits.Flower + Seeds.grains  , no data
trial.diseaseindex = train(Cohen_d ~ Virulence_factor + PathGenus + Necrotrophic  + Hemibiotrophic + Biotrophic + Soil + Water + Leaves.stem + Root.Vascular, data=diseaseindex, method="glmnet", tuneGrid=expand.grid(alpha=0,lambda=seq(0,10,by=0.01)),trControl=trainControl(method="repeatedcv",n=10,repeats=100))
trial.diseaseindex$bestTune$lambda #the best lambda, need to  plug into next lines of code
diseaseindex.out = linearRidge(Cohen_d ~ Virulence_factor + PathGenus + Necrotrophic  + Hemibiotrophic + Biotrophic + Soil + Water + Leaves.stem + Root.Vascular, data=diseaseindex,lambda=trial.diseaseindex$bestTune$lambda)
summary(diseaseindex.out)


###categories with fewer data points
#colony diameter - boot 100
colonydiameter = filter.data[filter.data$Disease_Variable_Recat=="colony diameter",]
colonydiameter = colonydiameter[complete.cases(colonydiameter),]
#remove + Seeds.grains + Saprophytic, no data
trial = train(Cohen_d ~ Virulence_factor + PathGenus + Pathogen.lifecycle +  Necrotrophic  + Hemibiotrophic + Biotrophic +  Air + Soil + Water + Leaves.stem + Fruits.Flower + Root.Vascular, data=colonydiameter, method="glmnet", tuneGrid=expand.grid(alpha=0,lambda=seq(0,10,by=0.01)),trControl=trainControl(method="repeatedcv",n=10,repeats=100))
trial$bestTune$lambda #the best lambda, need to  plug into next lines of code
colonydiameter.out = linearRidge(Cohen_d ~ Virulence_factor + PathGenus + Pathogen.lifecycle +  Necrotrophic  + Hemibiotrophic + Biotrophic +  Air + Soil + Water + Leaves.stem + Fruits.Flower + Root.Vascular, data=colonydiameter,lambda=trial$bestTune$lambda)
summary(colonydiameter.out)


#"disease spikelets"   "fungal biomass"      "lesion size 2D"  "percentage survival"   


##histograms
#disease severity
dsdata = data[order(data$Disease_Variable_Recat), ]
dsdata$Disease_Variable_Recat = as.factor(dsdata$Disease_Variable_Recat)
dsdata = data[order(data$Disease_Variable_Recat), ]
disvar = unique(dsdata$Disease_Variable_Recat)
nbars = length(unique(dsdata$Disease_Variable_Recat))

disvarcols = c("darkorchid3")
plot(-100,-100, xlim=c(0.5,(nbars+0.5)), ylim=c(0,250), xlab="", ylab="", xaxt = "n")
abline(h=100, col="grey50", lty=2)
axis(side=1, tick=T, xlim=c(1,(nbars)), at=1:nbars, labels=disvar, las = 2, cex.axis=0.5)
barwidth=0.5
for(r in 1:nbars){
  t = dsdata[dsdata$Disease_Variable_Recat==as.character(disvar[r]),,drop=F]
  xpts = c((r-barwidth), (r+barwidth), (r+barwidth), (r-barwidth))
  ypts = c(0, 0, nrow(t), nrow(t))
  barcol = disvarcols[1]
  polygon(x=xpts, y=ypts, col=barcol)
}




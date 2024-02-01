library(corrplot)
library(beeswarm)
library(lme4)
library(readxl)
library(ggplot2)
library(reshape)
library(RColorBrewer)  

## elena e. giorgi
## fred hutchinson, seattle, wa
## egiorgi@fredhutch.org

setwd("______________")

mydata <- read.csv("BRAVE_under5_vax_dataset_080823.cleaned.csv", stringsAsFactors=F)

length(which(!is.na(mydata$year_pcr1)))
mydata[which(!is.na(mydata$year_pcr1)), c(1,3,5,10,13:21)]
## one id had covid twice, 8/2021 & 1/2022

summary(mydata$weight_kg)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  6.380   9.225  11.890  12.173  13.825  22.220 

summary(mydata$timing_weight)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -150.000  -16.500    0.000   -4.693   13.000  164.000 

aggregate(timing_weight ~ age_cat, mydata, summary)

### Heatmaps
mfivars <- 27:40
id50vars <- c(41, 43, 45)
id80vars <- c(42, 44, 46)

svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
svars[6] <- paste(svars[6], "PLUS", sep=".")

H <- as.data.frame(matrix(NA, ncol=2, nrow=length(mfivars)))
names(H) <- c("mRNA-1273", "BNT162b2")
rownames(H) <- svars

for(i in 1:length(mfivars)){ 
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	xxx <- 	aggregate(mfivar ~ vax_type, mydata, function(x) mean(log10(x)))
	H[i, ] <- xxx[ , 2] 
	}

H_melt <- melt(cbind(rownames(H), H))
names(H_melt) <- c("Variant", "Group", "GeoMean")
H_melt$Group <- factor(H_melt$Group, levels=c("mRNA-1273", "BNT162b2"))
H_melt$Variant <- factor(H_melt$Variant, levels=svars[14:1])

ggplot(H_melt, aes(Group, Variant)) +  geom_tile(aes(fill = GeoMean)) + scale_fill_distiller(palette = "Purples", direction = +1, limits = log10(c(5000, 20000))) + theme(text = element_text(size = 18))


H <- as.data.frame(matrix(NA, ncol=4, nrow=length(mfivars)))
names(H) <- c("mRNA-1273 (-)", "BNT162b2 (-)", "mRNA-1273 (+)", "BNT162b2 (+)")
rownames(H) <- svars

for(i in 1:length(mfivars)){ 
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	xxx <- 	aggregate(mfivar ~ vax_type+prior_infection, mydata, function(x) mean(log10(x)))
	H[i, ] <- xxx[ , 3] 
	}

H_melt <- melt(cbind(rownames(H), H))
names(H_melt) <- c("Variant", "Group", "GeoMean")
H_melt$Group <- factor(H_melt$Group, levels=c("mRNA-1273 (-)", "mRNA-1273 (+)", "BNT162b2 (-)",  "BNT162b2 (+)"))
H_melt$Variant <- factor(H_melt$Variant, levels=svars[14:1])

ggplot(H_melt, aes(Group, Variant)) +  geom_tile(aes(fill = GeoMean)) + scale_fill_distiller(palette = "Purples", direction = +1, limits = log10(c(1000,25000)), breaks = log10(c(1,5000,10000,15000,20000))) + theme(text = element_text(size = 18))

# ggplot(Hmelt1, aes(Group, Variant)) +  geom_tile(aes(fill = GeoMean)) + scale_fill_gradient(low="white", high="blue", limits = c(0,25000), breaks =  c(0,10000,15000, 20000, 25000)) + geom_text(aes(label = round(GeoMean, 1))) + theme(text = element_text(size = 18))

### by age
H <- as.data.frame(matrix(NA, ncol=4, nrow=length(mfivars)))
names(H) <- c("mRNA-1273 (24-59 mo)", "BNT162b2 (24-59 mo)", "mRNA-1273 (6-23 mo)", "BNT162b2 (6-23 mo)")
rownames(H) <- svars

for(i in 1:length(mfivars)){ 
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	xxx <- 	aggregate(mfivar ~ vax_type+age_cat, mydata, function(x) mean(log10(x)))
	H[i, ] <- xxx[ , 3] 
	}

H_melt <- melt(cbind(rownames(H), H))
names(H_melt) <- c("Variant", "Group", "GeoMean")
H_melt$Group <- factor(H_melt$Group, levels=c("mRNA-1273 (6-23 mo)", "mRNA-1273 (24-59 mo)", "BNT162b2 (6-23 mo)",  "BNT162b2 (24-59 mo)"))
H_melt$Variant <- factor(H_melt$Variant, levels=svars[14:1])

ggplot(H_melt, aes(Group, Variant)) +  geom_tile(aes(fill = GeoMean)) + scale_fill_distiller(palette = "Purples", direction = +1, limits = log10(c(3500,20000))) + theme(text = element_text(size = 18))



mfivars <- id50vars
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))

H <- as.data.frame(matrix(NA, ncol=2, nrow=length(mfivars)))
names(H) <- c("mRNA-1273", "BNT162b2")
rownames(H) <- svars

for(i in 1:length(mfivars)){ 
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	xxx <- 	aggregate(mfivar ~ vax_type, mydata, function(x) mean(log10(x), na.rm=T))
	H[i, ] <- xxx[ , 2] 
	}

H_melt <- melt(cbind(rownames(H), H))
names(H_melt) <- c("Variant", "Group", "GeoMean")
H_melt$Group <- factor(H_melt$Group, levels=c("mRNA-1273", "BNT162b2"))
H_melt$Variant <- factor(H_melt$Variant, levels=svars[3:1])

ggplot(H_melt, aes(Group, Variant)) +  geom_tile(aes(fill = GeoMean)) + scale_fill_distiller(palette = "YlOrRd", direction = +1, limits = log10(c(100,500))) + theme(text = element_text(size = 18))

H <- as.data.frame(matrix(NA, ncol=4, nrow=length(mfivars)))
names(H) <- c("mRNA-1273 (-)", "BNT162b2 (-)", "mRNA-1273 (+)", "BNT162b2 (+)")
rownames(H) <- toupper(svars)

for(i in 1:length(mfivars)){ 
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	xxx <- 	aggregate(mfivar ~ vax_type+prior_infection, mydata, function(x) mean(log10(x)))
	H[i, ] <- xxx[ , 3] 
	}

H_melt <- melt(cbind(rownames(H), H))
names(H_melt) <- c("Variant", "Group", "GeoMean")
H_melt$Group <- factor(H_melt$Group, levels=c("mRNA-1273 (-)", "mRNA-1273 (+)", "BNT162b2 (-)",  "BNT162b2 (+)"))
H_melt$Variant <- factor(H_melt$Variant, levels=svars[3:1])

ggplot(H_melt, aes(Group, Variant)) +  geom_tile(aes(fill = GeoMean)) + scale_fill_distiller(palette = "YlOrRd", direction = +1, limits = log10(c(50,1500))) + theme(text = element_text(size = 18))




H <- as.data.frame(matrix(NA, ncol=4, nrow=length(mfivars)))
names(H) <- c("mRNA-1273 (24-59 mo)", "BNT162b2 (24-59 mo)", "mRNA-1273 (6-23 mo)", "BNT162b2 (6-23 mo)")
rownames(H) <- toupper(svars)

for(i in 1:length(mfivars)){ 
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	xxx <- 	aggregate(mfivar ~ vax_type+age_cat, mydata, function(x) mean(log10(x)))
	H[i, ] <- xxx[ , 3] 
	}

H_melt <- melt(cbind(rownames(H), H))
names(H_melt) <- c("Variant", "Group", "GeoMean")
H_melt$Group <- factor(H_melt$Group, levels=c("mRNA-1273 (6-23 mo)", "mRNA-1273 (24-59 mo)", "BNT162b2 (6-23 mo)",  "BNT162b2 (24-59 mo)"))
H_melt$Variant <- factor(H_melt$Variant, levels=svars[3:1])

ggplot(H_melt, aes(Group, Variant)) +  geom_tile(aes(fill = GeoMean)) + scale_fill_distiller(palette = "YlOrRd", direction = +1, limits = log10(c(90,700))) + theme(text = element_text(size = 18))


### Beeswarm plots for Fig 1
mfivars <- c(27, 31, (37:41), 43, 45)
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))

 figure <- "MFIplots.TwoGroupsOnly.Fig1.pdf"
 pdf(file=figure, pointsize=14, width=21, height=21, useDingbats=F) 
 par(mfrow=c(3,1), mar=c(6,6,6,6), xpd=TRUE)
  	
 plot(1,1, xlim=c(0.5, 3.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[1:3], at=1:3, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 1:3){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ vax_type, data=mydata, at=c((i-0.2), (i+0.2)), pch=c(16, 17), add=T, cex=3)
		
	boxplot(log10(mydata$mfivar) ~ vax_type, data=mydata, 
		at=c((i-0.2), (i+0.2)), add = T, names = c("", ""), col="#0000ff22", boxwex=0.2, medcol="red", lwd=2)
	 	
 	}
 	
 plot(1,1, xlim=c(3.5, 6.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=3, cex.lab=4)
 axis(1, labels=svars[4:6], at=4:6, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 4:6){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ vax_type, data=mydata, at=c((i-0.2), (i+0.2)), pch=c(16, 17), add=T, cex=3)
		
	boxplot(log10(mydata$mfivar) ~ vax_type, data=mydata, 
		at=c((i-0.2), (i+0.2)), add = T, names = c("", ""), col="#0000ff22", boxwex=0.2, medcol="red", lwd=2)
	
 	}

 plot(1,1, xlim=c(6.5, 9.5), ylim=c(0.5, 4.5), xlab="", ylab="ID50", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[7:9], at=7:9, cex.axis=4)
 axis(2, labels=c(10, expression(10^2), expression(10^3), expression(10^4)), at=1:4, cex.axis=4, las=2)

  for(i in 7:9){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ vax_type, data=mydata, at=c((i-0.2), (i+0.2)), pch=c(16, 17), add=T, cex=3)
		
	boxplot(log10(mydata$mfivar) ~ vax_type, data=mydata, 
		at=c((i-0.2), (i+0.2)), add = T, names = c("", ""), col="#0000ff22", boxwex=0.2, medcol="red", lwd=2)
 	
 	}

dev.off()




### Figure S1 (optional version)
mfivars <- c(28:36, 38, 39)
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))

 figure <- "MFIplots.TwoGroupsOnly.FigS1.pdf"
 pdf(file=figure, pointsize=14, width=21, height=21, useDingbats=F) 
 par(mfrow=c(3,1), mar=c(6,6,6,6), xpd=TRUE)
  	
 plot(1,1, xlim=c(0.5, 4.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[1:4], at=1:4, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 1:4){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ vax_type, data=mydata, at=c((i-0.2), (i+0.2)), pch=c(16, 17), add=T, cex=3)
		
	boxplot(log10(mydata$mfivar) ~ vax_type, data=mydata, 
		at=c((i-0.2), (i+0.2)), add = T, names = c("", ""), col="#0000ff22", boxwex=0.2, medcol="red", lwd=2)
	 	
 	}
 legend("top", legend=c("mRNA-1273", "BNT162b2"), pch=c(16, 17), bty='n', cex=3, horiz=T)

plot(1,1, xlim=c(4.5, 8.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[5:8], at=5:8, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 5:8){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ vax_type, data=mydata, at=c((i-0.2), (i+0.2)), pch=c(16, 17), add=T, cex=3)
		
	boxplot(log10(mydata$mfivar) ~ vax_type, data=mydata, 
		at=c((i-0.2), (i+0.2)), add = T, names = c("", ""), col="#0000ff22", boxwex=0.2, medcol="red", lwd=2)
	 	
 	}
  	
 plot(1,1, xlim=c(8.5, 11.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[9:11], at=9:11, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 9:11){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ vax_type, data=mydata, at=c((i-0.2), (i+0.2)), pch=c(16, 17), add=T, cex=3)
		
	boxplot(log10(mydata$mfivar) ~ vax_type, data=mydata, 
		at=c((i-0.2), (i+0.2)), add = T, names = c("", ""), col="#0000ff22", boxwex=0.2, medcol="red", lwd=2)
	 	
 	}
 
dev.off()
 

### Wilcoxon
pvals <- rep(0, length(svars))

for(i in 1:length(svars)){
	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
 	w <- wilcox.test(mydata$mfivar[which(mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$vax_type=="Pfizer")])
	pvals[i] <- format(w$p.value, 4)
	
}

cbind(svars, pvals, p.adjust(pvals, "fdr"))

### Figure 2
### Beeswarm plots 
mfivars <- c(27, 31, (37:41), 43, 45)
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))

 figure <- "MFIplots.ByPriorInfection.Fig2.pdf"
 pdf(file=figure, pointsize=14, width=21, height=21, useDingbats=F) 
 par(mfrow=c(3,1), mar=c(6,6,6,6), xpd=TRUE)
  	
 plot(1,1, xlim=c(0.5, 3.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[1:3], at=1:3, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 1:3){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), 
		col=c("red", "orange", "blue", "dodgerblue"), pch=c(16,16,17,17), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, 
		at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="black", lwd=2)
	 	
 	}
 
 legend("top", legend=c("mRNA-1273 (-)", "mRNA-1273 (+)", "BNT162b2 (-)",  "BNT162b2 (+)"), fill=c("red", "orange", "blue", "dodgerblue"), horiz=T, bty='n', cex=2)
 	
 plot(1,1, xlim=c(3.5, 6.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=3, cex.lab=4)
 axis(1, labels=svars[4:6], at=4:6, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 4:6){
 	
  	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), 
		col=c("red", "orange", "blue", "dodgerblue"), pch=c(16,16,17,17), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, 
		at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="black", lwd=2)
	
 	}

 plot(1,1, xlim=c(6.5, 9.5), ylim=c(0.5, 4.5), xlab="", ylab="ID50", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[7:9], at=7:9, cex.axis=4)
 axis(2, labels=c(10, expression(10^2), expression(10^3), expression(10^4)), at=1:4, cex.axis=4, las=2)

  for(i in 7:9){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), 
		col=c("red", "orange", "blue", "dodgerblue"), pch=c(16,16,17,17), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, 
		at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="black", lwd=2)
 	
 	}

dev.off()



### Wilcoxon for figure
outfile <- "WilcoxResults.Pfizer.vs.Moderna.txt"
mfivars <- 27:40
id50vars <- c(41, 43, 45)
id80vars <- c(42, 44, 46)

# mfivars <- id50vars
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
svars[6] <- paste(svars[6], "PLUS", sep=".")

X <- as.data.frame(matrix(NA, ncol=2, nrow=length(svars)))
names(X) <- c("Variant", "PValue")
cnt <- 1

 for(i in 1:length(mfivars)){

	mydata$mfivar <- mydata[ ,mfivars[i]]
	X$Variant[i] <- svars[i]
 
	w1 <- wilcox.test(mydata$mfivar[which(mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$vax_type=="Pfizer")])

	X[i, 2] <- format(w1$p.value, 4)

}

Assay <- rep("MFI", dim(X)[1])
X <- cbind(Assay, X)


mfivars <- id50vars
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
X1 <- as.data.frame(matrix(NA, ncol=2, nrow=length(svars)))
names(X1) <- c("Variant", "PValue")
cnt <- 1

  for(i in 1:length(mfivars)){

	mydata$mfivar <- mydata[ ,mfivars[i]]
	X1$Variant[i] <- svars[i]
 
	w1 <- wilcox.test(mydata$mfivar[which(mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$vax_type=="Pfizer")])

	X1[i, 2] <- format(w1$p.value, 4)

}

Assay <- rep("ID50", dim(X1)[1])
X1 <- cbind(Assay, X1)
X <- rbind(X, X1)

write.table(X, file=outfile, quote=F, row.names=F, sep="\t")



outfile <- "WilcoxResults.byPrior.txt"
mfivars <- 27:40
id50vars <- c(41, 43, 45)
id80vars <- c(42, 44, 46)

# mfivars <- id50vars
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
svars[6] <- paste(svars[6], "PLUS", sep=".")
X <- as.data.frame(matrix(NA, ncol=4, nrow=length(svars)*4))
names(X) <- c("Variant", "Test1", "Test2", "PValue")
cnt <- 1

 for(i in 1:length(mfivars)){

	mydata$mfivar <- mydata[ ,mfivars[i]]
	X$Variant[cnt : (cnt+3)] <- rep(svars[i], 4)
 
	w1 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Pfizer")])
	w2 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Pfizer")])

	w3 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Moderna")])
	w4 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Pfizer")], 
			mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Pfizer")])

	X[cnt, 2:4] <- c("Y", "Mod vs. Pf", format(w1$p.value, 4))
	X[(cnt+1), 2:4] <- c("N", "Mod vs. Pf", format(w2$p.value, 4))
	X[(cnt+2), 2:4] <- c("Y vs N", "Moderna", format(w3$p.value, 4))
	X[(cnt+3), 2:4] <- c("Y vs N", "Pfizer", format(w4$p.value, 4))
	cnt <- cnt+4

}

Assay <- rep("MFI", dim(X)[1])
X <- cbind(Assay, X)


mfivars <- id50vars
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
X1 <- as.data.frame(matrix(NA, ncol=4, nrow=length(svars)*4))
names(X1) <- c("Variant", "Test1", "Test2", "PValue")
cnt <- 1

 for(i in 1:length(mfivars)){

	mydata$mfivar <- mydata[ ,mfivars[i]]
	X1$Variant[cnt : (cnt+3)] <- rep(svars[i], 4)
 
	w1 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Pfizer")])
	w2 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Pfizer")])

	w3 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Moderna")])
	w4 <- wilcox.test(mydata$mfivar[which(mydata$prior_infection=="Y" & mydata$vax_type=="Pfizer")], 
			mydata$mfivar[which(mydata$prior_infection=="N" & mydata$vax_type=="Pfizer")])

	X1[cnt, 2:4] <- c("Y", "Mod vs. Pf", format(w1$p.value, 4))
	X1[(cnt+1), 2:4] <- c("N", "Mod vs. Pf", format(w2$p.value, 4))
	X1[(cnt+2), 2:4] <- c("Y vs N", "Moderna", format(w3$p.value, 4))
	X1[(cnt+3), 2:4] <- c("Y vs N", "Pfizer", format(w4$p.value, 4))
	cnt <- cnt+4

}

Assay <- rep("ID50", dim(X1)[1])
X1 <- cbind(Assay, X1)
X <- rbind(X, X1)

write.table(X, file=outfile, quote=F, row.names=F, sep="\t")

### figure S2
#### supplemental figure 1
 mfivars <- c(28:30, 32:36)
 svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))

figure <- "SuppFig2.pdf"
 pdf(file=figure, pointsize=14, width=42, height=21, useDingbats=F) 
 par(mfrow=c(2,1), mar=c(6,6,6,6), xpd=TRUE)
  	
 plot(1,1, xlim=c(0.5, 4.5), ylim=c(3,5), xlab="", ylab="MFI", main="MFI", axes=F, cex.main=6, cex.lab=3)
 axis(1, labels=svars[1:4], at=1:4, cex.axis=4)
 # axis(2, labels=c("1000", "10,000", "100,000"), at=3:5, cex.axis=4, las=2)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 1:4){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), 
		col=c("red", "orange", "blue", "dodgerblue"), pch=c(16,16,17,17), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, 
		at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), add = T, names = c("","", "", ""), col="#0000ff22", boxwex=0.15, lwd=2)
	 	
 	}
 legend("top", legend=c("mRNA-1273 (-)", "mRNA-1273 (+)", "BNT162b2 (-)",  "BNT162b2 (+)"), 
 		fill=c("red", "orange", "blue", "dodgerblue"), bty='n', cex=3, horiz=T)


plot(1,1, xlim=c(4.5, 8.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=3, cex.lab=3)
 axis(1, labels=svars[5:8], at=5:8, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 5:8){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), 
		col=c("red", "orange", "blue", "dodgerblue"), pch=c(16,16,17,17), add=T, cex=2, lwd=2)
		
	boxplot(log10(mydata$mfivar) ~ prior_infection+vax_type, data=mydata, 
		at=c((i-0.3), (i-0.1), (i+0.1), (i+0.3)), add = T, names = c("","", "", ""), col="#0000ff22", boxwex=0.15, lwd=2)
	 	
 	}


dev.off()


### Figure 3, by age
### Beeswarm plots 
mfivars <- c(27, 31, (37:41), 43, 45)
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))

 figure <- "MFIplots.ByAge.Fig3.pdf"
 pdf(file=figure, pointsize=14, width=21, height=21, useDingbats=F) 
 par(mfrow=c(3,1), mar=c(6,6,6,6), xpd=TRUE)
  	
 plot(1,1, xlim=c(0.5, 3.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[1:3], at=1:3, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 1:3){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), 
		pch=c(16,1,17,2), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, 
		at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="red", lwd=2)
	 	
 	}
 
 legend("top", legend=c("mRNA-1273 (6-23 mos)", "mRNA-1273 (24-59 mos)", "BNT162b2 (6-23 mos)",  "BNT162b2 (24-59 mos)"), 
 		pch=c(1, 16, 2, 17), horiz=T, bty='n', cex=2.5)
 	
 plot(1,1, xlim=c(3.5, 6.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=3, cex.lab=4)
 axis(1, labels=svars[4:6], at=4:6, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 4:6){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), 
		pch=c(16,1,17,2), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, 
		at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="red", lwd=2)
	
 	}

 plot(1,1, xlim=c(6.5, 9.5), ylim=c(0.5, 4.5), xlab="", ylab="ID50", main=" ", axes=F, cex.main=6, cex.lab=4)
 axis(1, labels=svars[7:9], at=7:9, cex.axis=4)
 axis(2, labels=c(10, expression(10^2), expression(10^3), expression(10^4)), at=1:4, cex.axis=4, las=2)

  for(i in 7:9){
 	
  	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), 
		pch=c(16,1,17,2), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, 
		at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="red", lwd=2)
 	
 	}

dev.off()

#### Figure S3


### Wilcoxon for figure
outfile <- "WilcoxResults.byAge.txt"
mfivars <- 27:40
id50vars <- c(41, 43, 45)
id80vars <- c(42, 44, 46)

# mfivars <- id50vars
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
svars[6] <- paste(svars[6], "PLUS", sep=".")

X <- as.data.frame(matrix(NA, ncol=4, nrow=length(svars)*4))
names(X) <- c("Variant", "Test1", "Test2", "PValue")
cnt <- 1

 for(i in 1:length(mfivars)){

	mydata$mfivar <- mydata[ ,mfivars[i]]
	X$Variant[cnt : (cnt+3)] <- rep(svars[i], 4)
 
	w1 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Pfizer")])
	w2 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Pfizer")])

	w3 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Moderna")])
	w4 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Pfizer")], 
			mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Pfizer")])

	X[cnt, 2:4] <- c("24-59 months", "Mod vs. Pf", format(w1$p.value, 4))
	X[(cnt+1), 2:4] <- c("6-23 months", "Mod vs. Pf", format(w2$p.value, 4))
	X[(cnt+2), 2:4] <- c("24-59 months vs 6-23 months", "Moderna", format(w3$p.value, 4))
	X[(cnt+3), 2:4] <- c("24-59 months vs 6-23 months", "Pfizer", format(w4$p.value, 4))
	cnt <- cnt+4

}

Assay <- rep("MFI", dim(X)[1])
X <- cbind(Assay, X)


mfivars <- id50vars
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
X1 <- as.data.frame(matrix(NA, ncol=4, nrow=length(svars)*4))
names(X1) <- c("Variant", "Test1", "Test2", "PValue")
cnt <- 1

 for(i in 1:length(mfivars)){

	mydata$mfivar <- mydata[ ,mfivars[i]]
	X1$Variant[cnt : (cnt+3)] <- rep(svars[i], 4)
 
	w1 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Pfizer")])
	w2 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Pfizer")])

	w3 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Moderna")], 
			mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Moderna")])
	w4 <- wilcox.test(mydata$mfivar[which(mydata$age_cat=="24-59 months" & mydata$vax_type=="Pfizer")], 
			mydata$mfivar[which(mydata$age_cat=="6-23 months" & mydata$vax_type=="Pfizer")])

	X1[cnt, 2:4] <- c("24-59 months", "Mod vs. Pf", format(w1$p.value, 4))
	X1[(cnt+1), 2:4] <- c("6-23 months", "Mod vs. Pf", format(w2$p.value, 4))
	X1[(cnt+2), 2:4] <- c("24-59 months vs 6-23 months", "Moderna", format(w3$p.value, 4))
	X1[(cnt+3), 2:4] <- c("24-59 months vs 6-23 months", "Pfizer", format(w4$p.value, 4))
	cnt <- cnt+4

}

Assay <- rep("ID50", dim(X1)[1])
X1 <- cbind(Assay, X1)
X <- rbind(X, X1)

X$FDRq <- c(p.adjust(X$PValue[which(X$Assay=="MFI")], method="fdr"), p.adjust(X$PValue[which(X$Assay=="ID50")], method="fdr"), p.adjust(X$PValue[which(X$Assay=="ID80")], method="fdr"))

write.table(X, file=outfile, quote=F, row.names=F, sep="\t")

### figure S3
#### supplemental figure 3
 mfivars <- c(28:30, 32:36)
 svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))

figure <- "SuppFig3.pdf"
 pdf(file=figure, pointsize=14, width=42, height=21, useDingbats=F) 
 par(mfrow=c(2,1), mar=c(6,6,6,6), xpd=TRUE)
  	
 plot(1,1, xlim=c(0.5, 4.5), ylim=c(3,5), xlab="", ylab="MFI", main="MFI", axes=F, cex.main=6, cex.lab=3)
 axis(1, labels=svars[1:4], at=1:4, cex.axis=4)
 # axis(2, labels=c("1000", "10,000", "100,000"), at=3:5, cex.axis=4, las=2)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 1:4){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), 
		pch=c(16,1,17,2), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, 
		at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="red", lwd=2)
	 	
 	}

 legend("top", legend=c("mRNA-1273 (6-23 mos)", "mRNA-1273 (24-59 mos)", "BNT162b2 (6-23 mos)",  "BNT162b2 (24-59 mos)"), 
 		pch=c(1, 16, 2, 17), horiz=T, bty='n', cex=2.5)

plot(1,1, xlim=c(4.5, 8.5), ylim=c(3,5), xlab="", ylab="MFI", main=" ", axes=F, cex.main=3, cex.lab=3)
 axis(1, labels=svars[5:8], at=5:8, cex.axis=4)
 axis(2, labels=c(expression(10^3), expression(10^4), expression(10^5)), at=3:5, cex.axis=4, las=2)

  for(i in 5:8){
 	
	mydata$mfivar <- mydata[ ,mfivars[i]]
 	mydata$mfivar <- mydata[ ,mfivars[i]]
	beeswarm(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), 
		pch=c(16,1,17,2), add=T, cex=2)
		
	boxplot(log10(mydata$mfivar) ~ age_cat+vax_type, data=mydata, 
		at=c((i-0.1), (i-0.3), (i+0.3), (i+0.1)), add = T, names = c("","", "", ""), 
		col="#0000ff22", boxwex=0.15, medcol="red", lwd=2)
	 	
 	}


dev.off()

########### GLM #############
mfivars <- 27:46
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
svars[6] <- paste(svars[6], "PLUS", sep=".")

### Anova analysis, MFI
AnovaPvalues <- as.data.frame(matrix(NA, ncol=6, nrow=length(mfivars)))
names(AnovaPvalues) <- c("Variant", "AllInteraction", "Age", "Vax.Inf.Interaction", "PriorInf", "VaxType")

for(i in 1:length(mfivars)){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
 	AnovaPvalues[i,1] <- names(mydata)[mfivars[i]]
 	g0 <- lm(log10(mydata$mfivar) ~ I(prior_infection)*I(vax_type)*I(age_cat), data=mydata)
 	g1 <- lm(log10(mydata$mfivar) ~ I(prior_infection)*I(vax_type)+I(age_cat), data=mydata)
	AnovaPvalues[i,2] <- anova(g0,g1)$Pr[2]
	
 	g0 <- lm(log10(mydata$mfivar) ~ I(prior_infection)*I(vax_type)+I(age_cat), data=mydata)
 	g1 <- lm(log10(mydata$mfivar) ~ I(prior_infection)*I(vax_type), data=mydata)
	AnovaPvalues[i,3] <- anova(g0,g1)$Pr[2]

	g0 <- lm(log10(mydata$mfivar) ~ I(prior_infection)*I(vax_type), data=mydata)
 	g1 <- lm(log10(mydata$mfivar) ~ I(prior_infection)+I(vax_type), data=mydata)
	AnovaPvalues[i,4] <- anova(g0,g1)$Pr[2]

	g0 <- lm(log10(mydata$mfivar) ~ I(prior_infection)+I(vax_type), data=mydata)
 	g1 <- lm(log10(mydata$mfivar) ~ I(vax_type), data=mydata)
	AnovaPvalues[i,5] <- anova(g0,g1)$Pr[2]
	
	g0 <- lm(log10(mydata$mfivar) ~ I(prior_infection)+I(vax_type), data=mydata)
 	g1 <- lm(log10(mydata$mfivar) ~ I(prior_infection), data=mydata)
	AnovaPvalues[i,6] <- anova(g0,g1)$Pr[2]
	
}

### age not significant, also checked breaking it down in quartiles

AnovaPvalues$FDRq <- p.adjust(AnovaPvalues$PriorInf, "fdr")
AnovaPvalues
#          Variant AllInteraction        Age Vax.Inf.Interaction     PriorInf     VaxType         FDRq
#1       d614g_mfi     0.32562261 0.27661532          0.58267405 3.240400e-01 0.166277615 5.167496e-01
#2       alpha_mfi     0.06779922 0.25143378          0.61685641 6.111569e-01 0.196747970 6.111569e-01
#3        beta_mfi     0.18674621 0.42896902          0.75619807 3.617247e-01 0.299018290 5.167496e-01
#4       gamma_mfi     0.28992607 0.76899589          0.93631286 5.229831e-01 0.490118837 5.810924e-01
#5       delta_mfi     0.16520489 0.73307509          0.89861987 4.957090e-01 0.250270401 5.810924e-01
#6  delta_plus_mfi     0.03377193 0.33785546          0.82814503 5.827375e-01 0.168496965 6.111569e-01
#7         eta_mfi     0.26823728 0.59135397          0.81705899 3.475083e-01 0.324928425 5.167496e-01
#8        iota_mfi     0.10550494 0.41854672          0.82988081 3.574924e-01 0.220603405 5.167496e-01
#9       kappa_mfi     0.05218754 0.29120717          0.66979189 5.056337e-01 0.151380446 5.810924e-01
#10     lambda_mfi     0.14338605 0.26394638          0.65060905 3.107151e-01 0.123267248 5.167496e-01
#11        ba1_mfi     0.09329599 0.13632777          0.59280495 4.872021e-01 0.146444506 5.810924e-01
#12      ba1.1_mfi     0.20169759 0.36421585          0.05863007 2.827623e-05 0.017024615 8.078922e-05
#13        ba2_mfi     0.13852644 0.94712759          0.53935654 1.892838e-05 0.021644794 6.309459e-05
#14      ba4_5_mfi     0.19162348 0.98466801          0.45823258 4.187761e-10 0.046697269 1.675104e-09
#15     d614g_nt50     0.14580496 0.20801865          0.41663308 6.749626e-05 0.006922199 1.499917e-04
#16     d614g_nt80     0.10802182 0.13005225          0.30257128 3.657664e-05 0.020261408 9.144159e-05
#17       ba1_nt50     0.03611032 0.26780737          0.05597923 1.183641e-11 0.823807253 6.437279e-11
#18       ba1_nt80     0.13200053 0.52325687          0.13099230 2.310665e-12 0.427521752 2.310665e-11
#19      ba45_nt50     0.53345136 0.01987185          0.35810631 5.662385e-13 0.271762448 1.132477e-11
#20      ba45_nt80     0.56409227 0.13302993          0.54328934 1.287456e-11 0.037767000 6.437279e-11


##### GLM tables
mfivars <- 27:46
svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1])))
svars[6] <- paste(svars[6], "PLUS", sep=".")

SuppTable1 <- as.data.frame(matrix(NA, ncol=5, nrow=3*length(mfivars)))
names(SuppTable1) <- c("Variant", "Variable", "FoldChange", "CI.95", "QValue")
pvals <- rep(0, 3*length(mfivars))

mydata$vax_type <- factor(mydata$vax_type, levels=c("Pfizer", "Moderna"))
mydata$age_cat <- factor(mydata$age_cat, levels=c("6-23 months", "24-59 months"))
mydata$prior_infection <- factor(mydata$prior_infection, levels=c("N", "Y"))

cnt <- 1
for(i in 1:length(mfivars)){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
 	SuppTable1$Variant[cnt:(cnt+2)] <- rep(svars[i], 2)
 	
 	g0 <- lm(log10(mydata$mfivar) ~ prior_infection + vax_type + age_cat, data=mydata)
 	s0 <- summary(g0)
 	c0 <- confint(g0, level=0.90)

 	## Prior Inf compared to No Prior Inf
 	SuppTable1$Variable[cnt] <- "Y"
 	SuppTable1$FoldChange[cnt] <- paste(round((10^(s0$coef[2,1])-1)*100), "%", sep="")
 	SuppTable1$CI.95[cnt] <- paste(format(10^c0[2,1], digits=2), format(10^c0[2,2], digits=2), sep=", ")
 	pvals[cnt] <- s0$coef[2,4]
 	
 	## Moderna compared to Pfizer
 	SuppTable1$Variable[cnt+1] <- "mRNA-1273"
 	SuppTable1$FoldChange[cnt+1] <- paste(round((10^(s0$coef[3,1])-1)*100), "%", sep="")
	SuppTable1$CI.95[cnt+1] <- paste(format(10^c0[3,1], digits=2), format(10^c0[3,2], digits=2), sep=", ")
 	pvals[cnt+1] <- s0$coef[3,4]
 	
 	## child compared to infant
 	SuppTable1$Variable[cnt+2] <- "24-59 months"
 	SuppTable1$FoldChange[cnt+2] <- paste(round((10^(s0$coef[4,1])-1)*100), "%", sep="")
	SuppTable1$CI.95[cnt+2] <- paste(format(10^c0[4,1], digits=2), format(10^c0[4,2], digits=2), sep=", ")
 	pvals[cnt+2] <- s0$coef[4,4]

	cnt <- cnt+3
	
	}

SuppTable1$QValue <- c(format(p.adjust(pvals[1:42], method="fdr"), digits=3, scientific=F),
						format(p.adjust(pvals[c(43,44,45,49,50,51,55,56,57)], method="fdr"), digits=3, scientific=F),
						  format(p.adjust(pvals[c(46,47,48,52,53,54,58,59,60)], method="fdr"), digits=3, scientific=F))
SuppTable1
write.table(SuppTable1, file="SuppTable1.txt", sep="\t", quote=F, row.names=F)


### table with interaction terms
SuppTable3 <- as.data.frame(matrix(NA, ncol=5, nrow=4*length(mfivars)))
names(SuppTable3) <- c("Variant", "Variable", "Coefficient", "CI.95", "QValue")
pvals <- rep(0, 4*length(mfivars))

cnt <- 1
for(i in 1:length(mfivars)){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
 	SuppTable3$Variant[cnt:(cnt+3)] <- rep(svars[i], 4)
 	
 	g0 <- lm(log10(mydata$mfivar) ~ prior_infection + vax_type*age_cat, data=mydata)
 	s0 <- summary(g0)
 	c0 <- confint(g0, level=0.90)
 	
 	## Prior Inf compared to No Prior Inf
 	SuppTable3$Variable[cnt] <- "Y"
 	SuppTable3$Coefficient[cnt] <- format(10^(s0$coef[2,1]), digits=2)
 	SuppTable3$CI.95[cnt] <- paste(format(10^c0[2,1], digits=2), format(10^c0[2,2], digits=2), sep=", ")
 	pvals[cnt] <- s0$coef[2,4]
 	
 	## Moderna compared to Pfizer
 	SuppTable3$Variable[cnt+1] <- "Moderna"
 	SuppTable3$Coefficient[cnt+1] <- format(10^(s0$coef[3,1]), digits=2)
	SuppTable3$CI.95[cnt+1] <- paste(format(10^c0[3,1], digits=2), format(10^c0[3,2], digits=2), sep=", ")
 	pvals[cnt+1] <- s0$coef[3,4]

 	## Older Kids Compared to Younger
 	SuppTable3$Variable[cnt+2] <- "24-59 months"
 	SuppTable3$Coefficient[cnt+2] <- format(10^(s0$coef[4,1]), digits=2)
	SuppTable3$CI.95[cnt+2] <- paste(format(10^c0[4,1], digits=2), format(10^c0[4,2], digits=2), sep=", ")
 	pvals[cnt+2] <- s0$coef[4,4]

 	## Age/Vaccine Interaction Term
 	SuppTable3$Variable[cnt+3] <- "Age/Vaccine Interaction"
 	SuppTable3$Coefficient[cnt+3] <- format(10^(s0$coef[5,1]), digits=2)
	SuppTable3$CI.95[cnt+3] <- paste(format(10^c0[5,1], digits=2), format(10^c0[5,2], digits=2), sep=", ")
 	pvals[cnt+3] <- s0$coef[5,4]

	cnt <- cnt+4
	
	}

SuppTable3$QValue <- format(p.adjust(pvals, method="fdr"), digits=3, scientific=F)
SuppTable3
write.table(SuppTable3, file="SuppTable3.txt", sep="\t", quote=F, row.names=F)



##### Supplemental GLM tables for neut data
mfivars <- c(41,43, 45)
 svars <- toupper(unlist(lapply(strsplit(names(mydata)[mfivars], split="_"), function(x) x[1]))) 
svars

SuppTable2 <- as.data.frame(matrix(NA, ncol=6, nrow=4*length(mfivars)))
names(SuppTable2) <- c("Variant", "ID", "Variable", "FoldChange", "CI.95", "QValue")
pvals <- rep(0, 4*length(mfivars))

cnt <- 1
for(i in 1:length(mfivars)){
 	
 	mydata$mfivar <- mydata[ ,mfivars[i]]
 	SuppTable2$Variant[cnt:(cnt+3)] <- rep(svars[i], 4)
 	SuppTable2$ID[cnt:(cnt+3)] <- c(rep("ID50", 2), rep("ID80", 2))
 	
 	g0 <- lm(log10(mydata$mfivar) ~ prior_infection + vax_type + age_cat, data=mydata)
 	s0 <- summary(g0)
 	c0 <- confint(g0, level=0.90)
 	
 	## ID50, Prior Inf compared to No Prior Inf
 	SuppTable2$Variable[cnt] <- "Y"
 	SuppTable2$FoldChange[cnt] <- paste(round((10^(s0$coef[2,1])-1)*100), "%", sep="")
 	SuppTable2$CI.95[cnt] <- paste(format(10^c0[2,1], digits=2), format(10^c0[2,2], digits=2), sep=", ")
 	pvals[cnt] <- s0$coef[2,4]
 	
 	## ID50, Moderna compared to Pfizer
 	SuppTable2$Variable[cnt+1] <- "Moderna"
 	SuppTable2$FoldChange[cnt+1] <- paste(round((10^(s0$coef[3,1])-1)*100), "%", sep="")
	SuppTable2$CI.95[cnt+1] <- paste(format(10^c0[3,1], digits=2), format(10^c0[3,2], digits=2), sep=", ")
 	pvals[cnt+1] <- s0$coef[3,4]
 	
 	mydata$mfivar <- mydata[ ,(1+mfivars[i])]
 	
 	g0 <- lm(log10(mydata$mfivar) ~ prior_infection + vax_type + age_cat, data=mydata)
 	s0 <- summary(g0)
 	c0 <- confint(g0, level=0.90)
 	
 	## ID80, Prior Inf compared to No Prior Inf
 	SuppTable2$Variable[cnt+2] <- "Y"
 	SuppTable2$FoldChange[cnt+2] <- paste(round((10^(s0$coef[2,1])-1)*100), "%", sep="")
 	SuppTable2$CI.95[cnt+2] <- paste(format(10^c0[2,1], digits=2), format(10^c0[2,2], digits=2), sep=", ")
 	pvals[cnt+2] <- s0$coef[2,4]
 	
 	## ID80, Moderna compared to Pfizer
 	SuppTable2$Variable[cnt+3] <- "Moderna"
 	SuppTable2$FoldChange[cnt+3] <- paste(round((10^(s0$coef[3,1])-1)*100), "%", sep="")
	SuppTable2$CI.95[cnt+3] <- paste(format(10^c0[3,1], digits=2), format(10^c0[3,2], digits=2), sep=", ")
 	pvals[cnt+3] <- s0$coef[3,4]

	cnt <- cnt+4
	
	}

SuppTable2$QValue <- format(p.adjust(pvals, method="fdr"), digits=3, scientific=F)
SuppTable2
write.table(SuppTable2, file="SuppTable2.txt", 
				sep="\t", quote=F, row.names=F)

##############

##############################################################################
# Script to generate Figures 2 - 4 in the manuscript by Bien, Gaynanova, Lederer and Mueller
##########################################################################################

# Export results from large n and small n simulation settings from Matlab
library(R.matlab)

out_large_n = readMat("LargeNExample_09-Apr-2016_p100_n500_nRep21.mat")
out_small_n = readMat("SmallNExample_26-Jan-2016_p100 _500_n50_nRep21.mat")

# Identify the range of kappa and sigma values
kappa = out_large_n$kappaVec
nk = length(kappa)
sigma = out_large_n$sigVec
ns = length(sigma)

library(ggplot2)
library(plyr)

####################################
# Picture for run times (Figure 3)
###################################

# Large n case
meanRunTime = apply(out_large_n$runTimeMat,c(2,3,5), mean)
stdRunTime = apply(out_large_n$runTimeMat,c(2,3,5), sd)
p = out_large_n$pVec
meanRun1 <- data.frame(meanvalue = c(as.vector(meanRunTime[,,1]),as.vector(meanRunTime[,,2])),sdvalue = c(as.vector(stdRunTime[,,1]),as.vector(stdRunTime[,,2])), method = c(rep('qTREX', nk*ns), rep("cTREX", nk*ns)), kappa = rep(kappa, ns*2), sigma = rep(rep(sigma, each = nk),2) , p = rep(p, nk*ns*2), n = rep(out_large_n$n, nk*ns*2))


# Small n case
meanRunTime = apply(out_small_n$runTimeMat,c(1,2,3,5), mean)
stdRunTime = apply(out_small_n$runTimeMat,c(1,2,3,5), sd)
p = out_small_n$pVec
np = length(p)
meanRun2 <- data.frame(meanvalue = c(as.vector(meanRunTime[1,,,1]),as.vector(meanRunTime[1,,,2]),as.vector(meanRunTime[2,,,1]),as.vector(meanRunTime[2,,,2])),sdvalue =  c(as.vector(stdRunTime[1,,,1]),as.vector(stdRunTime[1,,,2]),as.vector(stdRunTime[2,,,1]),as.vector(stdRunTime[2,,,2])), method = rep(c(rep('qTREX', nk*ns), rep("cTREX", nk*ns)), np), kappa = rep(rep(kappa, ns*2), np), sigma = rep(rep(sigma, each = nk),2*np) , p = rep(p, each=nk*ns*2), n = rep(out_small_n$n, nk*ns*2*np))

# Combine both cases and create a figure
meanRun = rbind(meanRun1, meanRun2)
meanRun$setting = paste("n = ", meanRun$n,", p = ", meanRun$p, sep="")

nRep = out_large_n$numRep
meanRun$SE = meanRun$sdvalue/sqrt(nRep)

pdf(file="RuntimeFig_ggplot.pdf", width = 14, height = 7)
pd = position_dodge(0.35)
ggplot(meanRun[meanRun$kappa %in% c(0,0.9),], aes(x = as.factor(sigma), y = meanvalue, fill= factor(kappa), shape = factor(method), group = interaction(kappa, method))) + facet_grid(.~setting) + scale_shape_manual(values = c(21,22))+ scale_fill_grey()+guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "white"))) + ylab("Run time [s]") + scale_x_discrete(expand = waiver())+ xlab(expression(sigma)) + labs(fill = expression(kappa), shape = "Method")+theme_bw(base_size=18)+
  geom_errorbar(aes(ymin=meanvalue-2*SE, ymax=meanvalue+2*SE), width=.1, position=pd)+ geom_point(size = 3, position = pd)
dev.off()


#######################################
# Picture for estimation (Figure 4)
########################################

# Large n case
meanEst = apply(out_large_n$estMat,c(2,3,5), mean)
stdEst = apply(out_large_n$estMat,c(2,3,5), sd)
p = out_large_n$pVec
meanEst1 <- data.frame(meanvalue = c(as.vector(meanEst[,,1]),as.vector(meanEst[,,2])),sdvalue = c(as.vector(stdEst[,,1]),as.vector(stdEst[,,2])), method = c(rep('qTREX', nk*ns), rep("cTREX", nk*ns)), kappa = rep(kappa, ns*2), sigma = rep(rep(sigma, each = nk),2) , p = rep(p, nk*ns*2), n = rep(out_large_n$n, nk*ns*2))


# Small n case
meanEst = apply(out_small_n$estMat,c(1,2,3,5), mean)
stdEst = apply(out_small_n$estMat,c(1,2,3,5), sd)
p = out_small_n$pVec
np = length(p)
meanEst2 <- data.frame(meanvalue = c(as.vector(meanEst[1,,,1]),as.vector(meanEst[1,,,2]),as.vector(meanEst[2,,,1]),as.vector(meanEst[2,,,2])),sdvalue =  c(as.vector(stdEst[1,,,1]),as.vector(stdEst[1,,,2]),as.vector(stdEst[2,,,1]),as.vector(stdEst[2,,,2])), method = rep(c(rep('qTREX', nk*ns), rep("cTREX", nk*ns)), np), kappa = rep(rep(kappa, ns*2), np), sigma = rep(rep(sigma, each = nk),2*np) , p = rep(p, each=nk*ns*2), n = rep(out_small_n$n, nk*ns*2*np))

# Combine both cases and create a figure
meanEst = rbind(meanEst1, meanEst2)
meanEst$setting = paste("n = ", meanEst$n,", p = ", meanEst$p, sep="")

meanEst$SE = meanEst$sdvalue/sqrt(nRep)


pdf(file="EstFig_ggplot.pdf", width = 14, height = 7)
pd = position_dodge(0.35)
ggplot(meanEst[meanEst$kappa %in% c(0,0.9),], aes(x = as.factor(sigma), y = meanvalue, fill= factor(kappa), shape = factor(method), group = interaction(kappa, method)))  + facet_grid(.~setting) + ylab("Estimation error") + scale_shape_manual(values = c(21,22))+ scale_fill_grey()+guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "white")))+ xlab(expression(sigma)) + labs(fill = expression(kappa), shape = "Method")+theme_bw(base_size=18)+
  geom_errorbar(aes(ymin=meanvalue-2*SE, ymax=meanvalue+2*SE), width=.1, position=pd)+ geom_point(size = 3, position = pd)
dev.off()


######################################################
# Picture for random restarts probabilities (Figure 2)
######################################################

# Large n case
eps=1e-4
nl =out_large_n$nl 
nRand=out_large_n$nTREXRep-(nl-1)
# Store the number of agreements for each replication as a function of number of random restarts
agreePnrep = array(0, dim = c(length(out_large_n$pVec),nk,ns, nRand))
iter = 0
for (r in 1:out_large_n$numRep){
  for (s in 1:ns){
    for (k in 1:nk){
      for (pind in 1:length(out_large_n$pVec)){
        iter = iter + 1
        # Extract best solution across all 2p solutions
        objC = min(unlist(out_large_n$funCellC[[iter]]))

       # Extract best solution from first j non-lasso restarts
        for (j in 1:nRand){
          if (j==1)
            # Extract best q-trex solution from first j restarts
            objQ = unlist(out_large_n$funCellQ[[iter]])[1]
          else
            objQ = min(unlist(out_large_n$funCellQ[[iter]])[c(1,(nl+1):(nl+(j-1)))])
          # Compare objectives
          if (abs(objQ-objC)<eps)
            agreePnrep[pind,k,s,j]=agreePnrep[pind,k,s,j]+1
        }
      }
    }
  }
}

Freq1 <- as.data.frame.table(agreePnrep)
colnames(Freq1) <- c("p", "kappa","sigma", "j", "Freq")
Freq1$p <- mapvalues(Freq1$p, from = levels(Freq1$p), to = out_large_n$pVec)
Freq1$kappa <- mapvalues(Freq1$kappa, from = levels(Freq1$kappa), to = out_large_n$kappaVec)
Freq1$sigma <- mapvalues(Freq1$sigma, from = levels(Freq1$sigma), to = out_large_n$sigVec)
Freq1$j <- as.numeric(Freq1$j)
Freq1$n <- rep(out_large_n$n, nk*ns*nRand)

# Small n case
eps=1e-4
nl = out_large_n$nl 
nRand=out_large_n$nTREXRep-(nl-1)
agreePnrep = array(0, dim = c(length(out_small_n$pVec),nk,ns, nRand))
iter = 0
for (r in 1:out_small_n$numRep){
  for (s in 1:ns){
    for (k in 1:nk){
      for (pind in 1:length(out_small_n$pVec)){
        iter = iter + 1
        # Extract best solution across all 2p solutions
        objC = min(unlist(out_small_n$funCellC[[iter]]))
        
        # Extract best solution from first j non-lasso restarts
        for (j in 1:nRand){
          if (j==1)
            # Extract best q-trex solution from first j restarts
            objQ = unlist(out_small_n$funCellQ[[iter]])[1]
          else
            objQ = min(unlist(out_small_n$funCellQ[[iter]])[c(1,(nl+1):(nl+(j-1)))])
          # Compare objectives
          if (abs(objQ-objC)<eps)
            agreePnrep[pind,k,s,j]=agreePnrep[pind,k,s,j]+1
        }
      }
    }
  }
}

Freq2 <- as.data.frame.table(agreePnrep)
colnames(Freq2) <- c("p", "kappa","sigma", "j", "Freq")
Freq2$p <- mapvalues(Freq2$p, from = levels(Freq2$p), to = out_small_n$pVec)
Freq2$kappa <- mapvalues(Freq2$kappa, from = levels(Freq2$kappa), to = out_small_n$kappaVec)
Freq2$sigma <- mapvalues(Freq2$sigma, from = levels(Freq2$sigma), to = out_small_n$sigVec)
Freq2$j <- as.numeric(Freq2$j)
Freq2$n <- rep(out_small_n$n, nk*ns*np*nRand)

# Combine both cases and create a figure
Freq = rbind(Freq1, Freq2)
Freq$setting = paste("n = ", Freq$n,", p = ", Freq$p, sep="")
Freq$kappa = as.numeric(as.character(Freq$kappa))
 
# Function to generate colors as ggplot2 does by default, so that later the order of red and blue can be changed
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
colors = ggplotColours(3)[c(1,3,2)]

pdf(file="FreqFig_ggplot.pdf", width = 14, height = 14)
ggplot(Freq, aes(x = j, y = Freq/out_large_n$numRep, linetype = sigma, col = sigma)) + geom_line(lwd = 2) + facet_grid(kappa~setting, labeller = label_bquote(kappa == .(kappa))) + ylab("Fraction of agreement")+ xlab("Number of random restarts") + scale_linetype_manual(values = c("solid","dashed","dotted"))+labs(col = expression(sigma), linetype = expression(sigma))+theme_bw(base_size=18) + theme(legend.key.width = unit(1.5,"cm")) + scale_color_manual(values = colors)
dev.off()


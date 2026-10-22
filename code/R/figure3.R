# Estimating the mutation prevalence in the ageing human proteome
# S Pechmann (2026)
# Figure 3

library(ggplot2)
library(reshape2)
library(cowplot)
library(limma)

setwd("~/M2/missense")


# --------------------- Figure 3A ---------------------

geno <- read.table("data/processed/mutations_genetic2.txt", header=T)

t <- c(0:10)*10
te_data <- data.frame(time=t, 
                      d1_alpha=geno$AlphaMiss_full[3]*10e-9*(t*365/2),
                      d1_esm1b=geno$ESM1b_full[3]*10e-9*(t*365/2),
                      d1_score=(geno$AlphaMiss_full[3]+geno$ESM1b_full[3])/2 *10e-9*(t*365/2),
                      d2_alpha=geno$AlphaMiss_full[3]*10e-9*(t*365/5),
                      d2_esm1b=geno$ESM1b_full[3]*10e-9*(t*365/5),
                      d2_score=(geno$AlphaMiss_full[3]+geno$ESM1b_full[3])/2 *10e-9*(t*365/5),
                      d3_alpha=geno$AlphaMiss_full[3]*10e-9*(t*365/10),
                      d3_esm1b=geno$ESM1b_full[3]*10e-9*(t*365/10),
                      d3_score=(geno$AlphaMiss_full[3]+geno$ESM1b_full[3])/2 *10e-9*(t*365/10),
                      d4_alpha=geno$AlphaMiss_full[3]*10e-9*(t*365/15),
                      d4_esm1b=geno$ESM1b_full[3]*10e-9*(t*365/15),
                      d4_score=(geno$AlphaMiss_full[3]+geno$ESM1b_full[3])/2 *10e-9*(t*365/15)
                      )


svg("figures/figure3/timeevol_tissues_B4GALT7.svg", height=4, width=4)
#title="B4GALT7"
ggplot(te_data) + 
  geom_hline(yintercept=0.01, linetype='dashed', color="#222222") + 
  geom_ribbon((aes(x=t, ymin=d1_esm1b, ymax=d1_alpha)), alpha=0.5, fill="#A0C06F") + 
  geom_line(aes(x=t, y=d1_score), color="#A0C06F", size=2) + 
  geom_ribbon((aes(x=t, ymin=d2_esm1b, ymax=d2_alpha)), alpha=0.5, fill="#88B04B") + 
  geom_line(aes(x=t, y=d2_score), color="#88B04B", size=2) + 
  geom_ribbon((aes(x=t, ymin=d3_esm1b, ymax=d3_alpha)), alpha=0.5, fill="#5F7B35") + 
  geom_line(aes(x=t, y=d3_score), color="#5F7B35", size=2) + 
  geom_ribbon((aes(x=t, ymin=d4_esm1b, ymax=d4_alpha)), alpha=0.5, fill="#445826") + 
  geom_line(aes(x=t, y=d4_score), color="#445826", size=2) + 
  labs(x="Time [a]", y=expression(P[g]^d))+
  theme_classic() + 
  theme(
    text = element_text(size=20),
    title = element_text(size=12), 
    axis.title = element_text(size=20)
  )

dev.off()


# --------------------- Figure 3B ---------------------

scoreDF <- data.frame(gene=geno$Gene, score=(geno$AlphaMiss_full+geno$ESM1b_full)/2 *10e-9 )
scoreDF <- scoreDF[is.na(scoreDF$score)==F,]

fractionMAT <- matrix(NA, nrow=100, ncol=4)
for (i in 1:100){
  fractionMAT[i, 1] <- sum(  (scoreDF$score * (i*365/2)) > 0.01 )
  fractionMAT[i, 2] <- sum(  (scoreDF$score * (i*365/5)) > 0.01 )
  fractionMAT[i, 3] <- sum(  (scoreDF$score * (i*365/10)) > 0.01 )
  fractionMAT[i, 4] <- sum(  (scoreDF$score * (i*365/15)) > 0.01 )
  }
colnames(fractionMAT) <- c("2d", "5d", "10d", "15d")
fractionDF <- as.data.frame(fractionMAT)
fractionDF$idx <- c(1:100)
f.m <- melt(fractionDF, id='idx')

svg("figures/figure3/timeevol_geno_01.svg", height=4, width=4)

ggplot(f.m) + 
  geom_line(aes(x=idx, y=value, color=variable), size=1.8) + 
  labs(x="Time [a]", y='Number of genes') + 
  scale_color_manual(values=c("#A0C06F", "#88B04B", "#5F7B35", "#445826")) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title=element_blank(), 
    legend.position = c(0.2, 0.9)
  )

dev.off()


# --------------------- Figure 3C ---------------------

arv <- read.table("data/processed/all_reg_vars_log.txt", header=T, sep='\t') 
el_data <- arv[,c(12:52)]
eln <- normalizeBetweenArrays(el_data)


cor.simil <- matrix(NA, ncol(el_data), ncol(el_data))
for (i in 1:ncol(eln)){
  for (j in 1:ncol(eln)){
    cor.simil[i,j] <- cor(eln[,i], eln[,j], use='complete', method='spearman')
  }
}
colnames(cor.simil) <- colnames(eln)
row.names(cor.simil) <- colnames(eln)

elDF <- data.frame(eldiv=rowMeans(cor.simil), risk=colSums(arv$geno_full*10e-9*eln))
elDF$class <- rep("non", 41)
elDF$class[c(which(colnames(el_data)=="cerebral.cortex"), which(colnames(el_data)=="retina"))] <- "neuron"

ct <- cor.test(elDF$eldiv[elDF$class=="non"], elDF$risk[elDF$class=="non"])


svg("figures/figure3/eldiv.svg", height=4, width=4)

ggplot(elDF, aes(x=eldiv, y=risk)) + 
  geom_point(aes(color=class), size=2.5) + 
  geom_smooth(data=elDF[elDF$class=="non",], method='glm', color="#BBBBBB", alpha=0.2, linetype='dashed', size=0.5) + 
  geom_text(x=0.75, y=0.036, hjust=0, label=paste("R =", round(ct$estimate, 3), sep=' '), col="#222222", size=4, stat="unique", show.legend=F) + 
  geom_text(x=0.75, y=0.0357, hjust=0, label=paste("p =", round(ct$p.value, 5), sep=' '), col="#222222", size=4, stat="unique", show.legend=F) + 
  
  scale_color_manual(values=c("#9826bb", "#222222")) + 
  labs(x="Expression similarity", y="Total mutation risk") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.position = c(0.83, 0.9),
    legend.title= element_blank()
  )

dev.off()


# --------------------- Figure 3D ---------------------

geno <- read.table("data/processed/mutations_genetic2.txt", header=T)
pheno <- as.data.frame(read.table("data/processed/mutations_phenotypic.txt", header=T))
row.names(pheno) <- geno$Gene
pheno <- pheno[rowSums(is.na(geno))<12,]    # only entries also in geno (i.e. in database)


a <- 0.5/exp(80^0.5)    # coefficient so that at 100yr a 50% increase
a2 <-0.2/exp(80^0.5)
a3 <- 0.1/exp(80^0.5)
t <- c(1:100)
y <- 1 + a*exp(t^0.5)
y2 <- 1 + a2*exp(t^0.5)
y3 <- 1 + a3*exp(t^0.5)
plot(t,y, t='l')

phenoage <- data.frame(age=c(1:100), alpha=pheno$alpha_ribo_del[2]*y, esm1b=pheno$esm1b_ribo_del[2]*y, score=(pheno$alpha_ribo_del[2]+pheno$esm1b_ribo_del[2])/2*y, score20=(pheno$alpha_ribo_del[2]+pheno$esm1b_ribo_del[2])/2*y2, score10=(pheno$alpha_ribo_del[2]+pheno$esm1b_ribo_del[2])/2*y3)

svg("figures/figure3/phenoage.svg", height=4, width=4)

ggplot(phenoage) + 
  geom_hline(yintercept=phenoage$score[1], linetype='dashed', col="#222222") + 
  geom_hline(yintercept=phenoage$score[1]*1.5, linetype='dashed', col="#222222") + 
  geom_vline(xintercept=80, linetype='dashed', col="#222222") + 
  geom_line(aes(x=age, y=score20), size=1, color="#222222", alpha=0.5) +
  geom_line(aes(x=age, y=score10), size=1, color="#222222", alpha=0.5) +
  geom_line(aes(x=age, y=score), color="#A0C06F", size=2) + 
  geom_ribbon((aes(x=t, ymin=alpha, ymax=esm1b)), alpha=0.5, fill="#A0C06F") + 
  geom_text(x=4, y=phenoage$score[1]*1.5+0.04, hjust=0, label="50% increase", col="#222222", size=5, stat="unique", show.legend=F) + 
  labs(x="Time [a]", y=expression(P[ph]^d))+
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

dev.off()


# --------------------- Figure 3E ---------------------

agescore <- (pheno$alpha_ribo_del+pheno$esm1b_ribo_del)/2
agescoreDF <- data.frame(base=agescore, age10=agescore*1.1, age20=agescore*1.2, age50=agescore*1.5)
colnames(agescoreDF) <- c("0%", "10%", "20%", "50%")
agescore.m <- melt(agescoreDF)

svg("figures/figure3/phenoage_pct.svg", height=4, width=4)

ggplot(agescore.m) + 
  geom_hline(yintercept=1, linetype='dashed', col="#222222") + 
  geom_boxplot(aes(y=value, fill=variable)) + 
  scale_y_continuous(limits=c(0, 3)) + 
  scale_x_continuous(breaks=c(-0.281, -0.095, 0.095, 0.281), labels=c("0%", "10%", "20%", "50%")) + 
  scale_fill_manual(values=c("#CCCCCC", "#A0C06F", "#88B04B", "#5F7B35")) + 
  labs(x="", y=expression(P[ph]^d)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )

dev.off()


# --------------------- Figure 3F ---------------------

arv <- read.table("data/processed/all_reg_vars.txt", header=T, sep='\t')
el_data <- arv[,c(12:52)]
eln <- normalizeBetweenArrays(el_data)

elDF <- data.frame(eldiv=rowMeans(cor.simil), risk=colSums(arv$geno_full*10e-9*eln))
elDF$class <- rep("non", 41)
elDF$class[c(which(colnames(el_data)=="cerebral.cortex"), which(colnames(el_data)=="retina"))] <- "neuron"

phenorisk0  <- ifelse(arv$pheno_full < 1, arv$pheno_full, 1)
phenorisk10 <- ifelse(arv$pheno_full*1.1 < 1, arv$pheno_full*1.1, 1)
phenorisk20 <- ifelse(arv$pheno_full*1.2 < 1, arv$pheno_full*1.2, 1)
phenorisk50 <- ifelse(arv$pheno_full*1.5 < 1, arv$pheno_full*1.5, 1)

#eldiv=rowMeans(cor.simil), 
phenoEL <- data.frame(risk0=colSums(phenorisk0*eln, na.rm=T)/colSums(eln, na.rm=T), risk10=colSums(phenorisk10*eln, na.rm=T)/colSums(eln, na.rm=T), risk20=colSums(phenorisk20*eln, na.rm=T)/colSums(eln, na.rm=T), risk50=colSums(phenorisk50*eln, na.rm=T)/colSums(eln, na.rm=T))
colnames(phenoEL) <- c("0%", "10%", "20%", "50%") 
phenoEL.m <- melt(phenoEL)


svg("figures/figure3/phenoage_fraction.svg", height=4, width=4)

ggplot(phenoEL.m) + 
  geom_density(aes(x=value, fill=variable), alpha=0.7) +   # , position=position_dodge2()
  scale_fill_manual(values=c("#CCCCCC", "#A0C06F", "#88B04B", "#5F7B35")) + 
  labs(x="Fraction pathogenic", y="Density") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    legend.title = element_blank(),
    legend.position = 'top'
  )

dev.off()


# --------------------- Figure 3G ---------------------

geno <- read.table("data/processed/mutations_genetic2.txt", header=T)
pheno <- as.data.frame(read.table("data/processed/mutations_phenotypic.txt", header=T))

rawDF <- data.frame(gene=geno$Gene, geno_alpha=geno$AlphaMiss_full, geno_esm1b=geno$ESM1b_full, pheno_alpha=pheno$alpha_ribo_del, pheno_esm1b=pheno$esm1b_ribo_del)

arv <- read.table("data/processed/all_reg_vars.txt", header=T, sep='\t')
el_data <- arv[,c(12:52)]
eln <- normalizeBetweenArrays(el_data)


idx <- rep(NA, nrow(arv))
for (i in 1:nrow(arv)){
  current_gene <- as.character(arv$Gene[i])
  current_idx <- which(rawDF$gene==current_gene)
  idx[i] <- current_idx
}

arv_raw <- rawDF[idx,]

stomach <- 35
cervix <- 6

t <- c(1:100)
score_d6_alpha <- matrix(NA, nrow=nrow(eln), ncol=100)
score_d6 <- matrix(NA, nrow=nrow(eln), ncol=100)
score_d6_esm1b <- matrix(NA, nrow=nrow(eln), ncol=100)
for (i in 1:100){
  score_d6[,i] <-  arv$geno_full*10e-9*(i*365/6) * eln[,cervix] 
  score_d6_alpha[,i] <-  arv_raw$geno_alpha*10e-9*(i*365/6) * eln[,cervix] 
  score_d6_esm1b[,i] <-  arv_raw$geno_esm1b*10e-9*(i*365/6) * eln[,cervix] 
}

a <- 0.5/exp(80^0.5)    # coefficient so that at 100yr a 50% increase
y <- 1 + a*exp(t^0.5)

pheno_50 <- matrix(NA, nrow(eln), ncol=100)
pheno_50_alpha <- matrix(NA, nrow(eln), ncol=100)
pheno_50_esm1b <- matrix(NA, nrow(eln), ncol=100)

for (i in 1:100){
  y <- 1 + (0.5/exp(80^0.5)) * exp(i^0.5)
  #phenorisk <- arv$pheno_full * y 
  phenorisk  <- ifelse(arv$pheno_full * y < 1, arv$pheno_full * y, 1)
  phenorisk_alpha  <- ifelse(arv_raw$pheno_alpha * y < 1, arv_raw$pheno_alpha * y, 1)
  phenorisk_esm1b  <- ifelse(arv_raw$pheno_esm1b * y < 1, arv_raw$pheno_esm1b * y, 1)
  
  pheno_50[,i] <- phenorisk * eln[,cervix]
  pheno_50_alpha[,i] <- phenorisk_alpha * eln[,cervix]
  pheno_50_esm1b[,i] <- phenorisk_esm1b * eln[,cervix]
  
}


combinedDF <- data.frame(age=t, geno=colSums(score_d6)/sum(eln[,cervix]), geno_alpha=colSums(score_d6_alpha)/sum(eln[,cervix]), geno_esm1b=colSums(score_d6_esm1b)/sum(eln[,cervix]) , pheno=colSums(pheno_50, na.rm=T)/sum(eln[,cervix]), pheno_alpha=colSums(pheno_50_alpha, na.rm=T)/sum(eln[,cervix]), pheno_esm1b=colSums(pheno_50_esm1b, na.rm=T)/sum(eln[,cervix]) )

svg("figures/figure3/combined_fraction.svg", height=4, width=4)

ggplot(combinedDF) + 
  geom_ribbon((aes(x=t, ymin=geno_esm1b, ymax=geno_alpha)), alpha=0.5, fill="#445826") + 
  geom_line(aes(x=t, y=geno), color="#445826", size=1.5) + 
  geom_ribbon((aes(x=t, ymin=pheno_alpha, ymax=pheno_esm1b)), alpha=0.5, fill="#88B04B") + 
  geom_line(aes(x=t, y=pheno), color="#88B04B", size=1.5) + 
  geom_line(aes(x=t, y=geno+pheno), col=2, size=1, linetype='dashed') + 
  labs(x="Time [a]", y="Fraction pathogenic") + 
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

dev.off()
                         
                       
# --------------------- Figure 3H ---------------------

resi_geno <- read.table("data/processed/resilience_geno.txt", header=T)
resi_pheno <- read.table("data/processed/resilience_pheno.txt", header=T)
resi_combi <- read.table("data/processed/resilience_combi.txt", header=T)  

t_d2 <- (c(0:10)*10*365/2)
t_d5 <- (c(0:10)*10*365/5)

df_geno <- data.frame(t=c(), div=c(), R=c(), class=c())

for (i in 1:11){
  current_d2 <- t_d2[i]
  current_d5 <- t_d5[i]
  
  rr <- rep(NA, 3)
  r2 <- rep(NA, 3)
  for (j in 1:3){
    tmp <- resi_geno[resi_geno$divs <= current_d2 & resi_geno$rep==(j-1), ]
    tmp2 <- resi_geno[resi_geno$divs <= current_d5 & resi_geno$rep==(j-1), ]

    rr[j] <- tmp$R_geno[nrow(tmp)]
    r2[j] <- tmp2$R_geno[nrow(tmp2)]
    
  }
  m <- mean(rr)
  m2 <- mean(r2)
  current_df <- data.frame(t=(i-1)*10, div=current_d2, R=m, class="d2")
  current_df2 <- data.frame(t=(i-1)*10, div=current_d5, R=m2, class="d5")

  df_geno <- rbind(df_geno, current_df, current_df2)

}

df_pheno <- data.frame(t=resi_pheno$time, div=rep(NA, nrow(resi_pheno)), R=resi_pheno$R, class=rep("pheno", nrow(resi_pheno)) )
df_combi_d2 <- data.frame(t=c(0:10)*10 , div=rep(NA, sum(resi_combi$tdiv=='d2')/3) , R=colMeans(matrix(resi_combi$R_combi[resi_combi$tdiv=='d2'], ncol=11, nrow=3)) , class=rep("c_d2", sum(resi_combi$tdiv=='d2')/3) )
df_combi_d5 <- data.frame(t=c(0:10)*10 , div=rep(NA, sum(resi_combi$tdiv=='d5')/3) , R=colMeans(matrix(resi_combi$R_combi[resi_combi$tdiv=='d5'], ncol=11, nrow=3)) , class=rep("c_d5", sum(resi_combi$tdiv=='d5')/3) )


df <- rbind(df_geno, df_pheno, df_combi_d2, df_combi_d5)

tmp <- rep(NA, nrow(df))
tmp[df$class=="d2"] <- "geno (2d)"
tmp[df$class=="d5"] <- "geno (5d)"
tmp[df$class=="c_d5"] <- "geno (5d) + pheno"
tmp[df$class=="c_d2"] <- "geno (2d) + pheno"
tmp[df$class=="pheno"] <- "pheno"

df$class <- tmp
df$class <- factor(df$class, levels=c("geno (2d)", "geno (5d)", "pheno", "geno (2d) + pheno", "geno (5d) + pheno"))

svg("figures/figure3/network_resilience.svg", height=4, width=4)

ggplot(df) + 
  geom_hline(yintercept=df$R[df$class=="geno (2d)"][1], linetype='dashed', col="#222222") + 
  geom_hline(yintercept=df$R[df$class=="pheno"][1], linetype='dashed', col="#777777") + 
  geom_line(aes(x=t, y=R, col=class), size=1.5) + 
  scale_color_manual(values=c("#5F7B35","#88B04B", "#777777", "#6b4c3f", "#ad8574")) + 
  labs(x="Time [a]", y="Network resilience") + 
  scale_y_continuous(limits=c(0.29, 0.45)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title = element_blank(),
    legend.position = c(0.395, 0.25)
  )

dev.off()


# --------------------- Figure 3J ---------------------

r0 <- df$R[1]

diff_d2 <- df_geno$R[df_geno$class=="d2"]/r0
diff_d5 <- df_geno$R[df_geno$class=="d5"]/r0

diff_c_d2 <- df_combi_d2$R / df_pheno$R
diff_c_d5 <- df_combi_d5$R / df_pheno$R

diffdf <- data.frame(t=c(0:10)*10, geno_d2=diff_d2, geno_d5=diff_d5, pheno_d2=diff_c_d2, pheno_d5=diff_c_d5)
diffdf.m <- melt(diffdf, id='t')


svg("figures/figure3/network_diff.svg", height=4, width=4)

ggplot(diffdf.m) + 
  geom_line(aes(x=t, y=value, color=variable), size=1.5) + 
  labs(x="Time [a]", y="Relative decline") + 
  scale_color_manual(values=c("#5F7B35","#88B04B", "#6b4c3f", "#ad8574")) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title=element_blank(),
    legend.position = c(0.3, 0.2)
  )

dev.off()





# --------------------- Figure S3E ---------------------


resi_pheno <- read.table("data/processed/resilience_pheno.txt", header=T)
resi_pheno.c1 <- read.table("data/processed/resilience_pheno_const.txt", header=T)
resi_pheno.c2 <- read.table("data/processed/resilience_pheno_const2.txt", header=T)

pheno.si <- data.frame(t=resi_pheno$time, P_ph_d=resi_pheno$R_pheno/resi_pheno$R_pheno[1], C_17=resi_pheno.c2$R_pheno/resi_pheno.c2$R_pheno[1], C_50=resi_pheno.c1$R_pheno/resi_pheno.c1$R_pheno[1])
pheno.si.m <- melt(pheno.si, id='t')


svg("figures/figure3/suppl_pheno_constant.svg", height=4, width=4)

ggplot(pheno.si.m) + 
  geom_line(aes(x=t, y=value, color=variable), size=2) + 
  labs(x="Time [a]", y="Relative network resilience") + 
  scale_color_manual(values=c("#9826bb", "#739bd0", "#0f4c81")) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title = element_blank(), 
    legend.position = c(0.25, 0.25)
  )

dev.off()



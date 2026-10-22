# Estimating the mutation prevalence in the ageing human proteome
# S Pechmann (2026)
# Figure 4


library(ggplot2)
library(reshape2)
library(cowplot)
library(limma)
library(edgeR)

setwd("~/missense")

# --------------------- Figure 4A ---------------------
arv <- read.table("data/processed/all_reg_vars.txt", header=T, sep='\t')
ess <- as.data.frame(read.table("data/human/essential_genes.txt"))

l = c()
for (i in 1:nrow(ess)){
  current_idx <- which(arv$Gene == as.character(ess$V1[i]))
  l <- c(l, current_idx)
}

ess_genofull <- rbind( data.frame(score=arv$geno_full[l], class=rep("essential",length(arv$geno_full[l]))), 
                       data.frame(score=arv$geno_full[-l], class=rep("non", length(arv$geno_full[-l]))) )

p1 <- ggplot(ess_genofull) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#9826bb", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("essential", "non")) + 
  #scale_y_continuous(limits=c(0, 700) ) + 
  scale_y_continuous(breaks=c(0, 200, 400, 600), labels=c(0, expression(2*e^-7), expression(4*e^-7), expression(6*e^-7) ), limits=c(0, 700) ) + 
  labs(x="", y=expression(P[g]^d)) + 
  theme_classic() + 
  theme(
    text = element_text(size=18), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )
wilcox.test(ess_genofull$score[ess_genofull$class=="essential"], ess_genofull$score[ess_genofull$class=="non"])
  

#---

ess_genoavrg <- rbind( data.frame(score=arv$geno_avrg[l], class=rep("essential",length(arv$geno_avrg[l]))), 
                       data.frame(score=arv$geno_avrg[-l], class=rep("non", length(arv$geno_avrg[-l]))) )

p2 <- ggplot(ess_genoavrg) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#9826bb", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("essential", "non")) + 
  scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4), labels=c(0, expression(1*e^-10), expression(2*e^-10), expression(3*e^-10), expression(4*e^-10) )  ) + 
  labs(x="", y=expression(paste("avrg ",  P[g]^d))) + 
  theme_classic() +
  theme(
    text = element_text(size=18), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )
wilcox.test(ess_genoavrg$score[ess_genoavrg$class=="essential"], ess_genoavrg$score[ess_genoavrg$class=="non"])



#---

ess_genofullEBI <- rbind( data.frame(score=arv$geno_full[l]/arv$PPI_EBI[l], class=rep("essential",length(arv$geno_full[l]))), 
                       data.frame(score=arv$geno_full[-l]/arv$PPI_EBI[-l], class=rep("non", length(arv$geno_full[-l]))) )

p3 <- ggplot(ess_genofullEBI) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#9826bb", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("essential", "non")) + 
  scale_y_continuous(breaks=c(0, 10, 20, 30), labels=c(0, expression(1*e^-8), expression(2*e^-8), expression(3*e^-8) ), limits=c(0, 35)) + 
  labs(x="", y=expression(paste(P[g]^d, " / degree"))) + 
  theme_classic() +
  theme(
    text = element_text(size=18), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )
wilcox.test(ess_genofullEBI$score[ess_genofullEBI$class=="essential"], ess_genofullEBI$score[ess_genofullEBI$class=="non"])


#---

ess_genoavrgEBI <- rbind( data.frame(score=arv$geno_avrg[l]/arv$PPI_EBI[l], class=rep("essential",length(arv$geno_avrg[l]))), 
                          data.frame(score=arv$geno_avrg[-l]/arv$PPI_EBI[-l], class=rep("non", length(arv$geno_avrg[-l])))
)

p4 <- ggplot(ess_genoavrgEBI) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#9826bb", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("essential", "non")) + 
  scale_y_continuous(breaks=c(0, 0.01, 0.02, 0.03), labels=c(0, expression(1*e^-11), expression(2*e^-11), expression(3*e^-11) ) , limits=c(0, 0.03) ) + 
  labs(x="", y=expression(paste(paste("avrg ",P[g]^d), " / degree"))) + 
  theme_classic() +
  theme(
    text = element_text(size=18), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )
wilcox.test(ess_genoavrgEBI$score[ess_genoavrgEBI$class=="essential"], ess_genoavrgEBI$score[ess_genoavrgEBI$class=="non"])



svg("figures/figure4/essential_genes.svg", height=4, width=7)

plot_grid(p1, p2, p3, p4, ncol = 4, align = 'h')

dev.off()


# --------------------- Figure 4B ---------------------

hubs <- as.data.frame(read.table("data/processed/betweenness_centrality.txt", header=T))

gf <- rep(NA, nrow(hubs))
dg <- rep(NA, nrow(hubs))
for (i in 1:nrow(hubs)){
  current_idx <- which(arv$Gene == as.character(hubs$gene[i]))
  if (length(current_idx) > 0) {
    gf[i] <- arv$geno_full[current_idx]
    dg[i] <- arv$PPI_EBI[current_idx]
  }
}
hubs$geno_full <- gf
hubs$ppi_EBI <- dg
top500 <- sort(hubs$betweenness, index.return=T, decreasing=T)$ix[1:500]

boxplot(hubs$geno_full[top500], hubs$geno_full[-top500])
boxplot(hubs$ppi_EBI[top500], hubs$ppi_EBI[-top500])


hubs_genofull <- rbind( data.frame(score=hubs$geno_full[top500], class=rep("bottlenecks",length(hubs$geno_full[top500]))), 
                       data.frame(score=hubs$geno_full[-top500], class=rep("non", length(hubs$geno_full[-top500]))) )

p.geno <- ggplot(hubs_genofull) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#039bd0", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("bottleneck", "non")) + 
  scale_y_continuous(breaks=c(0, 300, 600, 900), labels=c(0, expression(3*e^-7), expression(6*e^-7), expression(9*e^-7) ) ,limits=c(0, 900)) + 
  labs(x="", y=expression(P[g]^d)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )
wilcox.test(hubs_genofull$score[hubs_genofull$class=="bottlenecks"], hubs_genofull$score[hubs_genofull$class=="non"])


#---

hubs_degree <- rbind( data.frame(score=hubs$geno_full[top500]/hubs$ppi_EBI[top500], class=rep("bottlenecks",length(hubs$ppi_EBI[top500]))), 
                       data.frame(score=hubs$geno_full[-top500]/hubs$ppi_EBI[-top500], class=rep("non", length(hubs$ppi_EBI[-top500]))) )

p.degree <- ggplot(hubs_degree ) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#039bd0", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("bottleneck", "non")) + 
  scale_y_continuous(breaks=c(0, 10, 20, 30), labels=c(0, expression(1*e^-8), expression(2*e^-8), expression(3*e^-8) ) ,limits=c(0, 30)) + 
  labs(x="", y=expression(paste(P[g]^d, " / degree"))) + 
  theme_classic() +
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )
wilcox.test(hubs_degree$score[hubs_degree$class=="bottlenecks"], hubs_degree$score[hubs_degree$class=="non"])


svg("figures/figure4/bottlenecks.svg", height=4, width=4)

plot_grid(p.geno, p.degree, ncol = 2, align = 'h')

dev.off()



# --------------------- Figure 4C ---------------------

resi_rand <- read.table("data/processed/resilience_rand_geno.txt", header=T)
tmp <- rep(NA, nrow(resi_rand))
tmp[resi_rand$class=="rand"] <- "random"
tmp[resi_rand$class=="prot"] <- "mut. robustness"
resi_rand$class <- tmp


svg("figures/figure4/network_random.svg", height=4, width=4)

ggplot(resi_rand) + 
  geom_line(aes(x=drop, y=R, color=class), size=2) + 
  scale_color_manual(values=c("#ad8574", "#777777")) + 
  labs(x="Number damaged genes", y="Network resilience") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title = element_blank(), 
    legend.position = c(0.65, 0.85)
  )

dev.off()



# --------------------- Figure 4D ---------------------

arv <- read.table("data/processed/all_reg_vars.txt", header=T, sep='\t') 
el_data <- arv[,c(12:52)]
eln <- normalizeBetweenArrays(el_data)

y <- eln[,c(32, 7, 35, 31,  6, 27, 16, 30)]
rownames(y) <- arv$Gene
d <- DGEList(counts=y, group=c(1,1,1,1, 2,2,2,2) )
d <- estimateCommonDisp(d)
de <- exactTest(d)
topTags(de)


plot(de$table$logFC, -log(de$table$PValue) )

padj <- p.adjust(de$table$PValue)
DEdf <- data.frame(fc=de$table$logFC, pv=-log(padj), gf=arv$geno_full)

svg("figures/figure4/volcano.svg", height=4, width=4)

ggplot(DEdf) +
  geom_hline(yintercept=-log(0.05), linetype='dashed', col="#222222") + 
  geom_point(aes(x=fc, y=pv, color=log(gf)), size=1.5) + 
  labs(x="logFC", y="-log(pval)") + 
  scale_color_gradient(low="#b2cd8a", high="#546d2e", na.value="grey50") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.position = 'none'
  )

dev.off()




# --------------------- Figure 4E ---------------------

DE2 <- rbind( DEdf[DEdf$fc < -2 & DEdf$pv > -log(0.05), ], 
              DEdf[DEdf$fc > 2 & DEdf$pv > -log(0.05), ]
              )
DE2$class <- ifelse(DE2$fc > 0, "up", "down")

svg("figures/figure4/DE_geno.svg", height=4, width=4)

ggplot(DE2) + 
  geom_boxplot(aes(y=gf, fill=class)) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("high EL \n in SL cells", "high EL \n in LL cells")) + 
  scale_y_continuous(breaks=c(0, 200, 400, 600, 800), labels=c(0, expression(2*e^-7), expression(4*e^-7), expression(6*e^-7), expression(8*e^-7) ), limits=c(0, 800)) + 
  scale_fill_manual(values=c("#a3c373", "#546d2e")) + 
  labs(x="", y=expression(P[g]^d)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

dev.off()

wilcox.test(DE2$gf[DE2$class=="down"], DE2$gf[DE2$class=="up"])



boxplot(arv$geno_full[DEdf$fc < -2 & DEdf$pv < 0.05], arv$geno_full[DEdf$fc > 2 & DEdf$pv < 0.05],  col=c(2, 3), ylim=c(0, 1000))
boxplot(rowMeans(y[tmp$fc > 2, c(1:4)]), rowMeans(y[tmp$fc > 2, c(5:9)]) , ylim=c(0, 100))
boxplot(rowMeans(y[DEdf$fc < -2 & DEdf$pv < 0.05, c(1:4)]), rowMeans(y[DEdf$fc < -2 & DEdf$pv < 0.05, c(5:9)]) , ylim=c(0, 100))
 

# --------------------- Figure 4F ---------------------


hl <- as.data.frame(read.table("data/processed/human_proteinHL.txt", header=T))

low25 <- which( (hl$Bcells <= quantile(hl$Bcells, 0.35, na.rm=T)) & 
  (hl$NKcells <= quantile(hl$NKcells, 0.35, na.rm=T)) & 
  (hl$Hepatocytes <= quantile(hl$Hepatocytes, 0.35, na.rm=T)) & 
  (hl$Monocytes <= quantile(hl$Monocytes, 0.35, na.rm=T)) )

high25 <- which( (hl$Bcells >= quantile(hl$Bcells, 0.65, na.rm=T)) & 
  (hl$NKcells >= quantile(hl$NKcells, 0.65, na.rm=T)) & 
  (hl$Hepatocytes >= quantile(hl$Hepatocytes, 0.65, na.rm=T)) & 
  (hl$Monocytes >= quantile(hl$Monocytes, 0.65, na.rm=T)) )


idxlow <- c()
for (i in 1:length(low25)){
  current_idx <- which(arv$Gene == as.character(hl$gene_name[low25][i]))
  if (length(current_idx) > 0) {
    idxlow <- c(idxlow, current_idx)
  }
}

idxhigh <- c()
for (i in 1:length(high25)){
  current_idx <- which(arv$Gene == as.character(hl$gene_name[high25][i]))
  if (length(current_idx) > 0) {
    idxhigh <- c(idxhigh, current_idx)
  }
}


pfDF <- rbind(data.frame(pheno_full=arv$pheno_full[idxlow], class=rep("SL", length(idxlow))),
              data.frame(pheno_full=arv$pheno_full[idxhigh], class=rep("LL", length(idxhigh))))

p.pf <- ggplot(pfDF) + 
  geom_boxplot(aes(y=pheno_full, fill=class)) + 
  scale_fill_manual(values=c("#a3c373", "#546d2e")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("SL", "LL")) + 
  #scale_y_continuous(breaks=c(0, 500, 1000, 1500), labels=c(0, expression(5*e^-7), expression(1*e^-6), expression(1.5*e^-6) ) ) + 
  scale_y_continuous(limits=c(0, 3) ) + 
  labs(x="", y=expression(P[ph]^d)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

wilcox.test(pfDF$pheno_full[pfDF$class=="SL"], pfDF$pheno_full[pfDF$class=="LL"])
  

paDF <- rbind(data.frame(pheno_full=arv$pheno_avrg[idxlow], class=rep("SL", length(idxlow))),
              data.frame(pheno_full=arv$pheno_avrg[idxhigh], class=rep("LL", length(idxhigh))))

p.pa <- ggplot(paDF) + 
  geom_boxplot(aes(y=pheno_full, fill=class)) + 
  scale_fill_manual(values=c("#a3c373", "#546d2e")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("SL", "LL")) + 
  #scale_y_continuous(breaks=c(0.1, 0.2, 0.3), labels=c(expression(1*e^-10), expression(2*e^-10), expression(3*e^-10) ) ) + 
  labs(x="", y=expression(paste("avrg ",  P[ph]^d))) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )
wilcox.test(paDF$pheno_full[paDF$class=="SL"], paDF$pheno_full[paDF$class=="LL"])



svg("figures/figure4/halflife.svg", height=4, width=4)

plot_grid(p.pf, p.pa, ncol = 2, align = 'h')

dev.off()


#---

elDF <- rbind(data.frame(el=log(rowMeans(eln[idxlow,])),  class=rep("SL", length(idxlow)) ),
              data.frame(el=log(rowMeans(eln[idxhigh,])), class=rep("LL", length(idxhigh)) )
              )

p.el <- ggplot(elDF) + 
  geom_boxplot(aes(y=el, fill=class)) +
  scale_fill_manual(values=c("#a3c373", "#546d2e")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("SL", "LL")) + 
  scale_y_continuous(limits=c(0, 8)) + 
  labs(x="", y="Expression") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )


lenDF <- rbind(data.frame(len=arv$length[idxlow],  class=rep("SL", length(idxlow)) ),
              data.frame(len=arv$length[idxhigh], class=rep("LL", length(idxhigh)) )
              )

p.len <- ggplot(lenDF) + 
  geom_boxplot(aes(y=len, fill=class)) +
  scale_fill_manual(values=c("#a3c373", "#546d2e")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("SL", "LL")) + 
  scale_y_continuous(limits=c(0, 7000)) + 
  labs(x="", y="Length") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

svg("figures/figure4/halflife_ellen.svg", height=4, width=2.5)

plot_grid(p.el, p.len, ncol = 1, align = 'v')

dev.off()


# --------------------- Figure 4G ---------------------

hpo_pheno <- as.data.frame(read.table("data/processed/hpo_phenotype.txt", header=T))
rr <- as.data.frame(read.table("data/processed/hpo_rand.txt", header=T))

phDF <- data.frame(n=hpo_pheno$n_score[hpo_pheno$n_score >= 25 & hpo_pheno$n_score <=300],
                   score=hpo_pheno$score[hpo_pheno$n_score >= 25 & hpo_pheno$n_score <=300])
ph.rand <- data.frame(n=rr$n[rr$n >=25 & rr$n <= 300],
                      median=rr$median[rr$n >=25 & rr$n <= 300],
                      q95=rr$q95[rr$n >=25 & rr$n <= 300],
                      q99=rr$q99[rr$n >=25 & rr$n <= 300])
ph.rand.m <- melt(ph.rand, id='n')


svg("figures/figure4/hpo_pheno.svg", height=4, width=5)

ggplot(phDF) + 
  geom_point(aes(x=n, y=score), color="#222222") + 
  geom_line(data=ph.rand.m, aes(x=n, y=value, color=variable), size=1.5 ) + 
  scale_color_manual(values=c("#546d2e","#88b04b","#a3c373")) + 
  scale_y_continuous(breaks=c(200, 400, 600, 800), labels=c(expression(2*e^-7), expression(4*e^-7), expression(6*e^-7), expression(8*e^-7)), limits=c(200, 850)) + 
  labs(x="Number of genes", y=expression( paste(paste("< ", P[g]^d), " >") ) ) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title=element_blank(), 
    legend.position = c(0.8, 0.85)
  )

dev.off()



# --------------------- Figure 4H ---------------------

hpo_disease <- as.data.frame(read.table("data/processed/hpo_disease.txt", header=T), sep='\t')

diDF <- data.frame(n=hpo_disease$n_score[hpo_disease$n_score >= 1 & hpo_disease$n_score <=25],
                   score=hpo_disease$score[hpo_disease$n_score >= 1 & hpo_disease$n_score <=25])
di.rand <- data.frame(n=rr$n[rr$n >=1 & rr$n <= 25],
                   median=rr$median[rr$n >=1 & rr$n <= 25],
                   q95=rr$q95[rr$n >=1 & rr$n <= 25],
                   q99=rr$q99[rr$n >=1 & rr$n <= 25])
di.rand.m <- melt(di.rand, id='n')




svg("figures/figure4/hpo_disease.svg", height=4, width=5)

ggplot(diDF) + 
  geom_point(aes(x=n, y=score), color="#222222") + 
  geom_line(data=di.rand.m, aes(x=n, y=value, color=variable), size=1.5 ) + 
  scale_color_manual(values=c("#546d2e","#88b04b","#a3c373")) + 
  scale_y_continuous(breaks=c(0, 1000, 2000, 3000, 4000), labels=c(0, expression(1*e^-6), expression(2*e^-6), expression(3*e^-6), expression(4*e^-6))) + 
  labs(x="Number of genes", y=expression( paste(paste("< ", P[g]^d), " >") ) )   + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title=element_blank(), 
    legend.position = c(0.8, 0.85)
  )

dev.off()


# --------------------- Figure 4J ---------------------

schematic






#-------------------------------------------------------
# 
# Supplement
#
# --------------------- Figure S4A ---------------------

essEL <- rbind(data.frame(el=rowMeans(eln[ess_genofull$class=="essential", ]), class=rep("essential", sum(ess_genofull$class=="essential") ) ),
              data.frame(el=rowMeans(eln[ess_genofull$class=="non", ]), class=rep("non", sum(ess_genofull$class=="non") )) )

p.el <- ggplot(essEL) + 
  geom_boxplot(aes(y=el, fill=class)) + 
  scale_fill_manual(values=c("#9826bb", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("essential", "non")) + 
  labs(x="", y="Expression") + 
  scale_y_continuous(limits=c(0, 80)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

wilcox.test(essEL$el[essEL$class=="essential"], essEL$el[essEL$class=="non"])


essL <- rbind(data.frame(len=arv$length[ess_genofull$class=="essential"], class=rep("essential", sum(ess_genofull$class=="essential") ) ),
               data.frame(len=arv$length[ess_genofull$class=="non"], class=rep("non", sum(ess_genofull$class=="non") )) )

p.l <- ggplot(essL) + 
  geom_boxplot(aes(y=len, fill=class)) + 
  scale_fill_manual(values=c("#9826bb", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("essential", "non")) + 
  labs(x="", y="Length") + 
  scale_y_continuous(limits=c(0, 4000)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

wilcox.test(essL$len[essL$class=="essential"], essL$len[essL$class=="non"])


svg("figures/figure4/essential_suppl.svg", height=4, width=4)

plot_grid(p.el, p.l, ncol = 2, align = 'h')

dev.off()


# --------------------- Figure S4B ---------------------

hubs <- as.data.frame(read.table("data/processed/betweenness_centrality.txt", header=T))
el_data <- arv[,c(12:52)]
eln <- normalizeBetweenArrays(el_data)


el <- rep(NA, nrow(hubs))
len <- rep(NA, nrow(hubs))
for (i in 1:nrow(hubs)){
  current_idx <- which(arv$Gene == as.character(hubs$gene[i]))
  if (length(current_idx) > 0) {
    el[i] <- mean(eln[current_idx, ], na.rm=T)
    len[i] <- arv$len[current_idx]
  }
}
hubs$el <- el
hubs$len <- len
top500 <- sort(hubs$betweenness, index.return=T, decreasing=T)$ix[1:500]

boxplot(hubs$ppi_EBI[top500], hubs$ppi_EBI[-top500])


hubs_el <- rbind( data.frame(score=hubs$el[top500], class=rep("bottlenecks",length(hubs$el[top500]))), 
                        data.frame(score=hubs$el[-top500], class=rep("non", length(hubs$el[-top500]))) )

p.el2 <- ggplot(hubs_el) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#039bd0", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("bottleneck", "non")) + 
  scale_y_continuous(limits=c(0, 80)) + 
  labs(x="", y="Expression") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

wilcox.test(hubs_el$score[hubs_el$class=="bottlenecks"], hubs_el$score[hubs_el$class=="non"])


#---

hubs_len <- rbind( data.frame(score=hubs$len[top500], class=rep("bottlenecks",length(hubs$len[top500]))), 
                      data.frame(score=hubs$len[-top500], class=rep("non", length(hubs$len[-top500]))) )

p.len2 <- ggplot(hubs_len ) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#039bd0", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("bottleneck", "non")) + 
  scale_y_continuous(limits=c(0, 4400)) + 
  labs(x="", y="Length") + 
  theme_classic() +
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

wilcox.test(hubs_len$score[hubs_len$class=="bottlenecks"], hubs_len$score[hubs_len$class=="non"])


svg("figures/figure4/bottlenecks_suppl.svg", height=4, width=4)

plot_grid(p.el2, p.len2, ncol = 2, align = 'h')

dev.off()


# --------------------- Figure S4C ---------------------

R_degree = cor(arv$PPI_EBI, arv$geno_full, use="complete")
R_degree2 = cor(arv$PPI_EBI, arv$geno_full, use="complete", method='spearman')


svg("figures/figure4/degree_genofull.svg", height=4, width=4)

ggplot(arv) + 
  geom_point(aes(x=PPI_EBI, y=geno_full), color="#ad8574") + 
  scale_x_continuous(limits=c(0, 1000)) + 
  scale_y_continuous(breaks=c(0, 1000, 2000, 3000), labels=c(0, expression(1*e^-6), expression(2*e^-6), expression(3*e^-6) ), limits=c(0, 3000)) + 
  geom_text(x=550, y=2900, hjust=0, label=paste("R =", round(R_degree, 2), sep=' '), col="#222222", size=5, stat="unique", show.legend=F) + 
  geom_text(x=550, y=2700, hjust=0, label=paste("R_sp =", round(R_degree2, 2), sep=' '), col="#222222", size=5, stat="unique", show.legend=F) + 
  labs(x="Degree", y=expression(P[g]^d)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20)
  )

dev.off()


# --------------------- Figure S4D ---------------------

# requires DGEList object from 3D

all.m <- melt(d$counts)
p.all <- ggplot(all.m) + 
  geom_boxplot(aes(y=value, fill=Var2)) + 
  scale_x_continuous(breaks=c(-0.33, -0.235, -0.1415, -0.047, 0.047, 0.1415, 0.235, 0.33), labels=colnames(d$counts)) + 
  scale_y_continuous(limits=c(0, 80)) + 
  scale_fill_manual(values=c("#A0C06F", "#88B04B", "#5F7B35", "#445826", "grey80", "grey70", "grey60", "grey50")) + 
  labs(x="", y="Expression", title="all") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    title = element_text(size=12),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1), 
    legend.position = 'none'
  )


lodf <- d$counts[de$table$logFC < -2 & de$table$PValue < 0.05, ]
lo.m <- melt(lodf)
p.lo <- ggplot(lo.m) + 
  geom_boxplot(aes(y=value, fill=Var2)) + 
  scale_x_continuous(breaks=c(-0.33, -0.235, -0.1415, -0.047, 0.047, 0.1415, 0.235, 0.33), labels=colnames(d$counts)) + 
  scale_y_continuous(limits=c(0, 80)) + 
  scale_fill_manual(values=c("#A0C06F", "#88B04B", "#5F7B35", "#445826", "grey80", "grey70", "grey60", "grey50")) + 
  labs(x="", y="Expression", title="down") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    title = element_text(size=12),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1), 
    legend.position = 'none'
  )

updf <- d$counts[de$table$logFC > 2 & de$table$PValue < 0.05, ]
up.m <- melt(updf)
p.up <- ggplot(up.m) + 
  geom_boxplot(aes(y=value, fill=Var2)) + 
  scale_x_continuous(breaks=c(-0.33, -0.235, -0.1415, -0.047, 0.047, 0.1415, 0.235, 0.33), labels=colnames(d$counts)) + 
  scale_y_continuous(limits=c(0, 80)) + 
  scale_fill_manual(values=c("#A0C06F", "#88B04B", "#5F7B35", "#445826", "grey80", "grey70", "grey60", "grey50")) + 
  labs(x="", y="Expression", title="up") + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    title = element_text(size=12),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1), 
    legend.position = 'none'
  )


svg("figures/figure4/diffexp_suppl.svg", height=8, width=6)

plot_grid(p.all, p.lo, p.up, ncol = 1, align = 'v')

dev.off()


# --------------------- Figure S4E ---------------------

hi.m <- melt(hl[high25,c(2:5)])
lo.m <- melt(hl[low25,c(2:5)])


p.lo <- ggplot(lo.m) + 
  geom_boxplot(aes(y=value, fill=variable)) +
  scale_x_continuous(breaks=c(-0.282,  -0.094, 0.094,  0.282), labels=colnames(hl[,c(2:5)])) + 
  scale_y_continuous(limits=c(0, 1300)) + 
  scale_fill_viridis_d(begin=0.3, end=0.6) + 
  labs(x="", y="Protein half-life [min]", title="low") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    title = element_text(size=12),
    axis.title.y = element_text(size=20),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1), 
    legend.position = 'none'
  )
  

p.hi <- ggplot(hi.m) + 
  geom_boxplot(aes(y=value, fill=variable)) +
  scale_x_continuous(breaks=c(-0.282,  -0.094, 0.094,  0.282), labels=colnames(hl[,c(2:5)])) + 
  scale_y_continuous(limits=c(0, 1300)) + 
  scale_fill_viridis_d(begin=0.3, end=0.6) + 
  labs(x="", y="Protein half-life [min]", title="high") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    title = element_text(size=12),
    axis.title.y = element_text(size=20),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1), 
    legend.position = 'none'
  )



svg("figures/figure4/halflife_suppl.svg", height=4, width=5)

plot_grid(p.lo, p.hi, ncol = 2, align = 'h')

dev.off()


# --------------------- Figure S4F ---------------------

genes.pheno <- as.data.frame(read.table("data/processed/hpo_phenogenes.txt", header=T))
genes.disease <- as.data.frame(read.table("data/processed/hpo_diseasegenes.txt", header=T))

idxpheno <- c()
for (i in 1:length(genes.pheno$gene)){
  current_idx <- which(arv$Gene == as.character(genes.pheno$gene[i]))
  if (length(current_idx) > 0) {
    idxpheno <- c(idxpheno, current_idx)
  }
}

idxdisease <- c()
for (i in 1:length(genes.disease$gene)){
  current_idx <- which(arv$Gene == as.character(genes.disease$gene[i]))
  if (length(current_idx) > 0) {
    idxdisease <- c(idxdisease, current_idx)
  }
}



gphdf <- rbind(data.frame(score=arv$geno_full[idxpheno], class=rep("disease", length(arv$geno_full[idxpheno]))), 
               data.frame(score=arv$geno_full[-idxpheno], class=rep("non", length(arv$geno_full[-idxpheno]))))



gdidf <- rbind(data.frame(score=arv$geno_full[idxdisease], class=rep("disease", length(arv$geno_full[idxdisease]))), 
               data.frame(score=arv$geno_full[-idxdisease], class=rep("non", length(arv$geno_full[-idxdisease]))))


svg("figures/figure4/genes_phenodisease.svg", height=4, width=2.8)

ggplot(gdidf) + 
  geom_boxplot(aes(y=score, fill=class)) + 
  scale_fill_manual(values=c("#039bd0", "#777777")) + 
  scale_x_continuous(breaks=c(-0.188, 0.188), labels=c("disease / \n phenotype", "non")) + 
  scale_y_continuous(breaks=c(0, 300, 600, 900), labels=c(0, expression(3*e^-7), expression(6*e^-7), expression(9*e^-7)), limits=c(0, 900)) + 
  labs(x="", y=expression(P[g]^d) ) + 
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    axis.text.x = element_text(angle=30, hjust=1),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(), 
    legend.position = 'none'
  )

dev.off()

wilcox.test(gdidf$score[gdidf$class=="disease"], gdidf$score[gdidf$class=="non"])



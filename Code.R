
### Map
library(dplyr)
library(ggplot2)
library(maps)
library(tidyverse)
library(readxl)

world.dat<-map_data("world")

ggplot() +
  geom_polygon(data=world.dat,aes(x=long,y=lat,group=group),
               fill="#dedede")+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(x=NULL,y=NULL)-> world.map

world.map

# Experimental Site
df<-read.csv(file.choose())

world.map+
  geom_point(data=df,aes(x=Longitude,y=Latitude,color=sufferedECE, size = 10), shape = 21)+
  scale_color_manual(values = c("Cold"="#448DCD"),
                     name="EWEs")+
  theme(legend.position = c(0.1,0.3))

# 16s Expermental site
df<-read.csv(file.choose())
world.map+
  geom_point(data=df,aes(x=Longitude,y=Latitude,color=sufferedECE, size = 10), shape = 21)+
  scale_color_manual(values = c("Cold"="black"),
                     name="EWEs")+
  theme(legend.position = c(0.1,0.3))



### Alpha diveristy

library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)

df <- read.table("otutab_rare.txt",header = T, row.names = 1, check.names = F)
Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)#spe.rich =sobs
index <- as.data.frame(cbind(Shannon, Simpson, Richness))
tdf <- t(df)
tdf<-ceiling(as.data.frame(t(df)))
obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[rownames(index),]
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df ==1) / colSums(df)
write.table(cbind(sample=c(rownames(index)),index),'diversity.index.txt', row.names = F, sep = '\t', quote = F)



### PCoA and statistics

library(ggplot2)
library(vegan)
library(patchwork)

otu <- read.csv('otutab.csv', row.names = 1)
otu <- data.frame(t(otu))
otu.distance <- vegdist(otu)
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p

group <- read.csv('design.csv')
colnames(group) <- c("samples","group")
df <- merge(pc12,group,by="samples")
color=c("#F7AF34","#448DCD")
p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=group,shape=group))+
  theme_bw()+
  geom_point(size=1, shape=19)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
  #guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"))+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#F7AF34","#448DCD"))+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
p1

p1 + stat_ellipse(data=df,geom = "polygon",level=0.9,linetype = 2,size=0.5,aes(fill=group),alpha=0.2,show.legend = T)

write.table(pc12, 'pc12.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')


dt <- read.csv('pc12.csv',header = T,row.names = 1)
head(dt)
pca <- dt
p <- ggplot(data = pca, aes(x = V1,y = V2))+
  geom_point(size = 0.5,
             aes(col = Group))
p

p1 <- p +
  stat_ellipse(aes(color = Group),fill = NA, geom = 'polygon')
p1

p2 <- p1 +
  scale_shape_manual(values=c(16,15)) + 
  scale_color_manual(values=c("#F7AF34","#448DCD")) + 
  labs(x = 'PC1 (10.97%)', y = 'PC2 (10.03%)') + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'right', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p2


p3 <- ggplot(data = pca) +
  geom_density(aes(x = V1, fill=Group), 
               color = 'black', alpha = 0.66, position = 'identity') +
  scale_linetype_manual(values = c('solid','dashed'))
p3


p4 <- p3 +
  scale_fill_manual(values=c("#F7AF34","#448DCD")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p4

p5 <- ggplot(data = pca) +
  geom_density(aes(x = V2, fill=Group),
               color = 'black', alpha = 0.66, position = 'identity') +
  scale_linetype_manual(values = c('solid','dashed')) +
  scale_fill_manual(values=c("#F7AF34","#448DCD")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  coord_flip() 
p5

pp <- p4 + plot_spacer() + p2 + p5 +
  plot_layout(guides = 'collect') &
  theme(legend.position='right') 
pp

pp + plot_layout(widths = c(1.5, 0.5), heights = c(0.5,1.5))

### statistics
otu <- read.delim('otutab.txt', row.names = 1)
otu <- data.frame(t(otu))
group_data = read.table(file="metadatanew.txt", sep="\t", header=T, row=1, check.names=F)
bc <- vegdist(otu)
pcoa <- cmdscale (bc,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pcoa_points =  as.data.frame(pc12)
pcoa_points = data.frame(pc12, group = group_data[row.names(pcoa_points),1])

# Anosim分析
df_anosim <- anosim(bc,pcoa_points$group,permutations = 999)
raw_p_value <- df_anosim$signif
adjusted_p_value <- p.adjust(raw_p_value, method = "BH")
df_anosim$bh_p_value <- adjusted_p_value
print(df_anosim)
write.table(df_anosim, 'df_anosim.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
output_path <- "df_anosim.txt"
file_conn <- file(output_path, "w")
cat("Call:\n", df_anosim$call, "\n", file = file_conn)
cat("ANOSIM statistic R:", df_anosim$statistic, "\n", file = file_conn)
cat("Significance:", df_anosim$signif, "\n", file = file_conn)
cat("Permutation:", df_anosim$permutation, "\n", file = file_conn)
cat("Number of permutations:", df_anosim$perm, "\n", file = file_conn)
close(file_conn)

# Adonis
adonis_result <- adonis2(bc~pcoa_points$group,data=otu,
                         distance = "bray",
                         permutations = 999)
adonis_result
raw_p_values <- adonis_result$`Pr(>F)`
adjusted_p_values <- p.adjust(raw_p_values, method = "BH")
adonis_result$bh_p_values <- adjusted_p_values
print(adonis_result)
write.table(adonis_result, 'adonis_result.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

# MRPP分析
otu <- read.delim('otutab.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
group <- read.delim('metadatanew.txt', sep = '\t', stringsAsFactors = FALSE)
mrpp_result <- mrpp(otu, group$Group, distance = 'bray', permutations = 999)
summary(mrpp_result)
mrpp_result$Pvalue
write.table(data.frame(Group = 'all', Distance = 'Bray-Curtis', A = mrpp_result$A,
                       Observe_delta = mrpp_result$delta, Expect_delta = mrpp_result$E.delta, P_value = mrpp_result$Pvalue),
            'MRPP.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE)
group_name <- unique(group$Group)
mrpp_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group, Group %in% c(group_name[i], group_name[j]))
    otu_ij <- otu[group_ij$names, ]
    mrpp_result_otu_ij <- mrpp(otu_ij, group_ij$Group, permutations = 999, distance = 'bray')	#Bray-Curtis ???????ȣ????? 999 ???û?
    mrpp_result_two <- rbind(mrpp_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', mrpp_result_otu_ij$A, mrpp_result_otu_ij$delta, mrpp_result_otu_ij$E.delta, mrpp_result_otu_ij$Pvalue))
  }
}
mrpp_result_two <- data.frame(mrpp_result_two, stringsAsFactors = FALSE)
names(mrpp_result_two) <- c('group', 'distance', 'A', 'Observe_delta', 'Expect_delta', 'P_value')
mrpp_result_two$P_value <- as.numeric(mrpp_result_two$P_value)
mrpp_result_two$P_adj_BH <- p.adjust(mrpp_result_two$P_value, method = 'BH')
write.table(mrpp_result_two, 'MRPP.result_two.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')




### Microbial inoculant application amount
library(lme4)
# Unit: CFU/kg soil
tbl.Bacterialbiomass <-read.csv(file.choose())
##Finalinoculationquantification
mdFinalinoculationquantification<-lmer(R~ Finalinoculationquantification + (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdFinalinoculationquantification)   ### F value, 0.5525 0.4647
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.003183279 0.3380941
##Log10
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)   ### F value, 0.4056 0.5321
r.squaredGLMM(mdLog10)
##  R2m       R2c
##  0.005335883 0.3466111

# Unit: spores / kg soil 
tbl.Bacterialbiomass <-read.csv(file.choose())
##Finalinoculationquantification
mdFinalinoculationquantification<-lmer(R~ Finalinoculationquantification + (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdFinalinoculationquantification)   ### F value, 3.3041 0.09343
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.05091611 0.2923991
##Log10
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)   ### F value, 2.0972 0.1757
r.squaredGLMM(mdLog10)
##  R2m       R2c
##  0.03847196 0.3015796

# Unit: cells / kg soil 
tbl.Bacterialbiomass <-read.csv(file.choose())
##Finalinoculationquantification
mdFinalinoculationquantification<-lmer(R~ Finalinoculationquantification + (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdFinalinoculationquantification)   ### F value, 0.9291 0.3522
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.02642898 0.4707813
##Log10
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)   ### F value, 0.8021 0.3804
r.squaredGLMM(mdLog10)
##  R2m       R2c
##  0.02416478 0.4674141

# Unit: cfu / m2 soil 
tbl.Bacterialbiomass <-read.csv(file.choose())
##Finalinoculationquantification
mdFinalinoculationquantification<-lmer(R~ Finalinoculationquantification + (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdFinalinoculationquantification)   ### F value, 1.0566   0.42
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.06000594 0.2797504
##Log10
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)   ### F value, 0.0847 0.7965
r.squaredGLMM(mdLog10)
##  R2m       R2c
##  0.008666368 0.3509686

# Unit: CFU/plant soil 
tbl.Bacterialbiomass <-read.csv(file.choose())
##Finalinoculationquantification
mdFinalinoculationquantification<-lmer(R~ Finalinoculationquantification + (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdFinalinoculationquantification)   ### F value, 0.0853 0.7978
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.003487135 0.2123026
##Log10
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)   ### F value, 1.379  0.475
r.squaredGLMM(mdLog10)
##  R2m       R2c
##  0.06230838 0.1505668

# Unit: CFU/seed soil 
tbl.Bacterialbiomass <-read.csv(file.choose())
##Finalinoculationquantification
mdFinalinoculationquantification<-lmer(R~ Finalinoculationquantification + (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdFinalinoculationquantification)   ### F value, 6.233 0.1467
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.1727989 0.6366796
##Log10
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)   ### F value, 6.233 0.1467
r.squaredGLMM(mdLog10)
##  R2m       R2c
##  0.1727989 0.6366797

# Unit: propagules/kgsoil
tbl.Bacterialbiomass <-read.csv(file.choose())
##Finalinoculationquantification
mdFinalinoculationquantification<-lmer(R~ Finalinoculationquantification + (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdFinalinoculationquantification)   ### F value, 0.4619 0.5911
r.squaredGLMM(mdFinalinoculationquantification)
##  R2m       R2c
##  0.03475909 0.3320802
##Log10
mdLog10<-lmer(R~ +Log10+ (1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdLog10)   ### F value, 0.6284 0.5358
r.squaredGLMM(mdLog10)
##  R2m       R2c
##  0.04427324 0.3104697


### Liner mixed effects model
library(lme4) #package for mixed effect model
library(MuMIn)
tbl.Bacterialbiomass <-read.csv(file.choose())
##Study
unique_Study_values <- unique(tbl.Bacterialbiomass$Study)
non_repeated_count <- length(unique_Study_values)  ### Study = 95

###background character
##time
mdtime<-lmer(R~ scale(Time)+Microbialtype +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdtime)   ### F value, 3.3049 0.0697385
r.squaredGLMM(mdtime)
##  R2m       R2c
##  0.04139493 0.3287656

##pH
mdpH<-lmer(R~ scale(Time)+Microbialtype+scale(pH) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdpH)   ### F value, 4.5241 0.040334 *
r.squaredGLMM(mdpH)
##  R2m       R2c
##  0.1733739 0.3559127

##SOC
mdSOC<-lmer(R~ scale(Time)+Microbialtype+ scale(SOC) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdSOC)  ### F value, 3.0147 0.09808
r.squaredGLMM(mdSOC) 
##  R2m       R2c
##  0.07957889 0.3565182

##TN
mdTN<-lmer(R~ scale(Time)+Microbialtype+scale(TN) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdTN)   ### F value, 2.1737 0.1625
r.squaredGLMM(mdTN) 
##  R2m       R2c
##  0.1092205 0.3652737

##NO3
mdNO3<-lmer(R~ scale(Time)+Microbialtype+ scale(NO3) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdNO3)   ### F value, 0.8693 0.3596
r.squaredGLMM(mdNO3)
##  R2m       R2c
## 0.1216774 0.4897478

##NH4
mdNH4<-lmer(R~ scale(Time)+Microbialtype+ scale(NH4) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdNH4)   ### F value, 0.0731 0.7890
r.squaredGLMM(mdNH4)
##  R2m       R2c
##  0.05927191 0.60936

##AP
mdAP<-lmer(R~ scale(Time)+Microbialtype+  scale(AP) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdAP)   ### F value, 0.0879 0.769316 
r.squaredGLMM(mdAP)
##  R2m       R2c
##  0.1318039 0.3547187

##AK
mdAK<-lmer(R~ scale(Time)+Microbialtype+  scale(AK) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdAK)   ### F value, 0.2796 0.6064
r.squaredGLMM(mdAK)
##  R2m       R2c
##  0.02536249 0.2283544

##AN
mdAN<-lmer(R~ scale(Time)+Microbialtype+ scale(AN) +(1|Study), weights=Wr, data=tbl.Bacterialbiomass)
anova(mdAN)   ### F value, 4.4938   0.04317 *
r.squaredGLMM(mdAN)
##  R2m       R2c
##  0.4461712 0.7384364





### Pearson's Correlation
library(tidyverse)
library(patchwork)
library(dplyr) #加载dplyr包
library(ggpmisc) #加载ggpmisc包
library(ggpubr)
library(ggplot2)
library(ggpmisc)
tbl.Bacterialbiomass <-read.csv(file.choose())
########
#lnpH
sum(!is.na(tbl.Bacterialbiomass$lnpH)) ## n = 77
p11 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnpH n = 7" , y="lnBacterialbiomass")

fit <- lm(R ~ lnpH, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnpH", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p11 <- p11 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p11
pdf("lnpH.pdf",width=8,height=8)
p11
dev.off() 


#lnSOC
sum(!is.na(tbl.Bacterialbiomass$lnSOC)) ## n = 78
p12 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnSOC n = 78" , y="lnBacterialbiomass")
fit <- lm(R ~ lnSOC, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnSOC", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p12 <- p12 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p12
pdf("lnSOC.pdf",width=8,height=8)
p12
dev.off() 

#lnTN
sum(!is.na(tbl.Bacterialbiomass$lnTN)) ## n = 57
p13 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnTN n = 57" , y="lnBacterialbiomass")
fit <- lm(R ~ lnTN, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnTN", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p13 <- p13 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p13
pdf("lnTN.pdf",width=8,height=8)
p13
dev.off() 



#lnNO3
sum(!is.na(tbl.Bacterialbiomass$lnNO3)) ## n = 32
p15 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95) + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                 size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnNO3")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnNO3 n = 32" , y="lnBacterialbiomass")
fit <- lm(R ~ lnNO3, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnNO3", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p15 <- p15 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p15
pdf("lnNO3.pdf",width=8,height=8)
p15
dev.off() 

#lnNH4
sum(!is.na(tbl.Bacterialbiomass$lnNH4)) ## n = 31
p16 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnNH4 n = 31" , y="lnBacterialbiomass")
fit <- lm(R ~ lnNH4, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnNH4", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p16 <- p16 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p16
pdf("lnNH4.pdf",width=8,height=8)
p16
dev.off() 

#lnAP
sum(!is.na(tbl.Bacterialbiomass$lnAP)) ## n = 58
p17 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnAP n = 58" , y="lnBacterialbiomass")
fit <- lm(R ~ lnAP, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnAP", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p17 <- p17 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p17
pdf("lnAP.pdf",width=8,height=8)
p17
dev.off() 

#lnAK
sum(!is.na(tbl.Bacterialbiomass$lnAK)) ## n = 44
p18 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnAK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnAK n = 44" , y="lnBacterialbiomass")
fit <- lm(R ~ lnAK, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnAK", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p18 <- p18 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p18
pdf("lnAK.pdf",width=8,height=8)
p18
dev.off() 

#lnAN
sum(!is.na(tbl.Bacterialbiomass$lnAN)) ## n = 33
p19 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95) + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                 size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnAN")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnAN n = 33" , y="lnBacterialbiomass")
fit <- lm(R ~ lnAN, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnAN", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p19 <- p19 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p19
pdf("lnAN.pdf",width=8,height=8)
p19
dev.off() 

#lnUrease
sum(!is.na(tbl.Bacterialbiomass$lnUrease)) ## n = 193
p20 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnUrease)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95) + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                 size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnUrease")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnUrease n = 193" , y="lnBacterialbiomass")
fit <- lm(R ~ lnUrease, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnUrease", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p20 <- p20 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p20
pdf("lnUrease.pdf",width=8,height=8)
p20
dev.off() 

#lnInvertase
sum(!is.na(tbl.Bacterialbiomass$lnInvertase)) ## n = 148
p21 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnInvertase)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95) + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                 size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnInvertase")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnInvertase n = 148" , y="lnBacterialbiomass")
fit <- lm(R ~ lnInvertase, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnInvertase", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p21 <- p21 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p21
pdf("lnInvertase.pdf",width=8,height=8)
p21
dev.off() 

#lncatalase
sum(!is.na(tbl.Bacterialbiomass$lncatalase)) ## n = 131
p22 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lncatalase)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lncatalase")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lncatalase n = 131" , y="lnBacterialbiomass")
fit <- lm(R ~ lncatalase, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lncatalase", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p22 <- p22 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p22
pdf("lncatalase.pdf",width=8,height=8)
p22
dev.off() 

#lnAlkp
sum(!is.na(tbl.Bacterialbiomass$lnAlkp)) ## n = 60
p23 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnAlkp)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnAlkp")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnAlkp n = 60" , y="lnBacterialbiomass")
fit <- lm(R ~ lnAlkp, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnAlkp", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p23 <- p23 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p23
pdf("lnAlkp.pdf",width=8,height=8)
p23
dev.off() 

#lnAcip
sum(!is.na(tbl.Bacterialbiomass$lnAcip)) ## n = 154
p24 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnAcip)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnAcip")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnAcip n = 154" , y="lnBacterialbiomass")
fit <- lm(R ~ lnAcip, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnAcip", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p24 <- p24 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p24
pdf("lnAcip.pdf",width=8,height=8)
p24
dev.off() 

#lnDH
sum(!is.na(tbl.Bacterialbiomass$lnDH)) ## n = 192
p25 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnDH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95) + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                 size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnDH")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnDH n = 192" , y="lnBacterialbiomass")
fit <- lm(R ~ lnDH, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnDH", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p25 <- p25 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p25
pdf("lnDH.pdf",width=8,height=8)
p25
dev.off()

#lnBG
sum(!is.na(tbl.Bacterialbiomass$lnBG)) ## n = 89
p26 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnBG)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnBG")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnBG n = 89" , y="lnBacterialbiomass")
fit <- lm(R ~ lnBG, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnBG", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p26 <- p26 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p26
pdf("lnBG.pdf",width=8,height=8)
p26
dev.off()

#lnFDA
sum(!is.na(tbl.Bacterialbiomass$lnFDA)) ## n = 79
p27 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnFDA)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnFDA")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnFDA n = 79" , y="lnBacterialbiomass")
fit <- lm(R ~ lnFDA, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnFDA", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p27 <- p27 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p27
pdf("lnFDA.pdf",width=8,height=8)
p27
dev.off()

#lnRDW
sum(!is.na(tbl.Bacterialbiomass$lnRDW)) ## n = 135
p28 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnRDW)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnRDW")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRDW n = 135" , y="lnBacterialbiomass")
fit <- lm(R ~ lnRDW, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnRDW", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p28 <- p28 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p28
pdf("lnRDW.pdf",width=8,height=8)
p28
dev.off()

#lnSDW
sum(!is.na(tbl.Bacterialbiomass$lnSDW)) ## n = 135
p29 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnSDW)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnSDW")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnSDW n = 135" , y="lnBacterialbiomass")
fit <- lm(R ~ lnSDW, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnSDW", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p29 <- p29 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p29
pdf("lnSDW.pdf",width=8,height=8)
p29
dev.off()

#lnPDW
sum(!is.na(tbl.Bacterialbiomass$lnPDW)) ## n = 114
p30 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnPDW)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnPDW")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnPDW n = 114" , y="lnBacterialbiomass")
fit <- lm(R ~ lnPDW, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnPDW", "Pr(>|t|)"]
# Benjamini-Hochberg adjusted
p_adjusted <- p.adjust(p_value, method = "BH")
p30 <- p30 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p30
pdf("lnPDW.pdf",width=8,height=8)
p30
dev.off()

#lnPlantheight
sum(!is.na(tbl.Bacterialbiomass$lnPlantheight)) ## n = 68
p31 <- ggplot(tbl.Bacterialbiomass, aes(y=R, x=lnPlantheight)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialbiomass", x="lnPlantheight")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnPlantheight n = 68" , y="lnBacterialbiomass")
fit <- lm(R ~ lnPlantheight, data = tbl.Bacterialbiomass)
p_value <- summary(fit)$coef["lnPlantheight", "Pr(>|t|)"]
p_adjusted <- p.adjust(p_value, method = "BH")
p31 <- p31 + annotate("text", x = 0.05,y = 0.8,label = sprintf("p (BH adjusted) = %.4f", p_adjusted),size = 4,color = "black")
p31
pdf("lnPlantheight.pdf",width=8,height=8)
p31
dev.off()


#######Cohesion
#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}
#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}
#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}
# Read in dataset
## Data should be in a matrix where each row is a sample. 
b <- read.csv("CK.csv",header = T, row.names = 1,sep = ",")
# # Read in custom correlation matrix, if desired. Must set "use.custom.cors" to TRUE
# if(use.custom.cors == T) {
#   custom.cor.mat <- read.csv("CK.csv", header = T, row.names = 1)
#   custom.cor.mat <- as.matrix(custom.cor.mat)
#   #Check that correlation matrix and abundance matrix have the same dimension
#   print(dim(b)[2] == dim(custom.cor.mat)[2])
# }


# Suggested steps to re-format data. At the end of these steps, the data should be in a matrix "c" where there are no empty samples or blank taxon columns. 
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

# Optionally re-order dataset to be in chronological order. Change date format for your data. 
#c <- c[order(as.Date(rownames(c), format = "%m/%d/%Y")), ]

# Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance, but not if matrix c is count data
rowsums.orig <- rowSums(c)

# Based on persistence cutoff, define a cutoff for the number of zeroes allowed in a taxon's distribution

zero.cutoff <- ceiling(0.1 * dim(c)[1])

# Remove taxa that are below the persistence cutoff
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# Remove any samples that no longer have any individuals, due to removing taxa
d <- d[rowSums(d) > 0, ]

# #If using custom correlation matrix, need to remove rows/columns corresponding to the taxa below persistence cutoff
# if(use.custom.cors == T){
#   custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) < (dim(c)[1]-zero.cutoff), apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
# }

# Create relative abundance matrix.  
rel.d <- d / rowsums.orig
# Optionally, check to see what proportion of the community is retained after cutting out taxa
hist(rowSums(rel.d))

# Create observed correlation matrix
cor.mat.true <- cor(rel.d)

# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()
# Run this loop for the null model to get expected pairwise correlations
# Bypass null model if the option to input custom correlation matrix is TRUE
use.custom.cors <- F
tax.shuffle <- T
iter <- 200
if(use.custom.cors == F) {
  ifelse(tax.shuffle, {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create empty matrix of same dimension as rel.d
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #For each otu
        for(j in 1:dim(rel.d)[2]){ 
          # Replace the original taxon vector with a permuted taxon vector
          perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
        }
        
        # Do not randomize focal column 
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
        
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  } , {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create duplicate matrix to shuffle abundances
        perm.rel.d <- rel.d 
        
        #For each taxon
        for(j in 1:dim(rel.d)[1]){ 
          which.replace <- which(rel.d[j, ] > 0 ) 
          # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          
          #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted 
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
        }
        
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  }
  )
}

# Save observed minus expected correlations. Use custom correlations if use.custom.cors = TRUE
ifelse(use.custom.cors == T, {
  obs.exp.cors.mat <- custom.cor.mat.sub}, {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
)

diag(obs.exp.cors.mat) <- 0
#### 
#### Produce desired vectors of connectedness and cohesion 

# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg
####
#### Combine vectors into one list and print 
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", "Negative Cohesion", "Positive Cohesion")

print(output)
write.csv(output, file = "output.csv")



###########Natural connectaity

library(ggraph)
library(vegan)
library(MCL)
library(tidyverse)
library(qgraph)

##Calculate matrix
##read otu table
otutab<-read.table("CK.txt",header = T,row.names=1,sep="\t")
otutab[is.na(otutab)]<-0
##keep 12 otus
counts<-rowSums(otutab>0)
otutab<-otutab[counts>=1,]

comm<-t(otutab)
sp.ra<-colMeans(comm)/1   #relative abundance of each species

###### there are two choices to get the correlation matrix #######
###### choice 1 (slow): calculate correlation matrix from OTU table
cormatrix=matrix(0,ncol(comm),ncol(comm))

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i]>0,comm[k,i],ifelse(comm[k,j]>0,0.01,NA))
    })
    speciesj<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,i]>0,0.01,NA))
    })
    corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormatrix[i,j]<-cormatrix[j,i]<-corij
    
  }}
# ###### end of the two choices of correlation matrix ########

row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.

cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched
cor.cutoff = 0.3
cor.adj = ifelse(abs(network.raw) >= cor.cutoff,1,0)

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(cor.adj, check.names = FALSE), 'CK.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##read otu table
otutab<-read.table("TR.txt",header = T,row.names=1,sep="\t")
otutab[is.na(otutab)]<-0
##keep 12 otus
counts<-rowSums(otutab>0)
otutab<-otutab[counts>=1,]

comm<-t(otutab)
sp.ra<-colMeans(comm)/1   #relative abundance of each species

###### there are two choices to get the correlation matrix #######
###### choice 1 (slow): calculate correlation matrix from OTU table
cormatrix=matrix(0,ncol(comm),ncol(comm))

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i]>0,comm[k,i],ifelse(comm[k,j]>0,0.01,NA))
    })
    speciesj<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,i]>0,0.01,NA))
    })
    corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormatrix[i,j]<-cormatrix[j,i]<-corij
    
  }}
# ###### end of the two choices of correlation matrix ########

row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.

cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched
cor.cutoff = 0.3
cor.adj = ifelse(abs(network.raw) >= cor.cutoff,1,0)
#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(cor.adj, check.names = FALSE), 'TR.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
library(igraph)

#定义函数，计算网络中节点的重要性，并根据重要性依次移出节点后，计算网络的自然连通度
nc <- function(g) {
  natcon <- function(g) {
    N   <- vcount(g)
    adj <- get.adjacency(g)
    evals <- eigen(adj)$value
    nc <- log(mean(exp(evals)))
    nc / (N - log(N))
  }
  nc.attack <- function(g) {
    hubord <- order(rank(betweenness(g)), rank(degree(g)), decreasing=TRUE)
    sapply(1:round(vcount(g)*.8), function(i) {
      ind <- hubord[1:i]
      tmp <- delete_vertices(g, V(g)$name[ind])
      natcon(tmp)
    })
  }
  nc <- nc.attack(g)
  nc
}

##来自 3 种环境的微生物网络鲁棒性评估
library(igraph)
library(ggplot2)

#读取网络邻接矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
adj1 <- read.delim('CK.matrix.txt', row.names = 1, sep = '\t')
adj2 <- read.delim('TR.matrix.txt', row.names = 1, sep = '\t')

g1 <- graph_from_adjacency_matrix(as.matrix(adj1), mode = 'undirected', diag = FALSE)
g2 <- graph_from_adjacency_matrix(as.matrix(adj2), mode = 'undirected', diag = FALSE)

#计算自然连通度
g1_nc <- nc(g1)
g2_nc <- nc(g2)
dat <- data.frame(
  network = c(rep('g1', length(g1_nc)), rep('g2', length(g2_nc))), 
  'Proportion of removes nodes' = c((1:length(g1_nc))/length(g1_nc), (1:length(g2_nc))/length(g2_nc)),
  'Natural Connectivity' = c(g1_nc, g2_nc), check.names = FALSE
)
dat <- subset(dat, `Proportion of removes nodes` <= 0.8)
#作图
ggplot(dat, aes(`Proportion of removes nodes`, `Natural Connectivity`, color = network)) +
  geom_point() +
  theme_bw() +
  labs(x = 'Proportion of removes nodes', y = 'Natural Connectivity')


#############Robust random
##read otu table
otutab<-read.table("CK.txt",header = T,row.names=1,sep="\t")
otutab[is.na(otutab)]<-0
##keep 12 otus
counts<-rowSums(otutab>0)
otutab<-otutab[counts>=1,]

comm<-t(otutab)
sp.ra<-colMeans(comm)/1   #relative abundance of each species

###### there are two choices to get the correlation matrix #######
###### choice 1 (slow): calculate correlation matrix from OTU table
cormatrix=matrix(0,ncol(comm),ncol(comm))

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i]>0,comm[k,i],ifelse(comm[k,j]>0,0.01,NA))
    })
    speciesj<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,i]>0,0.01,NA))
    })
    corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormatrix[i,j]<-cormatrix[j,i]<-corij
    
  }}

# ###### choice 2: directely read in correlation matrix downloaded from MENAP. MENAP downloaded correlation matrix is an upper triangle. Need to make it into a symetric matrix.
# cormatrix.input  <- matrix(0, 1596, 1596)
# cormatrix.input[row(cormatrix.input) >= col(cormatrix.input)] <- scan("CKtra.txt")
# cormatrix <- t(cormatrix.input)
# for (i in 1:nrow(cormatrix)){
#   for (j in 1:ncol(cormatrix)){
#     if (i>j){cormatrix[i,j]<-cormatrix[j,i]}
#   }
# }
# ###### end of the two choices of correlation matrix ########

row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.

cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

## robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep(2014,20),treat=rep("CK",20))

currentdat<-dat1

write.csv(currentdat,"random_removal_result_CK.csv")


#################### Robust target
##read otu table
otutab<-read.table("CK.txt",header = T,row.names=1,sep="\t")
otutab[is.na(otutab)]<-0
##keep 12 otus
counts<-rowSums(otutab>0)
otutab<-otutab[counts>=1,]

comm<-t(otutab)
sp.ra<-colMeans(comm)/1   #relative abundance of each species

###### there are two choices to get the correlation matrix #######
###### choice 1 (slow): calculate correlation matrix from OTU table
cormatrix=matrix(0,ncol(comm),ncol(comm))

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i]>0,comm[k,i],ifelse(comm[k,j]>0,0.01,NA))
    })
    speciesj<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,i]>0,0.01,NA))
    })
    corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormatrix[i,j]<-cormatrix[j,i]<-corij
    
  }}

# ###### choice 2: directely read in correlation matrix downloaded from MENAP. MENAP downloaded correlation matrix is an upper triangle. Need to make it into a symetric matrix.
# cormatrix.input  <- matrix(0, 1596, 1596)
# cormatrix.input[row(cormatrix.input) >= col(cormatrix.input)] <- scan("CKtra.txt")
# cormatrix <- t(cormatrix.input)
# for (i in 1:nrow(cormatrix)){
#   for (j in 1:ncol(cormatrix)){
#     if (i>j){cormatrix[i,j]<-cormatrix[j,i]}
#   }
# }
# ###### end of the two choices of correlation matrix ########

row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.

cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

## robustness simulation 
#input network matrix, number of removed keystone species, keystonespecies list,  and ra of all species
#return the proportion of species remained

#get the keystone species list
node.attri<-read.table("NodeAttribute.txt",header = T,sep="\t")
module.hub<-as.character(node.attri$Name[node.attri$Zi > 2.5 & node.attri$Pi <= 0.62])

#consider cascade effects: removed species will further influence the remaining nodes

rand.remov2.once<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    
    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0) 
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]
      
      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }
      
      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }
      
    }
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=length(module.hub)),
                 year=rep(2014,2*length(module.hub)),treat=rep("warmed",2*length(module.hub)))

currentdat = dat1
write.csv(currentdat,"targeted_CK.csv")






######### (based on Ubuntu 22.04.3 LTS)

### download NCBI data
# mkdir ~/SRAToolkit
# cd ~/SRAToolkit
# wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz
# tar xzvf sratoolkit.3.0.2-ubuntu64.tar.gz
# echo 'export PATH=$PATH:$HOME/SRAToolkit/sratoolkit.3.0.2-ubuntu64/bin ' >> ~/.bashrc
# source ~/.bashrc
# vdb-config --interactive
# cd ~/SRAToolkit
# prefetch --option-file SRR_Acc_List.txt
# fastq-dump --split-3 --gzip ./1/SRR11815940/SRR11815940.sra

### cutadapt
##install virtualenv
# sudo apt install python3-virtualenv
##Create a new virtual environment
# virtualenv cutadapt-venv
# cutadapt-venv/bin/pip --upgrade pip
# cutadapt-venv/bin/pip install cutadapt
##activate the virtual environment
# source cutadapt-venv/bin/activate
# cutadapt –version
##cutadapt
# cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o OF+P3-4_1.fastq.gz -p OF+P3-4_2.fastq.gz SRR20822028_1.fastq.gz SRR20822028_2.fastq.gz

### evolution Tree
#step1: sudo apt-get update -y
#step2: sudo apt-get install -y fasttree
#sudo apt-get update 
#cd /mnt/b/nifH/Process/Fasttree
#mafft 
#fasttree -gtr -nt otus_aligned.fas > otus.nwk













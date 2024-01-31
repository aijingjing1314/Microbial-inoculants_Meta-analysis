
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













#Ms landslide

#PCoA--------
library(SYNCSA)
library(vegan)

#PCoA - functional
#Traits per individuals
traits_ind  <-	read.csv("traits_ind.csv",	header=T,	row.names=1, sep = ";")	
head(traits_ind)
str(traits_ind)

#individual per plot (community)
comm  <-	read.csv("ind_plot.csv",	header = T, sep = ";", row.names = 1)	
head(comm)
comm  <-	t(comm)
head(comm)

library(vegan)
matriz  <-	matrix.x(comm, traits_ind) #package SYNCSA
summary(matriz)

matX  <-	matriz[[3]] # matrix.x has - matriz[[3]] - has the trait-weighted species composition matriz
head(matX)

dist.func  <-	vegdist(matX, method = "euclidean")
dist.func

func <- pcoa(dist.func)
func$values
func1 <- data.frame(func$vectors[,1:2])
biplot(func)

#data to use in plot-------
trat <- c(rep("Control",25), rep("Landslide",25)) #treatments
colvec <- c(rep("orange", 5), rep("#d7191c", 5), rep("orange", 5), rep("#2b83ba", 5), rep("orange", 5),
            rep("orange", 5), rep("#d7191c", 5), rep("orange", 5), rep("#2b83ba", 5), rep("orange", 5))

place <-  as.data.frame(c(rep("2-y", 5), rep("5-y", 5), rep("2-y", 5), rep("39-y", 5), rep("2-y", 5),
                          rep("2-y", 5), rep("5-y", 5), rep("2-y", 5), rep("39-y", 5), rep("2-y", 5)))
place[,1] <- factor(place[,1], levels = c("2-y", "5-y","39-y"))

pointsymbol <- c(rep(21,5), rep(17, 5), rep(21, 5), rep(15, 5), rep(21,5),
                 rep(21,5), rep(17, 5), rep(21, 5), rep(15, 5), rep(21,5)) 

colvec2 <- c(rep("orange", 5), rep("#d7191c", 5), rep("green", 5), rep("#2b83ba", 5), rep("hotpink", 5),
            rep("orange", 5), rep("#d7191c", 5), rep("green", 5), rep("#2b83ba", 5), rep("hotpink", 5))

local1 <- c(rep("MI",5), rep("SM", 5), rep("SA", 5), rep ("QU", 5), rep("CA", 5),
            rep("MI",5), rep("SM", 5), rep("SA", 5), rep ("QU", 5), rep("CA", 5))

pointsymbol2 <- c(rep(21,5), rep(22, 5), rep(23, 5), rep(24, 5), rep(25,5),
                 rep(21,5), rep(22, 5), rep(23, 5), rep(24, 5), rep(25,5))

#PCoA - taxon-----

mat_dado <- read.csv("taxonomic_matrix.csv", header = T, sep = ";", row.names = 3)
head(mat_dado)

matW <- decostand(mat_dado[,-c(1:4)], method = "hellinger")

dist.taxon  <-	vegdist(matW, method = "euclidean")
dist.taxon

taxon <- pcoa(dist.taxon)
taxon$values
taxon1 <- data.frame(taxon$vectors[,1:2])
biplot(taxon)

#data to use in plots----
colvec1 <- c(rep("orange", 10), rep("#d7191c", 10), rep("orange", 10), rep("#2b83ba", 10), rep("orange", 10))

pointsymbol1 <- c(rep(21,10), rep(17, 10), rep(21, 10), rep(15, 10), rep(21,10))

mat_dado$age <- factor(mat_dado$age, levels = c("2-y", "5-y","39-y"))

colvec3 <- c(rep("orange", 10), rep("#d7191c", 10), rep("green", 10), rep("#2b83ba", 10), rep("hotpink", 10))

pointsymbol3 <- c(rep(21,10), rep(22, 10), rep(23, 10), rep(24, 10), rep(25,10))


#two plots for ages-----
#png("PCoA_taxon_func.tiff", width=20, height=10, units="cm", res=400)
par(mfrow = c(1,2), oma = c(0, 0, 0, 0), mai = c(1.3, 0.8, 0.2, 0.2), xpd = NA)
plot(taxon1, type = "n", xlab = "PCoA1 (13%)", ylab = "PCoA2 (6.4%)", xaxt='n', yaxt = 'n')
ordispider(taxon1, groups= mat_dado$cond, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(taxon1, col = colvec1, bg = colvec1, pch = pointsymbol1)
axis(1, at = seq(-0.7,0.3, by = 0.3))
axis(2, at = seq(-0.6,0.5, by = 0.3))
#abline(v=0, col = "lightgray", lty = 5)
#abline(h=0, col = "lightgray", lty = 5)
#legend(-1.5, 0.7, legend = "a", cex = 1.5, bty = "n", xpd = T)

plot(func1, type = "n", xlab = "PCoA1 (75.7%)", ylab = "PCoA2 (9.3%)", xaxt='n', yaxt = 'n')
ordispider(func1, groups= trat, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(func1, col = colvec, bg = colvec, pch = pointsymbol)
axis(1, at = seq(-0.002,0.003), labels = c(-0.002, 0,0.001))
axis(2, at = seq(-0.002,0.001, by = 0.0010))
#abline(v=0, col = "lightgray", lty = 5)
#abline(h=0, col = "lightgray", lty = 5)
#legend(-0.0043, 0.00127, legend = "b", cex = 1.5, bty = "n", xpd = T)

legend(-0.005,-0.003, legend = c("2-y", "5-y", "39-y"), bty = "n",
       col = c("orange", "#d7191c", "#2b83ba"), pch = c(21,17,15), pt.bg = c("orange", "#d7191c", "#2b83ba"), horiz = TRUE, pt.cex = 1.5)
legend(-0.0055,-0.0033, legend = c("control", "landslide"), lty = c(1,2), col = c("black","darkgrey"), bty = "n",horiz = TRUE,
       pt.cex = 1.5)
#dev.off()


#two plots for sites-----

library(scales) #for transparency

#png("PCoA_taxon_func_sites.tiff", width=20, height=10, units="cm", res=400)
par(mfrow = c(1,2), oma = c(0, 0, 0, 0), mai = c(1.5, 0.8, 0.3, 0.2), xpd = NA)
plot(taxon1, type = "n", xlab = "PCoA1 (13%)", ylab = "PCoA2 (6.4%)")
ordispider(taxon1, groups= mat_dado$cond, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(taxon1, col = alpha(colvec3, 0.5), bg = alpha(colvec3, 0.5), pch = pointsymbol3)
#axis(1, at = seq(-0.7,0.3, by = 0.3))
#axis(2, at = seq(-0.6,0.5, by = 0.3))
#abline(v=0, col = "lightgray", lty = 5)
#abline(h=0, col = "lightgray", lty = 5)
#legend(-1.5, 0.7, legend = "a", cex = 1.5, bty = "n", xpd = T)

plot(func1, type = "n", xlab = "PCoA1 (75.7%)", ylab = "PCoA2 (9.3%)")
ordispider(func1, groups= trat, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(func1, col = alpha(colvec2, 0.5), bg = alpha(colvec2, 0.5), pch = pointsymbol2)
#axis(1, at = seq(-0.002,0.003), labels = c(-0.002, 0,0.001))
#axis(2, at = seq(-0.002,0.001, by = 0.0010))
#abline(v=0, col = "lightgray", lty = 5)
#abline(h=0, col = "lightgray", lty = 5)
#legend(-0.0043, 0.00127, legend = "b", cex = 1.5, bty = "n", xpd = T)

legend(-0.006,-0.0033, legend = c("MI", "SM", "SA", "QU", "CA"), bty = "n",
       col = alpha(c("orange", "#d7191c", "green","#2b83ba", "hotpink"), 0.5), pch = c(21,22,23,24,25), 
       pt.bg = alpha(c("orange", "#d7191c", "green","#2b83ba", "hotpink"), 0.5), horiz = TRUE, pt.cex = 1.5)
legend(-0.0058,-0.0037, legend = c("control", "landslide"), lty = c(1,2), col = c("black","darkgrey"), bty = "n",horiz = TRUE,
       pt.cex = 1.5)
#dev.off()


#paired T test for PCoA first score (by local) - Taxonomic-------
taxon1$local <- rownames(taxon1)
taxon1$local <- mat_dado$local
taxon1$cond <- mat_dado$cond
taxon1 <- taxon1[,-2]

taxon_control <- subset(taxon1, taxon1$cond == "control")
colnames(taxon_control)[1] <- "AxisC"

taxon_land <- subset(taxon1, taxon1$cond == "landslide")
colnames(taxon_land)[1] <- "AxisL"

taxon2 <- data.frame(cbind(taxon_control$AxisC, taxon_land$AxisL))
colnames(taxon2) <- c("Control", "Landslide")

taxon2$local <- taxon_land$local

taxonMI <- subset(taxon2, taxon2$local == "MI")
taxonCA <- subset(taxon2, taxon2$local == "CA")
taxonSA <- subset(taxon2, taxon2$local == "SA_")
taxonSM <- subset(taxon2, taxon2$local == "SM")
taxonQU <- subset(taxon2, taxon2$local == "QU")

# t test paired
#MI
t.test(taxonMI$Control, taxonMI$Landslide, paired = T)
#CA
t.test(taxonCA$Control, taxonCA$Landslide, paired = T)
#SA
t.test(taxonSA$Control, taxonSA$Landslide, paired = T)
#SM
t.test(taxonSM$Control, taxonSM$Landslide, paired = T)
#QU
t.test(taxonQU$Control, taxonQU$Landslide, paired = T)

taxon3 <- taxon1
levels(taxon3$local)[levels(taxon3$local) == "SA_"] <- "SA"
taxon3$age <- c(rep("2-y", 10), rep("5-y", 10), rep("2-y", 10), rep("39-y", 10), rep("2-y", 10))

taxon3$local <- factor(taxon3$local, levels = c("CA", "MI", "SA", "SM", "QU"))
taxon3$age <- factor(taxon3$age, levels = c("2-y", "5-y", "39-y"))

taxon3MI <- subset(taxon3, taxon3$local == "MI")

#ploting
library(ggplot2)
library(gridExtra)
library(cowplot)

taxon3$index <- c(1:4,1:5,5,1:5,1:5,1,1,2,2,3,3,4,4,5,5,1:5,1:5,1:5,1:5)

#png("TtestPCoA.tiff", width=20, height=10, units="cm", res=400)
ggplot(taxon3, aes(x = interaction(cond, local), y = Axis.1, colour = age, fill = cond)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(aes(fill = cond)) +
  geom_line(aes(group = interaction(index, local)),
            alpha = 0.5, colour = "black") +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "PCoA Axis1") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0), size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text.x = element_text(size = 10, colour = "black")) +
  scale_fill_manual(values = c("white", "grey")) +
  scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  facet_wrap(~local, nrow = 1, scales = "free", switch = "x") +
  scale_x_discrete(element_blank()) 
#dev.off()


#paired T test for PCoA first score (by local) - Functional--------
func1$local <- local1
func1$cond <- trat
func1 <- func1[,-2]

func_control <- subset(func1, func1$cond == "Control")
colnames(func_control)[1] <- "AxisC"

func_land <- subset(func1, func1$cond == "Landslide")
colnames(func_land)[1] <- "AxisL"

func2 <- data.frame(cbind(func_control$AxisC, func_land$AxisL))
colnames(func2) <- c("Control", "Landslide")

func2$local <- func_land$local

funcMI <- subset(func2, func2$local == "MI")
funcCA <- subset(func2, func2$local == "CA")
funcSA <- subset(func2, func2$local == "SA")
funcSM <- subset(func2, func2$local == "SM")
funcQU <- subset(func2, func2$local == "QU")

# t test paired
#MI
t.test(funcMI$Control, funcMI$Landslide, paired = T)
#CA
t.test(funcCA$Control, funcCA$Landslide, paired = T)
#SA
t.test(funcSA$Control, funcSA$Landslide, paired = T)
#SM
t.test(funcSM$Control, funcSM$Landslide, paired = T)
#QU
t.test(funcQU$Control, funcQU$Landslide, paired = T)

taxon3 <- taxon1
levels(taxon3$local)[levels(taxon3$local) == "SA_"] <- "SA"
taxon3$age <- c(rep("2-y", 10), rep("5-y", 10), rep("2-y", 10), rep("39-y", 10), rep("2-y", 10))

taxon3$local <- factor(taxon3$local, levels = c("CA", "MI", "SA", "SM", "QU"))
taxon3$age <- factor(taxon3$age, levels = c("2-y", "5-y", "39-y"))

taxon3MI <- subset(taxon3, taxon3$local == "MI")

#Traits boxplots-----
traits_plot <- read.csv("traits_cwm.csv", header = T, sep = ";")

levels(traits_plot$cond) <- (c("control", "landslide"))
levels(traits_plot$local)[levels(traits_plot$local) == "ca"] <- "CA"
levels(traits_plot$local)[levels(traits_plot$local) == "mi"] <- "MI"
levels(traits_plot$local)[levels(traits_plot$local) == "qu"] <- "QU"
levels(traits_plot$local)[levels(traits_plot$local) == "sa"] <- "SA"
levels(traits_plot$local)[levels(traits_plot$local) == "sm"] <- "SM"

traits_plot$local <- factor(traits_plot$local, levels = c("CA", "MI", "SA", "SM", "QU"))

traits_plot$age <- c(rep("2-y",10), rep("5-y",10), rep("2-y",10), rep("39-y",10), rep("2-y",10)) 
traits_plot$age <- factor(traits_plot$age, levels = c("2-y", "5-y", "39-y"))

#thricome
library(SYNCSA)

trait_sp <- read.csv("traits_species.csv", header = T, sep = ";")
levels(trait_sp$situacao) <- c("control", "landslide")

tricome <- as.data.frame(trait_sp[,c(1,15)])
row.names(tricome) <- tricome[,1] 
tricome[1] <- NULL 
ind_plot <- table(trait_sp$local_plot, trait_sp$ID)

colnames(ind_plot) == rownames(tricome)

dados_teste <- organize.syncsa(comm = ind_plot, traits = tricome)
ls(dados_teste)

#calculating traits averages at community level = matrix.t

cwm <- matrix.t(dados_teste$community, dados_teste$traits, scale = F)
ls(cwm)
cwm_matrix_t <- cwm$matrix.T
length(cwm_matrix_t)

cwm_matrix_t <- as.data.frame(cwm_matrix_t)
cwm_matrix_t$local_plot <- rownames(cwm_matrix_t)
colnames(cwm_matrix_t)[1] <- "tricome" 

local_cond <- unique(trait_sp[,c(2,4,5)])

library(dplyr)

cwm_tric <- dplyr::left_join(cwm_matrix_t, local_cond, by = "local_plot")

pubs <- aov(tricome ~ local*situacao, data = cwm_tric)
summary(pubs)
TukeyHSD(pubs)

levels(cwm_tric$local)[levels(cwm_tric$local) == "ca"] <- "CA"
levels(cwm_tric$local)[levels(cwm_tric$local) == "mi"] <- "MI"
levels(cwm_tric$local)[levels(cwm_tric$local) == "qu"] <- "QU"
levels(cwm_tric$local)[levels(cwm_tric$local) == "sa"] <- "SA"
levels(cwm_tric$local)[levels(cwm_tric$local) == "sm"] <- "SM"

cwm_tric$local <- factor(cwm_tric$local, levels = c("CA", "MI", "SA", "SM", "QU"))

library(ggplot2)
library(cowplot)
library(gridExtra)

#leaf area
la <- ggplot(traits_plot, aes(local, area, fill = cond)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = expression(paste("Leaf area (", mg^2, ")")),
                     #breaks = seq(10, 100, 25),
                     limits = c(10,100))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        legend.title=element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values = c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "green", "#2b83ba","hotpink")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
la

#leaf margin
leafm <- ggplot(traits_plot, aes(local, margem_0, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf margin",
                     #breaks = seq(0.30, 1.1, 0.20),
                     limits = c(0.20,1.1))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "green", "#2b83ba","hotpink")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafm

#leaf division
leafd <- ggplot(traits_plot, aes(local, divlimbo_0, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf division",
                     #breaks = seq(0.50, 1.1, 0.1),
                     limits = c(0.5,1.1))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafd

#perimeter:area ratio
leafp <- ggplot(traits_plot, aes(local, periarea, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Perimeter:area \n ratio",
                     #breaks = seq(0.04,1.4,0.2),
                     limits = c(0.035,1.5))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafp

#leaf thickness
leaft <- ggplot(traits_plot, aes(local, esp, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf thickness (mm)",
                     #breaks = seq(0.07,0.23,0.05),
                     limits = c(0.06,0.24))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leaft

#specific leaf area
sla <- ggplot(traits_plot, aes(local, sla, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = (expression(atop("Specific leaf", paste("area ")(mg%.%mm^-2)))),
                     #breaks = seq(200,520,50),
                     limits = c(190,520))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
sla

#Leaf dry matter content
leafdry <- ggplot(traits_plot, aes(local, ldmc, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf dry matter \n content (mg)",
                     #breaks = seq(190,415,50),
                     limits = c(189,420))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafdry

#Leaf pubescence
leafpu <- ggplot(cwm_tric, aes(local, tricome, fill = situacao)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf pubescence") +#,
                     #breaks = seq(0.08,1.1,0.2),
                     #limits = c(0.07,1.1))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafpu

#getting the legend from the first plot
legenda <- get_legend(la)

grafico <- plot_grid(leafd  + theme(legend.position="none",
                                    axis.title.x=element_blank(),
                                    axis.text.x=element_blank()),
                     leafp + theme(legend.position="none",
                                   axis.title.x=element_blank(),
                                   axis.text.x=element_blank()), 
                     leaft + theme(legend.position="none",
                                   axis.title.x=element_blank(),
                                   axis.text.x=element_blank()),
                     sla + theme(legend.position="none",
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank()),
                     leafdry + theme(legend.position="none"), 
                     leafpu + theme(legend.position="none"),
                     align = 'hv',
                     labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                     hjust = -0.2,
                     vjust = 1.1,  
                     nrow = 3,
                     scale = 1)
grafico

#png("boxplot_traits.tiff", width=20, height=25, units="cm", res=400)
plot_grid(grafico, legenda, ncol = 1, rel_heights = c(1, .1), rel_widths = c(1,1),
          label_size = 16)
#dev.off()


#two-way anova and tukey a posteriori-----
diver <- read.csv("diversity.csv", header = T, sep = ";")

#taxonomic alpha diversity (H')
taxalfa_anova <- aov(H ~ local*cond, data = diver) 
summary(taxalfa_anova)
TukeyHSD(taxalfa_anova)
#local
tax_alfa <- as.data.frame(TukeyHSD(taxalfa_anova)$local)
tax_alfa1 <- round(tax_alfa[4], digits = 3)
subset(tax_alfa1, tax_alfa1$`p adj` <0.05)

#local:cond
tax_alfa2 <- as.data.frame(TukeyHSD(taxalfa_anova)$`local:cond`)
tax_alfa3 <- round(tax_alfa2[4], digits = 3)
subset(tax_alfa3, tax_alfa3$`p adj` <0.05)


#taxonomic evenness (J)
jalfa_anova <- aov(J ~ local*cond, data = diver) 
summary(jalfa_anova)
TukeyHSD(jalfa_anova)
#local
jax_alfa <- as.data.frame(TukeyHSD(jalfa_anova)$local)
jax_alfa1 <- round(jax_alfa[4], digits = 3)
subset(jax_alfa1, jax_alfa1$`p adj` <0.05)

#local:cond
jax_alfa2 <- as.data.frame(TukeyHSD(jalfa_anova)$`local:cond`)
jax_alfa3 <- round(jax_alfa2[4], digits = 3)
subset(jax_alfa3, jax_alfa3$`p adj` <0.05)

#functional alpha diversity
falfa_anova <- aov(RaoQ ~ local*cond, data = diver) 
summary(falfa_anova)

#functional eveness
fevealfa_anova <- aov(Even ~ local*cond, data = diver) 
summary(fevealfa_anova)
TukeyHSD(fevealfa_anova)
fev_alfa <- as.data.frame(TukeyHSD(fevealfa_anova)$`local:cond`)
fev_alfa1 <- round(fev_alfa[4], digits = 3)
subset(fev_alfa1, fev_alfa1$`p adj` <0.05)

#alpha diversities boxplots-----
diver <- read.csv("diversity.csv", header = T, sep = ";")

diver$age <- factor(c(rep("2-y",5), rep("5-y",5), rep("2-y",5), rep("39-y",5), rep("2-y",5),
                      rep("2-y",5), rep("5-y",5), rep("2-y",5), rep("39-y",5), rep("2-y",5)))
diver$age <- factor(diver$age, levels = c("2-y", "5-y", "39-y"))

levels(diver$local)[levels(diver$local) == "ca"] <- "CA"
levels(diver$local)[levels(diver$local) == "mi"] <- "MI"
levels(diver$local)[levels(diver$local) == "qu"] <- "QU"
levels(diver$local)[levels(diver$local) == "sa "] <- "SA"
levels(diver$local)[levels(diver$local) == "sm"] <- "SM"

diver$local <- factor(diver$local, levels = c("CA", "MI", "SA", "SM", "QU"))

#taxonomic alfa diversity (H')-----------
#local:condition
#taxdiv_a <- ggplot(diver, aes(local, H, fill = cond, colour = age)) + 
#  stat_boxplot(geom ='errorbar') +
#  geom_boxplot() +
#  scale_x_discrete(name = "Site") +
#  scale_y_continuous(name = expression(paste("Taxonomic ", alpha, " diversity (H')"))) +
#  theme_bw() +
#  theme(text = element_text(size = 12),
#        axis.text.x = element_text(size = 12),
#        axis.text.y = element_text(size = 12),
#        legend.position = "bottom",
#        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
#        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
#        legend.text=element_text(size = 12),
#        legend.key.size = unit(1, "cm")) +
#  scale_fill_manual(values=c("white", "grey")) +
#  scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
#  labs(fill = "") +
#  background_grid(major = "", minor = "") +
#  annotate("text", x = 0.8, y = 3.2, label = "a", size = 4) +
#  annotate("text", x = 1.2, y = 2.9, label = "ac", size = 4) +
#  annotate("text", x = 1.8, y = 2.2, label = "b", size = 4) +
#  annotate("text", x = 2.2, y = 2.1, label = "b", size = 4) +
#  annotate("text", x = 2.8, y = 2.3, label = "ab", size = 4) +
#  annotate("text", x = 3.2, y = 1.95, label = "bc", size = 4) +
#  annotate("text", x = 3.8, y = 2.25, label = "bc", size = 4) +
#  annotate("text", x = 4.2, y = 3.2, label = "ab", size = 4) +
#  annotate("text", x = 4.8, y = 2.9, label = "ab", size = 4) +
#  annotate("text", x = 5.2, y = 3, label = "ab", size = 4)
#taxdiv_a

#Local---------------
taxdiv_alocal <- ggplot(diver, aes(local, H)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = expression(paste("Taxonomic ", alpha, " diversity (H')"))) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.3,0.2,0.1,0.1), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values = "white") +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 1, y = 3.3, label = "a", size = 4) +
  annotate("text", x = 2, y = 2.3, label = "b", size = 4) +
  annotate("text", x = 3, y = 2.4, label = "b", size = 4) +
  annotate("text", x = 4, y = 3.3, label = "ab", size = 4) +
  annotate("text", x = 5, y = 3.1, label = "a", size = 4)
taxdiv_alocal

#condition
taxdiv_acond <- ggplot(diver, aes(cond, H)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = expression(paste("Taxonomic ", alpha, " diversity (H')"))) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.3,0.2,0.1,0.1), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values = "white") +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 1, y = 3.3, label = "a", size = 4) +
  annotate("text", x = 2, y = 3.3, label = "b", size = 4)
taxdiv_acond

#########################################################################
#taxonomic evenness (J)
#local:condition
#jeve <- ggplot(diver, aes(local, J, fill = cond, colour = age)) + 
#  stat_boxplot(geom ='errorbar') +
#  geom_boxplot() +
#  scale_x_discrete(name = "Site") +
#  scale_y_continuous(name = "Taxonomic eveness (J)") +
#  theme_bw() +
#  theme(text = element_text(size = 12),
#        axis.text.x = element_text(size = 12),
#        axis.text.y = element_text(size = 12),
#        legend.position = "bottom",
#        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
#        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
#        legend.text=element_text(size = 12),
#        legend.key.size = unit(1, "cm")) +
#  annotate("text", x = 0.8, y = 0.98, label = "a", size = 4) +
#  annotate("text", x = 1.2, y = 1.01, label = "a", size = 4) +
#  annotate("text", x = 1.8, y = 1.02, label = "a", size = 4) +
#  annotate("text", x = 2.2, y = 1.02, label = "ab", size = 4) +
#  annotate("text", x = 2.8, y = 1.01, label = "a", size = 4) +
#  annotate("text", x = 3.2, y = 0.9, label = "b", size = 4) +
#  annotate("text", x = 3.8, y = 1.02, label = "a", size = 4) +
#  annotate("text", x = 4.2, y = 0.985, label = "ab", size = 4) +
#  annotate("text", x = 4.8, y = 1.015, label = "a", size = 4) +
#  annotate("text", x = 5.2, y = 1.01, label = "a", size = 4)
#  scale_fill_manual(values=c("white", "grey")) +
#  scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
#  labs(fill = "") +
#  background_grid(major = "", minor = "") +
#jeve

#local-----------
jeve_site <- ggplot(diver, aes(local, J)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Taxonomic eveness (J)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.3,0.2,0.1,0.1), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_color_manual(values = "white") +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 1, y = 1.03, label = "a", size = 4) +
  annotate("text", x = 2, y = 1.04, label = "a", size = 4) +
  annotate("text", x = 3, y = 1.03, label = "ab", size = 4) +
  annotate("text", x = 4, y = 1.04, label = "a", size = 4) +
  annotate("text", x = 5, y = 1.03, label = "ac", size = 4)
jeve_site

#cond
jeve_cond <- ggplot(diver, aes(cond, J)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Taxonomic eveness (J)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.3,0.2,0.1,0.1), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_color_manual(values = "white") +
  labs(fill = "") +
  background_grid(major = "", minor = "")  +
  annotate("text", x = 1, y = 1.04, label = "a", size = 4) +
  annotate("text", x = 2, y = 1.04, label = "b", size = 4)
jeve_cond

#######################################################################################
#functional alfa diversity (Raos)

#local:condition
#fundiv_a <- ggplot(diver, aes(local, RaoQ, fill = cond)) + 
#  stat_boxplot(geom ='errorbar') +
#  geom_boxplot() +
#  scale_x_discrete(name = "Site") +
#  scale_y_continuous(name = expression(paste("Functional ", alpha, " diversity (Rao's)"))) +
#  theme_bw() +
#  theme(text = element_text(size = 12),
#        axis.text.x = element_text(size = 12),
#        axis.text.y = element_text(size = 12),
#        legend.position = "bottom",
#        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
#        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
#        legend.text=element_text(size = 12),
#        legend.key.size = unit(1, "cm")) +
#  scale_fill_manual(values=c("white", "grey")) +
#  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
#  labs(fill = "") +
#  background_grid(major = "", minor = "")
#fundiv_a


#functional eveness------
#local:condition
feve <- ggplot(diver, aes(local, Even, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Functional eveness") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "right",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 11),
        legend.key.size = unit(0.9, "cm")) +
  scale_fill_manual(values=c("white", "grey")) +
  #scale_color_manual(values = c("orange", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") #+
  #annotate("text", x = 0.8, y = 0.82, label = "a", size = 4) +
  #annotate("text", x = 1.2, y = 0.83, label = "a", size = 4) +
  #annotate("text", x = 1.8, y = 0.9, label = "ab", size = 4) +
  #annotate("text", x = 2.2, y = 0.98, label = "ab", size = 4) +
  #annotate("text", x = 2.8, y = 0.95, label = "bc", size = 4) +
  #annotate("text", x = 3.2, y = 0.84, label = "bc", size = 4) +
  #annotate("text", x = 3.8, y = 0.95, label = "ab", size = 4) +
  #annotate("text", x = 4.2, y = 0.84, label = "ab", size = 4) +
  #annotate("text", x = 4.8, y = 0.87, label = "b", size = 4) +
  #annotate("text", x = 5.2, y = 0.82, label = "b", size = 4)
feve

#getting the legend from the first plot---------
legenda2 <- get_legend(feve)

grafico2 <- plot_grid(taxdiv_alocal + theme(legend.position="none",
                                       axis.title.x=element_blank(),
                                       axis.text.x=element_blank()), 
                      taxdiv_acond + theme(legend.position="none",
                                   axis.title.x=element_blank(),
                                   axis.text.x=element_blank()), 
                      jeve_site  + theme(legend.position="none"),
                      jeve_cond + theme(legend.position="none"),
                      feve + theme(legend.position = "none"),
                      align = 'hv',
                      #labels = c("(a)", "(b)", "(c)", "(d)", "(e)"),
                      hjust = 0,
                      vjust = 1.2,
                      rel_widths = c(1,0.7),
                      nrow = 3,
                      ncol = 2,
                      scale = 1)
grafico2

#png("boxplot_alpha_diver.tiff", width=18, height=20, units="cm", res=400)
plot_grid(grafico2, legenda2, ncol = 1, rel_heights = c(1, .1), rel_widths = c(1,1),
          label_size = 16)
#dev.off()


#Beta diversity-----

library(vegan)

#Functional beta diversity for conditions
fbetac <- betadisper(dist.func, trat) 
#fbetac <- betadisper(dist.func, local1) # Calculate multivariate dispersions
#plot(fbetac)
disfcond_control <- as.data.frame(fbetac$distances[1:25])
colnames(disfcond_control)[1] <- "dist" 
disfcond_landslide <- as.data.frame(fbetac$distances[26:50])
colnames(disfcond_landslide)[1] <- "dist"

trat1 <- c(rep("control",25), rep("landslide",25))

fbdc <- data.frame("cond" = trat1, "dist" = rbind(disfcond_control, disfcond_landslide))
fbdc$divers <- "Functional diversity"
fbdc$site <- c(rep("MI", 5), rep("SM", 5), rep("SA", 5), rep("QU", 5), rep("CA", 5),
               rep("MI", 5), rep("SM", 5), rep("SA", 5), rep("QU", 5), rep("CA", 5))

fbetacl_aov <- aov(dist~cond, data = fbdc)
summary(fbetacl_aov)

#condition
fbetacl <- ggplot(fbdc, aes(cond, dist, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Condition") +
  scale_y_continuous(name = "Centroid distance (PCoA)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
  scale_fill_manual(values=c("white", "white")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 1.2, y = 0.0035, label= "F = 6.67, P<0.05, DF = 1")
fbetacl

#Functional beta diversity for sites

site <- c(rep("CA",10), rep ("MI",10), rep("QU",10), rep("SA",10), rep("SM",10))

dist.func1 <- as.matrix(dist.func)
dist.func1 <- as.dist(dist.func1 [order(rownames(dist.func1 )),order(colnames(dist.func1))])

fbetas <- betadisper(dist.func1, site)  # Calculate multivariate dispersions
plot(fbetas)

disfCA <- as.data.frame(fbetas$distances[1:10])
colnames(disfCA)[1] <- "dist" 
disfMI <- as.data.frame(fbetas$distances[11:20])
colnames(disfMI)[1] <- "dist" 
disfQU <- as.data.frame(fbetas$distances[21:30])
colnames(disfQU)[1] <- "dist" 
disfSA <- as.data.frame(fbetas$distances[31:40])
colnames(disfSA)[1] <- "dist" 
disfSM <- as.data.frame(fbetas$distances[41:50])
colnames(disfSM)[1] <- "dist" 

fbds <- data.frame("site" = site, "dist" = rbind(disfCA, disfMI, disfQU, disfSA, disfSM))
fbds$site <- factor(fbds$site, levels = c("CA", "MI", "SA", "SM", "QU"))
fbds$age <- factor(c(rep("2-y",20), rep("39-y",10), rep("2-y",10), rep("5-y",10)))
fbds$age <- factor(fbds$age, levels = c("2-y", "5-y", "39-y"))
fbds$divers <- "Functional diversity"

fbds$cond <- c(rep("control", 5), rep("landslide", 5), rep("control", 5), rep("landslide", 5), rep("control", 5),
               rep("control", 5), rep("landslide", 5), rep("control", 5), rep("landslide", 5), rep("control", 5))

fbetasl_aov <- aov(dist~site, data = fbds)
summary(fbetasl_aov)
TukeyHSD(fbetasl_aov)
fbetasl2 <- as.data.frame(TukeyHSD(fbetasl_aov)$site)
fbetasl3 <- round(fbetasl2[4], digits = 3)
subset(fbetasl3, fbetasl3$`p adj` >=0.05)
subset(fbetasl3, fbetasl3$`p adj` <0.05)

library(multcompView)
multcompLetters(TukeyHSD(fbetasl_aov)$site[,4])

#local

fbetasl <- ggplot(fbds, aes(site, dist)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Centroid distance (PCoA)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank()) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 1.6, y = 0.0045, label= "F = 4.67, P<0.01, DF = 4") +
  annotate("text", x = 1, y = 0.0017, label= "bc") +
  annotate("text", x = 2, y = 0.0028, label= "ab") +
  annotate("text", x = 3, y = 0.0043, label= "a") +
  annotate("text", x = 4, y = 0.0024, label= "abc") +
  annotate("text", x = 5, y = 0.0011, label= "c")
fbetasl

#Taxonomic beta diversity for conditions

mat_dado1 <- mat_dado[order(mat_dado$cond),]

levels(mat_dado1$cond) <- c("Control", "Landslide")

matW1 <- decostand(mat_dado1[,-c(1:4)], method = "hellinger")

dist.taxon1  <-	vegdist(matW1, method = "euclidean")
dist.taxon1

tbetac <- betadisper(dist.taxon1, mat_dado1$cond)  # Calculate multivariate dispersions
plot(tbetac)

distccon <- as.data.frame(tbetac$distances[1:25])
colnames(distccon)[1] <- "dist" 
distclan <- as.data.frame(tbetac$distances[26:50])
colnames(distclan)[1] <- "dist" 

tbdc <- data.frame("cond" = mat_dado1$cond, "dist" = rbind(distccon,distclan))
tbdc$divers <- "Taxonomic diversity" 

tbdc$site <- c(rep("MI", 5), rep("SM", 5), rep("SA", 5), rep("QU", 5), rep("CA", 5),
               rep("MI", 5), rep("SM", 5), rep("SA", 5), rep("QU", 5), rep("CA", 5))

tbetasc_aov <- aov(dist~cond, data = tbdc)
summary(tbetasc_aov)

tbetasc <- ggplot(tbdc, aes(cond, dist, fill = cond)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Condition") +
  scale_y_continuous(name = "Centroid distance (PCoA)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
  scale_fill_manual(values=c("white", "white")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 1.1, y = 1.1, label= "F = 3.36, P = NS, DF = 1")
tbetasc

#Taxonomic beta diversity for sites

site <- c(rep("CA",10), rep ("MI",10), rep("QU",10), rep("SA",10), rep("SM",10))

dist.taxon2 <- as.matrix(dist.taxon1)
dist.taxon2 <- as.dist(dist.taxon2[order(rownames(dist.taxon2)),order(colnames(dist.taxon2))])

tbetac <- betadisper(dist.taxon2, site)  # Calculate multivariate dispersions
plot(tbetac)

distCA <- as.data.frame(tbetac$distances[1:10])
colnames(distCA)[1] <- "dist" 
distMI <- as.data.frame(tbetac$distances[11:20])
colnames(distMI)[1] <- "dist" 
distQU <- as.data.frame(tbetac$distances[21:30])
colnames(distQU)[1] <- "dist" 
distSA <- as.data.frame(tbetac$distances[31:40])
colnames(distSA)[1] <- "dist" 
distSM <- as.data.frame(tbetac$distances[41:50])
colnames(distSM)[1] <- "dist" 

tbds <- data.frame("site" = site, "dist" = rbind(distCA, distMI, distQU, distSA, distSM))
tbds$site <- factor(tbds$site, levels = c("CA", "MI", "SA", "SM", "QU"))
tbds$age <- factor(c(rep("2-y",20), rep("39-y",10), rep("2-y",10), rep("5-y",10)))
tbds$age <- factor(fbds$age, levels = c("2-y", "5-y", "39-y"))
tbds$divers <- "Taxonomic diversity"

tbds <- merge(tbds, mat_dado,  by=0, all=TRUE)
tbds <- tbds[, c(2,3,8)]

tbds_aov <- aov(dist~site, data = tbds)
summary(tbds_aov)

#local
tbetasl <- ggplot(tbds, aes(site, dist)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Centroid distance (PCoA)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank()) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 2, y = 1.06, label= "F = 0.67, P = NS, DF = 4")
tbetasl

#Plot grid
grafico3 <- plot_grid(#tbetasc +
                      #  theme(axis.title.x = element_blank()),
                      fbetacl,# +
                        #theme(axis.title.x = element_blank()), 
                      #tbetasl,
                      fbetasl +
                        theme(axis.title.y = element_blank()),
                      align = 'hv',
                      labels = c("(a)", "(b)"),
                      hjust = -0.2,
                      vjust = 1.2,  
                      nrow = 1,
                      ncol = 2,
                      scale = 1,
                      rel_widths = c(0.7,1))
grafico3

#png("boxplot_beta_diver.tiff", width=20, height=10, units="cm", res=400)
grafico3
#dev.off()


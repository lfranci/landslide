#Ms landslide
#Code by Luciana Franci

###########################################################################################
#PCoA
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

matriz  <-	matrix.x(comm, traits_ind) #package SYNCSA
summary(matriz)

matX  <-	matriz[[3]] # matrix.x has - matriz[[3]] - has the trait-weighted species composition matriz
head(matX)

dist.func  <-	vegdist(matX, method = "euclidean")
dist.func

#PCoA using vegan package
#func <- cmdscale(dist.func, eig=TRUE)
#func$points

func <- pcoa(dist.func)
func$values
func1 <- data.frame(func$vectors[,1:2])

#NMDS
#funcnmds <- metaMDS(matX, k = 2, distance = "euclidean")
#funcnmds$stress #great! 0.0677
#stressplot(funcnmds) #great

################################################################################################################
#data to use in plot
trat <- c(rep("Control",25), rep("Landslide",25)) #treatments
colvec <- c(rep("#fdae6b", 5), rep("#d7191c", 5), rep("#fdae61", 5), rep("#2b83ba", 5), rep("#fdae6b", 5),
            rep("#fdae6b", 5), rep("#d7191c", 5), rep("#fdae61", 5), rep("#2b83ba", 5), rep("#fdae6b", 5))

place <-  as.data.frame(c(rep("2-y", 5), rep("5-y", 5), rep("2-y", 5), rep("39-y", 5), rep("2-y", 5),
            rep("2-y", 5), rep("5-y", 5), rep("2-y", 5), rep("39-y", 5), rep("2-y", 5)))
place[,1] <- factor(place[,1], levels = c("2-y", "5-y","39-y"))

pointsymbol <- c(rep(21,5), rep(17, 5), rep(21, 5), rep(15, 5), rep(21,5),
                 rep(21,5), rep(17, 5), rep(21, 5), rep(15, 5), rep(21,5)) 

################################################################################################################
#polygon
#plot(funcnmds, type = "n")
#points(funcnmds, display = "sites", col = colvec, bg = colvec, pch = 21)
#ordihull(funcnmds, groups= trat, draw = "polygon", 
#         col = c("orange","lightgreen"), label = F)
################################################################################################################


#spider plot
par(mfrow = c(1,1))
par(xpd = T, mar=c(7, 4, 0, 0) + 0.1)
plot(func1, type = "n", xlab = "PCoA1 (75.7%)", ylab = "PCoA2 (9.3%)")
ordispider(func1, groups= trat, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(func1, col = colvec, bg = colvec, pch = pointsymbol)
legend(-0.0005, -0.0038, legend = c("2-y", "5-y", "39-y"), bty = "n",
                        col = c("#fdae6b", "#d7191c", "#2b83ba"), pch = 15, pt.bg = c("#fdae6b", "#d7191c", "#2b83ba"), horiz = TRUE)
legend(-0.001,-0.0042, legend = c("Control", "Landslide"), lty = c(1,2), col = c("black","darkgrey"), bty = "n",horiz = TRUE)

################################################################################################################
#library(ggplot2)

#data.scores <- as.data.frame(scores(funcnmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
#data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
#data.scores$trat <- trat  #  add the grp variable created earlier
#head(data.scores)  #look at the data

#species.scores <- as.data.frame(scores(funcnmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data

#t_cont <- data.scores[data.scores$trat == "Control", ][chull(data.scores[data.scores$trat == 
 #                                                                  "Control", c("NMDS1", "NMDS2")]), ]  # hull values for Control
#t_land <- data.scores[data.scores$trat == "Landslide", ][chull(data.scores[data.scores$trat == 
#                                                                   "Landslide", c("NMDS1", "NMDS2")]), ]  # hull values for Landslide
#hull.data <- rbind(t_cont, t_land)
#hull.data
#
#ggplot() + 
#  geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = trat, group = trat, linetype = trat), colour = "grey", alpha = 0.30) + # add the convex hulls
#  #geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species), alpha=0.5) +  # add the species labels
#  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, shape = trat, colour = place[,1]), size = 3) + # add the point markers
#  #geom_text(data = data.scores, aes(x = NMDS1, y = NMDS2, label = site), size = 8, vjust = 0, hjust = 0) +  # add the site labels
#  coord_equal() +
#  theme_bw() +
#  guides(colour = guide_legend(override.aes = list(shape = c(15,15,15), size = 7))) +
#  scale_fill_manual(values=c("Control" = "#f0f0f0", "Landslide" = "#bdbdbd")) +
#  scale_color_manual(values=c("2-y" = "#fdcc8a", "5-y" = "#d7191c", "39-y" = "#2b83ba")) +
#  theme(axis.text.x = element_blank(),  # remove x-axis text
#        axis.text.y = element_blank(), # remove y-axis text
#        axis.ticks = element_blank(),  # remove axis ticks
#        axis.title.x = element_text(size=18), # remove x-axis labels
#        axis.title.y = element_text(size=18), # remove y-axis labels
#        panel.background = element_blank(), 
#        panel.grid.major = element_blank(),  #remove major-grid labels
#        panel.grid.minor = element_blank(),  #remove minor-grid labels
#        plot.background = element_blank(),
#        plot.margin = unit(c(0,0.1,0,0), "cm"),
#        legend.title = element_blank(),
#        legend.position = "bottom",
#        legend.box = "vertical")



################################################################################################################
#PCoA - taxon

mat_dado <- read.csv("taxonomic_matrix.csv", header = T, sep = ";", row.names = 3)
head(mat_dado)

matW <- decostand(mat_dado[,-c(1:4)], method = "hellinger")

#matW <- table(mat_dado$local, mat_dado$sp)
#matw <- unclass(matW) #table to matriz

dist.taxon  <-	vegdist(matW, method = "euclidean")
dist.taxon

#taxonPCoA  <-	pcoa(dist.taxon)
#summary(taxonPCoA)
#biplot(taxonPCoA)

#https://estatisticarblog.wordpress.com/2016/05/05/analises-multivariadas/

#using vegan
#taxon <- cmdscale(dist.taxon, k = 3, eig=TRUE)
#taxon_points <- taxon$points[,1:2]

taxon <- pcoa(dist.taxon)
taxon$values
taxon1 <- data.frame(taxon$vectors[,1:2])

#NMDS
#taxon <- metaMDS(matw, k = 3, distance = "bray", binary = T)
#taxon$stress
#taxon_score <- data.frame(scores(taxon, choices = 1))
#stressplot(taxon) #good stress! 0.0944

colvec1 <- c(rep("#fdae6b", 10), rep("#d7191c", 10), rep("#fdae61", 10), rep("#2b83ba", 10), rep("#fdae6b", 10))
            #rep("#fdae6b", 5), rep("#d7191c", 5), rep("#fdae61", 5), rep("#2b83ba", 5), rep("#fdae6b", 4))

pointsymbol1 <- c(rep(21,10), rep(17, 10), rep(21, 10), rep(15, 10), rep(21,10))
                 #rep(21,5), rep(17, 5), rep(21, 5), rep(15, 5), rep(21,5)) 

mat_dado$age <- factor(mat_dado$age, levels = c("2-y", "5-y","39-y"))

#spider plot
par(xpd = T, mar=c(7, 4, 0, 0) + 0.1)
plot(taxon1, type = "n", xlab = "PCoA1 (13%)", ylab = "PCoA2 (6.4%)", asp = 1)
ordispider(taxon1, groups= mat_dado$cond, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(taxon1, col = colvec1, bg = colvec1, pch = pointsymbol1)
legend(-1, -2.4, legend = c("2-y", "5-y", "39-y"), bty = "n",
       col = c("#fdcc8a", "#d7191c", "#2b83ba"), pch = 15, pt.bg = c("#fdae61", "#d7191c", "#2b83ba"), horiz = TRUE)
legend(-1.4,-2.7, legend = c("Control", "Landslide"), lty = c(1,2), col = c("black","darkgrey"), bty = "n",horiz = TRUE)


################################################################################################################
#two plots
png("PCoA_taxon_func.tiff", width=23, height=13, units="cm", res=400)
par(mfrow = c(1,2), oma = c(0, 0, 0, 0), mai = c(1.5, 0.8, 0.1, 0.1), xpd = NA)
plot(taxon1, type = "n", xlab = "PCoA1 (13%)", ylab = "PCoA2 (6.4%)", asp = 1)
ordispider(taxon1, groups= mat_dado$cond, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(taxon1, col = colvec1, bg = colvec1, pch = pointsymbol1)
legend(-1.25, 0.67, legend = "a", cex = 1.5, bty = "n", xpd = T)

plot(func1, type = "n", xlab = "PCoA1 (75.7%)", ylab = "PCoA2 (9.3%)")
ordispider(func1, groups= trat, col = c("black","darkgrey"), label = F, lty = c(1,2))
points(func1, col = colvec, bg = colvec, pch = pointsymbol)
legend(-0.004, 0.00125, legend = "b", cex = 1.5, bty = "n", xpd = T)

legend(-0.0045,-0.0029, legend = c("2-y", "5-y", "39-y"), bty = "n",
       col = c("#fdcc8a", "#d7191c", "#2b83ba"), pch = c(21,17,15), pt.bg = c("#fdae61", "#d7191c", "#2b83ba"), horiz = TRUE, pt.cex = 1.5)
legend(-0.005,-0.0031, legend = c("Control", "Landslide"), lty = c(1,2), col = c("black","darkgrey"), bty = "n",horiz = TRUE,
       pt.cex = 1.5)

dev.off()

################################################################################################################

#paired T test for NMDS first score (by local)
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
t.test(taxonMI$axisc, taxonMI$axisl, paired = T)
#CA
t.test(taxonCA$axisc, taxonMI$axisl, paired = T)
#SA
t.test(taxonSA$axisc, taxonMI$axisl, paired = T)
#SM
t.test(taxonSM$axisc, taxonMI$axisl, paired = T)
#QU
t.test(taxonQU$axisc, taxonMI$axisl, paired = T)


#ploting
library(ggplot2)


MIt <- ggplot(taxon1, aes(local, Axis.1, fill = cond)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "PCoA Axis1") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        legend.title=element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  #scale_fill_manual(values = c("white", "grey")) +
  #scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  geom_point(aes(Axis.1), size = 5) +
  geom_line(aes(x = rep(c(-1, 1), each = 10), group = id))

MIt

#https://stackoverflow.com/questions/31102162/r-paired-dot-plot-and-box-plot-on-same-graph-is-there-a-template-in-ggplot2


################################################################################################################
#Traits boxplots
traits_plot <- read.csv("traits_cwm.csv", header = T, sep = ";")

levels(traits_plot$cond) <- (c("control", "landslide"))
levels(traits_plot$local)[levels(traits_plot$local) == "ca"] <- "CA"
levels(traits_plot$local)[levels(traits_plot$local) == "mi"] <- "MI"
levels(traits_plot$local)[levels(traits_plot$local) == "qu"] <- "QU"
levels(traits_plot$local)[levels(traits_plot$local) == "sa"] <- "SA"
levels(traits_plot$local)[levels(traits_plot$local) == "sm"] <- "SM"

traits_plot$local <- factor(traits_plot$local, levels = c("CA", "MI", "SA", "SM", "QU"))

library(ggplot2)
library(cowplot)
library(gridExtra)

################################################################################################################
#leaf area

traits_plot$age <- c(rep("2-y",10), rep("5-y",10), rep("2-y",10), rep("39-y",10), rep("2-y",10)) 
traits_plot$age <- factor(traits_plot$age, levels = c("2-y", "5-y", "39-y"))

la <- ggplot(traits_plot, aes(local, area, fill = cond, colour = age)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = expression(paste("Leaf area (", mg^2, ")")),
                     #breaks = seq(10, 100, 25),
                     limits = c(10,100))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        legend.title=element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values = c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
la

#leaf margin
leafm <- ggplot(traits_plot, aes(local, margem_0, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf margin",
                     #breaks = seq(0.30, 1.1, 0.20),
                     limits = c(0.20,1.1))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafm

#leaf division
leafd <- ggplot(traits_plot, aes(local, divlimbo_0, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf division",
                     #breaks = seq(0.50, 1.1, 0.1),
                     limits = c(0.5,1.1))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafd

#perimeter:area ratio
leafp <- ggplot(traits_plot, aes(local, periarea, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Perimeter:area \n ratio",
                     #breaks = seq(0.04,1.4,0.2),
                     limits = c(0.035,1.5))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafp

#leaf thickness
leaft <- ggplot(traits_plot, aes(local, esp, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf thickness (mm)",
                     #breaks = seq(0.07,0.23,0.05),
                     limits = c(0.06,0.24))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leaft

#specific leaf area
sla <- ggplot(traits_plot, aes(local, sla, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = (expression(atop("Specific leaf", paste("area ")(mg%.%mm^-2)))),
                     #breaks = seq(200,520,50),
                     limits = c(190,520))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
sla

#Leaf dry matter content
leafdry <- ggplot(traits_plot, aes(local, ldmc, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf dry matter \n content (mg)",
                     #breaks = seq(190,415,50),
                     limits = c(189,420))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafdry

#Leaf pubescence - unidade de medida
leafpu <- ggplot(traits_plot, aes(local, tricoma_0, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Leaf pubescence",
                     #breaks = seq(0.08,1.1,0.2),
                     limits = c(0.07,1.1))+
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0))) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
leafpu


#getting the legend from the first plot
legenda <- get_legend(la)

grafico <- plot_grid(la + theme(legend.position="none",
                                axis.title.x=element_blank(),
                                axis.text.x=element_blank()), 
             leafm + theme(legend.position="none",
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank()),
             leafd  + theme(legend.position="none",
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
             labels = letters[1:8],
             hjust = -1.5,
             vjust = 1.2,  
             nrow = 4,
             scale = 1)
grafico

png("boxplot_traits_color.tiff", width=20, height=25, units="cm", res=400)
plot_grid(grafico, legenda, ncol = 1, rel_heights = c(1, .1), rel_widths = c(1,1),
                       label_size = 16)
dev.off()


##############################################################################################
#alpha diversities boxplots

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

#taxonomic alfa diversity (H')
taxdiv_a <- ggplot(diver, aes(local, H, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = expression(paste("Taxonomic ", alpha, " diversity (H')"))) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 0.8, y = 3.2, label = "a", size = 4) +
  annotate("text", x = 1.2, y = 2.9, label = "ac", size = 4) +
  annotate("text", x = 1.8, y = 2.2, label = "b", size = 4) +
  annotate("text", x = 2.2, y = 2.1, label = "b", size = 4) +
  annotate("text", x = 2.8, y = 2.3, label = "ab", size = 4) +
  annotate("text", x = 3.2, y = 1.95, label = "bc", size = 4) +
  annotate("text", x = 3.8, y = 2.25, label = "bc", size = 4) +
  annotate("text", x = 4.2, y = 3.2, label = "ab", size = 4) +
  annotate("text", x = 4.8, y = 2.9, label = "ab", size = 4) +
  annotate("text", x = 5.2, y = 3, label = "ab", size = 4)
taxdiv_a


#taxonomic evenness (J)
jeve <- ggplot(diver, aes(local, J, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Taxonomic eveness (J)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 0.8, y = 0.98, label = "a", size = 4) +
  annotate("text", x = 1.2, y = 1.01, label = "a", size = 4) +
  annotate("text", x = 1.8, y = 1.02, label = "a", size = 4) +
  annotate("text", x = 2.2, y = 1.02, label = "ab", size = 4) +
  annotate("text", x = 2.8, y = 1.01, label = "a", size = 4) +
  annotate("text", x = 3.2, y = 0.9, label = "b", size = 4) +
  annotate("text", x = 3.8, y = 1.02, label = "a", size = 4) +
  annotate("text", x = 4.2, y = 0.985, label = "ab", size = 4) +
  annotate("text", x = 4.8, y = 1.015, label = "a", size = 4) +
  annotate("text", x = 5.2, y = 1.01, label = "a", size = 4)
jeve


#functional alfa diversity (Raos)
fundiv_a <- ggplot(diver, aes(local, RaoQ, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = expression(paste("Functional ", alpha, " diversity (Rao's)"))) +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "")
fundiv_a

#functional eveness
feve <- ggplot(diver, aes(local, Even, fill = cond, colour = age)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_x_discrete(name = "Site") +
  scale_y_continuous(name = "Functional eveness") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.2,0,0), "cm"),
        axis.title.y = element_text(vjust = 0, margin = margin(0,5,0,0)),
        legend.text=element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  scale_fill_manual(values=c("white", "grey")) +
  scale_color_manual(values = c("#fdae6b", "#d7191c", "#2b83ba")) +
  labs(fill = "") +
  background_grid(major = "", minor = "") +
  annotate("text", x = 0.8, y = 0.82, label = "a", size = 4) +
  annotate("text", x = 1.2, y = 0.83, label = "a", size = 4) +
  annotate("text", x = 1.8, y = 0.9, label = "ab", size = 4) +
  annotate("text", x = 2.2, y = 0.98, label = "ab", size = 4) +
  annotate("text", x = 2.8, y = 0.95, label = "bc", size = 4) +
  annotate("text", x = 3.2, y = 0.84, label = "bc", size = 4) +
  annotate("text", x = 3.8, y = 0.91, label = "ab", size = 4) +
  annotate("text", x = 4.2, y = 0.7, label = "ab", size = 4) +
  annotate("text", x = 4.8, y = 0.87, label = "b", size = 4) +
  annotate("text", x = 5.2, y = 0.82, label = "b", size = 4)
feve

#getting the legend from the first plot
legenda1 <- get_legend(taxdiv_a)

grafico1 <- plot_grid(taxdiv_a + theme(legend.position="none",
                                axis.title.x=element_blank(),
                                axis.text.x=element_blank()), 
                     jeve + theme(legend.position="none",
                                  axis.title.x=element_blank(),
                                  axis.text.x=element_blank()), 
                     fundiv_a  + theme(legend.position="none"),
                     feve + theme(legend.position="none"),
                     align = 'hv',
                     labels = letters[1:4],
                     hjust = -1,
                     vjust = 1.2,  
                     nrow = 2,
                     ncol = 2,
                     scale = 1)
grafico1

png("boxplot_alpha_diver.tiff", width=20, height=20, units="cm", res=400)
plot_grid(grafico1, legenda1, ncol = 1, rel_heights = c(1, .1), rel_widths = c(1,1),
          label_size = 16)
dev.off()

#
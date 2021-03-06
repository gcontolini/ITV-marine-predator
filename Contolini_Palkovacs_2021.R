# Intraspecific variation in a marine predator shapes community diversity by altering foundation species

# Gina M. Contolini and Eric. P. Palkovacs

# Department of Ecology and Evolutionary Biology, University of California Santa Cruz, 130 McAllister Way, Santa Cruz, CA 95060, USA

# Corresponding author: Gina Contolini. Email: gina@contolini.com

# Updated 26 Nov 2021

#### Load packages ####
library(car) # companion to applied regression
library(ggplot2) # plotting
library(ggpubr) # plot aesthetics
library(piecewiseSEM) # structural equation modeling
library(vegan) # community ecology package
library(tidyr) # data manipulation
library(dplyr) # data manipulation
#### Create standard error function ####
se = function(x){
   sd(x) / sqrt(length(x))
}
#### Load data ####
setwd(file.choose()) 
## mussel data
mussel.data <- read.csv('mussel.lengths.csv')
## dogwhelk data
dogwhelk.data <- read.csv('dogwhelk.lengths.csv')
## algae cover
algae.data <- read.csv('algae.cover.csv')

## model data
model.data <- read.csv('Contolini_Palkovacs_model_data.csv')
## community data
comm.data <-  read.csv('Contolini_Palkovacs_community_data.csv')
#### Create plot theme ####
plot.theme <- theme( 
   panel.background = element_blank(), 
   panel.grid.minor.x = element_blank(),
   panel.grid.major.x = element_blank(),
   panel.grid = element_blank(),
   axis.text = element_text(size=13,color='black'), #makes axis labels bigger
   axis.title = element_text(size=15), #makes axis title bigger and bold
   axis.title.y = element_text(margin=margin(0,10,0,5)), #pushes the y axis title away from the y axis labels
   axis.text.x = element_text(margin=margin(0,20,20,0)), #pushes the x axis title away from the x axis labels and pushes the bars down closer to the x axis labels
   legend.title = element_text(size=16),
   legend.text = element_text(size=16),
   legend.background = element_blank(),
   panel.border = element_rect(colour = "black", fill=NA, size=0.5), #border around the whole plot
   plot.title = element_text(size=15),
   legend.key=element_blank()
)
#### Calculate Shannon-Wiener diversity ####
# Copy comm to modify it
comm.for.model <- comm.data
# Make row names each cage ID
rownames(comm.for.model) <- comm.for.model$cage 
# Remove column of cage IDs and trtmnts because all cells must be numeric
comm.for.model$cage <- NULL
comm.for.model$trtmnt <- NULL
# Compute Shannon-Wiener diversity
shannon <- diversity(comm.for.model, index = 'shannon') 

# Add community data to treatment and mussel data
model.data$shannon <- shannon

#### Piecewise structural equation model ####
summary(psem(
     lm(drilled.mean ~ trtmnt,                         data = model.data)
   , lm(remain.90th ~ drilled.mean + trtmnt,           data = model.data)
   , lm(shannon ~ remain.90th + drilled.mean + trtmnt, data = model.data) )
) 

#### Similarity Percentages (SIMPER) ####
## first two cols are ID cols; exclude them in the simper analysis.
comm.trtmnts <- comm.data$trtmnt
summary(simper(comm.data[,3:42], comm.trtmnts))

#### Figure 2 ####
## Figure 2a differences in mussel bed structure.
drill.mean <- model.data %>%
   group_by(trtmnt) %>%
   summarise(trtmnt.mean = mean(drilled.mean), sd = sd(drilled.mean), se = se(drilled.mean), n = n(), .groups = 'drop') # note all values are rounded.

## order factor levels for trtmnt
drill.mean$trtmnt <- factor(drill.mean$trtmnt, levels = c('Control','HOP','SOB','LOM'))
drill.mean.se <- aes(ymin = trtmnt.mean - se, ymax = trtmnt.mean + se)

fig.2a <- ggplot(drill.mean, aes(x = trtmnt, y = trtmnt.mean, color = trtmnt, shape = trtmnt)) +
   geom_errorbar(drill.mean.se, width = 0) +
   geom_point(size = 5) +
   annotate(x = 0.6, y = 43, 'text', label = 'a', size = 9) +
   scale_y_continuous(limits = c(0,45), expand = c(0,0)) +
   scale_color_manual(values = c("Control" = "grey50", "LOM" = "goldenrod3", "SOB" = 'maroon', 'HOP' = 'dodgerblue3')) +
   scale_shape_manual(values = c(15,18,17,16)) +
   labs(y = 'Mean drilled mussel length (mm)', x = 'Treatment') +
   plot.theme + theme(legend.position = 'none')

## Figure 2b. differences in remaining mussel length (90 percentile)
remain.90 <- model.data %>%
   group_by(trtmnt) %>%
   summarise(mean.remain.90th = mean(remain.90th), sd = sd(remain.90th), se = se(remain.90th), n = n(), .groups = 'drop')
   
## order factor levels for trtmnt
remain.90$trtmnt <- factor(remain.90$trtmnt, levels = c('Control','HOP','SOB','LOM'))
remain.90.se <- aes(ymin = mean.remain.90th - se, ymax = mean.remain.90th + se)
fig.2b <- ggplot(remain.90, aes(x = trtmnt, y = mean.remain.90th, color = trtmnt, shape = trtmnt)) + 
   geom_errorbar(remain.90.se, width = 0) +
   geom_point(size = 5) +
   annotate(x = 0.6, y = 60, 'text', label = 'b', size = 9) +
   scale_color_manual(values = c("Control" = "grey50", "LOM" = "goldenrod3", "SOB" = 'maroon', 'HOP' = 'dodgerblue3')) +
   scale_shape_manual(values = c(15,18,17,16)) +
   labs(y = '90th percentile\nremaining mussel length (mm)', x = 'Treatment') +
   plot.theme + theme(legend.position = 'none', axis.title.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)))
   
## Figure 2c. NMDS plot
comm.mds1 <- metaMDS(subset(comm.data, select=-c(cage, trtmnt)), distance = 'bray')
### extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(comm.mds1))
### add metadata
data.scores$trtmnt <- model.data$trtmnt
### create hulls
trtmnt.Control <- data.scores[data.scores$trtmnt == "Control", ][chull(data.scores[data.scores$trtmnt == 
                                                                            "Control", c("NMDS1", "NMDS2")]), ]  # hull values for trtmnt Control
trtmnt.LOM <- data.scores[data.scores$trtmnt == "LOM", ][chull(data.scores[data.scores$trtmnt == 
                                                                              "LOM", c("NMDS1", "NMDS2")]), ]  # hull values for trtmnt LOM
trtmnt.SOB <- data.scores[data.scores$trtmnt == "SOB", ][chull(data.scores[data.scores$trtmnt == 
                                                                              "SOB", c("NMDS1", "NMDS2")]), ]  # hull values for trtmnt SOB
trtmnt.HOP <- data.scores[data.scores$trtmnt == "HOP", ][chull(data.scores[data.scores$trtmnt == 
                                                                              "HOP", c("NMDS1", "NMDS2")]), ]  # hull values for trtmnt HOP
### combine hulls
hull.data <- rbind(trtmnt.Control, trtmnt.LOM, trtmnt.SOB, trtmnt.HOP)  
### create mds plot
fig.2c <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color = trtmnt, shape = trtmnt)) + 
   geom_polygon(data = hull.data, fill = NA, aes(x = NMDS1, y = NMDS2, group = trtmnt), alpha = 0.30) + # add the convex hulls
   geom_point(size = 5) + 
   scale_color_manual(name = "Treatment",
                      labels = c('Control','LOM','SOB','HOP'),
                      values = c("Control" = "grey50", "LOM" = "goldenrod3", "SOB" = 'maroon', 'HOP' = 'dodgerblue3')) +
   scale_shape_manual(name = "Treatment",
                      labels = c('Control','LOM','SOB','HOP'),
                      values = c('Control' = 15, 'LOM' = 16, 'SOB' = 17, 'HOP' = 18)) +
   annotate(x = -0.9, y = 0.9, 'text', label = 'c', size = 9) +
   plot.theme + theme(legend.position = 'none', axis.title.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)))

## Arrange panels
Figure.2 <- ggarrange(fig.2a, fig.2b, fig.2c, nrow = 1, ncol = 3) 

## Save the high resolution figure to your working directory
ggsave('Figure.2.jpeg', device = 'jpeg', width = 12, height = 4, dpi = 300) 

#### Starting dogwhelk length ANOVA ####
data.noctrl <- model.data[-which(model.data$trtmnt == 'Control'),] # remove control group
par(mfrow=c(2,2))
plot(lm(data.noctrl$mean.dogwhelk.len ~ data.noctrl$trtmnt)) # check for normality and homoscedasticity
Anova(aov(mean.dogwhelk.len ~ trtmnt, data = data.noctrl), type = 2) # no differences in starting length among populations
#### Number of mussels drilled among treatments ANOVA ####
shapiro.test(model.data$N.drilled) # normal
Anova(aov(N.drilled ~ trtmnt, data = model.data)) # differences in the number of mussels drilled
TukeyHSD(aov(N.drilled ~ trtmnt, data = model.data))

#### Drilled mussel length ANOVA ####
Anova(aov(drilled.mean ~ trtmnt, data = model.data)) # differences in drilled mussel length
TukeyHSD(aov(drilled.mean ~ trtmnt, data = model.data))

#### Linear relationship between mean drilled mussel length and 90th remaining mussel length ####
summary(lm(remain.90th ~ drilled.mean, data = model.data))
ggplot(model.data, aes(x = drilled.mean, y = remain.90th)) +
   geom_smooth(method = 'lm', color = 'black') +
   geom_point() +
   labs(x = 'Mean drilled mussel length (mm)', y = '90th percentile remaining\nmussel length (mm)') +
   plot.theme
   
#### Linear relationship between remaining mussel length and Shannon diversity ###
summary(lm(remain.90th ~ shannon, data = model.data))
ggplot(model.data, aes(x = remain.90th, y = shannon)) +
   geom_smooth(method = 'lm', color = 'black') +
   geom_point() +
   labs(x = '90th percentile remaining mussel length (mm)', y = 'Shannon-Wiener diversity') +
   plot.theme

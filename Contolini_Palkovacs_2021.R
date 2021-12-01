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
# setwd(file.choose()) 

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

#### Read in raw data ####
# community data
comm_data <-  read.csv('Contolini_Palkovacs_community_data.csv')
# mussel data
mussel_data <- read.csv('Contolini_Palkovacs_mussel_data.csv')
# dogwhelk data
dogwhelk_data <- read.csv('Contolini_Palkovacs_dogwhelk_data.csv')
# algae cover
algae_data <- read.csv('Contolini_Palkovacs_algae_data.csv')

# model data # This will be deleted bc I will create it using the raw data.
model_data <- read.csv('Contolini_Palkovacs_model_data.csv') ###!!! mean dogwhelk length here is actually mean mussel length of ALL mussels combined! Good thing I'm not using that!

#### Summarize mussel data ####
mussel_data <- mussel_data  %>% 
   filter(!trtmnt == 'NCC') %>% # remove this trtmnt bc it is not used in this analysis
   mutate(trtmnt = ifelse(trtmnt == 'CC', 'Control', trtmnt)) # rename CC to Control

# Drilled mussels
drilled_mussel_sum <- mussel_data %>% 
   filter(drilled == TRUE) %>% # only drilled mussels (loose or not)
   filter(!is.na(mussel.length)) %>% #remove na values 
   group_by(plot, trtmnt) %>%
   summarise(drilled.mean = mean(mussel.length)
           , drilled.n    = n()
           , drilled.min  = min(mussel.length)
           , drilled.max  = max(mussel.length)
           , drilled.sd   = sd(mussel.length)
           , drilled.se   = se(mussel.length)
           , drilled.90th = quantile(mussel.length, prob = .9)
         , .groups = 'drop')

# Remaining mussels
remain_mussel_sum <- mussel_data %>% 
   filter(drilled == FALSE & loose == FALSE) %>% # only remaining mussels not drilled
   filter(!is.na(mussel.length)) %>% #remove na values 
   group_by(plot, trtmnt) %>%
   summarise(  remain.mean = mean(mussel.length)
             , remain.n    = n()
             , remain.min  = min(mussel.length)
             , remain.max  = max(mussel.length)
             , remain.sd   = sd(mussel.length)
             , remain.se   = se(mussel.length)
             , remain.90th = quantile(mussel.length, .9)
             , .groups = 'drop')
# Loose mussels
loose_mussel_sum <- mussel_data %>%
   filter(loose == TRUE) %>% # only loose mussels. could be drilled
   filter(!is.na(mussel.length)) %>% # remove na values
   group_by(plot, trtmnt) %>%
   summarise(loose.mean = mean(mussel.length)
           , loose.n    = n()
           , loose.sd   = sd(mussel.length)
           , loose.se   = sd(mussel.length)
           , .groups = 'drop')
          
#### Summarize dogwhelk data ####
dogwhelk_sum <- dogwhelk_data %>%
   filter(!notes == 'shell only') %>% # remove this value because it wasn't a snail
   group_by(cage, site) %>%
   rename(plot = cage, trtmnt = site) %>% # rename cols to match other dataframes
   summarise(dogwhelk.start.mean = mean(start.length, na.rm = T)
           , dogwhelk.end.mean = mean(end.length, na.rm = T)
           , .groups = 'drop')

#### Summarize algae cover data ####
algae_data <- algae_data %>%
   rename(algae.perc = Plot.Total...Algae, plot = Plot.ID) %>% # rename percent algae col
   subset(select = c(plot, Week, algae.perc)) %>% # select only relevant cols
   filter(Week == 40) %>% # only the last week; the final algal coverage
   filter(!is.na(algae.perc) & !plot %in% c('NC01','NC02','NC03','NC04','NC05','NC06')) # remove na values and unused NC treatment
algae_data$Week = NULL # remove week col

#### Calculate Shannon-Wiener diversity ####
# Copy comm to modify it
comm_data_model <- comm_data 
# Make row names each cage ID
rownames(comm_data_model) <- comm_data_model$cage 
# Remove column of cage IDs and trtmnts because all cells must be numeric
comm_data_model$cage <- NULL
comm_data_model$trtmnt <- NULL
# Compute Shannon-Wiener diversity
shannon <- diversity(comm_data_model, index = 'shannon') 

#### Build model dataset using above summarized data ####
model_data_2 <- drilled_mussel_sum %>% 
   full_join(remain_mussel_sum) %>%
   full_join(loose_mussel_sum) %>%
   full_join(dogwhelk_sum) %>%
   full_join(algae_data)
as.data.frame(model_data_2)
model_data$shannon <- shannon

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

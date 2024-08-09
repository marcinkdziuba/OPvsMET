library(readxl)
library(ggplot2)
library(survival)
library(survminer)
library(coxme)
library(Rmisc)
library(emmeans)
library(dplyr)
library(ggpubr)
library(car)
library(ARTool)
library(tidyr)
library(forecast)
library(zoo)
library(gridExtra)
###################### GUT PENETRATION ########################################
gut_penetration_mBio <- read_excel("~/gut penetration_mBio.xlsx") #load data
gut_penetration_mBio$spores<-as.numeric(gut_penetration_mBio$spores) # change the "spores" data into values

gut_penetration_mBio$sqrtspores<-sqrt(gut_penetration_mBio$spores) # creates new variable - square rooted spore number
sqrtM1<-lm(sqrtspores~parasite*food, data = gut_penetration_mBio) # creates linear model on square rooted values
anova(sqrtM1) # prints stats

op <- par(mfrow = c(2, 2)) # following three lines create diagnostic plots
plot(sqrtM1)             
par(op)

#posthoc
emmeans(sqrtM1, pairwise ~ parasite * food)

gptemp<-na.omit(gut_penetration_mBio) # removes all NAs from the database (necessary for summarySE)
gptempSE<-summarySE(data=gptemp,measurevar = 'sqrtspores', groupvars = c('food','parasite'))

# command below makes violin plot in ggplot

SVP<-ggplot(gut_penetration_mBio, aes(x=food, y=sqrtspores, fill=interaction(parasite,food))) + 
  geom_violin(position=position_dodge(.8))+
  scale_alpha_discrete(range = c(1, 0.5))+
  #guides(color=guide_legend("Parasite"), alpha = "none")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, position=position_dodge(.8), color = "white")+
  scale_fill_manual(values = c('black','darkorange1','grey60','darkgoldenrod1'), labels=c("MB High","OP+MB High","MB Low","OP+MB Low"), name = "Pathogen")+
  ylab("SQRT Attacking spores")+
  xlab("Food")+
  #scale_x_discrete(alpha=c(0.9, 0.2))+
  theme_classic()+
  scale_y_continuous(limits= c(0,5.05),expand = c(0,0))+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.position = c(0.85,0.85),
        legend.background = element_blank())
SVP

ggsave("SVP.png",SVP, width = 8, height = 8) # exports the plot to .png
ggsave("SVP.eps",SVP, width = 8, height = 8)

################################################################################

################# MORTALITY ####################################################

mortality_mBio <- read_excel("~/mortality mBio.xlsx")
mortality_mBio$lifespan<-as.numeric(mortality_mBio$lifespan)

mortality_mBio<-mortality_mBio[mortality_mBio$lifespan>19,] # removing individuals that didn't
# survive to the M. bicuspidata exposure
mortality_mBio<-na.omit(mortality_mBio)      # removing the rows with NA
mortality_mBio<-mortality_mBio[!(mortality_mBio$food=="H" & mortality_mBio$parasite=="MG_LATE_MET" & mortality_mBio$rep==16),]
mortality_mBio<-mortality_mBio[!(mortality_mBio$food=="H" & mortality_mBio$parasite=="MG_LATE_MET" & mortality_mBio$rep==29),]
#two lines above remove the two animals who were exposed to M. bicuspidata but did not die or get infected

coxM<- coxph(Surv(lifespan, status) ~ parasite+food, data = mortality_mBio) # setting up cox model
coxM
Anova(coxM) # prints stats

cox.zph(coxM) # checks for test assumptions


mortality_mBio_HF<-mortality_mBio[mortality_mBio$food=="H",] # makes new database (subests only high food)
mortality_mBio_LF<-mortality_mBio[mortality_mBio$food=="L",] # makes new database (subsets low food)

LM1<- # produces Kaplan-Meier plot for high food treatments
  ggsurvplot(survfit(Surv(lifespan,   status) ~ parasite, data = mortality_mBio_HF),
             conf.int=F, palette = c("black","darkorange1"),xlim = c(0, 42))
LM1$plot<-LM1$plot+
  ggplot2::annotate("text", 
                    x = 10, y = 0.2, # x and y coordinates of the text
                    label = "",size = 10)+
  scale_color_manual(values = c("black","darkorange1"), labels=c("MB","OP+MB"))+
  xlab("Day")+
  ggtitle("High food")+
  scale_y_continuous(limits=c(0,1.05),expand = c(0,0))+
  scale_x_continuous(limits = c(0,40), expand = c(0,0))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = c(0.2,0.2),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(color = "black",size = 20),
        axis.title.x = element_text(color = "black", size = 20),
        axis.title.y = element_text(colour = "black", size = 20))

LM1

coxHF<- coxph(Surv(lifespan, status) ~ parasite, data = mortality_mBio_HF) # cox model for high food
summary(coxHF)
Anova(coxHF)

cox.zph(coxHF)

LM2<- # K-M plot for low food 
  ggsurvplot(survfit(Surv(lifespan,   status) ~ parasite, data = mortality_mBio_LF),
             conf.int=F,title = "Low food", palette = c('grey60','darkgoldenrod1'),xlim = c(0, 42))
LM2$plot<-LM2$plot+
  ggplot2::annotate("text", 
                    x = 10, y = 0.2, # x and y coordinates of the text
                    label = "", size = 10)+
  scale_color_manual(values = c('grey60','darkgoldenrod1'), labels=c("MB","OP+MB"))+
  xlab("Day")+
  ggtitle("Low food")+
  scale_y_continuous(limits=c(0,1.05),expand = c(0,0))+
  scale_x_continuous(limits = c(0,40), expand = c(0,0))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = c(0.2,0.2),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(color = "black",size = 20),
        axis.title.x = element_text(color = "black", size = 20),
        axis.title.y = element_blank())

LM2

coxLF<- coxph(Surv(lifespan, status) ~ parasite, data = mortality_mBio_LF) # cox model for low food
summary(coxLF)
Anova(coxLF)

cox.zph(coxLF)


################################################################################

############################ REPRODUCTION ######################################

reproduction_mBio<- read_excel("~/reproduction mBio.xlsx")
reproduction_mBio$age<-as.numeric(reproduction_mBio$age) # turning age and NB into numeric
reproduction_mBio$NB<-as.numeric(reproduction_mBio$NB)
reproduction_mBio<-na.omit(reproduction_mBio) # removes units with NA
reproduction_mBio<-reproduction_mBio[!(reproduction_mBio$food=="H" & reproduction_mBio$parasite=="MG_LATE_MET" & reproduction_mBio$rep==29),]


####
totNB_mBio <- reproduction_mBio %>% group_by(food,parasite,rep) %>% # next four lines sum all the offspring
  summarise(tot_NB = sum(NB),                                       # and gives total lifetime reproduction
            .groups = 'drop')                                       # as tot_NB variable in new dataframe
totNB_mBio <- totNB_mBio %>% as.data.frame()                        # named totNB_mBio

totNB_mBio<-totNB_mBio[-c(2,5,8,10,11,19,62,66,69,73,94,104,105,106),]  # this removes all ind. who lived shorter than 20 days

totNB_mBio<-totNB_mBio[!(totNB_mBio$food=="H" & totNB_mBio$parasite=="MG_LATE_MET" & totNB_mBio$rep==16),]
totNB_mBio<-totNB_mBio[!(totNB_mBio$food=="H" & totNB_mBio$parasite=="MG_LATE_MET" & totNB_mBio$rep==29),]
# the two lines above remove animals who did not get infected or killed by M. bicuspidata


repM1<-lm(tot_NB~parasite*food, data = totNB_mBio) # linear model
anova(repM1)                                       # prints stats

op <- par(mfrow = c(2, 2))                         # diagnostic plots
plot(repM1)             
par(op)

totNB_mBioSE<-summarySE(data = totNB_mBio, measurevar = "tot_NB", groupvars = c("food","parasite"), na.rm=T)
# the line above summarizes data within treatments and calculates sd, se and ci

Rep<- # produces the plot
  ggplot(totNB_mBioSE, aes(x=food, y=tot_NB, fill=interaction(parasite,food), shape = interaction(parasite,food), 
                           #                           alpha = food
  )) + 
  geom_errorbar(aes(ymin=tot_NB-se, ymax=tot_NB+se), colour="black", width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2), size=6)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0))+
  ylab("Lifetime offspring production")+
  xlab("Food")+
  ggtitle(" ")+
  scale_x_discrete(labels=c("High","Low"))+
  scale_fill_manual(values = c('black','darkorange1','grey60','darkgoldenrod1'), name = "Pathogen",
                    labels=c("MB High","OP+MB High","MB Low","OP+MB Low")
  )+
  scale_shape_manual(values = c(21,22,21,22 ), name = "Pathogen", 
                     labels=c("MB High","OP+MB High","MB Low","OP+MB Low")
  )+
  
  #annotate("text", x = 1.75, y = 120, label = "Parasite F = 47.775 *** \n
  #        Food F = 488.206 *** \n
  #       Parasite x Food F = 0.153, p-value = 0.696")+
  theme_classic()+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.position = c(0.75,0.85))
Rep


mplot1<-LM1$plot # extracts the high food mortality plot into mplot1 object
mplot2<-LM2$plot # extracts the low food mortality plot into mplot2 object

Fig2<-grid.arrange(mplot1,mplot2, Rep, ncol = 3) # combines three plots 
ggsave("Fig2.png",Fig2, width = 16, height = 6 ) # exports the plot
ggsave("Fig2.eps",Fig2, width = 16, height = 6 ) # exports the plot

################################################################################

####################### PREVALENCE AND BURDEN ##################################

infection_mBio <- read_excel("~/infection mBio.xlsx") # loads data
infection_mBio$Mature<-as.numeric(infection_mBio$Mature) # changes n. of mature spores into numeric
infection_mBio$MetschInfectious<-as.numeric(infection_mBio$MetschInfectious) # changes effective infection score into numeric
infection_mBio$MetschAll<-as.numeric(infection_mBio$MetschAll) # changes infection score into numeric
infection_mBio$food<-as.factor(infection_mBio$food) # changes food into factor
infection_mBio$parasite<-as.factor(infection_mBio$parasite) # changes parasite into factor
infection_mBio<-na.omit(infection_mBio) # removes all NA cases

### INFECTION PROBABILITY ###

prevM1<-glm(MetschAll~parasite*food,family = binomial(link = "logit"),data = infection_mBio) # global binomial glm model with interaction
Anova(prevM1, type = 3) # prints stats

infection_HF<-infection_mBio[!(infection_mBio$food=="L"),] # makes new dataset with only high food (MB vs. OP+MB)
prevHFM1<-glm(MetschAll~parasite,family = binomial(link = "logit"),data = infection_HF) # binomial glm on MB high food vs. OP+MB high food
Anova(prevHFM1, type = 2) # stats

infection_MET<-infection_mBio[!(infection_mBio$parasite=="MG_LATE_MET"),] # makes new dataset with only MB treatment (in high and low food)
prevMetM1<-glm(MetschAll~food,family = binomial(link = "logit"),data = infection_MET) # binomial glm on MB high food vs. MB low food
Anova(prevMetM1, type = 2) # stats


prevSE<-summarySE(data = infection_mBio, measurevar = "MetschAll", groupvars = c("parasite", "food"), na.rm = T)
# the line above summarizes data within treatments and calculates sd, se and ci

S1<- # makes the plot
  ggplot(prevSE, aes(x=food, y=MetschAll, fill=interaction(parasite,food), shape = interaction(parasite,food),
                     #alpha = food
  )) + 
  geom_hline(yintercept = 0, linetype = 2)+
  geom_line(position=position_dodge(.2))+
  scale_alpha_discrete(range = c(1, 0.5))+
  #guides(fill=guide_legend("Parasite"), alpha = "none")+
  #expand_limits(y=c(0,4))+
  geom_errorbar(aes(ymin=MetschAll-se, ymax=MetschAll+se), width=.2, 
                position=position_dodge(.2)) + 
  scale_y_continuous(limits = c(0,1.05))+
  ylab("M. bicuspidata \n infection probability")+
  xlab("Food")+
  scale_x_discrete(labels=c("High","Low"))+
  geom_point(position=position_dodge(.2), size=6)+
  scale_fill_manual(values = c('black','darkorange1','grey60','darkgoldenrod1'), name = "Food", labels = c("MB High","OP+MB High","MB Low","OP+MB Low"))+
  scale_shape_manual(values = c(21,22,21,22), name = "Parasite", labels = c("MB High","OP+MB High","MB Low","OP+MB Low"))+
  theme_classic()+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.position = "none")

S1 # prints the plot


### EFFECTIVE INFECTION PROBABILITY ###

eff_prevM1<-glm(MetschInfectious~parasite*food,family = binomial(link = "logit"),data = infection_mBio ) # global binomial glm model with interaction
Anova(eff_prevM1, type = 3) # prints stats


eff_prevMetM1<-glm(MetschInfectious~food,family = binomial(link = "logit"),data = infection_MET) # binomial glm on MB high food vs. MB low food
Anova(eff_prevMetM1, type = 2) # stats

eff_prevHFM1<-glm(MetschInfectious~parasite,family = binomial(link = "logit"),data = infection_HF) # binomial glm on MB high food vs. OP+MB high food
Anova(eff_prevHFM1, type = 2) # stats


eff_prevSE<-summarySE(data = infection_mBio, measurevar = "MetschInfectious", groupvars = c("parasite", "food"), na.rm = T)
# the line above summarizes data within treatments and calculates sd, se and ci

S2<- # makes the plot
  ggplot(eff_prevSE, aes(x=food, y=MetschInfectious, fill=interaction(parasite,food), shape = interaction(parasite,food), 
                         #                        alpha = food
  )) + 
  geom_hline(yintercept = 0, linetype = 2)+
  geom_line(position=position_dodge(.2))+
  scale_alpha_discrete(range = c(1, 0.5))+
  #guides(fill=guide_legend("Parasite"), alpha = "none")+
  #expand_limits(y=c(0,4))+
  geom_errorbar(aes(ymin=MetschInfectious-se, ymax=MetschInfectious+se), width=.2, 
                position=position_dodge(.2)) + 
  scale_y_continuous(limits = c(0,1.05))+
  ylab("M. bicuspidata effective \n infection probability")+
  xlab("Food")+
  geom_point(position=position_dodge(.2), size=6)+
  scale_x_discrete(labels=c("High","Low"))+
  scale_fill_manual(values = c('black','darkorange1','grey60','darkgoldenrod1'), name = "Food", labels = c("MB High","OP+MB High","MB Low","OP+MB Low"))+
  scale_shape_manual(values = c(21,22,21,22), name = "Parasite", labels = c("MB High","OP+MB High","MB Low","OP+MB Low"))+
  theme_classic()+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.position = "none")

S2



### MATURE SPORES BURDEN ###

InfNoZero<-infection_mBio[infection_mBio$Mature > 0,] # creates a new dataframe with just effectively infected individuals
InfNoZero<-na.omit(InfNoZero) # removes cases with NA
InfNoZero$food<-as.factor(InfNoZero$food) # sets food into factor
InfNoZero$parasite<-as.factor(InfNoZero$parasite) # sets parasite into factor

art_bM1<-art(Mature~parasite*food,data=InfNoZero) # creates artAnova model (non-parametric anova) on number of mature spores
anova(art_bM1, type = 2) # prints stats


ILNOzeroSE<-summarySE(data = InfNoZero, measurevar = "Mature", groupvars = c("parasite","food"), na.rm=T)
# the line above summarizes data within treatments and calculates sd, se and ci


new_row = c(parasite = "MG_LATE_MET", food = "L", N = 1, Mature = 0, sd = NA, se = NA, ci = NA)
# the line above makes extra row for low food OP+MB that will indicate 0 spores on the plot
# it needs to be added because all the 0 spores data were filtered out before
ILNOzeroSE<-rbind(ILNOzeroSE, new_row) # attaches the new row to the dataframe

ILNOzeroSE$Mature<-as.numeric(ILNOzeroSE$Mature) # sets Mature as numeric
ILNOzeroSE$se<-as.numeric(ILNOzeroSE$se)         # sets se as numeric

S3<- # makes the plot
  ggplot(ILNOzeroSE, aes(x=food, y=Mature, fill=interaction(parasite,food), shape = interaction(parasite,food), 
                         #                         alpha = food
  )) + 
  geom_hline(yintercept = 0, linetype = 2)+
  geom_line(position=position_dodge(.2))+
  #  scale_alpha_discrete(range = c(1, 0.5))+
  guides(fill=guide_legend("Pathogen"), alpha = "none")+
  geom_errorbar(aes(ymin=Mature-se, ymax=Mature+se), width=.2, 
                position=position_dodge(.2)) + 
  ylab("M. bicuspidata \n mature spore burden x 1000")+
  xlab("Food")+
  scale_y_continuous(limits = c(0,180))+
  geom_point(position=position_dodge(.2), size=6)+
  scale_x_discrete(labels=c("High","Low"))+
  scale_fill_manual(values = c('black','darkorange1','grey60','darkgoldenrod1'), name = "Food", labels = c("MB High","OP+MB High","MB Low","OP+MB Low"))+
  scale_shape_manual(values = c(21,22,21,22), name = "Pathogen", labels = c("MB High","OP+MB High","MB Low","OP+MB Low"))+
  theme_classic()+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.position = c(0.7,0.85))

S3

sporeplot<-ggarrange(S1, S2, S3, ncol = 3, nrow = 1) # arranges the three plots into one

ggsave("sporeplot.png", sporeplot, height = 6, width = 14) # saves the plot
ggsave("sporeplot.eps", sporeplot, height = 6, width = 14)
################################################################################

######################### FIELD PREVALENCE #####################################

Field<-read.csv("~/Field.csv",header = T, sep = ",") # loads data

Field_long <- pivot_longer(Field, cols = c(micG.prev, metsch.prev), names_to = "parasite", values_to = "prevalence")
# the line above reshuffles the data, making a separate record for each parasite x date x lake combination
# and exports it into a new dataframe

micGandmetschprevplot<-ggplot(Field_long,aes(x=Julian, y=prevalence, group=parasite, color=parasite)) + # makes the plot
  geom_hline(yintercept=0, linetype = "dashed", color = "grey")+
  geom_line() +
  labs(color="parasite", x="Day of Year", y="Prevalence") +
  scale_color_manual(values = c("black","darkorange1"), labels = c("MB","OP"))+
  facet_grid(Year ~ Lake) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.15,0.9),
        axis.text = element_text(colour = "black"))

micGandmetschprevplot

ggsave("Filedplot.png", micGandmetschprevplot, height = 6, width = 6) # exports the plot
ggsave("Filedplot.eps", micGandmetschprevplot, height = 6, width = 6) 

W21<-Field[Field$Lake=="Walsh" & Field$Year==2021,]     # subsets the "Field" dataframe for one lake x one year
ccf(W21$metsch.prev,W21$micG.prev, main = "Walsh 2021") # plots correlation coefficients for consecutive lags
print(ccf(W21$metsch.prev,W21$micG.prev))               # prints the correlation coefficients

W22<-Field[Field$Lake=="Walsh" & Field$Year==2022,]
ccf(W22$metsch.prev,W22$micG.prev)
print(ccf(W22$metsch.prev,W22$micG.prev))

S21<-Field[Field$Lake=="Sullivan" & Field$Year==2021,]
ccf(S21$metsch.prev,S21$micG.prev)
print(ccf(S21$metsch.prev,S21$micG.prev))

S22<-Field[Field$Lake=="Sullivan" & Field$Year==2022,]
ccf(S22$metsch.prev,S22$micG.prev)
print(ccf(S22$metsch.prev,S22$micG.prev))




###############################################################################

################# FOOD PILOT ##################################################

Food_pilot <- read_excel("~/Food pilot.xlsx") # loads data

Food_pilot$food<-as.character(Food_pilot$food)  # changes food into "character"
Food_pilot$food <- factor(Food_pilot$food, levels=c("100% (High)", "0.5", "0.25","0.125","6.25% (Low)")) 
# the line above rearranges the order of food levels

food<-summarySE(data = Food_pilot, measurevar = "NB", groupvars = c("food"), na.rm = T)
# summarizes data, calculates sd, se and ci

Pilot<-
  ggplot(food, aes(x=food, y=NB)) +  # makes the plot
  geom_errorbar(aes(ymin=NB-se, ymax=NB+se), colour="black", width=.2, position=position_dodge(0.2)) +
  #geom_line(position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.2), size=4)+
  scale_y_continuous(limits = c(0,50), expand = c(0,0))+
  ylab("# offspring till clutch 4")+
  xlab("Food level")+
  scale_x_discrete(limit=rev(levels(food$food)),labels=c( "6.25% (Low)","12.5%","25%", "50%","100% (High)"))+
  ggtitle(" ")+
  theme_classic()+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.position = c(0.85,0.85))

Pilot

ggsave("Pilot.png", Pilot, height = 6, width = 8) # exports the plot

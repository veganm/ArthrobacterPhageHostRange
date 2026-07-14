# Script for host range data analysis
# Isolation host A. globiformis B-2979
# Emory University phage 2023-2024

pacman::p_load(tidyverse, ggpubr, patchwork, cowplot, GGally, readxl)
xTextSize<-14

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in cleaned data
HostRange24.all<-read.csv("GlobiHostRangeData_24_25.csv", header=TRUE)
#view(HostRange24.all)

#unique(HostRange24.all$Bacteria)

# heatmap EOPs
HostRange24.EOPhosts<-HostRange24.all |>
  dplyr::filter(Bacteria=="B2880" | Bacteria=="B24025" | Bacteria=="B24478" | Bacteria=="B24479") |>
  mutate(logEOP=log10(EOP)) 

# reorder for plotting
HostRange24.EOPhosts$Bacteria <-factor(HostRange24.EOPhosts$Bacteria, 
                                       levels=c('B2880', 'B24025', 'B24478', 'B24479'))
  
# Plot heatmap of log10(EOP). NA (no plaque formation recorded) is in grey
pEOPbyPhage.filtered.hm<-HostRange24.EOPhosts %>%
  ggplot(aes(y=Name, x=Bacteria, fill=logEOP)) + 
  geom_tile() +
  scale_fill_gradient2(low="#edf8b1", mid="#7fcdbb", high="navy", midpoint=-3, na.value="grey50") +
  #scale_fill_distiller(palette = "RdPu", direction=-1, na.value="grey50") +
  theme_classic()+
  theme(text=element_text(size=xTextSize),
        axis.text.x = element_text(angle=90, hjust=1),
        plot.title = element_text(hjust=0.5))+
  labs(y="Phage", x="Bacteria", title="Plating Efficiency")+
  facet_grid(cols=vars(Set))
pEOPbyPhage.filtered.hm


# filter only the entries with EOP values.
# Note that this removes Ultraviolet,
# which did not have EOP on any non-host isolate
HostRange24.EOP.all<-HostRange24.all[complete.cases(HostRange24.all),]
HostRange24.EOP.all$logEOP<-log10(HostRange24.EOP.all$EOP)
#view(HostRange24.EOP.all)
unique(HostRange24.EOP.all$PhageID)

# count number of non-target hosts 
# (primary data set only)
HostRangeCountByPhage<-HostRange24.EOP.all %>%
  group_by(PhageID) %>%
  dplyr::filter(LFW!="Y" & Set == 1) %>%
  summarise(nNonTarget=n())
view(HostRangeCountByPhage)
length(which(HostRangeCountByPhage$nNonTarget==1))
length(which(HostRangeCountByPhage$nNonTarget==2))
length(which(HostRangeCountByPhage$nNonTarget==3))
length(which(HostRangeCountByPhage$nNonTarget==4))

# plot raw logEOP by phage (another way of seeing these data)
HostRange24.EOP.all %>%
  ggplot(aes(x=logEOP, y=Name, pch=LFW))+
  geom_point()+
  theme_bw()+
  geom_vline(xintercept=0, color="red")+
  facet_wrap(~Bacteria)

cols<-c("N"="black", "Y"="blue")
# and without faceting
pEOPbyPhage.all<-HostRange24.EOP.all %>%
  ggplot(aes(y=logEOP, x=Name))+
  geom_boxplot()+
  #  geom_point(aes(pch=LFW, color=LFW), size=2)+
  geom_jitter(aes(pch=Bacteria, color=LFW), width=0.1, size=2)+
  theme_bw()+
  scale_color_manual(values=cols)+
  scale_shape_manual(values=c(17, 16, 15, 2))+
  geom_hline(yintercept=0, color="red")+
  theme(text=element_text(size=xTextSize),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(hjust=0.5))+
  labs(y=expression(log[10](EOP)), x="", title="EOP by phage")+
  facet_wrap(~Set, ncol=1)
pEOPbyPhage.all

# and by host across all phages
pEOPbyHost.all<-HostRange24.EOP.all %>%
  #dplyr::filter(Set == 1) %>%
  ggplot(aes(y=logEOP, x=Bacteria))+
  geom_boxplot()+
  #geom_violin()+
  geom_jitter(aes(pch=Bacteria, color=LFW), width=0.05, size=2)+
#  geom_jitter(aes(color=LFW), width=0.05, size=2)+
  theme_bw()+
  scale_color_manual(values=cols)+
  scale_shape_manual(values=c(17, 16, 15, 2))+
  theme(text=element_text(size=xTextSize),
        #axis.text.x = element_text(angle=90),
        plot.title = element_text(hjust=0.5))+
  geom_hline(yintercept=0, color="red")+
  labs(y=expression(log[10](EOP)), pch="Bacteria", title="EOP by host", x="") +
  #facet_wrap(~HostSpecies, scales="free_x", nrow = 1)+
  facet_grid(cols=vars(Set))
pEOPbyHost.all


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How do EOPs compare on the two non-isolation-host globi strains?
HostRange24.EOP.globi<-  HostRange24.EOP.all %>%
  dplyr::filter(LFW!="Y") %>%
  dplyr::select(Set, PhageID, Name, Bacteria, logEOP) %>%
  dplyr::filter(Bacteria=="B24025" | Bacteria=="B2880" | Bacteria=="B24478") %>%
  pivot_wider(names_from=Bacteria, values_from=logEOP)


### stat tests

# Wilcox
wilcox.test(HostRange24.EOP.globi$B2880, HostRange24.EOP.globi$B24478) #p=0.3576
wilcox.test(HostRange24.EOP.globi$B2880, HostRange24.EOP.globi$B24025) #p=0.002453 
wilcox.test(HostRange24.EOP.globi$B24478, HostRange24.EOP.globi$B24025) #p=0.02794

# Linear regressions
HostRange24.EOP.globi.filtered<-HostRange24.EOP.globi |>
  dplyr::filter(B2880 >=-4 & B24478>=-4)
lm.2880.24478<-lm(B2880~B24478, data=HostRange24.EOP.globi.filtered)
summary(lm.2880.24478)
# n.s.

HostRange24.EOP.globi.filtered<-HostRange24.EOP.globi |>
  dplyr::filter(B2880 >=-4 & B24025>=-4)
lm.2880.24025<-lm(B2880~B24025, data=HostRange24.EOP.globi.filtered)
summary(lm.2880.24025)
#m = 0.7004, b = -0.9731, p=0.000143 on F test
#R2adj = 0.7572

HostRange24.EOP.globi.filtered<-HostRange24.EOP.globi |>
  dplyr::filter(B24478 >=-4 & B24025>=-4)
lm.24478.24025<-lm(B24478~B24025, data=HostRange24.EOP.globi.filtered)
summary(lm.24478.24025)
# n.s.

# replace NAs with -8 for plotting
HostRange24.EOP.globi<-HostRange24.EOP.globi %>%
  dplyr::mutate(B24025= replace_na(B24025,-8)) %>%
  dplyr::mutate(B2880= replace_na(B2880,-8)) %>%
  dplyr::mutate(B24478= replace_na(B24478,-8))

#view(HostRange24.EOP.globi)

# medians?
median(HostRange24.EOP.globi$B2880[HostRange24.EOP.globi$B2880 > -8]) #-1.3098
median(HostRange24.EOP.globi$B24025[HostRange24.EOP.globi$B24025>-8]) #-0.32

#plot globi v globi
pHostRange24.EOP.globi<-HostRange24.EOP.globi %>%
  ggplot(aes(x=B24025, y=B2880))+
  geom_point(size=2)+
  xlim(-8, 1)+
  ylim(-8, 1)+
  #scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, color="grey20")+
  geom_vline(xintercept = 0, color="grey20")+
  geom_abline(slope=1, intercept=0, color="grey80")+
  geom_abline(slope=lm.2880.24025$coefficients[[2]], intercept = lm.2880.24025$coefficients[[1]], linetype="dashed", color="red")+
  theme_bw()+
#  scale_color_viridis_d(option="turbo")+
  theme(text=element_text(size=xTextSize),
        legend.position = "none",
        plot.title = element_text(hjust=0.5))+
  #labs(x="B-24025", y="B-2880", title="Globi EOPs")
  labs(x="B-24025", y="B-2880")
pHostRange24.EOP.globi

# globi v oryzae
pHostRange24.EOP.globi.oryzae1<-HostRange24.EOP.globi %>%
  ggplot(aes(x=B24025, y=B24478))+
  geom_point(size=2)+
  xlim(-8, 1)+
  ylim(-8, 1)+
  #scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, color="grey20")+
  geom_vline(xintercept = 0, color="grey20")+
  geom_abline(slope=1, intercept=0, color="grey80")+
  geom_abline(slope=lm.24478.24025$coefficients[[2]], intercept = lm.24478.24025$coefficients[[1]], linetype="dashed", color="red")+  
  theme_bw()+
  theme(text=element_text(size=xTextSize),
        legend.position = "none",
        plot.title = element_text(hjust=0.5))+
  #labs(x="B-24025", y="B-24478", title="24025 vs oryzae")
  labs(x="B-24025", y="B-24478")
pHostRange24.EOP.globi.oryzae1

pHostRange24.EOP.globi.oryzae2<-HostRange24.EOP.globi %>%
  ggplot(aes(x=B2880, y=B24478))+
  geom_point(size=2)+
  xlim(-8, 1)+
  ylim(-8, 1)+
#  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, color="grey20")+
  geom_vline(xintercept = 0, color="grey20")+
  geom_abline(slope=1, intercept=0, color="grey80")+
  geom_abline(slope=lm.2880.24478$coefficients[[2]], intercept = lm.2880.24478$coefficients[[1]], linetype="dashed", color="red")+
  theme_bw()+
  theme(text=element_text(size=xTextSize),
        legend.position = "none",
        plot.title = element_text(hjust=0.5))+
  #labs(x="B-2880", y="B-24478", title="2880 vs oryzae")
  labs(x="B-2880", y="B-24478")
pHostRange24.EOP.globi.oryzae2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### summary
# how many cases (fraction of phage) have plaque formation on each host?
HostRange24.phagesummary.all<-HostRange24.all %>%
  dplyr::filter(StrainLabel != "TARGET" & Set==1) %>%
  group_by(Host, HostSpecies, Bacteria) %>%
  summarize(nPhage=n(),
            countAny=sum(!is.na(EOP)),
            countLFW=sum(LFW=="Y"),
            countPlaque=countAny-countLFW,
            fAny=countAny/nPhage,
            fPlaque=countPlaque/nPhage)

pAnyPlaques_bar_all<-HostRange24.phagesummary.all %>%
  dplyr::filter(Bacteria != "B2884" & Bacteria != "B1814" & Bacteria != "B4425") %>%
  ggplot(aes(x=Bacteria, y=fAny, fill=factor(Bacteria)))+
  geom_col()+
  ylim(0,1)+
  theme_bw()+
  scale_fill_viridis_d(option="magma", end=0.85)+
  theme(text=element_text(size=xTextSize),
        #axis.text.x = element_text(angle=90),
        legend.position = "none",
        plot.title = element_text(hjust=0.5))+
  facet_wrap(~HostSpecies, scales="free_x", nrow=1)+
  labs(y="Fraction of Phage", x="", fill="Host \nStrain", title="Any Plaque Formation")
pAnyPlaques_bar_all
#ggsave("GlobiHR_AnyPlaques_bar.png", width=12, height=5, units="in", dpi=300)

# plot out for countable plaquing only
pIndivPlaques_bar_all<-HostRange24.phagesummary.all %>%
  dplyr::filter(Bacteria != "B24479" & Bacteria != "B2884" & Bacteria != "B1814" & Bacteria != "B4425") %>%
  ggplot(aes(x=Bacteria, y=fPlaque, fill=factor(Bacteria)))+
  geom_col()+
  theme_bw()+
  ylim(0,1)+
  scale_fill_viridis_d(option="magma", end=0.85)+
  theme(text=element_text(size=xTextSize),
        #axis.text.x = element_text(angle=90),
        plot.title = element_text(hjust=0.5),
        legend.position = "none")+
  facet_wrap(~HostSpecies, scales="free_x", nrow=1)+
  labs(y="Fraction of Phage", x="", fill="Host \nStrain", title="Individual Plaques")
pIndivPlaques_bar_all


# Ratios of EOPs
HostRange24.EOP.globi <- HostRange24.EOP.globi |>
  dplyr::mutate(B2880vB24025 = B2880-B24025)

view(HostRange24.EOP.globi)

#~~~~~~~~~~~~~~~~~~~
# composite plot
#pHostRange24.EOP.multipanel<-plot_grid(pHostRange24.EOP.globi, pHostRange24.EOP.globi.oryzae1, pHostRange24.EOP.globi.oryzae2,
#          ncol=1)

#(pAnyPlaques_bar_all | pEOPbyHost.all) / (pEOPbyPhage.filtered.hm | pHostRange24.EOP.multipanel)  + 
#  plot_annotation(tag_levels = 'A') #+
#  plot_layout(guides = 'collect')

pHostRange24.EOP.multipanel<-plot_grid(pHostRange24.EOP.globi, pHostRange24.EOP.globi.oryzae1, pHostRange24.EOP.globi.oryzae2,
                                          nrow=1)
                                       
(pAnyPlaques_bar_all |  pEOPbyHost.all) / (pEOPbyPhage.filtered.hm / pHostRange24.EOP.multipanel)  + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights=c(1,4))

ggsave("GlobiPhageHostRangeSummaryPlotTall.png", width=14, height=14, units="in", dpi=500)

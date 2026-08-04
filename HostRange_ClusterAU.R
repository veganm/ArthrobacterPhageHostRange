### Code for plotting host range data for cluster AU phage TrixiePhattel (AU4) and Anandi (AU6)
pacman::p_load(tidyverse, cowplot, ggpubr, RColorBrewer,
               readxl, patchwork)

xTextSize<-14 # default text size 

# read in the plating data from this experiment
HostRange_ClusterAU<-read_xlsx("HostRange_ClusterAU.xlsx")
glimpse(HostRange_ClusterAU)

# remove any whitespaces
#HostRange_ClusterAU$BactSpecies<-str_squish(HostRange_ClusterAU$BactSpecies)
#HostRange_ClusterAU$Bacteria<-str_squish(HostRange_ClusterAU$Bacteria)

# note that combinations that resulted in no plaques
# have been given a log10(EOP) of -8 for plotting

# reorder
HostRange_ClusterAU$Bacteria<-factor(HostRange_ClusterAU$Bacteria, levels=c('B24025', 'B2880', 'B2979',
                                                                            'B24478', 'B24479', 'B1814'))

pHostRange_Anandi<-HostRange_ClusterAU %>%
  dplyr::filter(Phage=="Anandi") %>%
  ggplot(aes(x=logEOP, y=Bacteria, color=factor(BactSpecies), pch = factor(Rep)))+
  geom_jitter(width=0.1, size=3)+
  geom_vline(xintercept=0, color="red")+
  theme_bw()+
  scale_color_viridis_d(option="plasma", end=0.85)+
  #scale_color_manual(values=c("#30123BFF" ,"#28BBECFF" ,"#31F299FF" ,"#D23105FF" ,"#7A0403FF")) +
  labs(color="Species", pch = "Replicate", x="log10EOP", y="Isolate")+
  theme(text=element_text(size=xTextSize))+
  facet_grid(cols=vars(Phage))

pHostRange_TP<-HostRange_ClusterAU %>%
  dplyr::filter(Phage=="TrixiePhattel") %>%
  ggplot(aes(x=logEOP, y=Bacteria, color=factor(BactSpecies), pch = factor(Rep)))+
  geom_jitter(width=0.1, size=3)+
  geom_vline(xintercept=0, color="red")+
  theme_bw()+
  scale_color_viridis_d(option="plasma", end=0.85)+
  #scale_color_manual(values=c("#30123BFF" ,"#28BBECFF" ,"#31F299FF" ,"#D23105FF" ,"#7A0403FF")) +
  labs(color="Species", pch = "Replicate", x="log10EOP", y="")+
  theme(text=element_text(size=xTextSize))+
  facet_grid(cols=vars(Phage))

pHostRange_Anandi + pHostRange_TP +
  plot_layout(guides="collect")

ggsave("pHostRange_ClusterAU.svg", width=110, height = 27, units="mm")

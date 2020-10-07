####analysis of plant and animal domesticates for EN in Adriatic and Danube basins####
####used for paper on dispersal and niche construction (Springer, Sept. 2020)####
####load libraries####
library(ggplot2) # data plotting
library(dplyr) # data manipulation
library(vegan) # data analysis
library(ca) # correspondence analysis
library(factoextra) # correspondence analysis tools
library(ggrepel) # plotting tool
library(ggtern) # ternary graphs

####load data####
en.animals<-read.csv("/home/marc/Dropbox/m_pers/mobility_20/data/en_animals.csv")
en.animals$bostaur<-as.integer(en.animals$bostaur)
en.animals<-mutate(en.animals,SampleID=row_number()) # add id number to each entry
en.animals$SampleID<-as.character(en.animals$SampleID)
en.plants<-read.csv("/home/marc/Dropbox/m_pers/mobility_20/data/en_plants.csv")
en.plants<-en.plants%>%mutate_if(is.integer,as.numeric)
en.plants<-mutate(en.plants,SampleID=row_number()) # add id number to each entry
en.plants$SampleID<-as.character(en.plants$SampleID)

####alpha and beta diversity for plants + CA####
#species richness for plants only (alpha diversity)
en.plants$richness<-as.numeric(rowSums(en.plants[,6:17]!="0"))
plot.richness.group<-ggplot(en.plants,aes(Group, richness))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_jitter(height=0,width=0.1)
plot.richness.group+
  ylab("Species richness")+
  scale_x_discrete(name="Geographical groups",labels=c("Adriatic","Danube"))+
  theme_bw()
# create regional datasets
plants.adriatic<-filter(en.plants,Group=="Coast") # plants for Adriatic only
plants.danube<-filter(en.plants,Group=="Inland") # plants for Danube basin only
#testing for statistical differences
shapiro.test(plants.adriatic$richness) # check for normality
shapiro.test(plants.danube$richness)
t.test(plants.adriatic$richness,plants.danube$richness)
#beta diversity plants
beta.plants<-vegdist(en.plants[,6:17],method="bray",binary=TRUE) # compute bray-sorensen pairwise comparisons across all data
mds.plants<-metaMDS(beta.plants) # ordination
mds.plants.df<-as.data.frame(mds.plants$points) #recreate as df for plotting
mds.plants.df$SampleID<-rownames(mds.plants.df) #recreate as df for plotting
mds.plants.df<-left_join(mds.plants.df,en.plants) #recreate as df for plotting
mds.plants.plot<-ggplot(mds.plants.df,aes(x=MDS1,y=MDS2,color=Group))+ # plot by group
  geom_point()
mds.plants.plot+
  labs(colour="Geographical\nGroups")+
  scale_colour_hue(labels=c("Adriatic","Danube"))+
  theme_bw()
#compare beta diversity values by regions
beta.plants.adriatic<-vegdist(plants.adriatic[,6:17],method="bray",binary=TRUE) # compute bray-sorensen pairwise comparisons for Adriatic
shapiro.test(beta.plants.adriatic) # check for normality
beta.plants.danube<-vegdist(plants.danube[,6:17],method="bray",binary=TRUE) # compute bray-sorensen pairwise comparisons for Inland
shapiro.test(beta.plants.danube) #check for normality
wilcox.test(beta.plants.adriatic, beta.plants.danube)
beta.plants.danube.2<-as.numeric(beta.plants.danube)
beta.plants.danube.2<-as_tibble(beta.plants.danube.2)
beta.plants.danube.2<-mutate(beta.plants.danube.2,group="Danube")
beta.plants.adriatic.2<-as.numeric(beta.plants.adriatic)
beta.plants.adriatic.2<-as_tibble(beta.plants.adriatic.2)
beta.plants.adriatic.2<-mutate(beta.plants.adriatic.2,group="Adriatic")
beta.plants<-bind_rows(beta.plants.adriatic.2,beta.plants.danube.2)
plot.beta.plant<-ggplot(beta.plants,aes(x=group,y=value))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  scale_x_discrete(name="Geographical groups")
plot.beta.plant+
  ylab("Sørensen")+
  theme_bw()
#correspondence analysis
ca.plants <- ca(en.plants[,c(6:9,11:12)])
plot(ca.plants)
print(ca.plants)
#extract percentage of variance and round-up to two decimals to use for automatic legend of CA plot
eigenvalues <- get_eigenvalue(ca.plants)
var.percent<-data.frame(eigenvalues$variance.percent)
Dim1.raw<-var.percent[1,1]
Dim1<-round(Dim1.raw, digits=2)
Dim2.raw<-var.percent[2,1]
Dim2<-round(Dim2.raw, digits=2)
#get column variables
col<-get_ca_col(ca.plants)
row<-get_ca_row(ca.plants)
# create data.frame combining CA results and necessary information for plotting (e.g. phase, group)
x1<-row$coord[,1]
x2<-row$coord[,2]
y1<-col$coord[,1]
y2<-col$coord[,2]
group<-data.frame(en.plants$Group)
bioregion<-data.frame(en.plants$Bioregion)
sites.scores <- data.frame(cbind(group,bioregion,x1, x2))
names(sites.scores) <- c("group","bioregion","x1","x2")
colnames <- data.frame(names(en.plants[,c(6:9,11:12)]))
colscores <- data.frame(cbind(y1,y2))
species.scores <- cbind(colnames,colscores)
names(species.scores) <- c("species", "y1", "y2")
species.scores$species<-c("barley", "emmer","einkorn","free,threshing.wheat",
                          "lentil", "pea")
sites.scores$group<-as.factor(sites.scores$group)
sites.scores$bioregion<-as.factor(sites.scores$bioregion)
#plot species
plot.ca.group <- ggplot(sites.scores, aes(x=x1, y=x2))+
  geom_blank()+
  xlab(label=paste("Dimension 1 (",Dim1,"%)"))+
  ylab(label=paste("Dimension 2 (",Dim2,"%)"))+  
  geom_point(data=sites.scores, aes(x=x1, y=x2,colour=group),size=2.8)+
  geom_point(data=species.scores, aes(x=y1, y=y2), shape=3, size=2)+
  geom_text(data=species.scores,aes(x=y1, y=y2), label=species.scores$species,vjust=0,nudge_y=-0.06)
plot.ca.group +
  labs(colour="Geographical\nGroups")+
  scale_colour_hue(labels=c("Adriatic","Danube"))+
  theme_bw()+
  coord_fixed()+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted")+
  theme(
    legend.position=c(0.02,0.98),
    legend.justification=c(0,1),
    legend.background=element_blank(),
    legend.box.background=element_rect(colour="black"))
####alpha and beta diversity for animals + ternary graphs####
en.animals$simpson<-diversity(en.animals[,10:12],index="simpson")
plot.richness.animals.group<-ggplot(en.animals,aes(Zone, simpson))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_jitter(height=0,width=0.1)+
  scale_x_discrete(name="Geographical groups", labels=c("Adriatic","Danube"))
plot.richness.animals.group+
  ylab("Simpson")+
  theme_bw()
#create regional datasets
animals.adriatic<-filter(en.animals,Zone=="Coast") # animals for Adriatic only
animals.danube<-filter(en.animals,Zone=="Inland") # plants for Danube basin only
shapiro.test(animals.adriatic$simpson)
shapiro.test(animals.danube$simpson)
wilcox.test(animals.adriatic$simpson,animals.danube$simpson)
#beta diversity animals
beta.animals<-vegdist(en.animals[,10:12],method="bray",binary=FALSE) # compute bray-sorensen pairwise comparisons across all data
mds.animals<-metaMDS(beta.animals) # ordination
mds.animals.df<-as.data.frame(mds.animals$points) #recreate as df for plotting
mds.animals.df$SampleID<-rownames(mds.animals.df) #recreate as df for plotting
mds.animals.df<-left_join(mds.animals.df,en.animals) #recreate as df for plotting
mds.animals.plot<-ggplot(mds.animals.df,aes(x=MDS1,y=MDS2,color=Zone))+ # plot by group
  geom_point()
mds.animals.plot+
  labs(colour="Geographical\nGroups")+
  scale_colour_hue(labels=c("Adriatic","Danube"))+
  theme_bw()
#compare beta diversity values by regions
beta.animals.adriatic<-vegdist(animals.adriatic[,10:12],method="bray",binary=FALSE) # compute bray-sorensen pairwise comparisons for Adriatic
shapiro.test(beta.animals.adriatic)
beta.animals.danube<-vegdist(animals.danube[,10:12],method="bray",binary=FALSE) # compute bray-sorensen pairwise comparisons for Inland
shapiro.test(beta.animals.danube)
wilcox.test(beta.animals.danube,beta.animals.adriatic)
beta.animals.danube.2<-as.numeric(beta.animals.danube)
beta.animals.danube.2<-as_tibble(beta.animals.danube.2)
beta.animals.danube.2<-mutate(beta.animals.danube.2,group="Danube")
beta.animals.adriatic.2<-mutate(beta.animals.adriatic.2,group="Adriatic")
beta.animals.adriatic.2<-as.numeric(beta.animals.adriatic)
beta.animals.adriatic.2<-as_tibble(beta.animals.adriatic.2)
beta.animals.adriatic.2<-mutate(beta.animals.adriatic.2,group="Adriatic")
beta.animals<-bind_rows(beta.animals.adriatic.2,beta.animals.danube.2)
plot.beta.animal<-ggplot(beta.animals,aes(x=group,y=value))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  scale_x_discrete(name="Geographical groups")
plot.beta.animal+
  ylab("Bray-Curtis")+
  theme_bw()
#ternary graphs
ptern <- ggtern(en.animals,aes(x=susscrd, y=bostaur, z=oviscap)) +
  geom_point(aes(color=Zone))
ptern +
  labs(colour="Geographical\nGroups")+
  scale_colour_hue(labels=c("Adriatic","Danube"))+
  theme_bw() +
  theme_showarrows() +
  Tlab("Cattle") + Rlab("O/C") + Llab("Pigs") + Wlab("%")

########################################################################################################################
# This script is written by Mohammad Tauqeer Alam. Here we analyse intracellular enzyme-metabolite activation network  #
# We use Yeast9 model as the basis and created the activation netwok. We have then mapped the network to KEGG pathways #
# to check the position of activated enzyems. All the figures used in the paper can be generated using this code.      #
# For using Yeast9 model please cite the original paper.                                                               #
# For using Any tools, cite their original paper.                                                                      #
# For using this code from our analysis, cite the current paper.                                                       #
########################################################################################################################
library(reshape2)
library(stringr)
library(pheatmap)
library(igraph)
library(RColorBrewer)
library(ggplot2)
library(poweRlaw)
mycolor = c(brewer.pal(12,"Set3"),brewer.pal(8,"Dark2"),brewer.pal(12,"Paired"))

#### Set the path for working directoy
setwd("/Users/mtalam/OneDrive - UAE University/UAEU/research/masters_project/sultana/")

#############################################################################
##### Loading activation network, Yeast9 model, expression data and more
#############################################################################
load("./yeast_GEM9/data_for_github/gene_expression_all_genome2.Rdata") #variable: exp.all2 all conditions used for average expression variation
genetic.int.costanzo = read.delim("./yeast_GEM9/data_for_github/sgadata_costanzo2009_stringentCutoff_101120.txt",header = F) #genetic interactions network from costanzo et al 2009
load("./yeast_GEM9/data_for_github/gene_expression_all_avgConditionWise_genome2.Rdata") # variable:exp.all.avgConditionWise. This was used to create correlation between genes
load("./yeast_GEM9/data_for_github/gene_co_expression_avgConditionWise_correlation.Rdata") # variable: exp.cor, correlation between genes

model.met = read.csv("./yeast_GEM9/data_for_github/model_met.csv") #metabolites of the Yeast9 model with details
ec.list = list.files("./EC_activator_Yeast9/") # list of enzymes
act.net = read.csv("./yeast_GEM9/data_for_github/activationNetwork.csv") # activation network
act.deg = read.csv("./yeast_GEM9/data_for_github/act_deg.csv") # degree of activators in the activation network. also contain essentiality short path length 
ec.details = read.csv("./yeast_GEM9/data_for_github/ec_details.csv") # details of enzymes of the activation network.
ec.class = read.csv("./yeast_GEM9/data_for_github/ec_class.csv") # enzyme classes

model_metabolic = readxl::read_excel("./yeast_GEM9/yeast-GEM-9.0.0/model/yeast-GEM.xlsx",sheet = "RXNS") # Yeast9 metabolic model

##### split the EC lists of a rxns and make separate rows for separate Enzyme
model_metabolic2 = c() # this will be used for EC comparisons later and getting rxns
for(i in 1:dim(model_metabolic)[1]) {
  ec_new = unlist(strsplit(model_metabolic$`EC-NUMBER`[i],";"))
  my.rxns = cbind(model_metabolic[i,],ec_new)
  model_metabolic2 = rbind(model_metabolic2,my.rxns)
}
model_metabolic2 = model_metabolic2[!is.na(model_metabolic2$ec_new),]

#############################################################################
##### gene, metabolites, enzyme knockout analysis, predicting essentiality
##### KOs were done using COBRAtoolbox in MATLAB by switching reactions off.
#############################################################################
############ Essentiality of genes
gene.ko.input = read.delim("./yeast_GEM9/data_for_github/KO_Input_GEM9_Gene.txt",header = F) #gene KO input for COBRAtoolbox
gene.ko.output = read.csv("./yeast_GEM9/data_for_github/Gene_KO_output_GEM9.txt",header = F,sep = "\t") #gene KO output from COBRAtoolbox
gene.ko.gr = gene.ko.output[gene.ko.output$V1 == "r_4041",]
gene.ko.gr = gene.ko.gr[-1]
gene.ko.gr.ratio = as.numeric(gene.ko.gr)/as.numeric(gene.ko.gr[1])
essential.gene = gene.ko.input$V1[which(gene.ko.gr.ratio<0.1)-1]
essentialGenes = essential.gene

############ Essentiality of enzymes
ec.ko.input = read.delim("./yeast_GEM9/data_for_github/KO_Input_GEM9_Enzyme.txt",header = F) #enzyme KO input for COBRAtoolbox
ec.ko.output = read.csv("./yeast_GEM9/data_for_github/Enzyme_KO_output_GEM9.txt",header = F,sep = "\t") #enzyme KO output from COBRAtoolbox
ec.ko.gr = ec.ko.output[ec.ko.output$V1 == "r_4041",]
ec.ko.gr = ec.ko.gr[-1]
ec.ko.gr.ratio = as.numeric(ec.ko.gr)/as.numeric(ec.ko.gr[1])
essential.ec = ec.ko.input$V1[which(ec.ko.gr.ratio<0.1)-1]

############ Essentiality of activators
met.ko.input = read.delim("./yeast_GEM9/data_for_github/KO_Input_GEM9_Metabolite.txt",header = F) #metabolite KO input for COBRAtoolbox
met.ko.output = read.delim("./yeast_GEM9/data_for_github/metabolite_KO_output_GEM9.txt",header = F,sep = "\t") #metabolite KO output from COBRAtoolbox
colnames(met.ko.output) = c("ID","WT",met.ko.input$V1)
rownames(met.ko.output) = met.ko.output$ID
met.ko.output = met.ko.output[,-1]
met.ko.gr = met.ko.output[which(rownames(met.ko.output)=='r_4041'),]
met.ko.gr = as.numeric(met.ko.gr)
wt.gr = met.ko.gr[1]
met.ko.gr = met.ko.gr[-1]
met.ko.gr.ratio = as.numeric(met.ko.gr)/as.numeric(wt.gr)
essential.met = met.ko.input$V1[which(met.ko.gr.ratio<0.1)]
essential.met.kegg = model.met$KEGG[model.met$ID %in% essential.met]
essential.met.kegg = essential.met.kegg[!essential.met.kegg==""]
essential.met.kegg = unique(essential.met.kegg)

#############################################################################
######### activation NETWORK
#############################################################################
g <- graph.data.frame(act.net, directed=T)
nodes  <- rownames(as.matrix(V(g)))
edge.list <- as.data.frame(get.edgelist(g))
degree <- as.matrix(igraph::degree(g),ncol=1) #print the number of edges per vertex
degree <- as.data.frame(degree)
colnames(degree) <- "Degree"
V(g)$label.cex = 0.5

pdf("./yeast_GEM9/data_for_github/figures/figure 1b activation network.pdf",width=6, height=6);
plot(g,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
     mark.border="#BCBCBC",pch=19,vertex.label=NA,
     vertex.color = ifelse(nodes %in% paste0("EC_",essential.ec), "#CC0000",
                           ifelse(nodes %in% essential.met.kegg,"#CC0000","#BCBCBC")),
     vertex.shape = ifelse(nodes %in% act.net$EC,"square","circle"),
     layout = layout_with_fr, main = "activation network")
legend("topleft",legend = c("essential","non-essential"),
       fill = c("#CC0000","#BCBCBC"),bty = "n")
dev.off()


#############################################################################
######## Enzymes and metabolites from the Yeast9 model and activation network
#############################################################################
############ proportion of intra- and extracellular activated and non-activated enzymes
ec.list.activated = c()
for(i in 1:length(ec.list)) {
  if(!is.null(tryCatch(read.table(paste0("./EC_activator_Yeast9/",ec.list[i])), error=function(e) NULL))) {
    ec.list.activated = c(ec.list.activated,ec.list[i])
  }
}

ec.list.activated = gsub("EC_","",gsub(".txt","",ec.list.activated))
ec.list.notactivated = gsub("EC_","",gsub(".txt","",ec.list))
ec.list.notactivated = ec.list.notactivated[!ec.list.notactivated %in% ec.list.activated]

pdf("./yeast_GEM9/data_for_github/figures/figure 1c Proportion of enzymes activation.pdf",width=6, height=6);
par(mar=c(10,4,4,2))
pie(c(length(unique(act.net$EC)),
      length(ec.list.activated)-length(unique(act.net$EC)),
      length(ec.list)-length(ec.list.activated)), 
    c("intracellular activated","extracellular activated","not activated"),col = c("#9FC5E8","#BCBCBC","white"),
    main = "Proportion of enzymes activation")
dev.off()

############ proportion of activator and non-activator metabolites of the Yeast9 model
all.met = length(unique(unlist(lapply(model.met$ID,function(x)
  paste(unlist(strsplit(x,"\\["))[-length(unlist(strsplit(x,"\\[")))],collapse = "\\[")))))

pdf("./yeast_GEM9/data_for_github/figures/figure 1d count activators and non activator metabolites of the model PIE.pdf",width=6, height=6);
pie(c(dim(act.deg)[1],
      all.met-dim(act.deg)[1]))
dev.off()

############ Metabolite classes of Activators and non activator metabolites of the model
act.met.class.count = as.data.frame(table(unlist(lapply(act.net$Activators, function(x)
  act.deg$class[act.deg$Activators==x]))))
colnames(act.met.class.count) = c("Met_class","act.class.count")

all.met.class.count = as.data.frame(table(unlist(lapply(unique(model.met$KEGG), function(x)
  model.met$class[which(model.met$KEGG==x)][1]))))
colnames(all.met.class.count) = c("Met_class","allMet.class.count")

cell.act = as.data.frame(table(act.deg$class))
colnames(cell.act) = c("Met_class","cell.act.class")

met.class.counts = merge(merge(cell.act,act.met.class.count,by = "Met_class"),all.met.class.count,by = "Met_class")

pdf("./yeast_GEM9/data_for_github/figures/figure 1e activators and non activator metabolites of the model.pdf",width=6, height=6);
par(mar=c(15,4,4,2))
barplot(rbind(met.class.counts$cell.act.class/sum(met.class.counts$cell.act.class)*100,
              met.class.counts$act.class.count/sum(met.class.counts$act.class.count)*100,
              met.class.counts$allMet.class.count/sum(met.class.counts$allMet.class.count)*100)
        ,beside = T,
        main = "activators and non activator metabolites of the model",
        ylab = "% of metabolites from the class",
        names.arg = met.class.counts$Met_class,las=2,col = c("#9FC5E8","black","white"))

legend("topleft",legend = c("activators","activation interactions","non-activators"),
       fill = c("#9FC5E8","black","white"),bty = "n")
dev.off()

############ Enzyme Classes of Activated and non-activated enzymes of the model
model.ec = unique(model_metabolic2$ec_new)
model.ec = model.ec[!model.ec==""]
non.activated.EC = model.ec[!paste0("EC_","",model.ec) %in% act.net$EC]

count.ec.class.interactions = table(gsub("EC_","",unlist(lapply(act.net$EC, function(x)
  unlist(strsplit(x,"\\."))[1]))))

count.activated.ec.class = table(gsub("EC_","",unlist(lapply(ec.details$EC, function(x)
  unlist(strsplit(x,"\\."))[1]))))

count.nonactivated.ec.class = table(unlist(lapply(non.activated.EC, function(x)
  unlist(strsplit(x,"\\."))[1])))

pdf("./yeast_GEM9/data_for_github/figures/figure 1f activated and non-activated enzymes of the model.pdf",width=6, height=6);
barplot((rbind(as.numeric(count.activated.ec.class)/sum(count.activated.ec.class)*100,
               as.numeric(count.ec.class.interactions)/sum(count.ec.class.interactions)*100,
               as.numeric(count.nonactivated.ec.class)/sum(count.nonactivated.ec.class)*100)),beside = T,
        main = "activated and non-activated enzymes of the model",
        names.arg = paste0("EC:",names(count.activated.ec.class),"..."),las=2,col = c("#9FC5E8","black","white"),
        ylab = "% of activation and non activation")
legend("topleft",legend = c("activated","activation interactions","non-activated"),
       fill = c("#9FC5E8","black","white"),bty = "n")
dev.off()


#############################################################################
######## Distribution of Activated Enzymes in Yeast9 model
#############################################################################
list.pathways = as.data.frame(cbind(SUBSYSTEM = unique(model_metabolic$SUBSYSTEM),
                                    rxn_count = unlist(lapply(unique(model_metabolic$SUBSYSTEM), function(x)
                                      sum(model_metabolic$SUBSYSTEM == x)))))

list.pathways = list.pathways[list.pathways$rxn_count>1,] # remove pathways with less than 2 rxns

############ ignore Transport pathways, signalling
list.pathways = list.pathways[-grep("Transport",list.pathways$SUBSYSTEM),]
list.pathways = list.pathways[-grep("signaling",list.pathways$SUBSYSTEM),]
list.pathways = list.pathways[-grep("Other",list.pathways$SUBSYSTEM),]

activated.pathways = as.data.frame(cbind(
  activated_all = unlist(lapply(list.pathways$SUBSYSTEM, function(x)
    length(unique(ec.list.activated[ec.list.activated %in% 
                                      unique(unlist(strsplit(model_metabolic$`EC-NUMBER`[model_metabolic$SUBSYSTEM == x],";")))
    ])))),
  
  activated_intracellular = unlist(lapply(list.pathways$SUBSYSTEM, function(x)
    length(unique(
      gsub("EC_","",unique(act.net$EC))[gsub("EC_","",unique(act.net$EC)) %in% 
                                          unique(unlist(strsplit(model_metabolic$`EC-NUMBER`[model_metabolic$SUBSYSTEM == x],";")))])))),
  
  not_activated = unlist(lapply(list.pathways$SUBSYSTEM, function(x)
    length(unique(ec.list.notactivated[ec.list.notactivated %in% unique(unlist(strsplit(model_metabolic$`EC-NUMBER`[model_metabolic$SUBSYSTEM == x],";")))]))))))

rownames(activated.pathways)=list.pathways$SUBSYSTEM

activated.pathways = activated.pathways[apply(activated.pathways, 1, function(x) any(x>0)),]

pdf("./yeast_GEM9/data_for_github/figures/figure 2a Proportion of pathways activation.pdf",width=6, height=6);
par(mar=c(10,4,4,2))
pie(c(sum(activated.pathways$activated_intracellular>0),
      sum(activated.pathways$activated_all>0 & activated.pathways$activated_intracellular == 0),
      dim(activated.pathways)[1]-sum(activated.pathways$activated_all>0)), 
    c("intracellular activated","extracellular activated","not activated"),col = c("#9FC5E8","#BCBCBC","white"),
    main = "Proportion of pathways activation")
dev.off()

#############################################################################
######## Position of activated enzymes in different pathway (KEGG)
#############################################################################
kegg.path.sce = read.delim2("./yeast_GEM9/data_for_github/sce00001.txt",row.names=NULL,header = F)
kegg.path.sce2 = kegg.path.sce[-grep("^B",kegg.path.sce$V1),]

kegg.path.sce2 = kegg.path.sce2[-c(grep("#",kegg.path.sce2$V1)[1]:length(kegg.path.sce2$V1)),] # remove non-metabolic

kegg.path.names = kegg.path.sce2$V1[kegg.path.sce2$V2==""]

kegg.path.names=unlist(lapply(kegg.path.names, function(x)
  paste(unlist(strsplit(x," "))[-c(1:5,length(unlist(strsplit(x," "))))],collapse = " ")))
kegg.path.names[17] = kegg.path.sce2$V1[kegg.path.sce2$V2==""][17]

kegg.path.index = which(kegg.path.sce2$V2=="")

kegg.path.ec = vector("list",length(kegg.path.names))
names(kegg.path.ec) = kegg.path.names
i=kegg.path.index[1]
x = 1

for(j in kegg.path.index[-1]) {
  my.path = kegg.path.sce2[i:j,]
  my.path.name = my.path$V1[1]
  my.path = my.path[-1,]
  my.path.ec = unlist(lapply(my.path$V2, function(x)
    unlist(strsplit(str_match(x, "\\[EC.*\\]")," "))))
  
  my.path.ec = my.path.ec[!duplicated(my.path.ec)]
  my.path.ec = my.path.ec[!is.na(my.path.ec)]
  my.path.ec = gsub("EC:|\\[|\\]","",my.path.ec)
  #### consider only those EC which are in the model
  my.path.ec = my.path.ec[my.path.ec %in% gsub("EC_|.txt","",ec.list)]
  kegg.path.ec[[x]] = my.path.ec
  i = j
  x = x+1
}

############ ignore pathways with less than 3 enzymes
kegg.path.ec = kegg.path.ec[lapply(kegg.path.ec, length)>2]
length(kegg.path.ec)

count.act.kegg.path.ec = unlist(lapply(kegg.path.ec, function(x)
  sum(gsub("EC_","",unique(act.net$EC)) %in% x)))

kegg.path.ec = kegg.path.ec[count.act.kegg.path.ec>0]
kegg.path.ec = kegg.path.ec[order(as.numeric(unlist(lapply(kegg.path.ec, length))))]

############ check where activated enzymes are in the pathway
activated.ec.position.path = vector("list",length(kegg.path.ec))
activated.ec.position.path.degree = vector("list",length(kegg.path.ec))
activated.ec.position.path.ess =  vector("list",length(kegg.path.ec))
activated.ec.in.path.count = vector("list",length(kegg.path.ec))
activated.ec.path.position.1.3 = matrix(0,nrow=length(kegg.path.ec),ncol=3)
activated.ec.essMet = vector("list",length(kegg.path.ec))
activated.ec.in.path = vector("list",length(kegg.path.ec))

names(activated.ec.position.path) = names(kegg.path.ec)
names(activated.ec.position.path.ess) = names(kegg.path.ec)

for(i in 1:length(kegg.path.ec)) {
  my.path.ec = kegg.path.ec[[i]]
  my.path.act.pos = which(my.path.ec %in% gsub("EC_","",unique(act.net$EC)))
  #### activated EC in pathways
  activated.ec.in.path[[i]] = my.path.ec[my.path.ec %in% gsub("EC_","",unique(act.net$EC))]
  
  activated.ec.path.position.1.3[i,1] = sum(1 %in% my.path.act.pos)
  activated.ec.path.position.1.3[i,2] = sum(2 %in% my.path.act.pos)
  activated.ec.path.position.1.3[i,3] = sum(3 %in% my.path.act.pos)
  
  activated.ec.position.path[[i]] = my.path.act.pos
  #### degree in the activation network
  activated.ec.position.path.degree[[i]] = unlist(lapply(my.path.ec[which(my.path.ec %in% gsub("EC_","",unique(act.net$EC)))], function(x)
    ec.details$degree_activated[gsub("EC_","",ec.details$EC)==x]))
  #### essential in the activation network
  activated.ec.position.path.ess[[i]] = unlist(lapply(my.path.ec[which(my.path.ec %in% gsub("EC_","",unique(act.net$EC)))], function(x)
    x %in% essential.ec))
  #### enzymes in number of pathways counted here
  activated.ec.in.path.count[[i]]=unlist(lapply(my.path.ec, function(x)
    sum(unlist(lapply(kegg.path.ec, function(y)
      sum(x%in%y))))))
}

pdf("./yeast_GEM9/data_for_github/figures/figure 2b activated enzymes at the begining of the pathway.pdf",width=6, height=6);
pie(c(sum(apply(activated.ec.path.position.1.3, 1, sum)>0),
      dim(activated.ec.path.position.1.3)[1]-sum(apply(activated.ec.path.position.1.3, 1, sum)>0)),
    c("first 3 EC activated","not activated"),col = c("black","white"),
    main = "activated enzymes at the begining of the pathway")
dev.off()

pdf("./yeast_GEM9/data_for_github/figures/figure 2c position of activated enzyme in kegg pathways SUPPLEMENT.pdf",width=14, height=8);
par(mar=c(20,4,4,2))
activated.ec.kegg.count = unlist(lapply(kegg.path.ec, length))
df.bar = barplot(activated.ec.kegg.count,col="#9FC5E8",las=2)
for(i in 1:length(kegg.path.ec)) {
  points(rep(df.bar[i],length(activated.ec.position.path[[i]])),activated.ec.position.path[[i]],
         pch=19,col=ifelse(activated.ec.position.path.ess[[i]],"red","blue"),
         cex = as.numeric(activated.ec.position.path.degree[[i]])/(max(as.numeric(activated.ec.position.path.degree[[i]])))*1.4) 
}
dev.off()

#############################################################################
######## Significant of one pathway activating others
#############################################################################
list.pathways2 = list.pathways
list.pathways2 = list.pathways2[!list.pathways2$SUBSYSTEM == "SLIME reaction",]
list.pathways2 = list.pathways2[!list.pathways2$SUBSYSTEM == "Exchange reaction",]
list.pathways2 = list.pathways2[!list.pathways2$SUBSYSTEM == "Growth",]

path.ec = c() # list of EC for the pathway
path.prod = c() # list of metabolites produced in the pathway
path.act.prod = c() # list of metabolie products which also act as activator for different enzymes
path.ec.act = c() # list of activators activating enzymes of the pathyway

for(i in 1:dim(list.pathways2)[1]) {
  my.path.ec = c() # list of EC for the pathway
  my.path.prod = c() # list of metabolites produced in the pathway
  my.path.act.prod = c() # list of metabolie products which also act as activator for different enzymes
  my.path.ec.act = c() # list of activators activating enzymes of the pathyway
  
  my.path = model_metabolic2[model_metabolic2$SUBSYSTEM==list.pathways2$SUBSYSTEM[i],]
  my.path.rev = my.path[grepl("<=>",my.path$EQUATION),]
  my.path.irr = my.path[!grepl("<=>",my.path$EQUATION),]
  my.path.prod = c(unlist(strsplit(my.path.rev$EQUATION,"<=>")),
                   unlist(strsplit(my.path.irr$EQUATION,"=>"))[2])
  my.path.prod = tolower(unique(trimws(unlist(strsplit(unlist(strsplit(my.path.prod,"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*")))))
  
  #### list of metabolites produced in the pathway
  my.path.prod = unique(model.met$KEGG[tolower(model.met$ID) %in% my.path.prod])
  
  #### list of metabolie products which also act as activator for different enzymes
  my.path.act.prod = my.path.prod[my.path.prod %in% unique(act.net$Activators)]
  
  #### list of EC for the pathway
  my.path.ec = unique(my.path$ec_new)
  
  #### list of all activators activating enzymes of this pathyway
  my.path.ec.act = unique(unlist(lapply(my.path.ec, function(x)
    unlist(strsplit(ec.details$Activators[gsub("EC_","",ec.details$EC) == x],"\\|")))))
  
  path.ec = rbind(path.ec,paste(my.path.ec,collapse = "|"))
  path.prod = rbind(path.prod,paste(my.path.prod,collapse = "|"))
  path.act.prod = rbind(path.act.prod,paste(my.path.act.prod,collapse = "|"))
  path.ec.act = rbind(path.ec.act,paste(my.path.ec.act,collapse = "|"))
  print(i)
}

list.pathways2$EC_list = path.ec
list.pathways2$met_produced = path.prod
list.pathways2$act_produced = path.act.prod
list.pathways2$EC_activators = path.ec.act

############ pathway pathway activation pvalue
pathway.activation.pv = matrix(0,nrow=length(list.pathways2$SUBSYSTEM),ncol=length(list.pathways2$SUBSYSTEM))
rownames(pathway.activation.pv)=list.pathways2$SUBSYSTEM
colnames(pathway.activation.pv)=list.pathways2$SUBSYSTEM

for(i in 1:dim(pathway.activation.pv)[1]) {
  met1 <- unlist(strsplit(list.pathways2$met_produced[i],'\\|'))
  for(j in 1:dim(pathway.activation.pv)[1]) {
    met2 <- unlist(strsplit(list.pathways2$met_produced[j],'\\|'))
    act2 <- unlist(strsplit(list.pathways2$EC_activators[j],'\\|'))
    
    if(i!=j) met1 <- met1[!(met1 %in% met2)] ## removes common metabolites
    m <- length(met1)
    n = length(unique(act.net$Activators)) - m
    k <- length(act2)
    x = length(met1[met1 %in% act2])
    pathway.activation.pv[i,j] <- phyper(x-1, m, n, k,lower.tail = F)
  } 
}
thr.pv = 0.05
pathway.activation.pv[pathway.activation.pv>thr.pv] = 0
pathway.activation.pv[pathway.activation.pv>0] = 1

############ pathway pathway activation network
path.path.net = c()
for(i in 1:dim(pathway.activation.pv)[1]) {
  for(j in 1:dim(pathway.activation.pv)[2]) {
    if(pathway.activation.pv[i,j]==1) {
      path.path.net = rbind(path.path.net,cbind(rownames(pathway.activation.pv)[i],colnames(pathway.activation.pv)[j]))
    }
  }
}
colnames(path.path.net) = c("activating","activated")

g.path.path.nw <- graph.data.frame(path.path.net, directed=T)
edge.list.path.path.nw <- as.data.frame(get.edgelist(g.path.path.nw))
nodes.path.path.nw  <- rownames(as.matrix(V(g.path.path.nw)))
degree <- as.matrix(igraph::degree(g.path.path.nw),ncol=1) #print the number of edges per vertex
degree <- as.data.frame(degree)
colnames(degree) <- "Degree"

pdf("./yeast_GEM9/data_for_github/figures/figure 2d pathway to pathway activation network.pdf",width=20, height=20);
plot(g.path.path.nw,edge.arrow.size=0.5,vertex.size=3,
     vertex.frame.color=NA,vertex.label.dist=0.5,)
dev.off()

#############################################################################
######## Degree of Essential & non-essential metabolites in the activation network
#############################################################################
nonActivator.met = model.met[!(model.met$KEGG %in% act.deg$Activators),]
nonActivator.met.gr = met.ko.output[which(rownames(met.ko.output)=='r_4041'),colnames(met.ko.output) %in% nonActivator.met$ID]
nonActivator.met.gr.ratio = as.numeric(nonActivator.met.gr)/as.numeric(wt.gr)

############ distribution of essential & nonessential enzymes and activators in the network
pdf("./yeast_GEM9/data_for_github/figures/figure 3a percentage of essential enzymes activated or not activated.pdf",width=8, height=8);
barplot(t(rbind(cbind(sum(act.deg$essentiality==1)/dim(act.deg)[1]*100,
                      sum(act.deg$essentiality==0)/dim(act.deg)[1]*100),
                cbind(sum(nonActivator.met.gr.ratio<0.1)/length(nonActivator.met.gr.ratio)*100,
                      sum(nonActivator.met.gr.ratio>0.1)/length(nonActivator.met.gr.ratio)*100),
                cbind(length(as.numeric(ec.details$degree_activated[gsub("EC_","",ec.details$EC) %in% essential.ec]))/dim(ec.details)[1]*100,
                      100 - length(as.numeric(ec.details$degree_activated[gsub("EC_","",ec.details$EC) %in% essential.ec]))/dim(ec.details)[1]*100),
                cbind(sum(non.activated.EC %in% essential.ec)/length(non.activated.EC)*100,
                      sum(!non.activated.EC %in% essential.ec)/length(non.activated.EC)*100))),
        col = c("#CC0000","#BCBCBC"),names.arg = c("act","nonAct","activated EC","nonactivated EC"),
        main = "percentage of essential enzymes activated or not activated")
dev.off()

############ degree distribution of activators
pdf("./yeast_GEM9/data_for_github/figures/figure 3b i degree distribution of activators.pdf",width=6, height=6);
plot(jitter(as.numeric(act.deg$degree_activating)), xlab="Activators",ylab = "degree of activation",
     main = "degree distribution of activators", pch = 19,
     col = ifelse(act.deg$essentiality==1,"#CC0000","#BCBCBC"))
abline(h = 15,lt=2)
legend("topright",legend = c("essential","non-essentail"),
       fill = c("#CC0000","#BCBCBC"),bty = "n")
dev.off()

pdf("./yeast_GEM9/data_for_github/figures/figure 3b ii hub metabolites degree ess and noness.pdf",width=10, height=4);
par(mfrow=c(1,3))
pie(c(sum(as.numeric(act.deg$degree_activating)>12 & act.deg$essentiality==1),
      sum(as.numeric(act.deg$degree_activating)>12 & act.deg$essentiality==0)),
    main = "degree=12")

pie(c(sum(as.numeric(act.deg$degree_activating)>15 & act.deg$essentiality==1),
      sum(as.numeric(act.deg$degree_activating)>15 & act.deg$essentiality==0)),
    main = "degree=15")

pie(c(sum(as.numeric(act.deg$degree_activating)>20 & act.deg$essentiality==1),
      sum(as.numeric(act.deg$degree_activating)>20 & act.deg$essentiality==0)),
    main = "degree=20")

dev.off()

pdf("./yeast_GEM9/data_for_github/figures/figure 3b iii degree distribution of activators.pdf",width=6, height=6);
boxplot(as.numeric(act.deg$degree_activating[act.deg$Activators %in% essential.met.kegg]),
        as.numeric(act.deg$degree_activating[!act.deg$Activators %in% essential.met.kegg]),
        outline = F,names = c("essential Act", "nonessential Act"),ylab = "degree of activation",
        col=c("#CC0000","#BCBCBC"),main = "degree distribution of activators")

legend("topleft",legend = paste("P value (essential vs nonessential = ",
                                paste("P=",format(t.test(as.numeric(act.deg$degree_activating[act.deg$Activators %in% essential.met.kegg]),
                                                         as.numeric(act.deg$degree_activating[!act.deg$Activators %in% essential.met.kegg]))$p.value,digits=3))),bty = "n")
dev.off()

#############################################################################
######## mapping dene expression data onto activation network 
#############################################################################
ec.details2 = ec.details[ec.details$genes != "",]

############ gene expression of activated enzymes
gene.exp.activated.ec <- vector(mode='list', length=dim(ec.details2)[1])
names(gene.exp.activated.ec) = ec.details2$EC
for(i in 1:dim(ec.details2)[1]) {
  my.genes = model_metabolic2$`GENE ASSOCIATION`[model_metabolic2$ec_new==gsub("EC_","",ec.details2$EC[i])]
  my.genes = gsub("\\(|\\)| ","",unlist(strsplit(my.genes," or ")))
  my.genes = my.genes[!is.na(my.genes)]
  if(sum(my.genes %in% rownames(exp.all2)) == 0) {
    gene.exp.activated.ec[[i]] = NA
  } else {
    my.exp = lapply(my.genes, function(x)
      (as.numeric(unlist(exp.all2[rownames(exp.all2) %in% unlist(strsplit(x,"and")),]))))
    my.exp = my.exp[!is.na(unlist(lapply(my.exp, mean)))]
    gene.exp.activated.ec[[i]] = my.exp[[which(unlist(lapply(my.exp, mean)) == max(unlist(lapply(my.exp, mean))))[1]]]
  }
  print(i)
}
gene.exp.activated.ec = gene.exp.activated.ec[!is.na(lapply(gene.exp.activated.ec, mean))]

############ gene expression of non-activated enzymes
gene.exp.nonactivated.ec <- vector(mode='list', length=length(non.activated.EC))
names(gene.exp.nonactivated.ec) = non.activated.EC
for(i in 1:length(non.activated.EC)) {
  if(is.na(non.activated.EC[i])) {
    gene.exp.nonactivated.ec[[i]] = NA
  } else {
    my.genes = model_metabolic2$`GENE ASSOCIATION`[model_metabolic2$ec_new==non.activated.EC[i]]
    my.genes = gsub("\\(|\\)| ","",unlist(strsplit(unlist(strsplit(my.genes," or ")),"and")))
    my.exp = lapply(my.genes, function(x)
      (as.numeric(unlist(exp.all2[rownames(exp.all2) %in% unlist(strsplit(x,"and")),]))))
    my.exp = my.exp[!is.na(unlist(lapply(my.exp, mean)))]
    gene.exp.nonactivated.ec[[i]] = my.exp[[which(unlist(lapply(my.exp, mean)) == max(unlist(lapply(my.exp, mean))))[1]]]
  }
  print(i)
}
gene.exp.nonactivated.ec = gene.exp.nonactivated.ec[!is.na(lapply(gene.exp.nonactivated.ec, mean))]

pdf("./yeast_GEM9/data_for_github/figures/figure 3c i gene expression variation and the degree of activation.pdf",width=12, height=6);
boxplot(gene.exp.activated.ec,outline=FALSE, col=ifelse(ec.details2$essentiality==1,"#CC0000","#BCBCBC"),
        main = "gene expression variation and the degree of activation",
        ylab="gene expression (log2)",las =2,ylim=c(0,30))
points(1:dim(ec.details2)[1],as.numeric(ec.details2$degree_activated),pch=19,
       col=ifelse(ec.details2$essentiality==1,"#CC0000","#BCBCBC"))
dev.off()

############ degree of activation for highly connected enzymes in the activation network
pdf("./yeast_GEM9/data_for_github/figures/figure 3c ii hub enzyme degree ess and noness.pdf",width=10, height=4);
par(mfrow=c(1,3))
pie(c(sum(as.numeric(ec.details$degree_activated)>12 & ec.details$essentiality==1),
      sum(as.numeric(ec.details$degree_activated)>12 & ec.details$essentiality==0)),
    main = "degree=12")
pie(c(sum(as.numeric(ec.details$degree_activated)>15 & ec.details$essentiality==1),
      sum(as.numeric(ec.details$degree_activated)>15 & ec.details$essentiality==0)),
    main = "degree=15")
pie(c(sum(as.numeric(ec.details$degree_activated)>20 & ec.details$essentiality==1),
      sum(as.numeric(ec.details$degree_activated)>20 & ec.details$essentiality==0)),
    main = "degree=20")
dev.off()

############ degree of enzyme activation for essential and non-essential enzymes
degree_activated.ess    = as.numeric(ec.details$degree_activated[gsub("EC_","",ec.details$EC) %in% essential.ec])
degree_activated.noness = as.numeric(ec.details$degree_activated[!gsub("EC_","",ec.details$EC) %in% essential.ec])
pdf("./yeast_GEM9/data_for_github/figures/figure 3c iii Enzyme activation degree boxplot.pdf",width=6, height=6);
boxplot(degree_activated.ess,
        degree_activated.noness,
        outline = F,names = c("essential enzyme", "nonessential enzyme"),
        ylab = "Number of time the enzyme activated",col = c("#CC0000","#BCBCBC"),
        main="Enzyme activation degree")
legend("topleft",legend = paste("P value (essential vs nonessential = ",
                                paste("P=",format(t.test(degree_activated.ess,degree_activated.noness)$p.value,digits=3))),bty = "n")
dev.off()

############ Gene expression variation for essential, non-essential activated and non-activated enzymes
pdf("./yeast_GEM9/data_for_github/figures/figure 3c iv gene expression variation of essential and nonessential activated enzymes.pdf",width=10, height=10);
boxplot(as.numeric(unlist(gene.exp.activated.ec[gsub("EC_","",names(gene.exp.activated.ec)) %in% essential.ec])),
        as.numeric(unlist(gene.exp.activated.ec[!gsub("EC_","",names(gene.exp.activated.ec)) %in% essential.ec])),
        as.numeric(unlist(gene.exp.nonactivated.ec[names(gene.exp.nonactivated.ec) %in% essential.ec])),
        as.numeric(unlist(gene.exp.nonactivated.ec[!names(gene.exp.nonactivated.ec) %in% essential.ec])),
        outline=FALSE, names = c("activated ess EC", "activated noness EC","nonactivated ess EC", "nonactivated noness EC"),
        col=c("#CC0000","#BCBCBC"),
        main = "gene expression variation",
        ylab = "log2 expression")
dev.off()


#############################################################################
######## Short Path Length of activators and non activators in Yeast9 model
#############################################################################
compartment = unique(model.met$COMPARTMENT)

############ creating metabolite-metabolite connection network from Yeast9
met.met.nw = c()
for(i in 1:length(model_metabolic$ID)) {
  my.rev = ifelse(grepl("<=>",model_metabolic$EQUATION[i]),TRUE,FALSE)
  my.subs.prod = unlist(ifelse(my.rev,strsplit(model_metabolic$EQUATION[i],"<=>"),
                               strsplit(model_metabolic$EQUATION[i],"=>")))
  
  if(my.rev) {
    my.subs = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod,"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
    my.prod = my.subs
  } else {
    my.subs = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod[1],"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
    my.prod = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod[2],"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
  }
  my.subs = my.subs[!my.subs==""]
  my.prod = my.prod[!my.prod==""]
  met.met.nw=rbind(met.met.nw,expand.grid(my.subs,my.prod))
}
met.met.nw = met.met.nw[!is.na(met.met.nw$Var2),]

############ creating graph for metabolite-metabolite connection network
g.met <- graph.data.frame(met.met.nw,directed = F)
g.met = simplify(g.met)
degree.g.met <- as.matrix(igraph::degree(g.met),ncol=1) #print the number of edges per vertex
degree.g.met <- as.data.frame(degree.g.met)
colnames(degree.g.met) <- "Degree"
degree.g.met = degree.g.met[order(degree.g.met$Degree,decreasing = T),,drop=F]

############ highly connected metabolites to be removed to avoid unnecessary connections
cofactors.rem = rownames(degree.g.met)[degree.g.met$Degree>50]
cofactors.rem = unique(tolower(unlist(lapply(cofactors.rem, function(x) unlist(strsplit(x,"\\["))[1]))))
cofactors.rem = unique(unlist(lapply(cofactors.rem, function(x) paste0(paste0(x,'[',compartment),']'))))
#### Final model after removing highly connected metabolites
met.met.nw.final = met.met.nw[!((tolower(met.met.nw[,1]) %in% cofactors.rem) | (tolower(met.met.nw[,2]) %in% cofactors.rem)),]

############ creating graph for metabolite-metabolite connection network after filtering highly connected metabolites
g.met.nw <- graph.data.frame(met.met.nw.final,directed = F)
g.met.nw = simplify(g.met.nw)
degree.g.met.nw <- as.matrix(igraph::degree(g.met.nw),ncol=1) #print the number of edges per vertex
degree.g.met.nw <- as.data.frame(degree.g.met.nw)
colnames(degree.g.met.nw) <- "Degree"
degree.g.met.nw = degree.g.met.nw[order(degree.g.met.nw$Degree,decreasing = T),,drop=F]
edge.list.met.met.nw <- as.data.frame(get.edgelist(g.met.nw))

############ shortest path length for each metabolite of Yeast9
shortest.paths =  shortest.paths(g.met.nw)
shortest.paths.glc = shortest.paths[which(rownames(shortest.paths) == "D-glucose[e]"),]
shortest.paths.glc = shortest.paths.glc[!is.infinite(shortest.paths.glc)]

############ shortest path length for activator and nonactivator from activation network
act.deg$ActShortPath = unlist(lapply(act.deg$Activators, function(x)
  min(shortest.paths.glc[names(shortest.paths.glc) %in% model.met$ID[model.met$KEGG ==x]])))

############ shortest path length for activator and nonactivator from Yeast9
model.met$ActShortPath = unlist(lapply(model.met$ID, function(x)
  min(shortest.paths.glc[names(shortest.paths.glc) == x])))

non.act.short.path = model.met$ActShortPath[!(model.met$ID %in% unique(model.met$ID[model.met$KEGG %in% act.deg$Activators]))]
non.act.short.path = non.act.short.path[!is.infinite(non.act.short.path)]

short.path.act = rbind(data.frame(variable = "act",shortPath=as.numeric(act.deg$ActShortPath[!is.infinite(act.deg$ActShortPath)])),
                       data.frame(variable = "nonAct",shortPath = as.numeric(non.act.short.path)))

p <- ggplot(short.path.act, aes(x=variable, y=shortPath,col=as.factor(variable))) + 
  geom_violin(width=1, adjust=1.4)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(title = paste0("P = ",format(t.test(short.path.act$shortPath[short.path.act$variable=="act"],
                                           short.path.act$shortPath[short.path.act$variable=="nonAct"])$p.value,digits=3)))
pdf("./yeast_GEM9/data_for_github/figures/figure 4a Short Path Length of activators and non activators.pdf",width=6, height=6);
p
dev.off()


############ density plot for short Path Length activators and nonactivators
den.act = density(as.numeric(act.deg$ActShortPath[!is.infinite(act.deg$ActShortPath)]),adjust = 2)
den.nonact = density(as.numeric(non.act.short.path),adjust = 2)

my.xlim = round(max(max(den.act$x),max(den.nonact$x)))
my.ylim = max(max(den.act$y),max(den.nonact$y))

pdf("./yeast_GEM9/data_for_github/figures/figure 4b density plot of Short Path Length activators and non activators.pdf",width=6, height=6);
plot(den.act,ylim=c(0,my.ylim),xlim=c(1,15),
     xlab = "Short Path Length", 
     main="Short Path Length activators and non activators")
polygon(den.act, col = rgb(1,0,0,0.2))
polygon(den.nonact, col = rgb(0,1,0,0.2))

legend("topleft",legend = c("activators","non activators"),
       fill=c(rgb(1,0,0,0.2),rgb(0,1,0,0.2)),bty = "n")

dev.off()

############ Short Path Length vs Degree of activation
pdf("./yeast_GEM9/data_for_github/figures/figure 4c Short Path Length vs Degree of activation for activators.pdf",width=6, height=6);
plot(as.numeric(act.deg$ActShortPath[!(is.infinite(act.deg$ActShortPath) | act.deg$ActShortPath==0)]),
     as.numeric(act.deg$degree_activating[!(is.infinite(act.deg$ActShortPath) | act.deg$ActShortPath==0)]),
     main = "activation Degree vs Short Path Length",
     xlab = "Short Path Length",ylab = "Degree")
legend("topleft",legend = paste0("P = ",
                                 format(cor.test(as.numeric(act.deg$ActShortPath[!(is.infinite(act.deg$ActShortPath) | act.deg$ActShortPath==0)]),
                                                 as.numeric(act.deg$degree_activating[!(is.infinite(act.deg$ActShortPath) | act.deg$ActShortPath==0)]))$p.value,digits=3),
                                 " & r = ",format(cor.test(as.numeric(act.deg$ActShortPath[!(is.infinite(act.deg$ActShortPath) | act.deg$ActShortPath==0)]),
                                                           as.numeric(act.deg$degree_activating[!(is.infinite(act.deg$ActShortPath) | act.deg$ActShortPath==0)]))$estimate,digits=3)),bty = "n")
dev.off()

#############################################################################
######## Short Path Length of activated and non-activated enzymes in Yeast9 model
#############################################################################
############ creating a network of reaction to reaction connection based on shared metabolites
rxn.met.net = data.frame()
for(i in 1:dim(model_metabolic)[1]) {
  my.rev = ifelse(grepl("<=>",model_metabolic$EQUATION[i]),TRUE,FALSE)
  
  my.subs.prod1 = tolower(unlist(ifelse(my.rev,strsplit(model_metabolic$EQUATION[i],"<=>"),
                                        strsplit(model_metabolic$EQUATION[i],"=>"))))
  if(my.rev) {
    my.subs1 = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod1,"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
    my.prod1 = my.subs1
  } else {
    my.subs1 = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod1[1],"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
    my.prod1 = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod1[2],"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
  }
  
  my.subs1 = my.subs1[!my.subs1==""]
  my.prod1 = my.prod1[!my.prod1==""]
  my.subs1 = my.subs1[!my.subs1 %in% cofactors.rem]
  my.prod1 = my.prod1[!my.prod1 %in% cofactors.rem]
  
  for(j in 1:dim(model_metabolic)[1]) {
    my.rev2 = ifelse(grepl("<=>",model_metabolic$EQUATION[j]),TRUE,FALSE)
    
    my.subs.prod2 = tolower(unlist(ifelse(my.rev2,strsplit(model_metabolic$EQUATION[j],"<=>"),
                                          strsplit(model_metabolic$EQUATION[j],"=>"))))
    if(my.rev2) {
      my.subs2 = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod2,"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
      my.prod2 = my.subs2
    } else {
      my.subs2 = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod2[1],"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
      my.prod2 = unique(trimws(unlist(strsplit(unlist(strsplit(my.subs.prod2[2],"\\s\\+\\s")),"^\\d*\\.*\\d*\\s|^\\d*\\.*\\d*e-\\d*"))))
    }
    my.subs2 = my.subs2[!my.subs2==""]
    my.prod2 = my.prod2[!my.prod2==""]
    my.subs2 = my.subs2[!my.subs2 %in% cofactors.rem]
    my.prod2 = my.prod2[!my.prod2 %in% cofactors.rem]
    
    if(sum(my.subs2 %in% my.prod1)>0 | sum(my.subs1 %in% my.prod2)>0) {
      rxn.met.net = rbind(rxn.met.net,cbind(model_metabolic$ID[i],model_metabolic$ID[j]))
    }
  }
  print(i)
}

############ remove duplicate connections
rxn.met.net2 = rxn.met.net[!duplicated(rxn.met.net), ]
rxn.met.net2 = rxn.met.net2[!(rxn.met.net2$V1 == rxn.met.net2$V2),]

############ after reactions network, replace reactions with EC number
rxn.met.net2$V3 = unlist(lapply(rxn.met.net2$V1, function(x)
  unique(ifelse(model_metabolic$`EC-NUMBER`[model_metabolic$ID==x]=="",x,
                model_metabolic$`EC-NUMBER`[model_metabolic$ID==x]))))

rxn.met.net2$V4 = unlist(lapply(rxn.met.net2$V2, function(x)
  unique(ifelse(model_metabolic$`EC-NUMBER`[model_metabolic$ID==x]=="",x,
                model_metabolic$`EC-NUMBER`[model_metabolic$ID==x]))))

rxn.met.net2$V3[is.na(rxn.met.net2$V3)] = rxn.met.net2$V1[is.na(rxn.met.net2$V3)]
rxn.met.net2$V4[is.na(rxn.met.net2$V4)] = rxn.met.net2$V2[is.na(rxn.met.net2$V4)]

############ expanding the network by EC list per reactions
#### split the EC list into an individual enzymes.
rxn.met.net3 = c()
for(i in 1:dim(rxn.met.net2)[1]) {
  print (i)
  rxn.met.net3=rbind(rxn.met.net3,
                     cbind(rxn.met.net2$V1[i],rxn.met.net2$V2[i],expand.grid(unlist(strsplit(rxn.met.net2$V3[i],"\\;")),
                                                                             unlist(strsplit(rxn.met.net2$V4[i],"\\;")))))
}

############ remove duplicate connections again
rxn.met.net.ec = rxn.met.net3[!duplicated(rxn.met.net3), ]

############ create graph and the calculated shortest path length
g.metRxn.nw <- graph.data.frame(rxn.met.net.ec[,3:4],directed = F)
shortest.paths.rxn =  shortest.paths(g.metRxn.nw)

############ shortest path length for each enzymes from glucose transport reaction "r_1166"
shortest.paths.ex.glc = shortest.paths.rxn[which(rownames(shortest.paths.rxn) == "r_1166"),]
shortest.paths.ex.glc = shortest.paths.ex.glc[!is.infinite(shortest.paths.ex.glc)]
shortest.paths.ex.glc.ec = shortest.paths.ex.glc[-grep("r_",names(shortest.paths.ex.glc))]

############ shortest path length for activated enzymes
ec.details$rxnShortPath = unlist(lapply(gsub("EC_","",ec.details$EC), function(x)
  ifelse(length(shortest.paths.ex.glc.ec[names(shortest.paths.ex.glc.ec) == x])>0,
         shortest.paths.ex.glc.ec[names(shortest.paths.ex.glc.ec) == x],NA)))

############ shortest path length for non-activated enzymes
non.act.short.path.rxn = unlist(lapply(non.activated.EC, function(x)
  ifelse(length(shortest.paths.ex.glc.ec[names(shortest.paths.ex.glc.ec) == x])>0,
         shortest.paths.ex.glc.ec[names(shortest.paths.ex.glc.ec) == x],NA)))
non.act.short.path.rxn = non.act.short.path.rxn[!is.infinite(non.act.short.path.rxn)]
non.act.short.path.rxn = non.act.short.path.rxn[!is.na(non.act.short.path.rxn)]

############ shortest path length for enzymes combined (activated and nonactivated)
C1 <- data.frame(ShortPath = as.numeric(ec.details$rxnShortPath[!is.infinite(as.numeric(ec.details$rxnShortPath))]), variable = "EC_activated")
C1 = C1[!is.na(C1$ShortPath),]
C2 <- data.frame(ShortPath = as.numeric(non.act.short.path.rxn), variable = "EC_nonactivated")
short.path.rxn <- rbind(C1, C2)

p <- ggplot(short.path.rxn, aes(x=variable, y=ShortPath,col=as.factor(variable))) + 
  geom_violin(width=1, adjust=1.6)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(title = paste0("P = ",format(t.test(as.numeric(ec.details$rxnShortPath[!is.infinite(as.numeric(ec.details$rxnShortPath))]),
                                           as.numeric(non.act.short.path.rxn))$p.value,digits=3)))
pdf("./yeast_GEM9/data_for_github/figures/figure 4d Short Path Length of activated and non activated enzymes.pdf",width=6, height=6);
p
dev.off()

############ shortest path length of activated f and non-activated enzyme
#### density plot for short Path Length activators and nonactivators
den.act.ec = density(as.numeric(C1$ShortPath),adjust = 2)
den.nonact.ec = density(as.numeric(non.act.short.path.rxn),adjust = 2)

my.xlim = round(max(max(den.act.ec$x),max(den.nonact.ec$x)))
my.ylim = max(max(den.act.ec$y),max(den.nonact.ec$y))

pdf("./yeast_GEM9/data_for_github/figures/figure 4e density plot of Short Path Length Enzymes.pdf",width=6, height=6);
plot(den.act.ec,ylim=c(0,my.ylim),xlim=c(1,my.xlim),
     xlab = "Short Path Length", 
     main="Short Path Length activated and nonactivated enzymes")
polygon(den.act.ec, col = rgb(1,0,0,0.2))
polygon(den.nonact.ec, col = rgb(0,1,0,0.2))

legend("topleft",legend = c("activated EC","non activated EC"),
       fill=c(rgb(1,0,0,0.2),rgb(0,1,0,0.2)),bty = "n")

dev.off()

############ degree of enzyme activation vs shortest path length from glucose uptake
pdf("./yeast_GEM9/data_for_github/figures/figure 4f Short Path Length vs Degree of activation for enzymes.pdf",width=6, height=6);
plot(as.numeric(ec.details$rxnShortPath[!is.na(ec.details$rxnShortPath)]),
     as.numeric(ec.details$degree_activated[!is.na(ec.details$rxnShortPath)]),
     main = "activation Degree vs Short Path Length EC",
     xlab = "Short Path Length",ylab = "Degree")

legend("topleft",legend = paste0("P = ",
                                 format(cor.test(as.numeric(ec.details$rxnShortPath[!is.infinite(ec.details$rxnShortPath)]),
                                                 as.numeric(ec.details$degree_activated[!is.infinite(ec.details$rxnShortPath)]))$p.value,digits=3),
                                 " & r = ",format(cor.test(as.numeric(ec.details$rxnShortPath[!is.infinite(ec.details$rxnShortPath)]),
                                                           as.numeric(ec.details$degree_activated[!is.infinite(ec.details$rxnShortPath)]))$estimate,digits=3)),bty = "n")
dev.off()

########################################################################
############ metabolite-enzyme activation + Gene co-expression network
########################################################################
############ Loading gene co-expression 
exp.cor.p = exp.cor$P
exp.cor.r = exp.cor$r
exp.cor.p[is.na(exp.cor.p)] = 1
exp.cor.padj = matrix(p.adjust(exp.cor.p,method = "BH"),nrow=nrow(exp.cor.p),ncol=ncol(exp.cor.p))
rownames(exp.cor.padj) = rownames(exp.cor.p)
colnames(exp.cor.padj) = colnames(exp.cor.p)

############ create gene co-expression network
gene.coexp.net = matrix(0, nrow=nrow(exp.cor.padj),ncol=ncol(exp.cor.padj))
rownames(gene.coexp.net) = rownames(exp.cor.padj)
colnames(gene.coexp.net) = colnames(exp.cor.padj)
gene.coexp.net[exp.cor.r>0.75 & exp.cor.padj<0.001]=1
gene.coexp.net = gene.coexp.net[apply(gene.coexp.net, 2, sum)>0,apply(gene.coexp.net, 2, sum)>0]
gene.coexp.net = as.matrix(gene.coexp.net)

g.gene.coexp.net <- graph_from_adjacency_matrix(gene.coexp.net,mode = "directed")
nodes.gene.coexp.net  <- rownames(as.matrix(V(g.gene.coexp.net)))
edge.listcoexp.net <- as.data.frame(get.edgelist(g.gene.coexp.net))
degree.coexp.net <- as.matrix(igraph::degree(g.gene.coexp.net),ncol=1) #print the number of edges per vertex
degree.coexp.net <- as.data.frame(degree.coexp.net)
colnames(degree.coexp.net) <- "Degree"
V(g.gene.coexp.net)$label.cex = 0.5

############ Activation network with associated genes of nodes co-expressed.
#### from the activation network, take an EC and find all genes, which is in ec.details object
#### Similarly take an activator and find enzyme which produces it in the Yeast9.  Then get the genes
#### check if any of genes associated with node1 [EC] is co-expressed with genes associated with node2[activator]
act.net2 = act.net
co.exp = c()
act.genetic.int = c()

for(i in 1:dim(act.net2)[1]) {
  node1 = act.net2$EC[i]
  node1.genes = unlist(strsplit(ec.details$genes[ec.details$EC == node1],"\\|"))
  
  node2 = act.net2$Activators[i]
  node2.ec = ec.details$EC[unlist(lapply(ec.details$prod, function(x)
    ifelse(sum(unlist(strsplit(x,"\\|")) == node2)==0,FALSE,TRUE)))]
  node2.genes = unlist(strsplit(ec.details$genes[ec.details$EC %in% node2.ec],"\\|"))
  
  if(is.null(node1.genes)|is.null(node2.genes)) {
    co.exp = rbind(co.exp,0)
    act.genetic.int = rbind(act.genetic.int,0)
  } else {
    node1.node2.gene = expand.grid(node1.genes,node2.genes)
    node1.node2.gene.coexp = (node1.node2.gene[,1] %in% edge.listcoexp.net$V1) & node1.node2.gene[,2] %in% edge.listcoexp.net$V2|
      node1.node2.gene[,1] %in% edge.listcoexp.net$V2 & node1.node2.gene[,2] %in% edge.listcoexp.net$V1
    co.exp = rbind(co.exp,ifelse(sum(node1.node2.gene.coexp)==0,0,1))
    
    node1.node2.gene.int = (node1.node2.gene[,1] %in% genetic.int.costanzo$V1) & node1.node2.gene[,2] %in% genetic.int.costanzo$V3|
      node1.node2.gene[,1] %in% genetic.int.costanzo$V3 & node1.node2.gene[,2] %in% genetic.int.costanzo$V1
    
    act.genetic.int = rbind(act.genetic.int,ifelse(sum(node1.node2.gene.int)==0,0,sum(node1.node2.gene.int)))
  }
  print(i)
}

act.net2$coExpressed = as.numeric(co.exp)
act.net2$act_genetic_int = as.numeric(act.genetic.int)
sum(act.net2$act_genetic_int>0)/dim(act.net2)[1]

########### remove common co-factors
rem.met = c("C00009","C00020","C00008","C00002")
act.net2 = act.net2[!act.net2$Activators %in% rem.met,] # Remove ATP, ADP, AMP, PI from the network
act.net2$RxnName = unlist(lapply(act.net2$EC, function(x)
  model_metabolic2$RxnName[model_metabolic2$ec_new==gsub("EC_","",x)][1]))
act.net2$MetName = unlist(lapply(act.net2$Activators, function(x)
  model.met$Name[model.met$KEGG==x][1]))

########### Creating a graph of gene co-expression network associated with activation network
g.coExp.nw <- graph.data.frame(act.net2[act.net2$coExpressed==1,1:2], directed=T)
V(g.coExp.nw)$label.cex = 0.5

nodes.coExp.nw  <- rownames(as.matrix(V(g.coExp.nw)))
edge.list.coExp.nw <- as.data.frame(get.edgelist(g.coExp.nw))
degree.coExp.nw <- as.matrix(igraph::degree(g.coExp.nw),ncol=1) #print the number of edges per vertex
degree.coExp.nw <- as.data.frame(degree.coExp.nw)
colnames(degree.coExp.nw) <- "Degree"
E(g.coExp.nw)$color = "#29C66A"

pdf("./yeast_GEM9/data_for_github/figures/figure 5a activation network and coexpressed.pdf",width=8, height=8);
plot(g.coExp.nw,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
     mark.border="#BCBCBC",pch=19,vertex.label=NA,
     vertex.shape = ifelse(names(V(g.coExp.nw)) %in% act.net$EC,"square","circle"),
     vertex.color = ifelse(names(V(g.coExp.nw)) %in% paste0("EC_",essential.ec), "#CC0000",
                           ifelse(names(V(g.coExp.nw)) %in% essential.met.kegg,"#CC0000","#BCBCBC")),
     layout = layout_with_fr,
     main="activation and coexpressed")
dev.off()

pdf("./yeast_GEM9/data_for_github/figures/figure 5b Activation betwork edge coexp, genetic interaction.pdf",width=8, height=8);
barplot(t(rbind(cbind(sum(act.net2$coExpressed==1)/dim(act.net2)[1]*100, 
                      sum(act.net2$coExpressed==0)/dim(act.net2)[1]*100),
                cbind(sum(act.net2$act_genetic_int>0)/dim(act.net2)[1]*100, 
                      sum(act.net2$act_genetic_int==0)/dim(act.net2)[1]*100))),
        col = c("#29C66A","#BCBCBC"),ylab = "% representation of edge coexp, genetic int",
        main = paste("both coexp, genetic.int % = ",round(sum(act.net2$coExpressed==1 & act.net2$act_genetic_int>1)/dim(act.net2)[1]*100),2))
legend("topleft",legend = c("coExpressed","not coExpressed"),
       fill = c("#29C66A","#BCBCBC"),bty = "n")
dev.off()

##############################################################################
#################################### END #####################################
##############################################################################

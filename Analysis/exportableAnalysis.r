
## -------------------------------------------------------------------------------------------------
#install.packages("tidyverse")
#install.packages("PVR")
#install.packages("ape")
#install.packages("phytools")
#install.packages("igraph")
#install.packages("magrittr")


library(dplyr)
library(ggplot2)
library(PVR)
library(ape)
library(phytools)
library(igraph)
library(magrittr)



## -------------------------------------------------------------------------------------------------
galltree <- ape::read.nexus(file="../Data.nosync/Nyman_Phylogeny.nex")
galltree <- phytools::force.ultrametric(galltree)

litter <- PVR::PVRdecomp(phy=galltree, type="nexus")

EigVec <- data.frame(litter@Eigen$vectors)
EigVal <- data.frame(litter@Eigen$values)
EigVal$litter.Eigen.values <- EigVal$litter.Eigen.values/(sum(EigVal$litter.Eigen.values))

phytools::plotTree.wBars(tree=galltree, setNames(EigVec$c1+abs(min(EigVec$c1)), galltree$tip.label))
phytools::plotTree.wBars(tree=galltree, setNames(EigVec$c2+abs(min(EigVec$c2)), galltree$tip.label))


## -------------------------------------------------------------------------------------------------
output <- NULL
nreps <- 1000

for(i in 1:nreps){
  #Choose a random set of nodes to remove
  numDrops <- sample(1:(length(galltree$tip.label)-round(length(galltree$tip.label)/10)), size=1)
  
  Intdrops <- sample(1:length(galltree$tip.label), size=numDrops, replace = FALSE)
  
  #Get the descendants of those nodes and remove them
  randPruned <- ape::drop.tip(galltree, galltree$tip.label[Intdrops])
  randGraph <- igraph::as.igraph(randPruned)
  
  #Run a phylodecomp on that galltree
  randDecomp <- PVR::PVRdecomp(phy=randPruned, type="nexus")
  
  randDF <- data.frame("rep"= i, 
             "vectorID"=names(((randDecomp@Eigen$values)/sum(randDecomp@Eigen$values))[1:10]), 
             "PropVar"=((randDecomp@Eigen$values)/sum(randDecomp@Eigen$values))[1:10],
             "Ntips"=length(randPruned$tip.label), 
             "infomapModularity"=igraph::infomap.community(randGraph)$modularity, 
             "tipsDiscarded"=toString(sort(Intdrops)),
             row.names = NULL
             )
  
  #Save outputs
  output <- rbind(output, randDF)
}


## -------------------------------------------------------------------------------------------------
sd(output$PropVar, na.rm=TRUE)

output %>% dplyr::filter(., vectorID=="c1") %>%
  ggplot(data=., aes(x=Ntips, y=PropVar))+geom_point()+
  xlab("Number of Species in Tree")+
  ylab("Proportion of Var Explained by 1st Eigenvector")

output %>% dplyr::filter(., vectorID=="c2") %>%
  ggplot(data=., aes(x=Ntips, y=PropVar))+geom_point()+
  xlab("Number of Species in Tree")+
  ylab("Proportion of Var Explained by 2nd Eigenvector")

output %>% group_by(., rep) %>% mutate(eigenSD=sd(PropVar)) %>%
  dplyr::select(., Ntips, eigenSD) %>% unique() %>%
  ggplot(data=., aes(x=Ntips, y=eigenSD))+
  geom_point()+
  xlab("Number of Species in Tree")+
  ylab("Standard Deviaiton of First 10 Eigenvalues' Proportions")+xlim(10, 75)


## -------------------------------------------------------------------------------------------------
output %>% dplyr::filter(., vectorID=="c1") %>%
  ggplot(data=., aes(x=infomapModularity, y=PropVar))+geom_point()+
  xlab("Modularity")+
  ylab("Proportion of Var Explained by 1st Eigenvector")

output %>% dplyr::filter(., vectorID=="c1") %>% dplyr::filter(., infomapModularity > 0) %>%
  ggplot(data=., aes(x=infomapModularity, y=PropVar))+geom_point()+
  xlab("Non-Zero Modularity")+
  ylab("Proportion of Var Explained by 1st Eigenvector")


## -------------------------------------------------------------------------------------------------
sharks_full <- ape::read.nexus(file="../Data.nosync/Stein2018SharkTree/output.nex")
shark_tree <- sharks_full[[1]]

litter <- PVR::PVRdecomp(phy=shark_tree, type="nexus")

EigVec <- data.frame(litter@Eigen$vectors)
EigVal <- data.frame(litter@Eigen$values)
EigVal$litter.Eigen.values <- EigVal$litter.Eigen.values/(sum(EigVal$litter.Eigen.values))

phytools::plotTree.wBars(tree=shark_tree, setNames(EigVec$c1+abs(min(EigVec$c1)), shark_tree$tip.label)) #This plot looks terrible but it validates our method


## -------------------------------------------------------------------------------------------------
mammals <- phytools::read.newick(file="../Data.nosync/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")

mammals$node.label <- NULL #Igraph doesn't like the node labels given is this default tree, so we just remove them


## -------------------------------------------------------------------------------------------------
spagTree <-  ape::read.nexus(file="../Data.nosync/Piatkowski&Shaw_2019_MCC_Tree.nex")


## -------------------------------------------------------------------------------------------------

analyze <- function(treesList, treeNames, nreps){
  output <- NULL
  ntrees <- length(treesList)
  for(i in 1:nreps){
      for(j in 1:ntrees){
        tree <- treesList[[j]]
        #Choose a random set of nodes to remove
        numDrops <- sample(1:(length(tree$tip.label)-10), size=1)
        Intdrops <- sample(1:length(tree$tip.label), size=numDrops, replace = FALSE)
        
        #Get the descendants of those nodes and remove them
        randPruned <- ape::drop.tip(tree, tree$tip.label[Intdrops])
        
        randGraph <- igraph::as.igraph(randPruned)
        
        #Run a phylodecomp on that tree
        randDecomp <- PVR::PVRdecomp(phy=randPruned, type="nexus")
        
        randDF <- data.frame("rep"= i, 
                   "vectorID"=names(((randDecomp@Eigen$values)/sum(randDecomp@Eigen$values))[1:10]), 
                   "PropVar"=((randDecomp@Eigen$values)/sum(randDecomp@Eigen$values))[1:10],
                   "Ntips"=length(randPruned$tip.label), 
                   "infomapModularity"=igraph::infomap.community(randGraph)$modularity, 
                   "tipsDiscarded"=toString(sort(Intdrops)),
                   "tree"=treeNames[j], 
                   row.names = NULL
                   )
        #Save outputs
      output <- rbind(output, randDF)
      }
  }
  return(output)
}


## -------------------------------------------------------------------------------------------------
#FastResults <- analyze(treesList = list(shark_tree, spagTree, galltree), treeNames = c("sharks", "spag", "gallers"), nreps=1000)
#save(FastResults, file="FastResults.RData")

slowResults <- analyze(treesList = list(mammals), treeNames = c("mammals"), nreps=1000)
save(slowResults, file="slowResults.RData")


## -------------------------------------------------------------------------------------------------
temp <- analyze(treesList = list(shark_tree, spagTree, galltree, mammals), treeNames = c("sharks", "spag", "gallers", "mammals"), nreps=1)

temp <- analyze(treesList = list(shark_tree), treeNames = c("sharks"), nreps=1)
temp <- analyze(treesList = list(spagTree), treeNames = c("spagTree"), nreps=1)
temp <- analyze(treesList = list(galltree), treeNames = c("galltree"), nreps=1)
temp <- analyze(treesList = list(mammals), treeNames = c("mammals"), nreps=1)


## -------------------------------------------------------------------------------------------------
temp <- dplyr::filter(output, rep==785)

cut <- as.numeric(strsplit(temp$tipsDiscarded[1], split=", ")[[1]])
randPruned <- ape::drop.tip(tree, tree$tip.label[cut])
plot(randPruned)


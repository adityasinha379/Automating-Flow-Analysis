# Loading required packages ####
require(ggplot2)
require(reshape2)
require(gridExtra)
require(flowCore)
require(Rtsne)
require(cluster)
require(gplots)
require(mixtools)
require(pheatmap)

# Read and explore data ####
d <- read.table(file="Data/ccRCC_652T_Tcells.txt", header=TRUE,sep="\t")[8:20]
temp<- colnames(d)
colnames(d) <- gsub(".*\\.","",temp)
remove(temp)

# pick out required columns, leave out CD45 and LD (already gated)
d_comp <- d[,-which(names(d) %in% c("CD45","LD"))]
m<-nrow(d_comp)
n<-ncol(d_comp)
remove(d)

# for plotting distributions of markers
p <- vector("list", length = n)
for (i in c(1:n))
{p[[i]] <- local({
  i <- i
  p1 <- ggplot(d_comp, aes(x = d_comp[[i]])) + geom_histogram() + xlab(colnames(d_comp)[i])
})}
grid.arrange(grobs = p, ncol = 4)
remove(p)

# Transformation ####
trans <- vector("list",n)
w <- rep(0.8,11)
for (i in 1:n){
  trans[[i]] <- local({
    i<-i
    trans1 <- logicleTransform(w = w[i], t = max(d_comp), m = 4.5, a = 0)
  })
}
d_trans <- as.data.frame(do.call(cbind, lapply(1:ncol(d_comp), function(i) trans[[i]](d_comp[, i]))))
colnames(d_trans)<-colnames(d_comp)

#check transform using violin plot
temp <- melt(d_trans,id.vars=NULL)
colnames(temp)[1] <- "marker"
ggplot(temp, aes(x=marker, y=value)) + geom_violin(fill="gray", col="gray") + guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(min(d_trans),max(d_trans))
remove(temp)

ggplot(d_trans, aes(x=TCRab,y=CD49a)) + geom_point(color='coral', size=0.1) + ggtitle("Logicle")


# Clustering and  Visualization ####
# tsne embedding
d_tsne <- Rtsne(d_trans,dims=2)
d_tsne <- as.data.frame(d_tsne[['Y']])
colnames(d_tsne)[1] <- 'bhSNE1'
colnames(d_tsne)[2] <- 'bhSNE2'

write.csv(d_trans, "temp.csv", row.names = FALSE)
# louvain clustering in python jupyter notebook Clustering.ipynb
out <- read.csv("temp.csv")$cluster

c <- max(out)
d_tsne$cluster <- factor(out)

ggplot(d_tsne, aes(x=bhSNE1, y=bhSNE2, col=cluster)) + geom_point(size=0.1)
ggplot(d_tsne, aes(x=bhSNE1, y=bhSNE2)) +
  geom_point(aes(col = d_trans$CD16), size=0.1) + scale_color_gradientn(colors=bluered(16))

# cluster percentages
cluster_perc <- as.data.frame(table(out)*100/m)
colnames(cluster_perc) <- c("cluster","percent")
ggplot(data=cluster_perc, aes(x=cluster, y=percent)) + geom_bar(stat="identity")

# Violin plots ####
# marker distributions for each cluster
p <- vector("list", length = c)
titles <- c(1:c)
for(i in c(1:c)){
  temp <- d_trans[which(out == i),]
  temp <- melt(temp,id.vars=NULL)
  colnames(temp)[1] <- "marker"
  temp2 <- melt(d_trans,id.vars=NULL)
  colnames(temp2)[1] <- "marker"
  p[[i]] <- local({
    i <- i
    p1 <- ggplot(temp2, aes(x=marker, y=value)) + geom_violin(fill="gray", col="gray") +
      geom_violin(data=temp, aes(x=marker, y=value, fill=marker)) +
      ggtitle(titles[i]) + guides(fill=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(min(d_trans),max(d_trans))
  })
}
#png(filename="Results/Lineage_Tumor_635_comp/1.2.2.violin-5.png",width=1500,height=898)
grid.arrange(grobs = p, ncol = floor(c/3)+1)
#dev.off()
remove(p,temp,temp2,titles)

# Subpopulation granularities ####
# explore distributions within clusters (populations)
out <- out7
c <- max(out)
select <- c(6)
idx <- which(out %in% select)
sub <- d_trans[idx,]

ggplot(d_tsne, aes(x=bhSNE1, y=bhSNE2)) + geom_point(aes(col = d_trans$CD16)) + scale_color_gradientn(colors=bluered(16))

# to see the distribtutions in the tsne for markers in clusters
for(i in c(1:c)){
  select <- c(i)
  idx <- which(out %in% select)
  sub <- d_trans[idx,]
  for(j in c(1:n)){
    temp <- paste(paste("Results/Lineage_Tumor_638/markerplots/",(i-1)*n+j,sep=""),".png",sep="")
    png(filename=temp,width=1500,height=898)
    print(ggplot() + geom_point(data=d_tsne, aes(y=bhSNE2, x=bhSNE1), col="gray", size=0.8) +
      geom_point(data=d_tsne[idx,], aes(y=bhSNE2, x=bhSNE1, col=sub[[j]]), size=0.8) +
      scale_color_gradientn(colors=bluered(16), name=paste("sub$",colnames(sub)[j],sep="")) +
      ggtitle(paste("Clusters",paste(select, collapse = ","))))
    dev.off()
  }
}

ggplot()+
  geom_point(data=d_trans,aes(x=CD103,y=CD49a),col='coral',size=0.6,alpha=0.3)+
  geom_point(data=sub,aes(x=CD103,y=CD49a),col='black',size=0.6,alpha=0.3)+
  xlim(-0.5,4.5)+ylim(-0.5,4.5)
# GMM model for subdividing populations (unused)
temp<-normalmixEM(sub$GzmA,k=2)
temp2<-normalmixEM(sub$GzmB,k=2)
length(which(sub$GzmA<2.2 & sub$GzmB<2.2))*100/m
ggplot(sub, aes(x=GzmA)) + geom_histogram(binwidth=0.05)
ggplot(sub, aes(x=GzmB)) + geom_histogram(binwidth=0.05)

remove(p,temp,temp2,sub,select,idx)

# Subclustering populations ####
out <- out7
# selected cluster (parent population)
select <- c(6)
idx <- which(out %in% select)
sub <- d_trans[idx,]
write.csv(sub,"temp.csv",row.names = FALSE)
# perform subclustering of selected cluster using Clustering.ipynb
out <- read.csv("temp.csv")$cluster

# tsne and violin plots for visualization
ggplot(d_tsne[idx,], aes(x=bhSNE1, y=bhSNE2, col=factor(out))) + geom_point(size=0.8)

c <- max(out)
p <- vector("list", length = c)
titles <- c(1:c)
for(i in c(1:c)){
  temp <- sub[which(out == i),]
  temp <- melt(temp,id.vars=NULL)
  colnames(temp)[1] <- "marker"
  temp2 <- melt(sub,id.vars=NULL)
  colnames(temp2)[1] <- "marker"
  p[[i]] <- local({
    i <- i
    p1 <- ggplot(temp2, aes(x=marker, y=value)) + geom_violin(fill="gray", col="gray") +
      geom_violin(data=temp, aes(x=marker, y=value, fill=marker)) +
      ggtitle(titles[i]) + guides(fill=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(-1.5,4.5)
  })
}
grid.arrange(grobs = p, ncol = 2)
remove(p,temp,temp2,titles,sub,select,idx)


# Cluster Memberships in Scatter Plot ####
select <- c(5)
idx1 <- which(d_tsne$cluster %in% select)
idx2 <- c(1:m)[!c(1:m) %in% idx1]
temp <- d_trans
temp$select <- as.character(d_tsne$cluster)
temp$select[idx1] <- "1"
temp$select[idx2] <- "0"
temp$select <- as.factor(temp$select)
ggplot(d_tsne, aes(y=bhSNE2, x=bhSNE1, col=temp$select)) + geom_point(size=0.8) + 
  scale_colour_manual(values=c("gray", "magenta")) + guides(col=FALSE) + 
  ggtitle(paste("Clusters",paste(select, collapse = ",")))
ggplot(temp, aes(y=CD56, x=CD103, col=select)) + geom_point(size=0.1) + 
  scale_colour_manual(values=c("gray", "magenta")) + guides(col=FALSE) + 
  ggtitle(paste("Clusters",paste(select, collapse = ",")))
remove(temp, select, idx1, idx2)
# Gating vs Clustering correlation (unused) ####
pop <- read.csv("temp.csv")
pop$Gating<-round(pop$Gating,2)
pop <- pop[70:81,]
pop[,1] <- pop[,1]*100/35.95
pop[,2] <- pop[,2]*100/38.8
ggplot(pop,aes(x=Gating,y=Flow)) + geom_point(col='orange')+
  geom_smooth(method="lm",alpha=0) + ggtitle('Gating vs Flow: CD8')+
  annotate("text",x=80,y=20,label=paste("Pearson: ",round(cor(pop, method='pearson')[2],3)))
remove(pop)

# Agglomerative Clustering Routine (unused) ####
temp <- aggregate(d_trans, list(out), mean)
temp <- agnes(temp, method = "ward")
pltree(temp, cex = 0.6, hang = -1, main = "Amalgamation of Clusters using AGNES")
for(i in c(1:5)){
  print(unlist(cut(as.dendrogram(temp),h=3.93)$lower[[i]]))
}
remove(temp)
# Cluster Mean Heatmap (unused) ####
clustermean <- matrix(nrow = n, ncol = c, dimnames=list(colnames(d_trans),c(1:c)))
for (i in c(1:c)){
  idx1 <- which(d_tsne$cluster == i)
  for (j in c(1:n)){
    clustermean[j,i] <- mean(d_trans[idx1,j])
  }
}
temp <- melt(clustermean)
colnames(temp) <- c("gene","cluster","mean")
ggplot(data = temp, aes(x=cluster, y=gene, fill=mean)) + geom_tile() + scale_x_discrete(limits=c(1:c)) + scale_fill_gradientn(colours=bluered(16)) + ggtitle("Cluster means")
pheatmap(clustermean,col=bluered(16))

temp <- d_trans[order(d_tsne$cluster),]
rownames(temp) <- NULL
temp <- as.matrix(temp)
temp <- melt(temp)
colnames(temp) <- c("cell","gene","value")
cluster_freqs <- vector("integer", length = c)
cluster_freqs[2:c] <- cumsum(table(d_tsne$cluster))[1:c-1]
cluster_freqs <- cluster_freqs + 1
ggplot(data = temp, aes(x=cell, y=gene, fill=value)) + geom_tile() + scale_fill_gradientn(colours=bluered(16)) + ggtitle("Fluorescence by cluster") +
  scale_x_continuous(breaks=c(cluster_freqs),labels=c(1:c)) + xlab("cluster") + theme(axis.text.x = element_text(face="bold"))
remove(temp)

# Boxplots for gene expression (unused) ####
idx1 <- which(d_tsne$cluster %in% c(2,4,9,10,12,14,16,15,18,19,22)) # GzmB+
idx2 <- c(1:m)[!c(1:m) %in% idx1]

p <- vector("list", length = 2)

temp <- d_trans[4]
temp$expression <- as.character(d_tsne$cluster)
temp$expression[idx1] <- "High"
temp$expression[idx2] <- "Low"
temp$expression <- as.factor(temp$expression)
temp <- melt(temp, id.vars =c("expression"))
colnames(temp)[2] <- "marker"
p[[1]] <- ggplot(data = temp, aes(x=marker, y=value)) + geom_boxplot(aes(fill=expression)) + ggtitle("Gene expression by cluster")

ggplot(d_tsne, aes(x=bhSNE1, y=bhSNE2, col=temp$expression)) +
  geom_point(size=1) + scale_x_continuous(limits = c(-40, 40)) + 
  scale_y_continuous(limits = c(-40, 40)) + labs(col = "expression") + ggtitle("GzmB expression")

idx1 <- which(d_tsne$cluster %in% c(8,12,15,16,18,23)) # CD8a+
idx2 <- c(1:m)[!c(1:m) %in% idx1]

temp <- d_trans[5]
temp$expression <- as.character(d_tsne$cluster)
temp$expression[idx1] <- "High"
temp$expression[idx2] <- "Low"
temp$expression <- as.factor(temp$expression)
temp <- melt(temp, id.vars =c("expression"))
colnames(temp)[2] <- "marker"
p[[2]] <- ggplot(data = temp, aes(x=marker, y=value)) + geom_boxplot(aes(fill=expression)) + ggtitle("Gene expression by cluster")
ggplot(d_tsne, aes(x=bhSNE1, y=bhSNE2, col=temp$expression)) +
  geom_point(size=1) + scale_x_continuous(limits = c(-40, 40)) + 
  scale_y_continuous(limits = c(-40, 40)) + labs(col = "expression") + ggtitle("CD8a expression")

grid.arrange(grobs = p, ncol = 2)
remove(p,temp)
# Phenograph stability analysis (unused) ####
# based on a similarity matrix between cluster labelings
get_simmatrix <- function(d,imax,k_arr){
  require(pracma)
  sim <- matrix(nrow = length(k_arr), ncol = imax, dimnames=list(k_arr,c(1:imax)))
  m <- nrow(d)
  
  sink("/dev/null")
  for (j in c(1:length(k_arr))){
    for (i in c(1:imax)){
      require(Rphenograph)
      samp1 <- sample(m,as.integer(0.8*m))
      samp1 <- samp1[order(samp1)]
      samp2 <- sample(m,as.integer(0.8*m))
      samp2 <- samp2[order(samp2)]
      common <- intersect(samp1,samp2)
      s <- length(common)
      
      d_n <- d[samp1,]
      lab1 <- as.matrix(membership(Rphenograph(as.matrix(d_n), k = k_arr[j])[[2]]))
      rownames(lab1) <- samp1
      lab1 <- subset(lab1, rownames(lab1) %in% common)
      d_n <- d[samp2,]
      lab2 <- as.matrix(membership(Rphenograph(as.matrix(d_n), k = k_arr[j])[[2]]))
      rownames(lab2) <- samp2
      lab2 <- subset(lab2, rownames(lab2) %in% common)
      
      c1 <- matrix(0,nrow=s,ncol=s)
      c2 <- matrix(0,nrow=s,ncol=s)
      for (l in c(1:length(lab1))){
        c1[l,] <- as.integer(lab1==lab1[l])
        c1[l,l] <- 0
        c2[l,] <- as.integer(lab2==lab2[l])
        c2[l,l] <- 0
      }
      sim[j,i] <- sum(dot(c1,c2))/sqrt(sum(dot(c1,c1))*sum(dot(c2,c2)))
      remove(c1,c2)
    }
  }
  sink()
  return(sim)
}

imax <- 20
k_arr <- c(45,100,300,1000)
sim <- get_simmatrix(d_trans, imax, k_arr)

p <- vector("list", length = length(k_arr))
for (i in c(1:length(k_arr)))
{p[[i]] <- local({
  i <- i
  p1 <- qplot(x=sim[i,]) + geom_histogram(binwidth=0.04) +
    ggtitle(paste("k = ", k_arr[i])) + xlim(0,1) + xlab("")
})}
grid.arrange(grobs = p, ncol = 3)
remove(p,imax,k_arr)

# Differential Expression: Mann-Whitney (unused) ####
# get pval for specific cluster
get_pvals <- function(d,idx1,idx2){
  n <- dim(d)[2]
  
  pvals <- vector(mode="double", length=n)
  log_fold_change <- vector(mode="double", length=n)
  for (j in c(1:n)){
    pvals[j] <- wilcox.test(x=d[idx1,j], y=d[idx2,j])$p.value
    log_fold_change[j] <- log(mean(d[idx1,j])/mean(d[idx2,j]))
  }
  return(list(pvals,log_fold_change))
}

pvals_raw <- matrix(nrow = n, ncol = c, dimnames=list(colnames(d_trans),c(1:c)))
pvals_corr <- matrix(nrow = n, ncol = c, dimnames=list(colnames(d_trans),c(1:c)))
log_fold_change <- matrix(nrow = n, ncol = c, dimnames=list(colnames(d_trans),c(1:c)))

# pval correction using benjamini-hochberg
for (i in c(1:c)){
  idx1 <- which(d_tsne$cluster == i)
  idx2 <- c(1:m)[!c(1:m) %in% idx1]
  temp <- get_pvals(d_trans,idx1,idx2)
  pvals_raw[,i] <- temp[[1]]
  log_fold_change[,i] <- temp[[2]]
  pvals_corr[,i] <- p.adjust(pvals_raw[,i], method = "hochberg", n = length(pvals_raw[,i]))
}
remove(temp)

temp <- melt(pvals_raw)
colnames(temp) <- c("gene","cluster","pval")
ggplot(data = temp, aes(x=cluster, y=gene, fill=pval)) + geom_tile() + scale_x_discrete(limits=c(1:c)) + scale_fill_gradientn(colours=bluered(16)) + ggtitle("Raw pvals")
temp <- melt(pvals_corr)
colnames(temp) <- c("gene","cluster","pval")
ggplot(data = temp, aes(x=cluster, y=gene, fill=pval)) + geom_tile() + scale_x_discrete(limits=c(1:c)) + scale_fill_gradientn(colours=bluered(16)) + ggtitle("Corrected pvals")
temp <- melt(log_fold_change)
colnames(temp) <- c("gene","cluster","log_fold_change")
ggplot(data = temp, aes(x=cluster, y=gene, fill=log_fold_change)) + geom_tile() + scale_x_discrete(limits=c(1:c)) + scale_fill_gradientn(colours=bluered(16)) + ggtitle("Log fold change")
remove(temp)


# Outlier Removal (unused)####
d_comp <- d[c(8:20)]
colnames(d_comp) <- c("CD45","KIR2D3D","GzmA","GzmB","CD8a","NKG2A","CD56","TCRab","LD","CD4","CD49a","CD103","CD16")
d_trans<-as.data.frame(apply(d_comp,2,trans))
temp<-boxplot(d_trans$LD, plot=FALSE)$out
temp<-temp[temp>mean(d_trans$LD)]
d_trans <- d_trans[which(!(d_trans$LD %in% temp)),]
temp<-boxplot(d_trans$CD45, plot=FALSE)$out
temp<-temp[temp<mean(d_trans$CD45)]
d_trans <- d_trans[which(!(d_trans$CD45 %in% temp)),c(2:8,10:13)]
d_comp <- d_comp[rownames(d_trans),c(2:8,10:13)]
remove(temp)

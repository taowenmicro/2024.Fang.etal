
#--代谢组学数据
#-非靶向代谢组学数据类似环境数据
#但是由于代谢物数量较多，所以需要单独分析

library(tidyverse)
library(phyloseq)
library(ggClusterNet)

psG = readRDS("./data.ps.meta.ming/ps_Metabolites.rds")
#---开始分析#--------
ps = psG
Top = 35

path.id = "compound.dry.CK"
map = ps %>%sample_data()
head(map)
# map$Group = paste(map$drought,map$Variety.class,sep = "_")
map$Group = map$drought
map$Group  %>% unique()
map$Group = gsub("MD","SD",map$Group)

sample_data(ps) = map

tax = ps %>% tax_table() %>% as.data.frame()
head(tax)
tax$KEGG.Compound.ID = tax$KEGG.Compound.ID %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(ps)  = tax_table(as.matrix(tax))

pst = ps %>% subset_taxa.wt("KEGG.Compound.ID",c("-"),T) %>%
  tax_glom_wt("KEGG.Compound.ID")


otu = pst %>% vegan_otu() %>% t() %>% as.data.frame()
head(otu)
pst

db = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/compound.description.txt")
head(db)
colnames(db)[1] = "id"

tax = pst %>% vegan_tax()%>% as.data.frame()
head(tax)

tax2 = tax %>% rownames_to_column("id") %>% left_join(db,by = "id")
row.names(tax2) = tax2$id
tax_table(pst) = tax_table(as.matrix(tax2))
tax$Metabolite
head(tax2)
saveRDS(pst,"./data.ps.meta.ming/ps.compounds.ko.compounds.rds")


# 
# dat = read.csv("E:/Shared_Folder/Function_local/R_function/micro/hmdb_all_class_221208.csv")
# head(dat)
# dat = dat %>% filter(KEGGID %in%tax$KEGG.Compound.ID )
# tax2 =  tax %>% left_join(dat,by = c("KEGG.Compound.ID" = "KEGGID")) %>% distinct(metab_id, .keep_all = TRUE)  
# head(tax2)
# row.names(tax2) = tax2$metab_id
# tax_table(ps) = tax_table(as.matrix(tax2))
# saveRDS(ps,".//ps_MS_upper.rds")

ps = pst
#--提取有多少个分组
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
# 设定排序顺序
axis_order = phyloseq::sample_data(ps)$Group %>% unique()

# axis_order = c("CK_NM2_Sensitive" , "CK_BM12_Sensitive" ,"CK_BF5_Sensitive","CK_DX40_Tolerant", "CK_LM36_Tolerant",  "CK_LM33_Tolerant", 
#                "MD_NM2_Sensitive",  "MD_BM12_Sensitive", "MD_BF5_Sensitive","MD_DX40_Tolerant" , "MD_LM36_Tolerant" , "MD_LM33_Tolerant" ,
#                "SD_NM2_Sensitive","SD_BM12_Sensitive","SD_BF5_Sensitive","SD_DX40_Tolerant", "SD_LM36_Tolerant" ,"SD_LM33_Tolerant")  
# axis_order = c("CK_Sensitive","CK_Tolerant","MD_Sensitive","MD_Tolerant","SD_Sensitive","SD_Tolerant")

#---代谢物过滤：去除QC样本和非代谢物样本
# ps <- subset_samples(ps,Group %in% c("QC));ps


#---主题颜色设置#-------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro//total_amplicon.R")
#---扩增子环境布置
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]


#  主题颜色更新-可改method可选anhei,出来暗黑风格的图片，Top就是堆叠柱状图展示的数量
res = theme.col(ps = psG,method = "no",Top =c(Top + 1),gnum = gnum)
mytheme1 = res$mytheme[2]#可选1-4
mytheme2 = res$mytheme[6]; #可选5-8
colset1 = res$col.group[[3]];#1-5可选
colset2 = res$col.bar[[2]];#1-3可选
colset3 = res$col.time[[3]];#1-6可选
colset4 = colset3


# 构建保存结果文件夹
result<- dir.amp(smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path
#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
repath = paste(res1path,"/",path.id,a,"/",sep = "");repath
dir.create(repath)
print(repath)


# 6.2 PLS-DA排序#---------
betapath = paste(otupath,"/beta_orda_total.compounds/",sep = "")
dir.create(betapath)


count = vegan_otu(ps)
map = as.data.frame(sample_data(ps))
map$Group

library(mixOmics)

#PLS-DA分析，这里也是选取2个主成分
plsda.datatm <-plsda(count, map$Group, ncomp = 2)
#PLS-DA without centroids
# mi=c("#1B9E77" ,"#D95F02")
plotIndiv(plsda.datatm , comp = c(1,2),
          group = map$Group, style = 'ggplot2' )
plotIndiv(plsda.datatm , comp = c(1,2),
          group = map$Group, style = 'ggplot2',ellipse = TRUE, 
          size.xlabel = 20, size.ylabel = 20, size.axis = 25, pch = 15, cex = 5)
require(grid)
require(gridExtra)
FileName <- paste(betapath,"/a2_",method,"orig_plot.pdf", sep = "")
ggsave(FileName, width = 6, height = 6)


#----提取数据作图
a = unclass(plsda.datatm)
#--提取坐标值
plotdata = as.data.frame(a$variates$X)
plotdata$SampleType = map$Group
#-提取解释度
eig = a$explained_variance$X
eig[1]
library(ggalt)
library(BiocManager)
# install("ggalt")
p = ggplot(data = plotdata,aes(x=comp1,y=comp2,group=SampleType,color=SampleType))+geom_point(size=5)+
  stat_ellipse(type = "t", linetype = 2)+
  geom_encircle(s_shape=1, expand=0) +
  labs(x=paste("X-variate 1 (", format(100 * eig[1]), "%)", sep=""),
       y=paste("X-variate 2 (", format(100 * eig[2] ), "%)", sep=""))+
  labs(title = "PLS-DA") 
p
# mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p=p+theme_bw()+scale_colour_manual(values = colset1)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
p

method = "pls-da"
FileName <- paste(betapath,"/a2_",method,"_plot.pdf", sep = "")
ggsave(FileName, p, width = 8, height = 8)
FileName1 <- paste(betapath,"/a2_",method,"_plot.jpg", sep = "")
ggsave(FileName1 , p, width = 8, height = 8)


#-----差异代谢物#---------
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")
alppath = paste(repath,"/All_different_metabolites/",sep = "")
dir.create(alppath)

#--非参数检验
result = statSuper(ps = ps,group  = "Group",artGroup = NULL,method = "wilcox")
head(result)
FileName <- paste(alppath,"/data_wlx_all.compounds.csv", sep = "")
write.csv(result,FileName,sep = "")


#--t检验检验--建议四个重复以上
result2 = statSuper(ps = ps,group  = "Group",artGroup = NULL,method = "ttext")
head(result2)
FileName <- paste(alppath,"/data_ttest_all.compounds.csv", sep = "")
write.csv(result2,FileName,sep = "")

id = result2 %>% filter(`CK_SD_fdr`<0.05) %>% arrange(desc(`CK_SD_log2_FC`)) %>%
  .$KEGG.Compound.ID



#----热图和气泡图展示--------

heatpath = paste(repath,"/heapmap_boplot.compound/",sep = "")
dir.create(heatpath)

#--注意map文件中一定要有ID列
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")
map = phyloseq::sample_data(ps)
map$ID = row.names(map)
phyloseq::sample_data(ps) = map


result <- Microheatmap(ps_rela = ps,id = id,col_cluster = FALSE)

p1 <- result[[1]] 
p1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p2 <- result[[2]]
p2
filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
ggsave(filename,p1,width = 14,height = 20)

filename = paste(heatpath,"Topggbubble.pdf",sep = "")
ggsave(filename,p2,width = 14,height = 20)


#----11代谢网络#-------
netpath = paste(repath,"/network_metabolites/",sep = "")
dir.create(netpath)
library(ggClusterNet)
library(igraph)
library(sna)

tax = ps %>% vegan_tax() %>% as.data.frame()
head(tax)
tax$conpo = "conpound"
tax_table(ps) = tax_table(as.matrix(tax))
# source("E:\\Shared_Folder\\Function_local\\R_function/micro/lizi_network.R")
result = network.2(ps = ps,
                 N = 500,
                 big = TRUE,
                 select_layout = F,
                 layout_net = "model_maptree",
                 r.threshold=0.6,
                 p.threshold=0.05,
                 label = TRUE,
                 ncol = gnum,
                 path = netpath,
                 fill = "conpo",
                 zipi = FALSE)
# 全部样本的网络比对
p = result[[1]]
p
# 全部样本网络参数比对
data = result[[2]]

plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 16*gnum,height = 16,limitsize = FALSE)

tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)



# --随机森林特征探索分析#---------
matpath = paste(repath,"/Machine_learing.kegg.compounds/",sep = "")
dir.create(matpath )


source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\MicroMachine_learning.R")

library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
mapping = as.data.frame(phyloseq::sample_data(ps))

head(tax)
optimal = 50
#--随机森林全套-如果圈图尚未显示前面几个，就设定max大一点
rank.names(ps)
result = MicroRF(ps = ps ,
                 group  = "Group",
                 optimal = 50,
                 rfcv = F,
                 nrfcvnum = 5,
                 min = -1,
                 max = 10,
                 fill = NULL,
                 lab = "Metabolite"
)
#火柴图展示前二十个重要的OTU
p <- result[[1]] + 
  mytheme1
p

filename = paste(matpath,"/randonforest_loading.pdf",sep = "")
ggsave(filename,p,width = 16,height = optimal/2)
filename = paste(matpath,"/randonforest_loading.jpg",sep = "")
ggsave(filename,p,width = 16,height = optimal/2)

# 圈图展示
p <- result[[2]]
p
filename = paste(matpath,"/randonforest_loading_circle.pdf",sep = "")
ggsave(filename,p,width = 8,height = 10)
filename = paste(matpath,"/randonforest_loading_circle.jpg",sep = "")
ggsave(filename,p,width = 8,height = 10)

p <- result[[6]]
p
filename = paste(matpath,"/Show_model.pdf",sep = "")
ggsave(filename,p,width = 8,height = 4)
filename = paste(matpath,"/Show_model.jpg",sep = "")
ggsave(filename,p,width = 8,height = 4)


data <- result[[5]]
filename = paste(matpath,"/randomforest_data.csv",sep = "")
write.csv(data,filename,quote = F)
dat1 = data
dim(data)



data1 <- result[[7]]
filename = paste(matpath,"/randomforest_data.all.csv",sep = "")
write.csv(data1,filename,quote = F)


head(data1)



# head(otu_table)
otutab = as.data.frame(t(vegan_otu(ps)))
colnames(otutab) <- gsub("-","_",colnames(otutab))
colnames(otutab) <- gsub("[/]","_",colnames(otutab))
colnames(otutab) <- gsub("[(]","_",colnames(otutab))
colnames(otutab) <- gsub("[)]","_",colnames(otutab))
colnames(otutab) <- gsub("[:]","_",colnames(otutab))
colnames(otutab) <- gsub("[[]","",colnames(otutab))
colnames(otutab) <- gsub("[]]","_",colnames(otutab))
colnames(otutab) <- gsub("[#]","_",colnames(otutab))
colnames(otutab) <- gsub("[+]","_",colnames(otutab))
colnames(otutab) <- gsub(" ","_",colnames(otutab))
colnames(otutab) <- gsub("[,]","_",colnames(otutab))
design = as.data.frame(sample_data(ps))
## 计算相对丰度，计算每个物种丰度均值，按照均值排序
OTU = as.matrix(otutab)
norm = t(t(OTU)/colSums(OTU,na=TRUE)) * 100 # normalization to total 100
norma = norm %>% 
  t() %>% as.data.frame()
#数据分组计算平均值
iris.split <- split(norma,as.factor(design$Group))

iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
# 组合结果
norm2 <- do.call(rbind,iris.apply)%>% 
  t() 
norm2 = as.data.frame(norm2) %>% rownames_to_column("ID")
head(norm2)

data1 = data1 %>% rownames_to_column("ID")

data2 = data1 %>% left_join(norm2,by = "ID")
head(data2)

filename = paste(matpath,"/randomforest_data.all.csv",sep = "")
write.csv(data1,filename,quote = F)

filename = paste(matpath,"/randomforest_data.abundance.csv",sep = "")
write.csv(norm2,filename,quote = F)




# 分组样本展示

library("phyloseq")
library(tidyverse)

mapping = as.data.frame(sample_data(ps))
head(mapping)
mapping$ID = row.names(mapping) 
sample_data(ps) = mapping

dim(data)
head(data)


ps1 <- ps %>%
  subset_taxa(
    row.names(tax_table(ps)) %in% as.character(data$org.id)
  )
ps1


otu_table = as.data.frame(t(vegan_otu(ps1)))
head(otu_table)

design = as.data.frame(sample_data(ps1))
## 计算相对丰度，计算每个物种丰度均值，按照均值排序
OTU = as.matrix(otu_table)
norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
norma = norm %>% 
  t() %>% as.data.frame()
#数据分组计算平均值
iris.split <- split(norma,as.factor(design$Group))

iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
# 组合结果
norm2 <- do.call(rbind,iris.apply)%>% 
  t() 
norm2 = as.data.frame(norm2)
norm2$mean=apply(norm2,1,mean)
norm2$ID = row.names(norm2)
colnames(norm2)
##按照mean、列进行排序desc设置从大到小排序
norm3<- arrange(norm2, desc(mean))
row.names(norm3) = norm3$ID
norm3$ID = NULL
### 提取前30个属
wt = norm3
wt$mean = NULL

### 添加OTU注释信息
head(data)
tax_table = as.data.frame(vegan_tax(ps1))
head(tax_table)
wt_tax = merge(wt,tax_table,by = "row.names",all = F)
head(wt_tax)
row.names(wt_tax) = wt_tax$Row.names

res = wt_tax

res$ID = row.names(wt_tax)

head(data)
head(res)
row.names(res) = res$ID
wt = res[,c(colnames(wt))]
head(wt)
head(data)


# row.names(res) = paste0(res$DESCRPTION,".",res$ID)
wt = wt[match(as.character(data$org.id),row.names(wt)),]


color = colorRampPalette(c("navy", "white", "firebrick3"))(60)

wt2<-sqrt(wt)
wt2[wt2>0.5]<-0.5
wt2<-sqrt(wt2)

head(wt)

library(pheatmap)
p = pheatmap(wt2,fontsize=6,cellwidth = 12, cellheight =6,cluster_rows = FALSE,cluster_cols = FALSE,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
             display_numbers = FALSE,labels_row  = row.names(wt))


p

FileName2 <- paste(matpath,"/","mean_laading",".pdf", sep = "")
ggsave(FileName2, p, width = 6, height =12)
# ggsave(FileName2, p, width = 6, height =12, device = cairo_pdf, family = "Times New Roman" )



#--按照样本展示
head(data)

mapping = as.data.frame(sample_data(ps1))
head(mapping)
mapping$ID = row.names(mapping) 
sample_data(ps) = mapping





ps1 <- ps %>%
  subset_taxa(
    row.names(tax_table(ps)) %in% as.character(data$org.id)
  )
ps1


otu_table = as.data.frame(t(vegan_otu(ps1)))
head(otu_table)

design = as.data.frame(sample_data(ps1))
## 计算相对丰度，计算每个物种丰度均值，按照均值排序
OTU = as.matrix(otu_table)
norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
norma = norm %>% 
  as.data.frame()
head(norma)
### 提取前30个属
wt = norma


### 添加OTU注释信息

tax_table = as.data.frame(vegan_tax(ps1))
head(tax_table)
wt_tax = merge(wt,tax_table,by = "row.names",all = F)
head(wt_tax)
row.names(wt_tax) = wt_tax$Row.names

res = wt_tax

res$ID = row.names(wt_tax)

row.names(res) = res$ID
wt = res[,c(colnames(wt))]
head(wt)

wt = wt[match(as.character(data$org.id),row.names(wt)),]


color = colorRampPalette(c("navy", "white", "firebrick3"))(60)

wt2<-sqrt(wt)
wt2[wt2>0.5]<-0.5
wt2<-sqrt(wt2)
head(wt)


library(pheatmap)
p = pheatmap(wt2,fontsize=6,cellwidth = 12, 
             cellheight =6,cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
             display_numbers = FALSE,
             labels_row  = row.names(wt),labels_col  = colnames(wt))

p

FileName2 <- paste(matpath,"/","Sample_laading",".pdf", sep = "")
ggsave(FileName2, p, width = 15, height =12)
#ggsave(FileName2, p, width = 6, height =12, device = cairo_pdf, family = "Times New Roman" )
#---相对丰度标准化展示

ps1 <- ps %>%
  subset_taxa(
    row.names(tax_table(ps)) %in%as.character(data$org.id)
  )
ps1


otu_table = as.data.frame(t(vegan_otu(ps1)))
head(otu_table)

design = as.data.frame(sample_data(ps1))
## 计算相对丰度，计算每个物种丰度均值，按照均值排序
OTU = as.matrix(otu_table)
norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
norma = norm %>% 
  t() %>% as.data.frame()
#数据分组计算平均值
iris.split <- split(norma,as.factor(design$Group))

iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
# 组合结果
norm2 <- do.call(rbind,iris.apply)%>% 
  t() 
norm2 = as.data.frame(norm2)
norm2$mean=apply(norm2,1,mean)
norm2$ID = row.names(norm2)
colnames(norm2)
abun = norm2
head(abun)
abun$mean = NULL
abun_a = reshape2::melt(abun,
                        id.var = c("ID"),
                        variable.name = "id", 
                        value.name = "count")

head(abun_a)
abun_a$iid = rep(paste(1:(length(abun_a$id)/2)),2)

library(plyr)
abun_a1 = ddply(abun_a,"iid",transform,percount = count/sum(count)*100)
head(abun_a1)
mi = c( "firebrick3","navy")

head(abun_a1)
abun_a1$ID = factor(abun_a1$ID,levels = data$org.id)
p = abun_a1  %>%
  ggplot(aes(x = ID, y = percount, fill =id, group =id)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = mi)+
  # theme_classic()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1,size = 4)) 
p

plotname <- paste(matpath,"/a4_random_forst_loaing_abun",optimal,".pdf",sep = "")
ggsave(plotname, p, width = 15, height = 6)




# 机器学习特征网络

library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)

ps
map= sample_data(ps)
ids = map$Group %>% unique()

for (group in ids) {
  pst = ps %>% subset_samples.wt("Group",group) %>%
    subset_taxa.wt("OTU",data$org.id)
  
  
  #----------计算相关
  result = cor_Big_micro2 (ps = pst,
                           N = 0,
                           met.scale = "TMM",
                           r.threshold=0.6,
                           p.threshold=0.05,
                           method = "spearman"
  )
  #--提取相关矩阵
  cor = result[[1]]
  
  
  
  #-网络中包含的OTU的phyloseq文件提取
  ps_net = pst %>% subset_taxa.wt("OTU",colnames(cor))
  otu = ps_net %>% 
    vegan_otu() %>%
    t() %>%
    as.data.frame() %>% rowSums() %>% as.data.frame() %>%rownames_to_column("elements")
  
  netClu  = modulGroup( cor = cor,cut = NULL,method = "cluster_fast_greedy" )
  head(netClu)
  result2 = model_maptree_group(cor = cor,
                                nodeGroup = netClu,
  )
  node = result2[[1]]
  head(node)
  
  #-----计算边
  edge = edgeBuild(cor = cor,node = node) %>% as.data.frame()
  head(edge)
  head(node)
  
  nod = data.frame(row.names = node$elements,v.name = node$elements,ID =  node$elements)
  head(nod)
  
  colnames(otu)[2] = "MeanAbundance"
  nod = nod %>% left_join(otu,by = c("v.name" = "elements"))
  
  # Add label angle
  number_of_bar<-nrow(nod)
  nod$id = seq(1, nrow(nod))
  angle= 360 * (nod$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  nod$hjust<-ifelse(angle>180, 1.05, -0.05)
  nod$angle<-ifelse(angle>180, 90-angle+180, 90-angle)
  
  mygraph <- graph_from_data_frame( data.frame(from = edge$OTU_2,to = edge$OTU_1,wei_label = edge$cor), 
                                    vertices = nod, 
                                    directed = FALSE )
  
  # Find community
  com <- walktrap.community(mygraph)
  # mycolor = c("#377EB8","#E41A1C" )
  tem = edge$cor %>% unique() %>% length()
  if (tem == 1) {
    mycolor = c("#E41A1C" )
    
  } else if(tem == 2){
    mycolor = c("#377EB8","#E41A1C" )
  }
  
  igraph = make_igraph(cor )
  net.pro = net_properties(igraph)
  filename = paste(matpath,"/",group,"network.proterties.csv",sep = "")
  write.csv(net.pro,filename,quote = F)
  
  
  
  #(b)曲线链接
  library(ggraph)
  p = ggraph(mygraph,layout = 'linear', circular = TRUE) +
    geom_edge_arc(aes(edge_colour=as.factor(wei_label)), edge_alpha=0.5, edge_width=0.3) +
    geom_node_point( aes(size = MeanAbundance),shape=21,color='black',alpha=0.9,fill = "#FF7F00") +
    scale_size_continuous(range=c(0.5,10)) +
    scale_fill_manual(values=mycolor ) +
    geom_node_text(aes(x = x*1.06, y=y*1.06, label=ID, angle=angle,hjust=hjust),size=5,color = "black") +
    scale_color_manual(values=mycolor) +
    scale_edge_color_manual(values=mycolor) +
    expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))+
    coord_fixed()+
    theme_minimal() +
    theme(
      legend.position="none",
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks =element_blank(),
      axis.text =element_blank(),
      axis.title = element_blank(),
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) 
  
  p
  
  
  filename = paste(matpath,"/",group,"special.network.HIGH.pdf",sep = "")
  ggsave(filename,p,width = 8,height = 7)
  filename = paste(matpath,"/",group,"special.network.HIGH.jpg",sep = "")
  ggsave(filename,p,width = 8,height = 7)
  
  filename = paste(matpath,"/",group,"special.network.data.node.High.csv",sep = "")
  write.csv(nod,filename,quote = F)
  
  filename = paste(matpath,"/",group,"special.network.data.edge.High.csv",sep = "")
  write.csv(edge,filename,quote = F)
  
}


#----11代谢网络#-------
#  全部的更新都可以使用

#  新网络分析-2023年末更新#--------

netpath = paste(repath,"/network.new.compounds/",sep = "")
dir.create(netpath)

rank.names(ps)
library(ggrepel)
library(igraph)
# detach("package:MicrobiotaProcess")

# 8.1 网络分析主函数#--------
rank.names(ps)


tax = ps %>% tax_table() %>% as.data.frame()
head(tax)
tax$fill ="Metabolites"
tax$id = NULL
tax$ID = NULL
tax_table(ps)  = tax_table(as.matrix(tax))

tab.r = network.pip(
  ps = ps,
  N = 500,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  method = "spearman",
  label = TRUE,
  lab = "elements",
  group = "Group",
  fill = "fill",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)


#  建议保存一下输出结果为R对象，方便之后不进行相关矩阵的运算，节约时间
saveRDS(tab.r,paste0(netpath,"network.pip.sparcc.rds"))
tab.r = readRDS(paste0(netpath,"network.pip.sparcc.rds"))



dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab
# 大型相关矩阵跑出来不容易，建议保存，方便各种网络性质的计算
saveRDS(cortab,paste0(netpath,"cor.matrix.all.group.rds"))
cor = readRDS(paste0(netpath,"cor.matrix.all.group.rds"))

#-提取全部图片的存储对象
plot = tab.r[[1]]
# 提取网络图可视化结果
p0 = plot[[1]]

ggsave(paste0(netpath,"plot.network.pdf"),p0,width = 12,height = 5)
ggsave(paste0(netpath,"plot.network2.pdf"),p0,width = 16,height = 10)


p1 = plot[[3]]

ggsave(paste0(netpath,"plot.network.compare.pdf"),p1,width = 12,height = 5)
ggsave(paste0(netpath,"plot.network2.compare.pdf"),p1,width = 16,height = 10)


# 8.2 网络属性计算#---------
i = 1
id = names(cor)
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  dat = net_properties.4(igraph,n.hub = F)
  head(dat,n = 16)
  colnames(dat) = id[i]
  
  if (i == 1) {
    dat2 = dat
  } else{
    dat2 = cbind(dat2,dat)
  }
}

head(dat2)

FileName <- paste(netpath,"net.network.attribute.data.csv", sep = "")
write.csv(dat2,FileName,quote = F)

# 8.3 单个样本的网络属性#-------

for (i in 1:length(id)) {
  pst = ps %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(pst = pst,cor = cor[[id[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}

head(dat.f2)

FileName <- paste(netpath,"net.network.attribute.data.sample.csv", sep = "")
write.csv(dat.f2,FileName,quote = F)

map= sample_data(ps)
map$ID = row.names(map)
map = map %>% as.tibble()

dat3 = dat.f2 %>% rownames_to_column("ID") %>% inner_join(map,by = "ID")
FileName <- paste(netpath,"net.network.attribute.data.sample.add.group.info.csv", sep = "")
write.csv(dat3,FileName,quote = F)

# 8.4 计算节点属性#---------

for (i in 1:length(id)) {
  
  igraph= cor[[id[i]]] %>% make_igraph()
  nodepro = node_properties(igraph) %>% as.data.frame()
  nodepro$Group = id[i]
  head(nodepro)
  colnames(nodepro) = paste0(colnames(nodepro),".",id[i])
  nodepro = nodepro %>%
    as.data.frame() %>%
    rownames_to_column("ASV.name")
  
  
  # head(dat.f)
  if (i == 1) {
    nodepro2 = nodepro
  } else{
    
    nodepro2 = nodepro2 %>% full_join(nodepro,by = "ASV.name")
    
  }
}

head(nodepro2)

FileName <- paste(netpath,"net.node.attribute.data.sample.csv", sep = "")
write_csv(nodepro2,FileName)




# 8.5 定制输出图片#----

dat = tab.r[[2]]
node = dat$net.cor.matrix$node
edge = dat$net.cor.matrix$edge
head(edge)
head(node)

node2  = add.id.facet(node,"Group")
head(node2)

p <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor
),
data = edge, size = 0.03,alpha = 0.1) +
  geom_point(aes(X1, X2,
                 fill = fill,
                 size = igraph.degree),
             pch = 21, data = node,color = "gray40") +
  facet_wrap(.~ label,scales="free_y",nrow = 1) +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  geom_text(aes(X1, X2,label = Metabolite),pch = 21, data = node2) +
  scale_colour_manual(values = c("#6D98B5","#D48852")) +
  # scale_fill_hue()+
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p

ggsave(paste0(netpath,"plot.network.label.pdf"),p,width = 88,height = 50,limitsize = FALSE)



p <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor
),
data = edge, size = 0.03,alpha = 0.1) +
  geom_point(aes(X1, X2,
                 fill = fill,
                 size = igraph.degree),
             pch = 21, data = node,color = "gray40") +
  facet_wrap(.~ label,scales="free_y",nrow = 1) +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  geom_text_repel(aes(X1, X2,label = Metabolite),pch = 21, data = node2) +
  scale_colour_manual(values = c("#6D98B5","#D48852")) +
  # scale_fill_hue()+
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p

ggsave(paste0(netpath,"plot.network.label.repel.pdf"),p,width = 88,height = 50,limitsize = FALSE)

# g <- plot_segment_across_facets(p, from = 1, to = 2, 
#                                 point_idx1 = 11,
#                                 point_idx2 = 1,
#                                 gp=gpar(lty=2, alpha=1))
# grid.newpage()
# grid.draw(g)
# 
# 
# 
# ps <- c(2,3,3,4,4,5,5,6,6,7)
# 
# while (length(ps) > 0){
#   g <- plot_segment_across_facets(g, from=ps[1], to=ps[2])
#   ps <- ps[-c(1:2)]
# }
# 
# grid.newpage()
# grid.draw(g)


#-8.6 zipi可视化-定制#-----

dat.z = dat$zipi.data
head(dat.z)
x1<- c(0, 0.62,0,0.62)
x2<- c( 0.62,1,0.62,1)
y1<- c(-Inf,2.5,2.5,-Inf)
y2 <- c(2.5,Inf,Inf,2.5)
lab <- c("peripheral",'Network hubs','Module hubs','Connectors')
roles.colors <- c("#E6E6FA","#DCDCDC","#F5FFFA", "#FAEBD7")
tab = data.frame(x1 = x1,y1 = y1,x2 = x2,y2 = y2,lab = lab)
tem = dat.z$group %>% unique() %>% length()
for ( i in 1:tem) {
  if (i == 1) {
    tab2 = tab
  } else{
    tab2 = rbind(tab2,tab)
  }
}


p <- ggplot() +
  geom_rect(data=tab2,
            mapping=aes(xmin=x1,
                        xmax=x2,
                        ymin=y1,
                        ymax=y2,
                        fill = lab))+
  guides(fill=guide_legend(title="Topological roles")) +
  scale_fill_manual(values = roles.colors)+
  geom_point(data=dat.z,aes(x=p, y=z,color=module)) + theme_bw()+
  guides(color= F) +
  ggrepel::geom_text_repel(data = dat.z,
                           aes(x = p, y = z,
                               color = module,label=label),size=4)+
  # facet_wrap(.~group) +
  facet_grid(.~ group, scale='free') +
  theme(strip.background = element_rect(fill = "white"))+
  xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
p

ggsave(paste0(netpath,"plot.network2.zipi.pdf"),p,width = 16,height = 10)



# 8.7 随机网络，幂率分布#-------
dat.r = dat$random.net.data

p3 <- ggplot(dat.r) +
  geom_point(aes(x = ID,y = network,
                 group =group,fill = group),pch = 21,size = 2) +
  geom_smooth(aes(x = ID,y = network,group =group,color = group))+
  facet_grid(.~g,scales = "free") +
  theme_bw() + theme(
    plot.margin=unit(c(0,0,0,0), "cm")
  )
p3




#-8.8 网络显著性比较#-----
dat = module.compare.net.pip(
  ps = NULL,
  corg = cor,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 3)
res = dat[[1]]
head(res)


FileName <- paste(netpath,"net.compare.diff.sig.csv", sep = "")
write.csv(res,FileName,quote = F)



# 8.91  网络稳定性-模块相似性#--------
library(tidyfst)

res1 = module.compare.m(
  ps = NULL,
  corg = cor,
  zipi = FALSE,
  zoom = 0.2,
  padj = F,
  n = 3)

#不同分组使用一个圆圈展示，圆圈内一个点代表一个模块，相连接的模块代表了相似的模块。
p1 = res1[[1]]
p1
#--提取模块的OTU，分组等的对应信息
dat1 = res1[[2]]
head(dat1)
#模块相似度结果表格
dat2 = res1[[3]]
head(dat2)


dat2$m1 = dat2$module1 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$m2 = dat2$module2 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$cross = paste(dat2$m1,dat2$m2,sep = "_Vs_")
# head(dat2)
dat2 = dat2 %>% filter(module1 != "none")

p2 = ggplot(dat2) + geom_bar(aes(x = cross,fill = cross)) +
  labs(x = "",
       y = "numbers.of.similar.modules"
  )+ theme_classic()

p2


#--发现分组1和分组3网络更相似一些
FileName <- paste(netpath,"module.compare.groups.pdf", sep = "")
ggsave(FileName, p1, width = 10, height = 10)

FileName <- paste(netpath,"numbers.of.similar.modules.pdf", sep = "")
ggsave(FileName, p2, width = 8, height = 8)

FileName <- paste(netpath,"module.otu.csv", sep = "")
write.csv(dat1,FileName, quote = F)

FileName <- paste(netpath,"module.compare.groups.csv", sep = "")
write.csv(dat2,FileName, quote = F)



#-8.92 网络稳定性-去除关键节点-网络鲁棒性#-----
# 鲁棒性计算需要物种丰富，所以即使计算好了相关矩阵，也需要输入ps对象
library(patchwork)
res2= Robustness.Targeted.removal(ps = ps,
                                  corg = cor,
                                  degree = TRUE,
                                  zipi = FALSE
)
p3 = res2[[1]]
p3
#提取数据
dat4 = res2[[2]]


# dir.create("./Robustness_Random_removal/")
path = paste(netpath,"/Robustness_Random_removal/",sep = "")
fs::dir_create(path)
write.csv(dat4,
          paste(path,"random_removal_network.csv",sep = ""))
ggsave(paste(path,"random_removal_network.pdf",sep = ""),  p3,width = 8,height = 4)



#-8.93 网络稳定性-随即取出任意比例节点-网络鲁棒性#---------
res3 = Robustness.Random.removal(ps = ps,
                                 corg = cortab,
                                 Top = 0
)
p4 = res3[[1]]
p4
#提取数据
dat5 = res3[[2]]
# head(dat5)
path = paste(netpath,"/Robustness_Targeted_removal/",sep = "")
fs::dir_create(path)
write.csv(dat5,
          paste(path,"Robustness_Targeted_removal_network.csv",sep = ""))
ggsave(paste(path,"Robustness_Targeted_removal_network.pdf",sep = ""),  p4,width = 8,height = 4)



#-8.94 网络稳定性-计算负相关的比例#----
res4 = negative.correlation.ratio(ps = ps,
                                  corg = cortab,
                                  # Top = 500,
                                  degree = TRUE,
                                  zipi = FALSE)

p5 = res4[[1]]
p5
dat6 = res4[[2]]
#-负相关比例数据
# head(dat6)
path = paste(netpath,"/Vulnerability/",sep = "")
fs::dir_create(path)

write.csv(dat6,
          paste(path,"Vnegative.correlation.ratio_network.csv",sep = ""))
ggsave(paste(path,"negative.correlation.ratio_network.pdf",sep = ""),  p5,width = 4,height = 4)




#8.96 网络稳定性-网络抗毁性#------
library("pulsar")
res6 = natural.con.microp (
  ps = ps,
  corg = cor,
  norm = TRUE,
  end = 150,# 小于网络包含的节点数量
  start = 0
)
p7 = res6[[1]]
p7
dat8  = res6[[2]]
path = paste(netpath,"/Natural_connectivity/",sep = "")
fs::dir_create(path)
write.csv(dat8,
          paste(path,"/Natural_connectivity.csv",sep = ""))
ggsave(paste(path,"/Natural_connectivity.pdf",sep = ""),  p7,width = 5,height = 4)





#  差异基因和差异代谢物网络#----------
psr = base::readRDS("./data.ps.meta.ming/ps_transp.rds")

map = sample_data(psr)
head(map)
# map$Group = gsub("[?]","",map$treatment)
# map$map$Group  %>% unique()
head(map)
map$Group = paste(map$drought,sep = ".")
map$Group = gsub("MD","SD",map$Group)

map$Group  %>% unique()
sample_data(psr) = map


tax = as.data.frame(vegan_tax(psr))
head(tax)
tax$KO_id = tax$KO_id  %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(psr)  = tax_table(as.matrix(tax))
psr2 <- psr %>% subset_taxa.wt("KO_id" ,"" ,T) %>% 
  tax_glom_wt("KO_id")
psr2

diffpath.2 = "./tem"
dir.create(diffpath.2)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
res = EdgerSuper(ps = psr2,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)
head(res)
id = res %>% filter(`CK-SDlevel` != "nosig") %>% arrange(desc(`CK-SDlogFC`)) %>% .$species
psr3 = psr2 %>% subset_taxa.wt("OTU",id)


psG = readRDS("./data.ps.meta.ming/ps_Metabolites.rds")
#---开始分析#--------
ps = psG
Top = 35

path.id = "compound.dry.CK"
map = ps %>%sample_data()
head(map)
# map$Group = paste(map$drought,map$Variety.class,sep = "_")
map$Group = map$drought
map$Group  %>% unique()
map$Group = gsub("MD","SD",map$Group)

sample_data(ps) = map

tax = ps %>% tax_table() %>% as.data.frame()
head(tax)
tax$KEGG.Compound.ID = tax$KEGG.Compound.ID %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(ps)  = tax_table(as.matrix(tax))

psG2 = ps %>% subset_taxa.wt("KEGG.Compound.ID",c("-"),T) %>%
  tax_glom_wt("KEGG.Compound.ID")

#--t检验检验--建议四个重复以上
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")
result2 = statSuper(ps = psG2,group  = "Group",artGroup = NULL,method = "ttext")
head(result2)

id2  = result2 %>% filter(`CK_SD_fdr`<0.05) %>% arrange(desc(`CK_SD_log2_FC`)) %>%
  .$KEGG.Compound.ID
psG3 = psG2 %>% subset_taxa.wt("OTU",id2)




#  网络图中的标签我们要进行很好的设定


ps.all = merge.ps(ps1 = psr3,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "rna",
                  dat2.lab = "compounds")
ps.all 

# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.network/")
dir.create(path)
path
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 60,height = 30,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 90,height = 40,limitsize = FALSE)


map = ps.all %>% sample_data()
map$Group = "one"
sample_data(ps.all) = map


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.8,
  p.threshold=0.01,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.network.all.togather/")
dir.create(path)
path
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 40,height = 40,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 50,height = 50,limitsize = FALSE)



#  差异Mkegg和差异代谢物网络#----------
psr = base::readRDS("./data.ps.meta.ming/ps_transp.rds")

map = sample_data(psr)
head(map)
# map$Group = gsub("[?]","",map$treatment)
# map$map$Group  %>% unique()
head(map)
map$Group = paste(map$drought,sep = ".")
map$Group = gsub("MD","SD",map$Group)

map$Group  %>% unique()
sample_data(psr) = map


tax = as.data.frame(vegan_tax(psr))
head(tax)
tax$KO_id = tax$KO_id  %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(psr)  = tax_table(as.matrix(tax))
psr2 <- psr %>% subset_taxa.wt("KO_id" ,"" ,T) %>% 
  tax_glom_wt("KO_id")
psr2

getOption("clusterProfiler.download.method")
R.utils::setOption( "clusterProfiler.download.method",'auto' )
#--基于Mkegg进行分析
Mkegg = clusterProfiler::download_KEGG('ko',keggType = "MKEGG")
PATH2ID <- Mkegg$KEGGPATHID2EXTID
PATH2NAME <- Mkegg$KEGGPATHID2NAME
head(PATH2NAME)
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')
head(PATH_ID_NAME )


ko = psr2 %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)

tem2 = ko %>% left_join(PATH_ID_NAME,by = c("KO" = "KO")) %>%
  tidyfst::filter_dt(KEGGID != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("KEGGID","DESCRPTION")) %>%
  as.data.frame()
head(tem2)
colnames(tem2)

otu = tem2[,sample_names(psr2)]
tax = data.frame(MDESCRPTION = tem2$DESCRPTION,MDESCRPTION2 = tem2$DESCRPTION,row.names = tem2$KEGGID)


row.names(otu) = tem2$KEGGID
head(otu)
psr4 = phyloseq(
  otu_table(as.matrix(otu),taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(ps)
)


otupath = "./"
diffpath.2 = paste(otupath,"/EDgeR.Mkegg/",sep = "")
dir.create(diffpath.2)

res = EdgerSuper(ps = ps2,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)

head(res)


id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`CK-SDlevel` != "nosig") %>%
  arrange(desc(`CK-SDlogFC`))
id = id$ID
psr5 = psr4 %>% subset_taxa.wt("OTU",id)



psG = readRDS("./data.ps.meta.ming/ps_Metabolites.rds")
#---开始分析
ps = psG
Top = 35

path.id = "compound.dry.CK"
map = ps %>%sample_data()
head(map)
# map$Group = paste(map$drought,map$Variety.class,sep = "_")
map$Group = map$drought
map$Group  %>% unique()
map$Group = gsub("MD","SD",map$Group)

sample_data(ps) = map

tax = ps %>% tax_table() %>% as.data.frame()
head(tax)
tax$KEGG.Compound.ID = tax$KEGG.Compound.ID %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(ps)  = tax_table(as.matrix(tax))

psG2 = ps %>% subset_taxa.wt("KEGG.Compound.ID",c("-"),T) %>%
  tax_glom_wt("KEGG.Compound.ID")

#--t检验检验--建议四个重复以上
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")
result2 = statSuper(ps = psG2,group  = "Group",artGroup = NULL,method = "ttext")
head(result2)

id2  = result2 %>% filter(`CK_SD_fdr`<0.05) %>% arrange(desc(`CK_SD_log2_FC`)) %>%
  .$KEGG.Compound.ID
psG3 = psG2 %>% subset_taxa.wt("OTU",id2)


#  网络图中的标签我们要进行很好的设定




ps.all = merge.ps(ps1 = psr5,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "rna",
                  dat2.lab = "compounds")


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.network.mkegg/")
dir.create(path)
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 60,height = 30,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 90,height = 40,limitsize = FALSE)


map = ps.all %>% sample_data()
map$Group = "one"
sample_data(ps.all) = map


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.network.all.togather/")
dir.create(path)
path
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 40,height = 40,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 50,height = 50,limitsize = FALSE)



#  差异rection和差异代谢物网络#----------
psr = base::readRDS("./data.ps.meta.ming/ps_transp.rds")

map = sample_data(psr)
head(map)
# map$Group = gsub("[?]","",map$treatment)
# map$map$Group  %>% unique()
head(map)
map$Group = paste(map$drought,sep = ".")
map$Group = gsub("MD","SD",map$Group)

map$Group  %>% unique()
sample_data(psr) = map


tax = as.data.frame(vegan_tax(psr))
head(tax)
tax$KO_id = tax$KO_id  %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(psr)  = tax_table(as.matrix(tax))
psr2 <- psr %>% subset_taxa.wt("KO_id" ,"" ,T) %>% 
  tax_glom_wt("KO_id")
psr2

dat = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/ko-reaction.txt",
                 header = F
)

colnames(dat) = c("ID","reaction")
head(dat)



# kegg <- clusterProfiler::download_KEGG('ko')
# PATH2ID <- kegg $KEGGPATHID2EXTID
# PATH2NAME <- kegg$KEGGPATHID2NAME
# head(PATH2NAME)
# PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
# colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')
# head(PATH_ID_NAME )

ko = psr2 %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)
dat$ID = gsub("ko:","",dat$ID)
dat$reaction = gsub("rn:","",dat$reaction)


tem2 = ko %>% left_join(dat,by = c("KO" = "ID")) %>%
  tidyfst::filter_dt(reaction != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("reaction")) %>%
  as.data.frame()

colnames(tem2)

otu = tem2[,sample_names(psr2)]
row.names(otu) = tem2$reaction
head(otu)


dat2 = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/reaction.txt",
                  header = F
)
head(dat2)
colnames(dat2) = c("reaction","DESCRPTION")

tax = data.frame(row.names = dat2$reaction,DESCRPTION = dat2$DESCRPTION,DESCRPTION2 = dat2$DESCRPTION)
head(tax)

psr6 = phyloseq(
  otu_table(as.matrix(otu),taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(ps)
)


otupath = "./"
diffpath.2 = paste(otupath,"/EDgeR.reaction/",sep = "")
dir.create(diffpath.2)

res = EdgerSuper(ps = psr6,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)

head(res)


id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`CK-SDlevel` != "nosig") %>%
  arrange(desc(`CK-SDlogFC`))
id = id$ID
psr7 = psr6 %>% subset_taxa.wt("OTU",id)



psG = readRDS("./data.ps.meta.ming/ps_Metabolites.rds")
#---开始分析
ps = psG
Top = 35

path.id = "compound.dry.CK"
map = ps %>%sample_data()
head(map)
# map$Group = paste(map$drought,map$Variety.class,sep = "_")
map$Group = map$drought
map$Group  %>% unique()
map$Group = gsub("MD","SD",map$Group)

sample_data(ps) = map

tax = ps %>% tax_table() %>% as.data.frame()
head(tax)
tax$KEGG.Compound.ID = tax$KEGG.Compound.ID %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(ps)  = tax_table(as.matrix(tax))

psG2 = ps %>% subset_taxa.wt("KEGG.Compound.ID",c("-"),T) %>%
  tax_glom_wt("KEGG.Compound.ID")

#--t检验检验--建议四个重复以上
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")
result2 = statSuper(ps = psG2,group  = "Group",artGroup = NULL,method = "ttext")
head(result2)

id2  = result2 %>% filter(`CK_SD_fdr`<0.05) %>% arrange(desc(`CK_SD_log2_FC`)) %>%
  .$KEGG.Compound.ID
psG3 = psG2 %>% subset_taxa.wt("OTU",id2)


#  网络图中的标签我们要进行很好的设定




ps.all = merge.ps(ps1 = psr7,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "rna",
                  dat2.lab = "compounds")


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.network.reaction/")
dir.create(path)
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 60,height = 30,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 90,height = 40,limitsize = FALSE)


map = ps.all %>% sample_data()
map$Group = "one"
sample_data(ps.all) = map


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.network.all.togather.reaction/")
dir.create(path)
path
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 40,height = 40,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 50,height = 50,limitsize = FALSE)



merge.ps <- function(ps1 ,
                     ps2,
                     N1 = 100,
                     N2 = 100,
                     scale = TRUE,
                     onlygroup = FALSE,#不进行列合并，只用于区分不同域
                     dat1.lab = "bac",
                     dat2.lab = "fun") {
  
  if (scale == TRUE) {
    if (!is.null(ps16s)) {
      ps1  = phyloseq::transform_sample_counts(ps1, function(x) x / sum(x) )
    }
    if (!is.null(psITS)) {
      ps2  = phyloseq::transform_sample_counts(ps2, function(x) x / sum(x) )
    }
  }
  if (!is.null(ps1)) {
    # ps_16s = phyloseq::filter_taxa(ps16s, function(x) mean(x) > N16s, TRUE)#select OTUs according to  relative abundance
    ps_16s  =  filter_OTU_ps(ps = ps1,Top = N1)
    ###
    otu_table_16s = as.data.frame(t(vegan_otu(ps_16s)))
    row.names(otu_table_16s) = paste(dat1.lab,row.names(otu_table_16s),sep = "_")
    ## change the OTU name of bac and fungi OTU table
    tax_table_16s = as.data.frame(vegan_tax(ps_16s))
    row.names(tax_table_16s) = paste(dat1.lab,row.names(tax_table_16s),sep = "_")
    #-- add a col marked the bac and fungi
    tax_table_16s$filed = rep(dat1.lab,length(row.names(tax_table_16s)))
  }
  if (!is.null(ps2)) {
    # ps_ITS = phyloseq::filter_taxa(psITS, function(x) mean(x) > NITS , TRUE)#select OTUs according to  relative abundance
    ps_ITS = filter_OTU_ps(ps = ps2,Top = N2)
    otu_table_ITS = as.data.frame(t(vegan_otu(ps_ITS)))
    row.names(otu_table_ITS) = paste(dat2.lab,row.names(otu_table_ITS ),sep = "_")
    tax_table_ITS = as.data.frame(vegan_tax(ps_ITS))
    row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
    tax_table_ITS$filed = rep(dat2.lab,length(row.names(tax_table_ITS)))
    
  }
  
  
  if (!is.null(ps2) & !is.null(ps1) ) {
    ## merge OTU table of bac and fungi
    
    
    
    otu_table = rbind(otu_table_16s[,intersect(names(otu_table_ITS),names(otu_table_16s))],
                      
                      otu_table_ITS[,intersect(names(otu_table_ITS),names(otu_table_16s))])
    
    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s,tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table
    
    mapping = as.data.frame( phyloseq::sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))
    
    
  } else if(is.null(psITS) & !is.null(ps16s) ) {
    otu_table = rbind(otu_table_16s)
    
    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,row.names = row.names(otu_table),id = row.names(otu_table)))
    }
    #on of map table as final map table
    mapping = as.data.frame(sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))
    
    
  } else if (!is.null(ps2) & is.null(ps1)){
    otu_table = rbind(otu_table_ITS)
    
    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table
    mapping = as.data.frame( phyloseq::sample_data(psITS))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = T),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))
    
  }
  
  tax = pallps %>% vegan_tax() %>%
    as.data.frame() %>% dplyr::select(filed,everything())
  phyloseq::tax_table(pallps) = as.matrix(tax)
  
  
  return(pallps)
}


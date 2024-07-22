

# 宏基因组#------




#--宏基因组KEGG注释结果分析#--------
library(tidyverse)
library(phyloseq)
library(ggClusterNet)
library(fs)


psko = readRDS("./data.ps.meta.ming/ps_KEGG.rds")


map = psko %>%sample_data()
head(map)
# map$Group = paste(map$drought,map$Variety.class,sep = "_")
map$Group = map$drought
map$Group  %>% unique()
map$Group = gsub("MD","SD",map$Group)
sample_data(psko) = map
tax = psko %>% tax_table() %>% as.data.frame()
head(tax)
#  去除空缺值
otu = psko %>% vegan_otu() %>% t()
otu[is.na(otu)] =0
otu_table(psko) = otu_table(otu,taxa_are_rows = TRUE)

psko = psko %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

path.id = "meta.kegg"
map = psko %>%sample_data()
head(map)
# map$Group  %>% unique()
# head(map)
# map$Group = paste(map$drought,sep = ".")
# map$Group = gsub("MD","SD",map$Group)

map$Group  %>% unique()


sample_data(psko) = map


ps = psko


#--提取有多少个分组#-----------
gnum = sample_data(ps)$Group %>%unique() %>% length()

axis_order = sample_data(ps)$Group %>%unique()

# 设定排序顺序--ID
map = sample_data(ps)
map$ID = row.names(map)
map = map %>% 
  as.tibble() %>%
  dplyr::arrange(desc(Group))
axis_order.s = map$ID

# ps = ps %>% filter_OTU_ps(8000)
#--物种多样性分析#------
# source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/total_amplicon_Meta.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
#  主题颜色更新-可改method可选anhei,出来暗黑风格的图片，Top就是堆叠柱状图展示的数量
Top = 20
res = theme.col(ps = ps,method = "no",Top =c(Top + 1),gnum = gnum)
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
otupath = paste(res1path,"/",path.id,a,"/",sep = "");otupath
dir.create(otupath)
print(otupath)


#-----差异基因#---------

diffpath = paste(otupath,"/Different.gene.ko/",sep = "")
dir.create(diffpath)
md = c("edgr","t","desep2","wlx")
# 准备脚本
source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/EdgerSuper_Meta.R")
# source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/Plot.CompareWithCK.R",encoding = "utf-8")
supath = paste0(diffpath,"/EdgeR/") 
dir.create(supath)
res = EdgerSuper(ps = psko,group  = "Group",artGroup =NULL,
                 path = supath
)
head(res)

id = res %>% 
  rownames_to_column("ID") %>% 
  filter(`CK-SDlevel` != "nosig") %>% arrange(desc(`CK-SDlogFC`)) %>% .$ID
id


#----热图和气泡图展示基因--------

heatpath = paste(otupath,"/heapmap_boplot.gene.meta/",sep = "")
dir.create(heatpath)

#--注意map文件中一定要有ID列
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")
map = phyloseq::sample_data(psko)
map$ID = row.names(map)
phyloseq::sample_data(psko) = map


result <- Microheatmap(ps_rela = psko,id = id,col_cluster = FALSE)

p1 <- result[[1]] 
p1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p2 <- result[[2]]
p2
filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
ggsave(filename,p1,width = 14,height = 40)

filename = paste(heatpath,"Topggbubble.pdf",sep = "")
ggsave(filename,p2,width = 14,height = 40)



#-----差异Mkegg#---------

diffpath = paste(otupath,"/Different.gene.Mkegg/",sep = "")
dir.create(diffpath)
md = c("edgr","t","desep2","wlx")
# 准备脚本
source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/EdgerSuper_Meta.R")
# source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/Plot.CompareWithCK.R",encoding = "utf-8")
supath = paste0(diffpath,"/EdgeR/") 
dir.create(supath)
res = EdgerSuper(ps = psko3,group  = "Group",artGroup =NULL,
                 path = supath
)
head(res)

id = res %>% 
  rownames_to_column("ID") %>% 
  filter(`CK-SDlevel` != "nosig") %>% arrange(desc(`CK-SDlogFC`)) %>% .$ID
id


#----热图和气泡图展示基因--------

heatpath = paste(otupath,"/heapmap_boplot.Mkegg.meta/",sep = "")
dir.create(heatpath)

#--注意map文件中一定要有ID列
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")
map = phyloseq::sample_data(psko3)
map$ID = row.names(map)
phyloseq::sample_data(psko3) = map


result <- Microheatmap(ps_rela = psko3,id = id,col_cluster = FALSE)

p1 <- result[[1]] 
p1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p2 <- result[[2]]
p2
filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
ggsave(filename,p1,width = 14,height = 15)

filename = paste(heatpath,"Topggbubble.pdf",sep = "")
ggsave(filename,p2,width = 14,height = 15)



#-----差异rection#---------

diffpath = paste(otupath,"/Different.gene.reaction/",sep = "")
dir.create(diffpath)
md = c("edgr","t","desep2","wlx")
# 准备脚本
source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/EdgerSuper_Meta.R")
# source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/Plot.CompareWithCK.R",encoding = "utf-8")
supath = paste0(diffpath,"/EdgeR/") 
dir.create(supath)
res = EdgerSuper(ps = psko5,group  = "Group",artGroup =NULL,
                 path = supath
)
head(res)

id = res %>% 
  rownames_to_column("ID") %>% 
  filter(`CK-SDlevel` != "nosig") %>% arrange(desc(`CK-SDlogFC`)) %>% .$ID
id


#----热图和气泡图展示基因--------

heatpath = paste(otupath,"/heapmap_boplot.reaction.meta/",sep = "")
dir.create(heatpath)

#--注意map文件中一定要有ID列
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")
map = phyloseq::sample_data(psko5)
map$ID = row.names(map)
phyloseq::sample_data(psko5) = map


result <- Microheatmap(ps_rela = psko5,id = id,col_cluster = FALSE)

p1 <- result[[1]] 
p1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p2 <- result[[2]]
p2
filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
ggsave(filename,p1,width = 14,height = 35)

filename = paste(heatpath,"Topggbubble.pdf",sep = "")
ggsave(filename,p2,width = 14,height = 35)



#--2.1_KEGG基因整体差异排序分析#------
betapath = paste(otupath,"/ko_ordinate/",sep = "")
dir.create(betapath)


source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/BetaDiv__Meta.R")
source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/MicroTest_Meta.R")
source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function//pairMicroTest_Meta.R")


# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA
library(pacman)
library(MicrobiotaProcess)
p_unload(MicrobiotaProcess)
methodlist = c("PCoA", "PCA")

for (method in methodlist) {
  result = BetaDiv(ps = psko, group = "Group", dist = "bray", method = method, Micromet = "MRPP", 
                   pvalue.cutoff = 0.05,
                   pair = F
  )
  p3_1 = result[[1]] + 
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) +
    mytheme1 
  # theme(legend.position = c(0.2,0.2))
  p3_1
  #带标签图形出图
  p3_2 = result[[3]] +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 
  # theme(legend.position = c(0.2,0.2))
  p3_2
  
  FileName <- paste(betapath,"/a2_",method,"bray.pdf", sep = "")
  ggsave(FileName, p3_1, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"",method,"bray.jpg", sep = "")
  ggsave(FileName1 , p3_1, width = 12, height = 12)
  
  FileName <- paste(betapath,"/a2_",method,"bray_label.pdf", sep = "")
  ggsave(FileName, p3_2, width = 12, height = 12)
  FileName1 <- paste(betapath,"/a2_",method,"bray_label.jpg", sep = "")
  ggsave(FileName1 , p3_2, width = 12, height = 12)
  
  # 提取出图数据
  plotdata = result[[2]]
  FileName <-  paste(betapath,"/a2_",method,"bray.csv", sep = "")
  write.csv(plotdata,FileName)
  #---------排序-精修图
  plotdata =result[[2]]
  head(plotdata)
  # 求均值
  cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
  cent
  # 合并到样本坐标数据中
  segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
                by = 'Group', sort = FALSE)
  
  # p2$layers[[2]] = NULL
  # library(ggcor)
  library(ggsci)
  p3_3 = p3_1 +geom_segment(data = segs,
                            mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
    geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 
  # theme(legend.position = c(0.2,0.2))
  p3_3
  
  FileName <- paste(betapath,"/a2_",method,"bray_star.pdf", sep = "")
  ggsave(FileName, p3_3, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"bray_star.jpg", sep = "")
  ggsave(FileName1 , p3_3, width = 8, height = 8)
  
}

map
#提取总体比较
TResult =result[[5]]
head(TResult)

# 提取两两检测结果
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_adonis.csv", sep = "")
write.csv(pair,FileName)
FileName <- paste(betapath,"Total_adonis.csv", sep = "")
write.csv(TResult,FileName)


#  功能基因驱动微生物#---------


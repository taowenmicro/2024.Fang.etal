

# 根际细菌#--------
# 16s数据
ps01 = readRDS("data/ps_16s.rds") %>% subset_samples.wt("zone","rhi")
map = ps01 %>%sample_data()
head(map)
map$id2 = gsub("RS_","",row.names(otu1))
map$Group = map$drought
map$Group  %>% unique()
map$Group = gsub("MD","SD",map$Group)

sample_data(ps01) = map
ps01 = ps01 %>% changeSamplenames("id2")
ps01 = ps01 %>% remove.zero()

ps02 = ps01 %>% tax_glom_wt("Genus")
ps = ps01
path.id = "16S.dry.rhi"
#--提取有多少个分组
Top = 35
gnum = sample_data(ps)$Group %>%unique() %>% length()
gnum
# 设定排序顺序
map$ID

axis_order = sample_data(ps)$Group %>%unique()


map = sample_data(ps)
map$ID = row.names(map)
# map = map %>% 
#   as.tibble() %>%
#   dplyr::arrange(desc(Group))

axis_order.s = map$ID

source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]


#  主题颜色更新-可改method可选anhei,出来暗黑风格的图片，Top就是堆叠柱状图展示的数量
res = theme.col(ps = ps,method = "no",Top =c(Top + 1),gnum = gnum)
# res = theme.col(ps,method = "no",Top = 10)
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




#--alpha多样性#---------
alppath = paste(otupath,"/alpha/",sep = "")
dir.create(alppath)


# source("../micro/alpha-diversity.R")

#---多种指标alpha多样性分析加出图-标记显著性
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")
index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )

#--多种组合alpha分析和差异分析出图
alp = alpha(ps = ps,inde="Shannon",group = "Group",Plot = TRUE )
index= alp
head(index)

#--提取三个代表指标作图
all.alpha = T

if (all.alpha) {
  sel = c(match("Inv_Simpson",colnames(index)),
          match("Pielou_evenness",colnames(index)),
          match("Simpson_evenness",colnames(index)),
          match("Richness",colnames(index)),
          match("Chao1",colnames(index)),
          match("ACE",colnames(index)),
          match("Shannon",colnames(index))
          
  )
  h = 3
} else{
  sel = c(match("Shannon",colnames(index)),match("Richness",colnames(index)),
          match("Pielou_evenness",colnames(index)))
  h = 1
}

n = length(sel) + 3


data = cbind(data.frame(ID = 1:length(index$Group),group = index$Group),index[sel])
head(data)


result = EasyStat::MuiKwWlx2(data = data,num = c(3:(n -1)))

FileName <- paste(alppath,"/alpha_diversity_different_label.csv", sep = "")
write.csv(result,FileName,sep = "")
FileName <- paste(alppath,"/alpha_diversity_index.csv", sep = "")
write.csv(index,FileName,sep = "")


result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:(n -1)),
                                          result = result,
                                          sig_show ="abc",ncol = 3 )
p1_1 = result1[[1]] + 
  ggplot2::scale_x_discrete(limits = axis_order) + 
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)
p1_1

#如何升级展示-提取数据用小提琴图展示
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) + 
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 3) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat)) +
  ggplot2::scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(color=guide_legend(title = NULL),
         shape=guide_legend(title = NULL),
         fill = guide_legend(title = NULL)
  ) +
  ggplot2::scale_fill_manual(values = colset1)
p1_0


res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(n -1)),result = result,sig_show ="abc",ncol = 3)
p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) +
  mytheme1+ 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
# p1_2

res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(n -1)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + 
  mytheme1 + 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
# p1_3
w = 1
h = 3


FileName <- paste(alppath,"Alpha_Facet_vl", ".pdf", sep = "")
ggsave(FileName, p1_0, width = ((1 + gnum) *w), height =4 *h,limitsize = FALSE)
FileName <- paste(alppath,"Alpha_Facet_vl", ".png", sep = "")
ggsave(FileName, p1_0, width = ((1 + gnum) *w), height =4 *h,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = (1 + gnum)*w, height =4*h,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = (1 + gnum)*w , height = 4*h,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = (1 + gnum) *w, height = 4*h,limitsize = FALSE)


#--总体差异检测alpha多样性
krusk1 = ggpubr::compare_means( Shannon ~ group, data=data, method = "kruskal.test")
krusk2 = ggpubr::compare_means( Richness ~ group, data=data, method = "kruskal.test")
krusk3 = ggpubr::compare_means( Pielou_evenness ~ group, data=data, method = "kruskal.test")

dat = rbind(krusk1,krusk2,krusk3) %>% as.data.frame()
FileName <- paste(alppath,"/alpha_diversity_index_all_p_Kruskal-Wallis.csv", sep = "")
write_csv(dat,FileName)



#---排序分析beta-diversity#----
betapath = paste(otupath,"/beta/",sep = "")
dir.create(betapath)

# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA tsne 

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

methodlist = c("NMDS","PCoA", "PCA")
map= sample_data(ps)
map$Group %>% unique()
pst = ps
# methodlist = c("LDA")
for (method in methodlist) {
  result = BetaDiv(ps = pst, group = "Group", dist = "bray",
                   method = method, Micromet = "anosim", pvalue.cutoff = 0.05,
                   pair = F)
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
  ggsave(FileName, p3_1, width = 8, height = 7)
  FileName1 <- paste(betapath,"/a2_",method,"",method,"bray.jpg", sep = "")
  ggsave(FileName1 , p3_1, width = 12, height = 11)
  
  FileName <- paste(betapath,"/a2_",method,"bray_label.pdf", sep = "")
  ggsave(FileName, p3_2, width = 12, height = 12)
  FileName1 <- paste(betapath,"/a2_",method,"bray_label.jpg", sep = "")
  ggsave(FileName1 , p3_2, width = 12, height = 11)
  
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
  ggsave(FileName, p3_3, width = 8, height = 7)
  FileName1 <- paste(betapath,"/a2_",method,"bray_star.jpg", sep = "")
  ggsave(FileName1 , p3_3, width = 8, height = 7)
  
}

map

#提取总体比较
TResult =result[[5]]
head(TResult)

# 提取两两检测结果
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
write.csv(pair,FileName)
FileName <- paste(betapath,"Total_anosim.csv", sep = "")
write.csv(TResult,FileName)

#--换用adonis差异分析
title1 = MicroTest(ps = ps, Micromet = "adonis", dist = "bray")
title1
FileName <- paste(betapath,"Total_adonis.csv", sep = "")
write.csv(title1,FileName)

#---物种组成展示#---------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")
barpath = paste(otupath,"/Microbial_composition/",sep = "")
dir.create(barpath)

phyloseq::rank_names(ps)

if (gnum < 12) {
  mytheme2 = mytheme1
}
w = 0.5
for (j in c("Phylum" , "Class" ,  "Order"  , "Family" , "Genus","Species")) {
  result = barMainplot(ps = ps01,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = Top)
  p4_1 <- result[[1]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset2) +
    scale_x_discrete(limits = axis_order) +
    mytheme2
  p4_1
  
  p4_2  <- result[[3]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset2) +
    scale_x_discrete(limits = axis_order) + 
    mytheme2
  p4_2
  
  databar <- result[[2]] %>% group_by(Group,aa) %>%
    dplyr::summarise(sum(Abundance)) %>% as.data.frame()
  head(databar)
  colnames(databar) = c("Group",j,"Abundance(%)")
  
  
  FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
  ggsave(FileName1, p4_2, width = (5+ gnum)*w, height =15 ,limitsize = FALSE)
  FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
  ggsave(FileName2, p4_2, width = (5+ gnum)*w, height =15,limitsize = FALSE)
  
  FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
  ggsave(FileName1, p4_1, width = (5+ gnum)*w, height =15 ,limitsize = FALSE)
  FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
  ggsave(FileName2, p4_1, width = (5+ gnum)*w, height =15 ,limitsize = FALSE)
  
  FileName <- paste(barpath,"/a2_",j,"_bar_data",".csv", sep = "")
  write.csv(databar,FileName,quote = F)
}

detach("package:ggalluvial")


#---共有微生物特有微生物
#---flower plot#-----
flowpath = paste(otupath,"/flowplot/",sep = "")
dir.create(flowpath)


source("E:\\Shared_Folder\\Function_local\\R_function\\micro/ggflowerplot.R")
p0_1 <- ggflower(ps = ps ,
                 # rep = 1,
                 group = "ID",
                 start = 1, # 风车效果
                 m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                 a = 0.2, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "yellow",
                 N = 0.5
) 

p0_1 

# p + scale_fill_brewer(palette = "Paired")
FileName1 <- paste(flowpath,"ggflowerID.pdf", sep = "")
ggsave(FileName1, p0_1, width = 8, height = 8)
FileName2 <- paste(flowpath,"ggflowerID.jpg", sep = "")
ggsave(FileName2, p0_1, width = 8, height = 8 )



p0_2 <- ggflower(ps = ps,
                 # rep = 1,
                 group = "Group",
                 start = 1, # 风车效果
                 m1 = 1.8, # 花瓣形状，方形到圆形到棱形，数值逐渐减少
                 a = 0.3, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "yellow",
                 N = 0.1
) + scale_fill_manual(values = colset1) 
p0_2

FileName1 <- paste(flowpath,"ggflowerGroup.pdf", sep = "")
ggsave(FileName1, p0_2, width = 14, height = 14)
FileName2 <- paste(flowpath,"ggflowerGroup.jpg", sep = "")
ggsave(FileName2, p0_2, width = 14, height = 14 )


# ggplot升级版本韦恩图和Upset#-------
library(grid)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven.Upset.gg.R")

if (gnum < 6) {
  Venpath = paste(otupath,"/Ven_Upset_super/",sep = "")
  dir.create(Venpath)
  
  library(ggVennDiagram)
  res = Ven.Upset(ps =  ps,
                  group = "Group",
                  N = 0.5,
                  size = 3)
  
  p1 = res[[1]]
  p2 = res[[2]]
  
  filename3 <- paste(Venpath,"Ven_gg.pdf", sep = "")
  ggsave(filename3, p1, width = 8, height = 8)
  filename3 <- paste(Venpath,"Upset_gg.pdf", sep = "")
  ggsave(filename3, p2, width = 8, height = 8)
}

#---Ven-Upset#----------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/VenSeper.R")
source("E:/Shared_Folder/Function_local/R_function/micro/barMainplot.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven-Upset.R")
# j = "Genus"
group = "Group"
ps_Ven = ps
# BiocManager::install("VennDiagram")
# otutab = as.data.frame(otu_table(ps))

map = as.data.frame(phyloseq::sample_data(ps_Ven))
gnumven <- map[,group] %>% unique() %>% dim()

if (gnumven[1] < 6) {
  
  
  Venpath = paste(otupath,"/Ven_Upset_super/",sep = "")
  dir.create(Venpath)
  
  
  result = VenUpset(ps = ps_Ven,
                    group = group,
                    path = Venpath
  )
  
  #---每个部分
  result = VenSeper(ps = ps_Ven,
                    path = Venpath,
                    group = group,
                    j = j,
                    Top = 10
                    
  )
  # 提取韦恩图中全部部分的otu极其丰度做门类柱状图
  p7_1 <- result[[1]]
  #每个部分序列的数量占比，并作差异
  p8 <- result[[2]]
  # 每部分的otu门类冲积图
  p7_2 <- result[[3]]
  
  
  FileName <- paste(Venpath,j,"count_Facet_ven", ".pdf", sep = "")
  ggsave(FileName, p7_1, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"diff_count_box", ".pdf", sep = "")
  ggsave(FileName, p8, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"count_Facet_ven_flow", ".pdf", sep = "")
  ggsave(FileName, p7_2, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"count_Facet_ven", ".jpg", sep = "")
  ggsave(FileName, p7_1, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"diff_count_box", ".jpg", sep = "")
  ggsave(FileName, p8, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"count_Facet_ven_flow", ".jpg", sep = "")
  ggsave(FileName, p7_2, width = 15, height = 12)
  
}




#--维恩网络#-------
library(ggClusterNet)
library(phyloseq)
biospath = paste(otupath,"/biospr_network_Ven/",sep = "")
dir.create(biospath)


N = 0.5
result = ggClusterNet::div_network(ps)
edge = result[[1]]
data = result[[3]]

result <- ggClusterNet::div_culculate(table = result[[3]],distance = 1.1,distance2 = 1.5,distance3 = 1.3,order = FALSE)
# result <- div_culculate(table = result[[3]],distance = 1,distance2 = 1.2,distance3 = 1.1,order = FALSE)
edge = result[[1]]

plotdata = result[[2]]

#--这部分数据是样本点数据
groupdata <- result[[3]]
# table(plotdata$elements)
node =  plotdata[plotdata$elements == unique(plotdata$elements), ]

otu_table = as.data.frame(t(ggClusterNet::vegan_otu(ps)))
tax_table = as.data.frame(ggClusterNet::vegan_tax(ps))
res = merge(node,tax_table,by = "row.names",all = F)
row.names(res) = res$Row.names
res$Row.names = NULL
plotcord = res

xx = data.frame(mean  =rowMeans(otu_table))

plotcord = merge(plotcord,xx,by = "row.names",all = FALSE)
head(plotcord)
# plotcord$Phylum
row.names(plotcord) = plotcord$Row.names
plotcord$Row.names = NULL
head(plotcord)
library(ggrepel)
head(plotcord)
plotcord$Species[!plotcord$Species %in%c("Desulfatiglans_parachlorophenolica",
                                         "Sulfuricella_denitrificans",
                                         "Thiobacillus_thiophilus"
)] = ""


head(groupdata)
p = ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
                            data = edge, size = 0.3,color = "yellow") +
  geom_point(aes(X1, X2,fill = Species,size =mean ),pch = 21, data = plotcord) +
  geom_point(aes(X1, X2),pch = 21, data = groupdata,size = 5,fill = "blue",color = "black") +
  geom_text(aes(X1, X2,label = elements ), data = groupdata,hjust = 1,vjust = -1) +
  theme_void()

p

filename = paste(biospath,"/","biostr_Ven_network.species.several.pdf",sep = "")
ggsave(filename,p,width = (15),height = (12))
filename = paste(biospath,"/","biostr_Ven_network.jpg",sep = "")
ggsave(filename,p,width = (15),height = (12))

detach("package:ggClusterNet")
detach("package:phyloseq")


# 差异微生物#--------

source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")

res = EdgerSuper(ps = ps01,
                 group  = "Group",
                 artGroup = NULL,
                 j = "Genus",
                 path = diffpath.2
)
head(res)

id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`CK-SDlevel` != "nosig") %>%
  arrange(desc(`CK-SDlogFC`))
id = id$ID




source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\MicroMachine_learning.R")

library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
library(ggClusterNet)
library(phyloseq)
ps = ps02
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)
i= 1
for (i in 1:length(aaa[1,])) {
  matpath = paste(otupath,"/Machine_learing.",paste(aaa[,i][1],aaa[,i][2],sep = "."),sep = "")
  dir.create(matpath )
  
  pst = ps %>% subset_samples.wt("Group",group) %>%
    filter_taxa(function(x) sum(x ) > 10, TRUE)
  randomforest.wt(
    pst = pst,
    ROC = F,
    rfcv = F,
    optimal = 50,
    matpath=matpath)
}







#----热图和气泡图展示基因--------

heatpath = paste(otupath,"/heapmap_boplot.micro/",sep = "")
dir.create(heatpath)

#--注意map文件中一定要有ID列
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")
map = phyloseq::sample_data( ps02)
map$ID = row.names(map)
phyloseq::sample_data( ps02) = map


result <- Microheatmap(ps_rela = ps02,id = id,col_cluster = FALSE)

p1 <- result[[1]] 
p1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p2 <- result[[2]]
p2
filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
ggsave(filename,p1,width = 24,height = 15)

filename = paste(heatpath,"Topggbubble.pdf",sep = "")
ggsave(filename,p2,width = 14,height = 25)

#  全部的更新都可以使用

#  新网络分析-2023年末更新#--------



netpath = paste(otupath,"/network.new/",sep = "")
dir.create(netpath)

rank.names(ps)
library(ggrepel)
library(igraph)
detach("package:MicrobiotaProcess")

# 8.1 网络分析主函数#--------
tab.r = network.pip(
  ps = ps,
  N = 980,
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
  fill = "Phylum",
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

ggsave(paste0(netpath,"plot.network.pdf"),p0,width = 22,height = 5,limitsize = FALSE)
ggsave(paste0(netpath,"plot.network2.pdf"),p0,width = 36,height = 10,limitsize = FALSE)

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
                 fill = Phylum,
                 size = igraph.degree),
             pch = 21, data = node,color = "gray40") +
  facet_wrap(.~ label,scales="free_y",nrow = 1) +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
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



#-8.94 网络稳定性-计算负相关的比例#-----
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




#-8.95 网络稳定性-群落稳定性-只有pair样本使用#-----
treat = ps %>% sample_data()
treat$pair = paste( "A",c(rep(1:6,3)),sep = "")
# head(treat)
sample_data(ps) = treat
#一般性的处理，没有时间梯度的，这里设定time为F，意味着每两个群落的组合都进行比较
res5 = community.stability( ps = ps,
                            corg = cor,
                            time = FALSE)
p6 = res5[[1]]
p6
dat7 = res5[[2]]

path = paste(netpath,"/community.stability/",sep = "")
fs::dir_create(path)

write.csv(dat7,
          paste(path,"community.stability.data.csv",sep = ""))
ggsave(paste(path,"community.stability..boxplot.pdf",sep = ""),  p6,width = 4,height = 4)



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


# 参数文件; Tao Wen, 2023.8.28

#--输入和设定文件
# ps0 = base::readRDS("./data/dataNEW/ps.rds")
ps = base::readRDS("./data.ps.meta.ming/ps_transp.rds")
# ps = base::readRDS("./ps_transp.Majie.rds")

path.id = "mRNA.wheat"
# path.id = "mRNA.sesame"

library(tidyverse)
library(phyloseq)
library(ggClusterNet)



map = sample_data(ps)
head(map)
# map$Group = gsub("[?]","",map$treatment)
# map$map$Group  %>% unique()
head(map)
map$Group = paste(map$drought,sep = ".")
map$Group = gsub("MD","SD",map$Group)

map$Group  %>% unique()



sample_data(ps) = map

#--最终确定的phyloseq对象定义为ps

# 
# ps = ps %>%
#   tax_glom_wt.3(ranks = "gene_id") %>%
#   filter_taxa(function(x) sum(x ) > 10, TRUE)
# ps

tax = as.data.frame(vegan_tax(ps))
head(tax)
tax$KO_id = tax$KO_id  %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(ps)  = tax_table(as.matrix(tax))
ps1 <- ps %>% subset_taxa.wt("KO_id" ,"" ,T) %>% 
  tax_glom_wt("KO_id")
ps1
#  合并kegg#--------
otu = ps1 %>% vegan_otu() %>% t() %>%
  as.data.frame()
head(otu)
saveRDS(ps1,"ps.trans.fj.kegg.rds")
#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
Top = 20
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c("CK","T.guizhouense","F.oxysporum")


# 设定排序顺序--ID
map = sample_data(ps)
map$ID = row.names(map)
map = map %>% 
  as.tibble() %>%
  dplyr::arrange(desc(Group))
axis_order.s = map$ID
sam.n = axis_order.s %>% length()
#-主题--颜色等
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




# 构建保存结果文件夹#---------
result<- dir.amp(smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path


#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/",path.id,a,"/",sep = "");otupath
dir.create(otupath)
print(otupath)


#-1.1排序分析beta-diversity#----
betapath = paste(otupath,"/ordination/",sep = "")
dir.create(betapath)

# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA tsne 

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

methodlist = c("PCoA")

# methodlist = c("LDA")
for (method in methodlist) {
  result = BetaDiv(ps = ps1, group = "Group",
                   dist = "bray",
                   method = method, 
                   Micromet = "anosim", 
                   pvalue.cutoff = 0.05,
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
# pairResult = pairMicroTest(ps = ps, Micromet = "adonis")
# FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
# write.csv(pair,FileName)


# 3.1 EdgeR-Desep2差异分析#-------
diffpath = paste(otupath,"/diff_tax/",sep = "")
dir.create(diffpath)



diffpath.1 = paste(diffpath,"/DEsep2/",sep = "")
dir.create(diffpath.1)
# 准备脚本
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\DESep2_micro.R")
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Plot.CompareWithCK.R",encoding = "utf-8")

res = DESep2_Meta2(ps = ps,
                   group  = "Group",
                   artGroup =NULL,
                   j = "gene",
                   path = diffpath.1
)

head(res)
tail(res)
filename = paste(diffpath.1,"/","_","OTU","_","DESep2_all.csv",sep = "")
write.csv(res,filename,quote = F)


diffpath.2 = paste(diffpath,"/EDgeR/",sep = "")
dir.create(diffpath.2)

res = EdgerSuper(ps = ps1,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)
head(res)
filename = paste(diffpath.2,"/","_","OTU","_","edger_all.csv",sep = "")
write.csv(res,filename)



#--edger--曼哈顿图绘制#-------
source("./coding.mrna/differ.mahadun.R")

diffpath = paste(otupath,"/diff_Manhattan/",sep = "")
dir.create(diffpath)

edge_Manhattan(
  ps = ps1,
  pvalue = 0.05,
  lfc = 0,
  diffpath = diffpath 
)

head(res)
id = res %>% filter(`CK-SDlevel` != "nosig") %>% arrange(desc(`CK-SDlogFC`)) %>% .$species



#----热图和气泡图展示变异系数大的微生物丰度--------

heatpath = paste(otupath,"/heapmap_boplot/",sep = "")
dir.create(heatpath)

#--注意map文件中一定要有ID列
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")



map = phyloseq::sample_data(ps1)
map$ID = row.names(map)
phyloseq::sample_data(ps1) = map

  
  result <- Microheatmap(ps_rela = ps1,id = id,col_cluster = FALSE)
  
  p1 <- result[[1]] 
  p1
  # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
  p2 <- result[[2]]
  p2
  filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
  ggsave(filename,p1,width = 14,height = 35)
  
  filename = paste(heatpath,"Topggbubble.pdf",sep = "")
  ggsave(filename,p2,width = 14,height = 35)


  #6.3 基于 Mkegg 通路合并#-------
  #匹配kegg基因名称
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
  
  
  ko = ps1 %>% vegan_otu() %>% t() %>%
    as.data.frame() %>% rownames_to_column("KO")
  head(ko)
  
  tem3 = ko %>% left_join(PATH_ID_NAME,by = c("KO" = "KO")) %>%
    tidyfst::filter_dt(KEGGID != "")
  head(tem3)
  dat = read.delim(
    "E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/ko00001.tsv",
                   header = F
  )
  
  head(dat)
  dat = dat[-1,]
  colnames(dat)  = c("SuperClass","Class","Subclass","kegg")
  
  dat$gene = dat$kegg %>%strsplit("[ ]") %>% 
    sapply(`[`, 1)

  tem4 = tem3 %>%left_join(dat, by = c("KO" = "gene"))
  head(tem4)
  
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
  head(norm2)
  
  tem5 = tem4 %>% left_join(norm2,by = c("KO" = "ID"))
  head(tem5)
  dim(tem5)
  filename = paste(otupath,"kegg.database.all.info.csv",sep = "")
  write_csv(tem5,filename)
  
  
  
  
  tem2 = ko %>% left_join(PATH_ID_NAME,by = c("KO" = "KO")) %>%
    tidyfst::filter_dt(KEGGID != "") %>%
    tidyfst::summarise_vars(is.numeric,sum,by =c("KEGGID","DESCRPTION")) %>%
    as.data.frame()
  head(tem2)
  colnames(tem2)
  
  otu = tem2[,sample_names(ps)]
  tax = data.frame(MDESCRPTION = tem2$DESCRPTION,MDESCRPTION2 = tem2$DESCRPTION,row.names = tem2$KEGGID)
  # tem2$DESCRPTION = gsub("[,]",".",tem2$DESCRPTION)
  # tem2$DESCRPTION = gsub("[-]",".",tem2$DESCRPTION)
  # tem2$DESCRPTION = gsub("[+]",".",tem2$DESCRPTION)
  # tem2$DESCRPTION = gsub("[/]",".",tem2$DESCRPTION)
  # tem2$DESCRPTION = gsub("[(]",".",tem2$DESCRPTION)
  # tem2$DESCRPTION = gsub("[)]",".",tem2$DESCRPTION)
  # tem2$DESCRPTION = gsub("[=>]",".",tem2$DESCRPTION)
  # tem2$DESCRPTION = gsub("[ ]","_",tem2$DESCRPTION)
  # tem2$DESCRPTION %>% unique()
  
  
  
  row.names(otu) = tem2$KEGGID
  head(otu)
  
  ps2 = phyloseq(
    otu_table(as.matrix(otu),taxa_are_rows = TRUE),
    tax_table(as.matrix(tax)),
    sample_data(ps)
  )
  saveRDS(ps2,"./data.ps.meta.ming//ps_kegg_Mfunction.rds")

#  差异Mkegg通路分析#-------
  
  diffpath.2 = paste(otupath,"/EDgeR.Mkegg/",sep = "")
  dir.create(diffpath.2)
  
  res = EdgerSuper(ps = ps2,group  = "Group",artGroup = NULL,
                   j = "gene",
                   path = diffpath.2
  )
  
  head(res)
  
  filename = paste(diffpath.2,"/","_","OTU","_","edger_all.csv",sep = "")
  write.csv(res,filename)
  
  id = res %>% filter(`CK-SDlevel` != "nosig") %>%
    arrange(desc(`CK-SDlogFC`))
  id = row.names(id)
  
  
  #----热图和气泡图展示变异系数大的微生物丰度--------
  
  heatpath = paste(otupath,"/heapmap_boplot.Mkegg/",sep = "")
  dir.create(heatpath)
  
  #--注意map文件中一定要有ID列
  source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")
  map = phyloseq::sample_data(ps2)
  map$ID = row.names(map)
  phyloseq::sample_data(ps2) = map
  
  
  result <- Microheatmap(ps_rela = ps2,id = id,col_cluster = FALSE)
  
  p1 <- result[[1]] 
  p1
  # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
  p2 <- result[[2]]
  p2
  filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
  ggsave(filename,p1,width = 14,height = 25)
  
  filename = paste(heatpath,"Topggbubble.pdf",sep = "")
  ggsave(filename,p2,width = 14,height = 25)
  
  
  
  #  7.1 按照 reaction 合并基因#------
  
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
  
  ko = ps1 %>% vegan_otu() %>% t() %>%
    as.data.frame() %>% rownames_to_column("KO")
  head(ko)
  dat$ID = gsub("ko:","",dat$ID)
  dat$reaction = gsub("rn:","",dat$reaction)
  
  
  tem2 = ko %>% left_join(dat,by = c("KO" = "ID")) %>%
    tidyfst::filter_dt(reaction != "") %>%
    tidyfst::summarise_vars(is.numeric,sum,by =c("reaction")) %>%
    as.data.frame()
  
  colnames(tem2)
  
  otu = tem2[,sample_names(ps)]
  row.names(otu) = tem2$reaction
  head(otu)
  
  
  dat2 = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/reaction.txt",
                    header = F
  )
  head(dat2)
  colnames(dat2) = c("reaction","DESCRPTION")
  
  tax = data.frame(row.names = dat2$reaction,DESCRPTION = dat2$DESCRPTION,DESCRPTION2 = dat2$DESCRPTION)
  head(tax)
  
  ps3 = phyloseq(
    otu_table(as.matrix(otu),taxa_are_rows = TRUE),
    tax_table(as.matrix(tax)),
    sample_data(ps)
  )
  
  fs::dir_create("./data.ps.meta.ming/")
  saveRDS(ps3,"./data.ps.meta.ming//ps_kegg_function.reaction.rds")

  
  
  #  差异reaction通路分析#-------
  
  diffpath.2 = paste(otupath,"/EDgeR.reaction/",sep = "")
  dir.create(diffpath.2)
  
  res = EdgerSuper(ps = ps3,group  = "Group",artGroup = NULL,
                   j = "gene",
                   path = diffpath.2
  )
  
  head(res)
  
  filename = paste(diffpath.2,"/","_","OTU","_","edger_all.csv",sep = "")
  write.csv(res,filename)
  
  id = res %>% filter(`CK-SDlevel` != "nosig") %>%
    arrange(desc(`CK-SDlogFC`))
  id = row.names(id)
  
  
  #----热图和气泡图展示变异系数大的微生物丰度--------
  
  heatpath = paste(otupath,"/heapmap_boplot.rection/",sep = "")
  dir.create(heatpath)
  
  #--注意map文件中一定要有ID列
  source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")
  map = phyloseq::sample_data(ps2)
  map$ID = row.names(map)
  phyloseq::sample_data(ps2) = map
  
  
  result <- Microheatmap(ps_rela = ps3,id = id,col_cluster = FALSE)
  
  p1 <- result[[1]] 
  p1
  # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
  p2 <- result[[2]]
  p2
  filename = paste(heatpath,"/","Topggheatmap.pdf",sep = "")
  ggsave(filename,p1,width = 14,height = 30)
  
  filename = paste(heatpath,"Topggbubble.pdf",sep = "")
  ggsave(filename,p2,width = 14,height = 30)
  
  
  
  
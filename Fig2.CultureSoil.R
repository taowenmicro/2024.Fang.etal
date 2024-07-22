
# 图2.4 养土微生物结果绘制#--------
library(tidyverse)
library(ggClusterNet)
library(phyloseq)

dat = readxl::read_xlsx("./data/randomforest_data.all.abundance(1).xlsx",sheet = 1)
head(dat)
colnames(dat)
dat = dat %>% arrange(`富集标签`)
dat$富集标签

ps01 = readRDS("./data/Cullture.soil/ps_16s.gg.rds")

map= sample_data(ps01)
head(map)
map$Compounds %>% unique()

data = data.frame(ID =map$Compounds %>% unique(),Group = "enrich")
write_csv(data,"./data/fig2.culture.compound.order.csv")
map$Group = map$Compounds
sample_data(ps01) = map




# ps %>% tax_glom_wt(6)


#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

Top = 120

source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

#  升级颜色挑选方式
set.seed(122)
res = theme.col(ps = ps,method = "no",Top =c(Top + 1),gnum = gnum)
mytheme1 = res$mytheme[sample(c(1:4),size=1)]#可选1-4
mytheme2 = res$mytheme[sample(c(5:8),size=1)]; #可选5-8
colset1 = res$col.group[[sample(c(1:5),size=1)]];#1-5可选
colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
colset3 = res$col.time[[sample(c(1:6),size=1)]];#1-6可选
colset4 = colset3



#  排序设定#-------
 
library(ggalluvial)
map = sample_data(ps)
head(map)
map$Group %>% unique()


#  属水平全部物质展示#-------
#--最终确定的phyloseq对象定义为ps
ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
  subset_taxa.wt("Genus"," g__",TRUE)
tax = ps %>% vegan_tax() %>%
  as.data.frame()
head(tax)
tax$Genus
ps = ps #%>% subset_taxa.wt("Genus"," g__")
j = "Genus"

label = FALSE
sd = FALSE
Top = Top
group  = "Group"
axis_ord = ord
label = TRUE 
sd = FALSE
tran = TRUE


  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  

    psdata = psdata%>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x,na.rm = TRUE)} )
 
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()

  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }

  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  

  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[2]];#1-3可选
  colset2 = res$col.bar[[3]];#1-3可选
  colset2 = res$col.bar[[1]];#1-3可选
  
  
  
  
  p4 <- ggplot(Taxonomies_x , aes(x =  group, y = Abundance, fill = aa, order = aa)) +
    geom_bar(stat = "identity",width = 0.8) +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() +
    guides(fill = guide_legend( ncol = 1, byrow = TRUE))

  # p4
  p4 = p4 + scale_x_discrete(limits = ord) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) +
    scale_fill_manual(values = colset2)
  
  dir.create("./paper/result2/Culture.soil.bar/")
  ggsave("./paper/result2/Culture.soil.bar/89.bar.order.remv.no.120.muicol.6col.pdf",p4,
         width = 36,
         height = 46
  )

#  展示特征菌的丰度--89种物质#-------
  
  map= sample_data(ps)
  head(map)
  
  tem = data.frame(row.names = ord,Group = ord)
  map2 = data.frame(row.names = map$ID,ID= map$ID,Concentration  =map$Concentration,
                    Compounds = map$Compounds,
                    Group = map$Compounds
  )
  head(map2)
  
  map3 = tem %>% full_join(map2,by = "Group")
  map3$Group
  head(map3)
  map3$Group2  = map3$Group
  row.names(map3) = map3$ID
  map3$Group = paste0(map3$Group,map3$Concentration)
  
  
  ord2 = map3$Group %>% unique()
  sample_data(ps) = map3
  
  
  
  ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    subset_taxa.wt("Genus"," g__",TRUE)
  
  
  tax = ps %>% vegan_tax() %>%
    as.data.frame() %>% arrange(Genus)
  head(tax)
  tax$Genus %>% unique()
  
  ps = ps %>% subset_taxa.wt("Genus"," g__Bradyrhizobium")
  j = j
  label = FALSE
  sd = FALSE
  Top = Top
  group  = "Group"
  axis_ord = ord
  label = TRUE 
  sd = FALSE
  tran = FALSE
  
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  
  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  
  
  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
  
  
  
  head(Taxonomies_x)
  dat3 = Taxonomies_x %>% left_join(ord.t,by = c("group" = "ID"))
  
  p4 <- ggplot(dat3 , aes(x =  group, y = Abundance, fill = group.y, order = aa)) +
    geom_bar(stat = "identity",width = 0.8,color = "black") +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() 
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))
  
  # p4
  p4 = p4 + scale_x_discrete(limits = ord) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) 
  
  ggsave("./paper/result2/Culture.soil.bar/89.bar.order.慢生根瘤菌2.pdf",p4,
         width = 10,
         height = 12
  )
  
  
  

  
  #  展示特征菌的丰度--89种物质加上浓度#-------
  ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    subset_taxa.wt("Genus"," g__",TRUE)
  map= sample_data(ps)
  head(map)
  
  tem = data.frame(row.names = ord,Group = ord,ord.t$group)
  map2 = data.frame(row.names = map$ID,ID= map$ID,Concentration  =map$Concentration,
                    Compounds = map$Compounds,
                    Group = map$Compounds
  )
  head(map2)
  map2$Group %>% unique()
  map3 = tem %>% full_join(map2,by = "Group")
  map3$Group
  map3$Group2 = map3$Group
  row.names(map3) = map3$ID
  map3$Group = paste0(map3$Group,map3$Concentration)
  map3$Group %>% unique()
  
  ord2 = map3$Group %>% unique()
  sample_data(ps) = map3
  

  tax = ps %>% vegan_tax() %>%
    as.data.frame() %>% arrange(Genus)
  head(tax)
  tax$Genus %>% unique()
  
  ps = ps #%>% subset_taxa.wt("Genus"," g__Bradyrhizobium")
  j = j
  label = FALSE
  sd = FALSE
  Top = Top
  group  = "Group"
  axis_ord = ord
  label = TRUE 
  sd = FALSE
  tran = TRUE
  
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  
  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  
  
  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
  
  
  
  head(Taxonomies_x)
  head(map3)
  dat3 = Taxonomies_x %>% left_join(map3,by = c("group" = "Group"))
  head(dat3)
  colnames(dat3)
  p4 <- ggplot(dat3 , aes(x =  group, y = Abundance, fill = ord.t.group, order = aa)) +
    geom_bar(stat = "identity",width = 0.8) +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() 
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))
  
  # p4
  p4 = p4 + 
    scale_x_discrete(limits = ord2) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) 
  
  ggsave("./paper/result2/Culture.soil.bar210126//89.bar.89.bar.order.ra.sample2.pdf.pdf",p4,
         width = 10,
         height = 32,
         limitsize = FALSE
  )
  
  

  #  展示发病属的丰度--89种物质加上浓度#-------
  ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    subset_taxa.wt("Genus"," g__",TRUE)
  tax = ps %>% vegan_tax() %>%
    as.data.frame() %>% arrange(Genus)
  head(tax)
  tax$Genus %>% unique()
  

  
  map= sample_data(ps)
  head(map)
  
  tem = data.frame(row.names = ord,Group = ord,ord.t$group)
  map2 = data.frame(row.names = map$ID,ID= map$ID,Concentration  =map$Concentration,
                    Compounds = map$Compounds,
                    Group = map$Compounds
  )
  head(map2)
  map2$Group %>% unique()
  map3 = tem %>% full_join(map2,by = "Group")
  map3$Group
  map3$Group2 = map3$Group
  row.names(map3) = map3$ID
  map3$Group = paste0(map3$Group,map3$Concentration)
  map3$Group %>% unique()
  
  ord2 = map3$Group %>% unique()
  sample_data(ps) = map3
  
  
id.D = c(" g__Bradyrhizobium"," g_Candidatus_Koribacter",
    " g__Mesorhizobium"," g__Rhodanobacter"," g__Anaeromyxobacter",
    " g__Candidatus_Solibacter",
    " g__Shinella",
    " g__Ralstonia"
  )
  
  ps = ps %>% 
    # scale_micro() %>%
    subset_taxa.wt("Genus",
                             id.D)
  j = j
  label = FALSE
  sd = FALSE
  Top = Top
  group  = "Group"
  axis_ord = ord
  label = TRUE 
  sd = FALSE
  tran = FALSE
  
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  
  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  
  
  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
  
  
  
  head(Taxonomies_x)
  head(map3)
  dat3 = Taxonomies_x %>% left_join(map3,by = c("group" = "Group"))
  head(dat3)
  colnames(dat3)
  p4 <- ggplot(dat3 , aes(x =  group, y = Abundance, fill = ord.t.group, order = aa)) +
    geom_bar(stat = "identity",width = 0.8) +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() 
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))
  
  # p4
  p4 = p4 + 
    scale_x_discrete(limits = ord2) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) 
  
  ggsave("./paper/result2/Culture.soil.bar/89.bar.order.D.all.pdf",p4,
         width = 10,
         height = 32,
         limitsize = FALSE
  )
  
  
  
  
  
  
  
  
  
  
  #  展示健康属的丰度--89种物质加上浓度#-------
  ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    subset_taxa.wt("Genus"," g__",TRUE)
  tax = ps %>% vegan_tax() %>%
    as.data.frame() %>% arrange(Genus)
  head(tax)
  tax$Genus %>% unique()
  
  map= sample_data(ps)
  head(map)
  
  tem = data.frame(row.names = ord,Group = ord,ord.t$group)
  map2 = data.frame(row.names = map$ID,ID= map$ID,Concentration  =map$Concentration,
                    Compounds = map$Compounds,
                    Group = map$Compounds
  )
  head(map2)
  map2$Group %>% unique()
  map3 = tem %>% full_join(map2,by = "Group")
  map3$Group
  map3$Group2 = map3$Group
  row.names(map3) = map3$ID
  map3$Group = paste0(map3$Group,map3$Concentration)
  map3$Group %>% unique()
  
  ord2 = map3$Group %>% unique()
  sample_data(ps) = map3
  
  
  id.H =c(" g__Candidatus_Nitrososphaera",
          " g__Paenisporosarcina",
          " g__Aeromicrobium",
          " g__Pirellula","g__Anaeromyxobacter",
          " g__Planomicrobium",
          " g__Aquicella",
          " g__Phycicoccus"
          )
  ps = ps %>% 
    # scale_micro() %>%
    subset_taxa.wt("Genus",
                   id.H)
  j = j
  label = FALSE
  sd = FALSE
  Top = Top
  group  = "Group"
  axis_ord = ord
  label = TRUE 
  sd = FALSE
  tran = FALSE
  
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  
  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  
  
  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
  
  
  
  head(Taxonomies_x)
  head(map3)
  dat3 = Taxonomies_x %>% left_join(map3,by = c("group" = "Group"))
  head(dat3)
  colnames(dat3)
  p4 <- ggplot(dat3 , aes(x =  group, y = Abundance, fill = ord.t.group, order = aa)) +
    geom_bar(stat = "identity",width = 0.8) +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() 
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))
  
  # p4
  p4 = p4 + 
    scale_x_discrete(limits = ord2) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) 
  
  ggsave("./paper/result2/Culture.soil.bar/89.bar.order.H.all.pdf",p4,
         width = 10,
         height = 32,
         limitsize = FALSE
  )
  
  
  
  
  
  
  
  
  
  
  #  展示氨氧化属的丰度--89种物质加上浓度#-------
  ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    subset_taxa.wt("Genus"," g__",TRUE)
  tax = ps %>% vegan_tax() %>%
    as.data.frame() %>% arrange(Genus)
  head(tax)
  tax$Genus %>% unique()
  
  map= sample_data(ps)
  head(map)
  
  tem = data.frame(row.names = ord,Group = ord,ord.t$group)
  map2 = data.frame(row.names = map$ID,ID= map$ID,Concentration  =map$Concentration,
                    Compounds = map$Compounds,
                    Group = map$Compounds
  )
  head(map2)
  map2$Group %>% unique()
  map3 = tem %>% full_join(map2,by = "Group")
  map3$Group
  map3$Group2 = map3$Group
  row.names(map3) = map3$ID
  map3$Group = paste0(map3$Group,map3$Concentration)
  map3$Group %>% unique()
  
  ord2 = map3$Group %>% unique()
  sample_data(ps) = map3
  
  id.N = c(" g__Nitratireductor" ," g__Nitrobacteria",
           " g__Nitrosomonas" , 
           " g__Nitrosovibrio"," g__Nitrospira" )
  # id.N = c(" g__Nitrospira" )
  ps = ps %>% 
    # scale_micro() %>%
    subset_taxa.wt("Genus",
                   id.N)
  j = j
  label = FALSE
  sd = FALSE
  Top = Top
  group  = "Group"
  axis_ord = ord
  label = TRUE 
  sd = FALSE
  tran = FALSE
  
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  
  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  
  
  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
  
  
  
  head(Taxonomies_x)
  head(map3)
  dat3 = Taxonomies_x %>% left_join(map3,by = c("group" = "Group"))
  head(dat3)
  colnames(dat3)
  p4 <- ggplot(dat3 , aes(x =  group, y = Abundance, fill = ord.t.group, order = aa)) +
    geom_bar(stat = "identity",width = 0.8) +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() 
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))
  
  # p4
  p4 = p4 + 
    scale_x_discrete(limits = ord2) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) 
  
  ggsave("./paper/result2/Culture.soil.bar/89.bar.order.N.2.pdf",p4,
         width = 10,
         height = 32,
         limitsize = FALSE
  )
  
  
  
  
#  全部8 ML特征菌病原菌的丰度展示#--------
  
  dat = readxl::read_xlsx("./data/Bigdata/最后一遍。总体.xlsx")
  head(dat)
  
  dat2 = dat %>% as.matrix() %>% as.vector() %>% table() %>% as.data.frame() %>%
    
    arrange(desc(Freq)) %>% 
    filter(Freq > 0) 
  head(dat2)
  psF = readRDS("./data/Bigdata/ps_filter_FBwilt_rhi.rds")
  psF
  
  
  psF = psF %>% tax_glom_wt("Genus")
  map = psF %>% sample_data()
  map$Group = map$Group3
  map$Group = map$Group %>%  strsplit( "_") %>% 
    sapply(`[`, 1)
  
  
  map$Group %>% table()
  
  sample_data(psF) = map
  id <- psF %>% 
    ggClusterNet::vegan_otu() %>%
    as.data.frame() %>%
    names()
  ord.x = id[length(id):1]
  
  ps1 <- psF %>%
    subset_taxa(
      row.names(tax_table(psF)) %in% as.character(id)
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
  library(reshape2)
  abun_a = melt(abun,
                id.var = c("ID"),
                variable.name = "id", 
                value.name = "count")
  
  head(abun_a)

  abun_a$iid = rep(paste(1:(length(abun_a$id)/2)),2)
  library(plyr)
  abun_a1 = ddply(abun_a,"iid",transform,percount = count/sum(count)*100)
  head(abun_a1)
  
  id.D = abun_a1 %>% filter(id =="Disease",percount > 50) %>% .$ID
  
  
  ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    subset_taxa.wt("Genus"," g__",TRUE)
  tax = ps %>% vegan_tax() %>%
    as.data.frame() %>% arrange(Genus)
  head(tax)
  tax$Genus %>% unique()
  
  
  
  map= sample_data(ps)
  head(map)
  
  tem = data.frame(row.names = ord,Group = ord,ord.t$group)
  map2 = data.frame(row.names = map$ID,ID= map$ID,Concentration  =map$Concentration,
                    Compounds = map$Compounds,
                    Group = map$Compounds
  )
  head(map2)
  map2$Group %>% unique()
  map3 = tem %>% full_join(map2,by = "Group")
  map3$Group
  map3$Group2 = map3$Group
  row.names(map3) = map3$ID
  map3$Group = paste0(map3$Group,map3$Concentration)
  map3$Group %>% unique()
  
  ord2 = map3$Group %>% unique()
  sample_data(ps) = map3
  
  
  id.HD = paste0(" g__", as.character(id.D))
  
  ps = ps %>% 
    # scale_micro() %>%
    subset_taxa.wt("Genus",
                   id.HD)
  j = j
  label = FALSE
  sd = FALSE
  Top = Top
  group  = "Group"
  axis_ord = ord
  label = TRUE 
  sd = FALSE
  tran = FALSE
  
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  
  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  
  
  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
  
  
  
  head(Taxonomies_x)
  head(map3)
  dat3 = Taxonomies_x %>% left_join(map3,by = c("group" = "Group"))
  head(dat3)
  colnames(dat3)
  p4 <- ggplot(dat3 , aes(x =  group, y = Abundance, fill = ord.t.group, order = aa)) +
    geom_bar(stat = "identity",width = 0.8) +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() 
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))
  
  # p4
  p4 = p4 + 
    scale_x_discrete(limits = ord2) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) 
  
  ggsave("./paper/result2/Culture.soil.bar/89.bar.order.D.all.eight.ML.pdf",p4,
         width = 10,
         height = 32,
         limitsize = FALSE
  )
  
  
  
  
  
  
  
  
  
  #  全部8 ML特征菌健康菌的丰度展示#--------
  
  dat = readxl::read_xlsx("./data/Bigdata/最后一遍。总体.xlsx")
  head(dat)
  
  dat2 = dat %>% as.matrix() %>% as.vector() %>% table() %>% as.data.frame() %>%
    
    arrange(desc(Freq)) %>% 
    filter(Freq > 0) 
  head(dat2)
  psF = readRDS("./data/Bigdata/ps_filter_FBwilt_rhi.rds")
  psF
  
  
  psF = psF %>% tax_glom_wt("Genus")
  map = psF %>% sample_data()
  map$Group = map$Group3
  map$Group = map$Group %>%  strsplit( "_") %>% 
    sapply(`[`, 1)
  
  
  map$Group %>% table()
  
  sample_data(psF) = map
  id <- psF %>% 
    ggClusterNet::vegan_otu() %>%
    as.data.frame() %>%
    names()
  ord.x = id[length(id):1]
  
  ps1 <- psF %>%
    subset_taxa(
      row.names(tax_table(psF)) %in% as.character(id)
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
  library(reshape2)
  abun_a = melt(abun,
                id.var = c("ID"),
                variable.name = "id", 
                value.name = "count")
  
  head(abun_a)
  
  abun_a$iid = rep(paste(1:(length(abun_a$id)/2)),2)
  library(plyr)
  abun_a1 = ddply(abun_a,"iid",transform,percount = count/sum(count)*100)
  head(abun_a1)
  
  id.D = abun_a1 %>% filter(id =="Disease",percount < 50) %>% .$ID
  
  
  ps = ps01 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    subset_taxa.wt("Genus"," g__",TRUE)
  tax = ps %>% vegan_tax() %>%
    as.data.frame() %>% arrange(Genus)
  head(tax)
  tax$Genus %>% unique()
  
  
  
  map= sample_data(ps)
  head(map)
  
  tem = data.frame(row.names = ord,Group = ord,ord.t$group)
  map2 = data.frame(row.names = map$ID,ID= map$ID,Concentration  =map$Concentration,
                    Compounds = map$Compounds,
                    Group = map$Compounds
  )
  head(map2)
  map2$Group %>% unique()
  map3 = tem %>% full_join(map2,by = "Group")
  map3$Group
  map3$Group2 = map3$Group
  row.names(map3) = map3$ID
  map3$Group = paste0(map3$Group,map3$Concentration)
  map3$Group %>% unique()
  
  ord2 = map3$Group %>% unique()
  sample_data(ps) = map3
  
  
  id.HD = paste0(" g__", as.character(id.D))
  
  ps = ps %>% 
    # scale_micro() %>%
    subset_taxa.wt("Genus",
                   id.HD)
  j = j
  label = FALSE
  sd = FALSE
  Top = Top
  group  = "Group"
  axis_ord = ord
  label = TRUE 
  sd = FALSE
  tran = FALSE
  
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- tax_glom_wt(ps = ps
                        ,ranks = j)
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  # otu [1:20,1:20]
  dim(otu)
  otu[is.na(otu)] = 0
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  
  
  if (tran == TRUE) {
    Taxonomies$Abundance = Taxonomies$Abundance * 100
  }
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  
  colnames(Taxonomies)[match(j,colnames(Taxonomies))] <- "aa"
  data = c()
  
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    head(Taxonomies)
    colnames(Taxonomies)
    c <- Taxonomies %>% 
      dplyr::filter(Group == as.character(a))
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  
  Taxonomies = table
  head(Taxonomies)
  
  #按照分组求均值
  
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), stats::sd(Abundance,na.rm = TRUE))
  head(zhnagxu2)
  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  # zhnagxu3 = zhnagxu2
  zhnagxu3 = zhnagxu2 %>% left_join(cc) %>% arrange(allsum)
  
  dat = plyr::ddply(zhnagxu3,'group',
                    dplyr::summarize,label_sd = cumsum(Abundance),
                    label_y = cumsum(Abundance) - 0.5*Abundance)
  
  
  
  dat$aa = zhnagxu3$aa %>% unique()
  head(dat)
  # Taxonomies_x$group = NULL
  # Taxonomies_x = Taxonomies_x %>% arrange(desc(aa))
  # Taxonomies_x$aa = NULL
  # Taxonomies_x$label_y =
  
  Taxonomies_x = zhnagxu3 %>% inner_join(dat)
  
  
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }
  
  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  colours
  
  colset2 = colorRampPalette(c("#FFFFD9" ,"#41B6C4", "#1D91C0" , "#081D58"))(130) #=c("#FFFFCC","#C7E9B4","#7FCDBB", "#41B6C4" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0" ,"#1D91C0")
  
  # colset2 =  colorRampPalette(res$col.time[[sample(c(1:6),size=1)]])(130)
  colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
  
  
  
  head(Taxonomies_x)
  head(map3)
  dat3 = Taxonomies_x %>% left_join(map3,by = c("group" = "Group"))
  head(dat3)
  colnames(dat3)
  p4 <- ggplot(dat3 , aes(x =  group, y = Abundance, fill = ord.t.group, order = aa)) +
    geom_bar(stat = "identity",width = 0.8) +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "") + coord_flip() +
    theme_classic() 
  # guides(fill = guide_legend( ncol = 1, byrow = TRUE))
  
  # p4
  p4 = p4 + 
    scale_x_discrete(limits = ord2) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) 
  
  ggsave("./paper/result2/Culture.soil.bar/89.bar.order.H.all.eight.ML.pdf",p4,
         width = 10,
         height = 32,
         limitsize = FALSE
  )


#  混合健康代谢物和发病代谢物调控
# }
  
  







  
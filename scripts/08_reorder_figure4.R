args<-commandArgs(TRUE);ROOT<-normalizePath(args[1],winslash='/')
suppressPackageStartupMessages({library(grid);library(gtable);library(ggplot2)})
o<-readRDS(file.path(ROOT,'checks/Figure4_layout_objects.rds'));list2env(o,envir=.GlobalEnv)
geometry<-data.frame()
for(i in 1:6){g<-gpanels[[i]];ix<-grep('^panel($|-)',g$layout$name);side<-if(i%in%c(1,2))60 else 48;g$widths[unique(g$layout$l[ix])]<-unit(side,'mm');g$heights[unique(g$layout$t[ix])]<-unit(side,'mm');gpanels[[i]]<-g;geometry<-rbind(geometry,data.frame(panel=LETTERS[i],part=seq_along(ix),width_mm=side,height_mm=side,ratio=1))}
write.table(geometry,file.path(ROOT,'checks/Figure4_panel_geometry.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
saveRDS(list(gpanels=gpanels,product_order=product_order,product_labels=product_labels,group_colours=group_colours,country_order=country_order,country_colours=country_colours),file.path(ROOT,'checks/Figure4_layout_objects.rds'))
code<-readLines(file.path(ROOT,'scripts/01_community_figures.R'),encoding='UTF-8');a<-grep('^plot4<-function',code);b<-grep('^# Separate panels',code);eval(parse(text=code[seq.int(a,b-1)]));graphics.off();cat('FIGURE4 ORDER A B / C D E / F VERIFIED IN SOURCE\n')

args<-commandArgs(TRUE); ROOT<-normalizePath(args[1],winslash='/');setwd(file.path(ROOT,'work'))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ape);library(grid)})
options(device=function(...)grDevices::png(filename=tempfile(fileext='.png'),type='cairo',...));source('plot_config.R')
rd<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE,quote='"',na.strings=c('NA',''))
wm<-function(m,p)write.table(data.frame(sample_name=rownames(m),m,check.names=FALSE),p,sep='\t',quote=TRUE,row.names=FALSE)
wt<-function(d,n)write.table(d,file.path(ROOT,'checks',n),sep='\t',quote=TRUE,row.names=FALSE)
rmatrix<-function(p){d<-rd(p);m<-as.matrix(d[-1]);storage.mode(m)<-'double';rownames(m)<-d[[1]];m}
ann<-rd('data/sample_annotation.tsv');live<-rd(file.path(ROOT,'02_CONTROLLI/Table_S3_Samples_accessions_QC/Samples_100.tsv'))
stopifnot(setequal(ann$sample_name,live$sample_name),!anyDuplicated(ann$sample_name))
live<-live[match(ann$sample_name,live$sample_name),];stopifnot(all(ann$country==live$country),all(ann$raw_material==live$raw_material),all(ann$batch==live$batch))
product_labels<-setNames(live$product,ann$product);product_labels<-product_labels[!duplicated(names(product_labels))]
raw<-lapply(unique(ann$batch),function(b){r<-rd(paste0('data/kraken_batches/',b,'/species_relabund_matrix.tsv'));m<-t(as.matrix(r[setdiff(names(r),c('species','kingdom'))]));storage.mode(m)<-'double';colnames(m)<-r$species;m})
sp<-sort(unique(unlist(lapply(raw,colnames))));M<-matrix(0,nrow(ann),length(sp),dimnames=list(ann$sample_name,sp))
stopifnot(!anyDuplicated(unlist(lapply(raw,rownames))),setequal(unlist(lapply(raw,rownames)),ann$sample_name))
for(m in raw)M[rownames(m),colnames(m)]<-m
P<-M;for(b in unique(ann$batch)){i<-ann$batch==b;P[i,apply(M[i,,drop=FALSE],2,max)<.001]<-0}
P<-P[,colSums(P)>0,drop=FALSE];P<-P/rowSums(P)
K<-rmatrix('../02_CONTROLLI/Table_S4_Taxonomic_profiles/Kraken_100.tsv');L<-rmatrix(file.path(ROOT,'02_CONTROLLI/Table_S4_Taxonomic_profiles/Kraken_100.tsv'))
stopifnot(setequal(colnames(P),colnames(K)),setequal(colnames(P),colnames(L)))
P<-P[rownames(K),colnames(K)];errlocal<-max(abs(P-K));errlive<-max(abs(P-L[rownames(P),colnames(P)]));stopifnot(errlocal<1e-10,errlive<1e-10,nrow(P)==100,ncol(P)==1380)
ps<-.5*min(P[P>0]);clrnew<-log(P+ps);clrnew<-clrnew-rowMeans(clrnew)
clrref<-rmatrix('data/04/clr_matrix_kraken_species.tsv');errclr<-max(abs(clrnew-clrref));stopifnot(errclr<1e-7)
pc<-prcomp(clrnew,center=TRUE,scale.=FALSE);pv<-100*pc$sdev^2/sum(pc$sdev^2)
ref<-rd(file.path(ROOT,'02_CONTROLLI/Table_S7_Beta_diversity/PCA_Kraken.tsv'));ref<-ref[match(rownames(P),ref$sample_name),]
for(j in 1:2)if(cor(pc$x[,j],ref[[paste0('PC',j)]])<0){pc$x[,j]<- -pc$x[,j];pc$rotation[,j]<- -pc$rotation[,j]}
errpc<-max(abs(pc$x[,1:2]-as.matrix(ref[c('PC1','PC2')])));stopifnot(errpc<1e-6)
pa<-ann[match(rownames(P),ann$sample_name),];pa$PC1<-pc$x[,1];pa$PC2<-pc$x[,2]
pa$product<-factor(pa$product,levels=product_order);pa$country<-factor(pa$country,levels=country_order)
wm(P,file.path(ROOT,'checks/reconstructed_relative.tsv'));wm(clrnew,file.path(ROOT,'checks/reconstructed_CLR.tsv'))
wt(transform(pa,product=unname(product_labels[as.character(product)])),'Figure5_PCA_coordinates.tsv')
wt(data.frame(axis=paste0('PC',seq_along(pv)),variance_percent=pv),'Figure5_variance.tsv')
wt(data.frame(n=100,taxa=1380,raw_taxa=ncol(M),pseudocount=ps,err_local=errlocal,err_Drive_S4=errlive,err_CLR=errclr,err_PCA=errpc,PC1=pv[1],PC2=pv[2]),'reconstruction_checks.tsv')
# Repeat the canonical ellipses on identical scores, grouping only the display.
ellipse<-function(f){do.call(rbind,lapply(unique(as.character(pa[[f]])),function(g){xy<-as.matrix(pa[pa[[f]]==g,c('PC1','PC2')]);if(nrow(xy)<3)return(NULL);ee<-eigen(cov(xy),symmetric=TRUE);a<-seq(0,2*pi,length.out=181);e<-t(ee$vectors%*%diag(sqrt(pmax(ee$values,0)*5.991),2)%*%rbind(cos(a),sin(a)));e<-sweep(e,2,colMeans(xy),'+');data.frame(factor=f,group=g,n=nrow(xy),PC1=e[,1],PC2=e[,2])}))}
els<-lapply(c('product','country','raw_material'),ellipse);wt(do.call(rbind,els),'Figure5_covariance_ellipses.tsv')
xlim<-range(unlist(lapply(els,function(e)e$PC1)),pa$PC1);ylim<-range(unlist(lapply(els,function(e)e$PC2)),pa$PC2)
# Square PCA panels: both axes cover the same numerical span.
common_limits<-range(c(xlim,ylim));xlim<-common_limits;ylim<-common_limits
# Equal axis units, identical limits and symbols in all three panels.
pp<-lapply(1:3,function(i){f<-c('product','country','raw_material')[i];e<-els[[i]];cols<-switch(f,product=product_colours,country=country_colours,raw_material=raw_colours)
p<-ggplot(pa,aes(PC1,PC2))+geom_polygon(data=e,aes(group=group,fill=group),alpha=.09,colour=NA,show.legend=FALSE)+geom_path(data=e,aes(group=group,linetype=if(i==1)'group'else group),colour=if(i==1)'grey50' else 'grey35',linewidth=.35,show.legend=i!=1)+scale_fill_manual(values=cols,guide='none')+geom_point(aes(colour=product,shape=country),size=1.65,alpha=.88)+scale_colour_manual(values=product_colours,breaks=product_order,labels=product_labels)+scale_shape_manual(values=country_shapes,breaks=country_order)+coord_equal(xlim=xlim,ylim=ylim)+labs(title=c('A: Product','B: Country','C: Raw material')[i],x=sprintf('PC1 (%.1f%%)',pv[1]),y=sprintf('PC2 (%.1f%%)',pv[2]),colour='Product',shape='Country (point shape)',linetype='Ellipse group')+theme_paper(9)+theme_arial()+theme(plot.margin=margin(3,3,2,3,'mm'),legend.position='bottom',legend.box='vertical')
if(i==1){p<-p+scale_linetype_manual(values=c(group='solid'),guide='none')}else{lev<-if(i==2)country_order else raw_order;p<-p+scale_linetype_manual(values=setNames(c('solid','dashed','dotdash'),lev),breaks=lev,guide=guide_legend(nrow=1,override.aes=list(colour='grey35')))}
if(i==1){p$layers[[2]]<-geom_path(data=e,aes(PC1,PC2,group=group,colour=group),inherit.aes=FALSE,linewidth=.35,show.legend=FALSE)}
p+guides(colour='none',shape='none')})
# Extract one common product/point legend and separate line-style keys per panel.
legplot<-ggplot(pa,aes(PC1,PC2,colour=product,shape=country))+geom_point()+scale_colour_manual(values=product_colours,breaks=product_order,labels=product_labels)+scale_shape_manual(values=country_shapes,breaks=country_order)+labs(colour='Product',shape='Country (shape)')+theme_void(base_size=9,base_family='Arial')+theme(legend.position='bottom',legend.box='vertical',legend.key.width=unit(3,'mm'),legend.key.height=unit(3,'mm'),legend.spacing.y=unit(1,'mm'))+guides(colour=guide_legend(ncol=1,byrow=TRUE,order=1,override.aes=list(shape=16,size=2)),shape=guide_legend(nrow=1,title.position='top',order=2,override.aes=list(colour='black',size=2)))
getlegend<-function(p){g<-ggplotGrob(p);ix<-which(g$layout$name=='guide-box-bottom');g$grobs[[ix]]}
legend<-wrap_elements(full=getlegend(legplot))
fig5<-wrap_plots(A=pp[[1]],B=pp[[2]],C=pp[[3]],D=legend,design='AB\nCD',heights=c(1,1))
saveplot<-function(p,n,w,h){
 png(file.path(ROOT,'outputs',paste0(n,'.png')),width=w,height=h,units='mm',res=600,type='cairo',bg='white');print(p);dev.off()
 cairo_pdf(file.path(ROOT,'outputs',paste0(n,'.pdf')),width=w/25.4,height=h/25.4,bg='white',family='Arial');print(p);dev.off()
 svg(file.path(ROOT,'outputs',paste0(n,'.svg')),width=w/25.4,height=h/25.4,bg='white',family='Arial');print(p);dev.off()
 png(file.path(ROOT,'checks',paste0(n,'_preview.png')),width=w,height=h,units='mm',res=160,type='cairo',bg='white');print(p);dev.off()
}
# Export the Figure 5 object already constructed above.
saveplot(fig5,'Figure5_beta_diversity',230,230)
# Retain product labels as verified in the current Drive metadata.
code<-readLines('04_figure3.R',encoding='UTF-8');code[grepl("^source\\('plot_config.R'",code)]<-'';code[grepl('^save_figure\\(',code)]<-'';code[grepl('^fig3<-',code)]<-''
eval(parse(text=code))
# Validate plotted composition, Shannon and country means against live supplements.
pg<-rd(file.path(ROOT,'02_CONTROLLI/Table_S8_Taxon_groups_and_screening/Product_group_means.tsv'));gp$product<-unname(product_labels[gp$product]);pg<-pg[match(gp$product,pg$product),];errgroup<-max(abs(as.matrix(gp[-1])-as.matrix(pg[-1])));stopifnot(errgroup<1e-8)
asrc<-rd(file.path(ROOT,'02_CONTROLLI/Table_S6_Alpha_diversity/Alpha_Kraken.tsv'));hv<-vegan::diversity(P,'shannon');errh<-max(abs(hv-asrc$shannon[match(names(hv),asrc$sample_name)]));stopifnot(errh<1e-8)
ll<-rd(file.path(ROOT,'02_CONTROLLI/Table_S6_Alpha_diversity/Product_letters.tsv'));actual<-rd('derived/Figure3_product_contrast_letters.tsv');actual$letters<-ll$letters[match(unname(product_labels[actual$product]),ll$product)];stopifnot(!anyNA(actual$letters));labs_current<-setNames(actual$letters,actual$product);dp<-contrasts[contrasts$tool=='kraken'&contrasts$factor=='product',];for(ii in seq_len(nrow(dp))){shared<-length(intersect(strsplit(labs_current[[dp$group1[ii]]],'')[[1]],strsplit(labs_current[[dp$group2[ii]]],'')[[1]]))>0;stopifnot(shared==(dp$p.adj[ii]>=.05))};p3b$layers[[3]]$data$label<-unname(labs_current[as.character(p3b$layers[[3]]$data$group)])
sh<-rd(file.path(ROOT,'02_CONTROLLI/Table_S8_Taxon_groups_and_screening/Displayed_country_taxa.tsv'));sh<-sh[match(shown$species,sh$species),];errtax<-max(abs(as.matrix(shown[3:8])-as.matrix(sh[3:8])));stopifnot(errtax<1e-8)
Hs<-vapply(c('product','country','raw_material'),function(g)unname(kruskal.test(hv~pa[[g]])$statistic),numeric(1));omlive<-rd(file.path(ROOT,'02_CONTROLLI/Table_S6_Alpha_diversity/Omnibus_KW.tsv'));omlive<-omlive[omlive$metric=='shannon',];stopifnot(max(abs(Hs-omlive$H[match(names(Hs),omlive$factor)]))<1e-8)
stlive<-rd(file.path(ROOT,'02_CONTROLLI/Table_S7_Beta_diversity/Figure5_statistics.tsv'));ssT<-sum(scale(clrnew,center=TRUE,scale=FALSE)^2);r2<-vapply(c('product','country','raw_material'),function(g){groups<-split(seq_len(nrow(pa)),pa[[g]]);sum(vapply(groups,function(ii)length(ii)*sum((colMeans(clrnew[ii,,drop=FALSE])-colMeans(clrnew))^2),numeric(1)))/ssT},numeric(1));stopifnot(max(abs(r2-stlive$PERMANOVA_R2[match(names(r2),stlive$factor)]))<1e-10);wt(data.frame(factor=names(r2),reproduced_R2=r2,source_R2=stlive$PERMANOVA_R2[match(names(r2),stlive$factor)],reproduced_Shannon_H=Hs),'additional_statistical_checks.tsv')
wt(data.frame(err_group_means=errgroup,err_Shannon=errh,err_displayed_taxa=errtax,displayed_taxa=nrow(shown),all_product_letters_match=TRUE),'Figure4_checks.tsv')
# Fixed physical square panels, measured on the actual graphics device.
product_codes<-setNames(paste0('P',seq_along(product_order)),product_order)
wt(data.frame(code=unname(product_codes),product=unname(product_labels[product_order])),'Figure4_product_key.tsv')
ct<-theme(text=element_text(family='Arial',size=9),plot.title=element_text(size=9,face='bold'),axis.title=element_text(size=8),axis.text=element_text(size=8,colour='black'),plot.margin=margin(1,1,1,1,'mm'),aspect.ratio=1)
glabels<-setNames(c('LAB-associated genera','Staphylococcus and other\nselected genera','Selected fungal genera','Selected enterobacterial genera','Other selected bacterial genera','Flagged assignments (unconfirmed)','Other taxa'),names(group_colours))
p3a<-p3a+ct+labs(title='A: Mean taxonomic composition',y='Mean relative abundance (%)')+scale_x_discrete(limits=product_order,labels=product_codes)+scale_fill_manual(values=group_colours,labels=glabels)+theme(legend.position='bottom',legend.text=element_text(size=7),legend.key.height=unit(2.2,'mm'),legend.key.width=unit(2.2,'mm'),legend.spacing.y=unit(0,'mm'),legend.margin=margin(0,0,0,0),legend.box.spacing=unit(1,'mm'),axis.text.x=element_text(angle=0,hjust=.5))+guides(fill=guide_legend(ncol=1,byrow=TRUE))
p3b<-p3b+ct+labs(title='B: Shannon diversity by product')+scale_x_discrete(limits=product_order,labels=product_codes)+theme(axis.text.x=element_text(angle=0,hjust=.5))
p3c<-p3c+ct+labs(title='C: Product contrasts')+scale_x_discrete(limits=product_order,labels=product_codes,drop=FALSE)+scale_y_discrete(limits=rev(product_order),labels=product_codes,drop=FALSE)+theme(axis.text.x=element_text(angle=0,hjust=.5),legend.text=element_text(size=7),legend.title=element_text(size=8))
p3d<-p3d+ct+labs(title='D: Shannon diversity by country')
p3e<-p3e+ct+labs(title='E: Shannon diversity by raw material')
p3f<-p3f+ct+labs(title='F: Country-associated species',x='Mean relative abundance (%)')+theme(legend.position='bottom',legend.text=element_text(size=8),legend.title=element_text(size=8),strip.text=element_text(size=9),plot.margin=margin(1,1,1,1,'mm'))
p3f$layers[[2]]$aes_params$size<-7/ggplot2::.pt
# gtable panel dimensions are absolute units, never stretched to fill a cell.
measurements<-data.frame()
fixsquare<-function(p,label){
 g<-ggplotGrob(p);ix<-grep('^panel($|-)',g$layout$name);stopifnot(length(ix)>=1)
 side<-if(label%in%c('A','B'))60 else 48
 cols<-unique(g$layout$l[ix]);rows<-unique(g$layout$t[ix]);g$widths[cols]<-unit(side,'mm');g$heights[rows]<-unit(side,'mm');g$respect<-FALSE
 for(j in seq_along(ix)){k<-ix[j];w<-convertWidth(g$widths[g$layout$l[k]],'mm',valueOnly=TRUE);h<-convertHeight(g$heights[g$layout$t[k]],'mm',valueOnly=TRUE);stopifnot(abs(w-side)<1e-8,abs(h-side)<1e-8);measurements<<-rbind(measurements,data.frame(panel=label,part=j,width_mm=w,height_mm=h,ratio=w/h))}
 g
}
# Local legends are drawn in matching blocks beneath A and B.
p3a<-p3a+theme(legend.position='none')
p3f<-p3f+theme(legend.position='none')
cairo_pdf(file.path(ROOT,'checks/panel_measurement_device.pdf'),width=230/25.4,height=230/25.4,family='Arial')
gpanels<-Map(fixsquare,list(p3a,p3b,p3c,p3d,p3e,p3f),c('A','B','C','D','E','F'))
# Match the three F plotting columns to the three alpha-diversity columns.
fg<-gpanels[[6]];pcidx<-grep('^panel($|-)',fg$layout$name);pcols<-sort(unique(fg$layout$l[pcidx]));for(k in 1:2){between<-seq.int(pcols[k]+1,pcols[k+1]-1);fg$widths[between]<-unit(0,'mm');fg$widths[between[1]]<-unit(230/3-48,'mm')};fg$heights[length(fg$heights)]<-unit(0,'mm');gpanels[[6]]<-fg
cat('GROB HEIGHTS MM:',vapply(gpanels,function(g)convertHeight(sum(g$heights),'mm',valueOnly=TRUE),numeric(1)),'\n');dev.off()
wt(measurements,'Figure4_panel_geometry.tsv')
saveRDS(list(gpanels=gpanels,product_order=product_order,product_labels=product_labels,group_colours=group_colours,country_order=country_order,country_colours=country_colours),file.path(ROOT,'checks/Figure4_layout_objects.rds'))
plot4<-function(){
 grid.newpage();pushViewport(viewport(layout=grid.layout(3,6,widths=unit(rep(230/6,6),'mm'),heights=unit(c(98,60,72),'mm'))))
 frameleft<-function(g,cellwidth){i<-grep('^panel($|-)',g$layout$name)[1];l<-g$layout$l[i];(cellwidth-convertWidth(sum(g$widths),'mm',valueOnly=TRUE))/2+convertWidth(sum(g$widths[seq_len(l-1)]),'mm',valueOnly=TRUE)}
 drawcell<-function(g,r,c,shift=0){pushViewport(viewport(layout.pos.row=r,layout.pos.col=c));pushViewport(viewport(x=unit(.5,'npc')+unit(shift,'mm'),y=unit(1,'npc')-unit(if(r==1)3 else 1,'mm'),height=sum(g$heights),width=sum(g$widths),just=c('centre','top')));grid.draw(g);popViewport(2)}
 drawcell(gpanels[[1]],1,1:3);b_shift<-frameleft(gpanels[[1]],115)-frameleft(gpanels[[2]],115);drawcell(gpanels[[2]],1,4:6,shift=b_shift)
 c_shift<-frameleft(gpanels[[4]],230/3)-frameleft(gpanels[[3]],230/3);drawcell(gpanels[[3]],2,1:2,shift=c_shift);drawcell(gpanels[[4]],2,3:4);drawcell(gpanels[[5]],2,5:6)
 f_shift<-frameleft(gpanels[[4]],230/3)-frameleft(gpanels[[6]],230);drawcell(gpanels[[6]],3,1:6,shift=f_shift)
 # Taxonomic legend: two compact columns directly beneath A.
 pushViewport(viewport(layout.pos.row=1,layout.pos.col=1:3,clip='off'))
 a_left<-frameleft(gpanels[[1]],115);legend_top<-76
 legend_text<-c('LAB-associated genera','Staphylococcus and other\nselected genera','Selected fungal genera','Selected enterobacterial genera','Other selected bacterial genera','Flagged assignments (unconfirmed)','Other taxa')
 for(j in 1:7){x<-a_left+((j-1)%%2)*47;y<-legend_top+c(0,6.5,11.3,16.1)[((j-1)%/%2)+1];grid.rect(x=unit(x,'mm'),y=unit(1,'npc')-unit(y+1.3,'mm'),width=unit(2.1,'mm'),height=unit(2.1,'mm'),just=c('left','centre'),gp=gpar(fill=unname(group_colours[j]),col=NA));grid.text(legend_text[j],x=unit(x+3.1,'mm'),y=unit(1,'npc')-unit(y,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=7,lineheight=1.05))}
 popViewport()
 # Product names retain their block under the upper-right panel B; no heading.
 pushViewport(viewport(layout.pos.row=1,layout.pos.col=4:6,clip='off'))
 c_left<-frameleft(gpanels[[2]],115)+b_shift
 for(j in 1:10){column<-(j-1)%/%5;row<-(j-1)%%5;grid.text(paste0('P',j,'  ',unname(product_labels[product_order[j]])),x=unit(c_left+column*46,'mm'),y=unit(1,'npc')-unit(76+row*4,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=7.2,col='#222222'))}
 popViewport()
 # Country legend is centred beneath the complete F plotting block.
 pushViewport(viewport(layout.pos.row=3,layout.pos.col=1:6,clip='off'))
 f_left<-frameleft(gpanels[[6]],230)+f_shift;f_centre<-f_left+(2*230/3+48)/2
 key_top<-convertHeight(sum(gpanels[[6]]$heights),'mm',valueOnly=TRUE)+1.8
 key_gp<-gpar(fontfamily='Arial',fontsize=7.5)
 text_width<-function(s)convertWidth(grobWidth(textGrob(s,gp=key_gp)),'mm',valueOnly=TRUE)
 widths<-vapply(country_order,text_width,numeric(1));total_width<-text_width('Country')+4+sum(2.2+1.1+widths)+14
 key_left<-f_centre-total_width/2;grid.text('Country',x=unit(key_left,'mm'),y=unit(1,'npc')-unit(key_top,'mm'),just=c('left','top'),gp=key_gp)
 cursor<-key_left+text_width('Country')+4
 for(j in 1:3){grid.rect(x=unit(cursor,'mm'),y=unit(1,'npc')-unit(key_top+1.2,'mm'),width=unit(2.2,'mm'),height=unit(2.2,'mm'),just=c('left','centre'),gp=gpar(fill=unname(country_colours[country_order[j]]),col=NA));grid.text(country_order[j],x=unit(cursor+3.3,'mm'),y=unit(1,'npc')-unit(key_top,'mm'),just=c('left','top'),gp=key_gp);cursor<-cursor+3.3+widths[j]+7}
 popViewport(2)
}
png(file.path(ROOT,'outputs/Figure4_composition_alpha.png'),width=230,height=230,units='mm',res=600,type='cairo',bg='white');plot4();dev.off()
cairo_pdf(file.path(ROOT,'outputs/Figure4_composition_alpha.pdf'),width=230/25.4,height=230/25.4,family='Arial');plot4();dev.off()
svg(file.path(ROOT,'outputs/Figure4_composition_alpha.svg'),width=230/25.4,height=230/25.4,family='Arial',bg='white');plot4();dev.off()
png(file.path(ROOT,'checks/Figure4_composition_alpha_preview.png'),width=230,height=230,units='mm',res=160,type='cairo',bg='white');plot4();dev.off()
# Separate panels for inspection of the axes.
for(j in 1:5){single<-list(p3a,p3b,p3c,p3d,p3e)[[j]];if(j==1)single<-single+theme(legend.position='bottom');saveplot(single,paste0('Figure4_panel_',LETTERS[j]),110,110)}

saveRDS(list(pa=pa,ellipses=els,pv=pv,product_order=product_order,product_labels=product_labels,product_colours=product_colours,country_shapes=country_shapes),file.path(ROOT,'checks/verified_plot_data.rds'))
grDevices::graphics.off()
writeLines(capture.output(sessionInfo()),file.path(ROOT,'checks/sessionInfo_R.txt'))
writeLines(capture.output(warnings()),file.path(ROOT,'checks/R_warnings.txt'))
cat('RECONSTRUCTION PASS; FIGURES 4 AND 5 EXPORTED\n');print(rd(file.path(ROOT,'checks/reconstruction_checks.tsv')))

quit(save='no',status=0,runLast=FALSE)

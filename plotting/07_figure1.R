# Plot literature summaries from the 282 curated records.
# The Figure 1 entry point supplies metadata, matrix and plot configuration.
stopifnot(nrow(meta)==282,ncol(mat)==282)
world<-read_tsv('data/review/world_map_maps.tsv')
original<-readLines('original_scripts/review/FS1_geo_temporal.R',encoding='UTF-8')
i<-grep('^centro <- list',original);j<-grep('^# ---- per-country',original)
stopifnot(length(i)==1,length(j)==1,j>i)
eval(parse(text=original[i:(j-1)]))
ct<-trimws(meta$country);ct[grepl('turk|rkiye',ct,ignore.case=TRUE)]<-'Turkey';joint<-ct=='India/Nepal'
lev<-names(review_product_colours)
tab<-table(ct[!joint],factor(meta$product[!joint],levels=lev));agg<-data.frame(country=rownames(tab),as.data.frame.matrix(tab),check.names=FALSE)
agg$total<-rowSums(tab)
stopifnot(all(agg$country%in%names(centro)),sum(agg$total)==279)
agg$lon<-vapply(agg$country,function(x)centro[[x]][1],numeric(1));agg$lat<-vapply(agg$country,function(x)centro[[x]][2],numeric(1))
# Preserve the original minimum radius and report it explicitly in the caption.
agg$r<-pmax(1.3,sqrt(agg$total)*.95)
write_tsv(agg,'Figure1_country_counts.tsv')
# Polygon wedges reproduce scatterpie geometry without adding its compiled dependencies.
wedges<-function(d,scale=1){do.call(rbind,lapply(seq_len(nrow(d)),function(i){
  counts<-as.numeric(d[i,lev]);ends<-c(0,cumsum(counts/sum(counts))*2*pi)
  do.call(rbind,lapply(which(counts>0),function(k){angle<-seq(ends[k],ends[k+1],length.out=max(8,ceiling(100*counts[k]/sum(counts))));data.frame(long=c(d$lon[i],d$lon[i]+d$r[i]*scale*sin(angle),d$lon[i]),lat=c(d$lat[i],d$lat[i]+d$r[i]*scale*cos(angle),d$lat[i]),product=lev[k],slice=paste(i,k))}))
}))}
map_panel<-function(d,xlim,ylim,scale=1,ratio=1.25){ggplot()+geom_polygon(data=world,aes(long,lat,group=group),fill='#F2F2F2',colour='grey75',linewidth=.12)+geom_polygon(data=wedges(d,scale),aes(long,lat,group=slice,fill=product),colour='grey30',linewidth=.1,alpha=.92)+scale_fill_manual(values=review_product_colours,limits=lev,name='Product category')+coord_fixed(ratio,xlim=xlim,ylim=ylim,expand=FALSE)+theme_void(base_size=9)+theme(legend.position='bottom',plot.title=element_text(face='bold',size=10),plot.tag=element_text(face='bold',size=12))}
# Map cropped to the sampled latitudes; product and size keys placed in empty ocean areas.
p1b<-map_panel(agg,c(-128,150),c(-50,66),ratio=1)+labs(title='B: Geographical distribution of corpus records')+
  theme(legend.position='inside',legend.position.inside=c(.006,.03),legend.justification=c(0,0),legend.direction='vertical',
        legend.key.size=unit(3.2,'mm'),legend.title=element_text(face='bold',size=8),legend.text=element_text(size=8),
        legend.background=element_rect(fill='white',colour='grey80',linewidth=.2),legend.margin=margin(3,4,3,4),
        plot.title=element_text(face='bold',size=10,margin=margin(b=8)),plot.margin=margin(14,2,2,2))+
  guides(fill=guide_legend(ncol=1))
size_key<-data.frame(n=c(1,10,50),x=c(72,84,101),y=-35)
key_paths<-do.call(rbind,lapply(seq_len(nrow(size_key)),function(i){angle<-seq(0,2*pi,length.out=101);r<-pmax(1.3,sqrt(size_key$n[i])*.95);data.frame(long=size_key$x[i]+r*sin(angle),lat=size_key$y[i]+r*cos(angle),group=i)}))
p1b<-p1b+geom_path(data=key_paths,aes(long,lat,group=group),linewidth=.25)+geom_text(data=size_key,aes(x,y-11,label=n),size=2.8)+annotate('text',x=88,y=-21,label='Records per country',size=2.8,fontface='bold')
yc<-as.data.frame(table(year=as.integer(meta$year),technique=factor(meta$technique,levels=names(technique_colours))));names(yc)[3]<-'n';yt<-aggregate(n~year,yc,sum)
p1a<-ggplot(yc,aes(year,n,fill=technique))+geom_col(width=.82,colour='white',linewidth=.1)+geom_text(data=yt,aes(year,n,label=n),inherit.aes=FALSE,vjust=-.45,size=2.7)+
  scale_fill_manual(values=technique_colours,name='Technique')+scale_y_continuous(expand=expansion(mult=c(0,.24)))+
  labs(x='Publication year',y='Corpus records',title='A: Corpus records by publication year and technique')+theme_paper()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7.5),legend.position='inside',legend.position.inside=c(.015,.96),legend.justification=c(0,1),
        legend.direction='horizontal',legend.background=element_rect(fill='white',colour=NA),legend.title=element_text(size=8))+guides(fill=guide_legend(nrow=1))
pt<-as.data.frame(table(product=factor(meta$base_product,levels=lev),raw=factor(meta$raw_material_group,levels=names(review_raw_colours))));names(pt)[3]<-'n'
tot<-aggregate(n~product,pt,sum);sm<-tapply(meta$smoked,factor(meta$base_product,levels=lev),sum)
tot$smoked<-sm[as.character(tot$product)];stopifnot(sum(tot$n)==282,sum(tot$smoked)==51)
write_tsv(pt,'Figure1_product_raw_counts.tsv');write_tsv(tot,'Figure1_smoking_counts.tsv')
# Horizontal bars keep product names readable; smoked counts follow each name in brackets.
rt<-as.data.frame(table(factor(meta$raw_material_group,levels=names(review_raw_colours))));names(rt)<-c('raw','n')
raw_labels<-setNames(paste0(rt$raw,' (',rt$n,')'),rt$raw);stopifnot(sum(rt$n)==282)
product_smoke_labels<-setNames(paste0(lev,' [',tot$smoked,']'),lev)
p1c<-ggplot(pt,aes(n,product,fill=raw))+geom_col(width=.42,colour='white',linewidth=.2,position=position_stack(reverse=TRUE))+
  geom_text(aes(label=ifelse(n>=30,n,'')),position=position_stack(vjust=.5,reverse=TRUE),size=2.3,colour='white')+
  geom_text(data=tot,aes(n,product,label=n),inherit.aes=FALSE,hjust=-.25,size=2.8,fontface='bold')+
  geom_text(data=tot,aes(0,product,label=product_smoke_labels[as.character(product)]),inherit.aes=FALSE,hjust=0,nudge_x=1,nudge_y=.37,size=2.8)+
  scale_fill_manual(values=review_raw_colours,labels=raw_labels,name='Raw material (total)')+scale_y_discrete(limits=rev(lev),expand=expansion(add=c(.45,.85)))+
  scale_x_continuous(expand=expansion(mult=c(0,.16)))+labs(x='Corpus records',y=NULL,title='C: Records per product category',subtitle='Segments: raw material; [ ] smoked records')+
  theme_paper()+theme(panel.grid.major.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position='inside',legend.position.inside=c(.99,.03),legend.justification=c(1,0),
        legend.background=element_rect(fill='white',colour='grey80',linewidth=.2),legend.title=element_text(size=8),legend.margin=margin(2,3,2,3),plot.subtitle=element_text(size=8))+
  guides(fill=guide_legend(ncol=1))
dj<-vegan::vegdist(t(mat),method='jaccard',binary=TRUE);pco<-cmdscale(dj,k=2,eig=TRUE);ev<-pmax(pco$eig,0);ve<-100*ev/sum(ev)
stopifnot(abs(ve[1]-7.8)<.1,abs(ve[2]-4.5)<.1)
ord<-data.frame(code=meta$code,PC1=pco$points[,1],PC2=pco$points[,2],technique=meta$technique);write_tsv(ord,'review_PCoA_technique.tsv')
p1e<-ggplot(ord,aes(PC1,PC2,colour=technique))+geom_point(size=1.3,alpha=.75)+stat_ellipse(level=.90,linewidth=.45,linetype=2,show.legend=FALSE)+
  scale_colour_manual(values=technique_colours)+labs(x=sprintf('PCo1 (%.1f%%)',ve[1]),y=sprintf('PCo2 (%.1f%%)',ve[2]),colour='Technique',
  title='D: Reported labels by technique',subtitle='Binary Jaccard PCoA; dashed: 90% ellipses')+
  theme_paper()+theme(aspect.ratio=1,legend.position='inside',legend.position.inside=c(.99,.99),legend.justification=c(1,1),legend.background=element_rect(fill='white',colour='grey80',linewidth=.2),legend.margin=margin(2,3,2,3),legend.title=element_text(size=8),plot.subtitle=element_text(size=8))+
  guides(colour=guide_legend(ncol=1,override.aes=list(size=2.2,alpha=1)))
bottom<-p1c|p1e
bottom<-bottom+plot_layout(widths=c(1.05,1))
save_figure((free(p1a)/free(p1b)/bottom+plot_layout(heights=c(.8,1.85,1.55)))&theme(plot.title.position='plot'),'Figure_1_literature',180,235)
write_tsv(data.frame(metric=c('publications','corpus records','reported labels (973 taxonomic; 1 unspecified starter)','shotgun records','amplicon records','isolate records','smoked records','PCo1_pct','PCo2_pct'),value=c(nrow(inventory),ncol(mat),nrow(mat),22,139,121,51,ve[1:2])),'review_counts.tsv')

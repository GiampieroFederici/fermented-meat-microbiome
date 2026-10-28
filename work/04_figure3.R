# Canonical R community matrix; original-paper Conover/BH and CLR screening.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
source('plot_config.R')
read_tsv<-function(p)read.delim(p,check.names=FALSE,stringsAsFactors=FALSE,quote='"')
write_tsv<-function(x,n)write.table(x,file.path('derived',n),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
read_matrix<-function(p){x<-read_tsv(p);m<-as.matrix(x[-1]);storage.mode(m)<-'double';rownames(m)<-x[[1]];m}
kr<-read_matrix('../02_CONTROLLI/Table_S4_Taxonomic_profiles/Kraken_100.tsv')
clr<-read_matrix('data/04/clr_matrix_kraken_species.tsv')
ann<-read_tsv('data/sample_annotation.tsv');ann<-ann[match(rownames(kr),ann$sample_name),]
stopifnot(identical(rownames(kr),rownames(clr)),identical(colnames(kr),colnames(clr)),nrow(kr)==100)
alpha<-read_tsv('derived/samples_kraken.tsv');stats<-read_tsv('derived/community_statistics.tsv');contrasts<-read_tsv('derived/paper_conover_BH.tsv');omni<-read_tsv('derived/paper_alpha_omnibus.tsv')
conover_check<-read_tsv('derived/paper_conover_reproduction.tsv');stopifnot(conover_check$result=='PASS',conover_check$max_absolute_error<1e-10)
stopifnot(all(c('group1','group2','p','p.adj','n1','n2')%in%names(contrasts)))
for(tool in unique(contrasts$tool)) for(g in unique(contrasts$factor)) {
  d<-contrasts[contrasts$tool==tool&contrasts$factor==g,]
  stopifnot(max(abs(d$p.adj-p.adjust(d$p,'BH')))<1e-12)
}
stars<-function(q)ifelse(q<.001,'***',ifelse(q<.01,'**',ifelse(q<.05,'*','ns')))
# Insert-and-absorb display, as in the historical plotting script. It encodes
# the existing Conover/BH matrix, not a new statistical test.
contrast_letters<-function(lev,d){
  groups<-list(lev)
  sig<-d[d$p.adj<.05,]
  if(nrow(sig))for(i in seq_len(nrow(sig))){
    a<-sig$group1[i];b<-sig$group2[i];next_groups<-list()
    for(s in groups){
      if(all(c(a,b)%in%s))next_groups<-c(next_groups,list(setdiff(s,a),setdiff(s,b))) else next_groups<-c(next_groups,list(s))
    }
    groups<-next_groups[!vapply(seq_along(next_groups),function(j)any(vapply(seq_along(next_groups),function(k)j!=k&&all(next_groups[[j]]%in%next_groups[[k]])&&(!setequal(next_groups[[j]],next_groups[[k]])||j>k),logical(1))),logical(1))]
  }
  stopifnot(length(groups)<=26)
  used_letters<-letters[seq_along(groups)]
  lab<-setNames(vapply(lev,function(g)paste0(used_letters[vapply(groups,function(s)g%in%s,logical(1))],collapse=''),character(1)),lev)
  stopifnot(all(nchar(lab)<=length(groups)))
  for(i in seq_len(nrow(d))){shared<-length(intersect(strsplit(lab[[d$group1[i]]],'')[[1]],strsplit(lab[[d$group2[i]]],'')[[1]]))>0;stopifnot(shared==(d$p.adj[i]>=.05))}
  lab
}
add_contrasts<-function(p,g,lev){
  d<-contrasts[contrasts$tool=='kraken'&contrasts$factor==g,];stopifnot(nrow(d)==choose(length(lev),2))
  d$symbol<-stars(d$p.adj);write_tsv(d,paste0('Figure3_contrasts_',g,'.tsv'))
  ymax<-max(alpha$shannon);gap<-diff(range(alpha$shannon))*.16
  if(length(lev)>3){
    labels<-contrast_letters(lev,d)
    a<-data.frame(group=factor(lev,levels=lev),label=unname(labels),y=vapply(lev,function(v)max(alpha$shannon[alpha[[g]]==v]),numeric(1))+gap*.45)
    write_tsv(data.frame(product=lev,letters=unname(labels)),'Figure3_product_contrast_letters.tsv')
    p+geom_text(data=a,aes(group,y,label=label),inherit.aes=FALSE,size=text_mm,family='Arial')
  } else {
    d$x1<-match(d$group1,lev);d$x2<-match(d$group2,lev)
    d<-d[order(abs(d$x2-d$x1),pmin(d$x1,d$x2)),]
    # Place the first bracket above the highest observation.
    step<-.8;d$y<-ymax+.25+step*(seq_len(nrow(d))-1);d$label<-paste0(d$symbol,'  q=',formatC(d$p.adj,digits=2,format='g'))
    p+geom_segment(data=d,aes(x=x1,xend=x2,y=y,yend=y),inherit.aes=FALSE,linewidth=.3)+geom_segment(data=d,aes(x=x1,xend=x1,y=y,yend=y-.12),inherit.aes=FALSE,linewidth=.3)+geom_segment(data=d,aes(x=x2,xend=x2,y=y,yend=y-.12),inherit.aes=FALSE,linewidth=.3)+geom_text(data=d,aes(x=(x1+x2)/2,y=y+.07,label=label),inherit.aes=FALSE,size=text_mm,family='Arial',vjust=0)+scale_y_continuous(limits=c(NA,max(d$y)+.72),expand=expansion(mult=c(.03,0)))
  }
}

# Inherited species membership is preserved; labels avoid unsupported phenotypes.
mapping<-read.csv('data/species_to_group_original.csv',check.names=FALSE,stringsAsFactors=FALSE)
map_keys<-gsub('_',' ',mapping$species,fixed=TRUE);matrix_keys<-gsub('_',' ',colnames(kr),fixed=TRUE)
stopifnot(!anyDuplicated(map_keys),!anyDuplicated(matrix_keys),setequal(matrix_keys,map_keys))
mapping<-mapping[match(matrix_keys,map_keys),];mapping$historical_species_label<-mapping$species;mapping$species<-colnames(kr)
group_labels<-c('LAB'='LAB-associated genera','CNS/Staph'='Staphylococcus and other selected genera','Yeast/Fungi'='Selected fungal genera','Enterobacteria'='Selected enterobacterial genera','Spoilage/other Gram-neg'='Other selected bacterial genera','Pathogen (flagged)'='Flagged assignments (unconfirmed)','Other'='Other taxa')
mapping$group_label<-unname(group_labels[mapping$guild]);stopifnot(!anyNA(mapping$group_label))
write_tsv(mapping,'species_group_mapping.tsv')
group_colours<-setNames(c('#2E7D32','#F5A623','#8E44AD','#C1272D','#5D6D7E','#000000','#D5DBDB'),unname(group_labels))
gm<-sapply(names(group_labels),function(g)rowSums(kr[,mapping$guild==g,drop=FALSE])*100)
stopifnot(max(abs(rowSums(gm)-100))<1e-8)
gp<-aggregate(gm,by=list(product=ann$product),FUN=mean)
write_tsv(gp,'taxonomic_groups_product_means.tsv')
gpl<-do.call(rbind,lapply(names(group_labels),function(g)data.frame(product=gp$product,group=group_labels[[g]],mean_pct=gp[[g]])))
gpl$product<-factor(gpl$product,levels=product_order);gpl$group<-factor(gpl$group,levels=unname(group_labels))
group_n<-table(ann$product);group_axis<-setNames(paste0(product_labels[product_order],' (n=',group_n[product_order],')'),product_order)
# Use the same product order in panels A and B.
p3a<-ggplot(gpl,aes(product,mean_pct,fill=group))+geom_col(width=.8,position=position_stack(reverse=TRUE))+scale_fill_manual(values=group_colours)+scale_x_discrete(limits=product_order,labels=group_axis)+scale_y_continuous(expand=expansion(mult=c(0,0)))+coord_cartesian(ylim=c(0,100))+labs(x=NULL,y='Mean relative\nabundance (%)\nKraken2/Bracken',fill=NULL,title='A: Mean composition by taxonomic group in each product')+theme_paper(9)+theme(legend.position='right',axis.text.x=element_text(angle=30,hjust=1),plot.margin=margin(5,5,0,7,'mm'))+guides(fill=guide_legend(ncol=1))
alpha_box<-function(g,lev,colours,labels=lev){a<-alpha;a[[g]]<-factor(a[[g]],levels=lev);s<-omni[omni$metric=='shannon'&omni$factor==g,];stopifnot(nrow(s)==1,is.finite(s$H),is.finite(s$p),is.finite(s$q_BH));ggplot(a,aes(.data[[g]],shannon,fill=.data[[g]]))+geom_boxplot(outlier.shape=NA,width=.6,alpha=.6)+geom_point(aes(colour=.data[[g]]),position=position_jitter(width=.12,seed=123),size=.8,alpha=.8)+scale_fill_manual(values=colours,guide='none')+scale_colour_manual(values=colours,guide='none')+scale_x_discrete(labels=labels)+labs(x=NULL,y="Shannon H'",title=paste('Alpha diversity by',g),subtitle=sprintf('Kruskal-Wallis H=%.2f; p=%.2g; q=%.2g',s$H,s$p,s$q_BH))+theme_paper(8)}
p3b<-add_contrasts(alpha_box('product',product_order,product_colours,product_labels),'product',product_order)+theme(axis.text.x=element_text(angle=40,hjust=1))
dd<-contrasts[contrasts$tool=='kraken'&contrasts$factor=='product',];i<-match(dd$group1,product_order);j<-match(dd$group2,product_order)
dd$x<-factor(product_order[pmin(i,j)],levels=product_order);dd$y<-factor(product_order[pmax(i,j)],levels=rev(product_order));dd$score<--log10(dd$p.adj);dd$mark<-ifelse(dd$p.adj<.05,stars(dd$p.adj),'')
stopifnot(nrow(dd)==45,all(is.finite(dd$score)))
# Place the colour bar inside the empty upper-right triangle.
p3c<-ggplot(dd,aes(x,y,fill=score))+geom_tile(colour='white',linewidth=.2)+geom_text(aes(label=mark),size=text_mm,family='Arial')+scale_fill_gradient(low='#FFF7BC',high='#D7301F',name='-log10(q)')+scale_x_discrete(labels=product_labels,drop=FALSE)+scale_y_discrete(labels=product_labels,drop=FALSE)+labs(x=NULL,y=NULL,title='Product contrasts',subtitle='Conover-Iman post-hoc (BH-FDR)')+theme_paper(8)+theme(axis.text.x=element_text(angle=40,hjust=1),legend.position='inside',legend.position.inside=c(1,1),legend.justification=c(1,1),legend.direction='horizontal',legend.title.position='top',
  legend.key.width=grid::unit(5,'mm'),legend.key.height=grid::unit(2.5,'mm'),legend.background=element_blank())
p3d<-add_contrasts(alpha_box('country',country_order,country_colours),'country',country_order)
p3e<-add_contrasts(alpha_box('raw_material',raw_order,raw_colours),'raw_material',raw_order)+labs(title='Alpha diversity by raw material')

# Same prevalence and display rules as the historical panel, canonical CLR input.
keep<-colMeans(kr>0)>=.10
screen<-function(m){do.call(rbind,lapply(which(keep),function(j){test<-kruskal.test(m[,j]~ann$country);med<-tapply(m[,j],factor(ann$country,levels=country_order),median);data.frame(species=colnames(m)[j],H=unname(test$statistic),p=test$p.value,highest_median_CLR=names(which.max(med)),mean_abundance_pct=100*mean(kr[,j]),prevalence=mean(kr[,j]>0))}))}
da<-screen(clr);da$q<-p.adjust(da$p,'BH');da<-da[order(da$q,da$species),]
write_tsv(da,'country_taxon_screening.tsv')
# Sensitivity to the old multiplicative zero replacement; never used to select p-values.
delta<-1/ncol(kr)^2;replaced<-kr;zeros<-kr==0
replaced<-replaced*(1-rowSums(zeros)*delta);replaced[zeros]<-delta
oldclr<-log(replaced)-rowMeans(log(replaced));old<-screen(oldclr);old$q<-p.adjust(old$p,'BH')
sens<-merge(da[,c('species','q','highest_median_CLR')],old[,c('species','q','highest_median_CLR')],by='species',suffixes=c('_canonical','_historical_zero_rule'))
write_tsv(sens,'country_screening_zero_sensitivity.tsv')
# A CLR median can favour a country where the species was never detected.
# Display rule: among the 380 BH q<0.05 taxa (test unchanged), assign the country with the highest mean
# relative abundance, keep it only if detected (>0) in >=50% of that country's samples, show the top six per country.
sig<-da$species[da$q<.05]
cm<-sapply(country_order,function(c)100*colMeans(kr[ann$country==c,sig,drop=FALSE]))
cp<-sapply(country_order,function(c)colMeans(kr[ann$country==c,sig,drop=FALSE]>0))
assigned<-country_order[max.col(cm,ties.method='first')]
eligible<-cp[cbind(seq_along(sig),match(assigned,country_order))]>=.5
selected<-unlist(lapply(country_order,function(c){i<-which(eligible&assigned==c);head(sig[i][order(-cm[i,c])],6)}),use.names=FALSE)
stopifnot(length(sig)==380,sum(eligible)==140,length(selected)==18,!anyDuplicated(selected))
shown<-data.frame(species=selected,assigned_country=assigned[match(selected,sig)],mean_pct_Italy=cm[selected,'Italy'],mean_pct_Spain=cm[selected,'Spain'],mean_pct_China=cm[selected,'China'],
  prevalence_Italy=cp[selected,'Italy'],prevalence_Spain=cp[selected,'Spain'],prevalence_China=cp[selected,'China'])
write_tsv(shown,'Figure3_country_taxon_means.tsv')
means<-do.call(rbind,lapply(country_order,function(c)data.frame(country=c,species=selected,mean_pct=cm[selected,c])))
means$country<-factor(means$country,levels=country_order);means$species<-factor(means$species,levels=rev(selected))
# One box per assigned country (Italy | Spain | China), each with its six species and the
# three country bars; species order, colours and values unchanged. Names sit above their bars, coloured by the assigned
# country. Use one x scale (0-65) and a horizontal country legend.
means$box<-factor(shown$assigned_country[match(as.character(means$species),selected)],levels=country_order)
lab<-data.frame(species=factor(selected,levels=rev(selected)),country=shown$assigned_country,box=factor(shown$assigned_country,levels=country_order))
p3f<-ggplot(means,aes(mean_pct,species,fill=country))+geom_col(position=position_dodge(width=.42,reverse=TRUE),width=.4)+
  geom_text(data=lab,aes(0,species,label=species,colour=country),inherit.aes=FALSE,hjust=0,nudge_x=.4,nudge_y=.55,size=text_mm,family='Arial',fontface='italic',show.legend=FALSE)+
  facet_wrap(~box,nrow=1,scales='free_y')+
  scale_fill_manual(values=country_colours)+scale_colour_manual(values=country_colours)+scale_y_discrete(expand=expansion(add=c(.45,.95)))+
  scale_x_continuous(limits=c(0,65),breaks=seq(0,60,20),expand=expansion(mult=c(0,0)))+labs(x='Mean relative abundance (%) - Kraken2/Bracken',y=NULL,fill='Country',
  title='F: Country-associated species (top six per country)')+
  theme_paper(9)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.major.y=element_blank(),strip.background=element_blank(),
  legend.position='bottom',legend.direction='horizontal',legend.margin=margin(0,0,0,0),legend.box.spacing=grid::unit(1,'mm'))
# Arrange panels with A and F full width and the intermediate panels paired.
# Titles keep their texts with the new letters; the Kruskal-Wallis and screening subtitles move to the caption.
pos<-theme(plot.title.position='plot')
p3a<-p3a+pos
p3b<-p3b+labs(title='B: Shannon diversity differs\namong products',subtitle=NULL)+pos+theme(plot.margin=margin(1,5,2,6,'mm'))
p3c<-p3c+labs(title='C: Pairwise product differences\nin Shannon diversity',subtitle=NULL)+pos
p3d<-p3d+labs(title='D: Shannon diversity by country',subtitle=NULL)+pos
p3e<-p3e+labs(title='E: Shannon diversity by raw material',subtitle=NULL)+pos
p3f<-p3f+pos
# Keep the Shannon axis title next to its axis in panel B (free label).
fig3<-wrap_elements(full=p3a+theme_arial())/((free(p3b,type='label',side='l')|p3c)+plot_layout(widths=c(1.25,1)))/(p3d|p3e)/p3f+plot_layout(heights=c(62,41,35,46))&theme_arial()
save_figure(fig3,'Figure_3_alpha_composition',180,240,dpi=600,max_h=240)
write_tsv(data.frame(tested=nrow(da),significant=sum(da$q<.05),zero_rule_significance_agreement=sum((sens$q_canonical<.05)==(sens$q_historical_zero_rule<.05)),displayed_taxa=length(selected)),'Figure3_screening_summary.tsv')
writeLines(capture.output(sessionInfo()),'derived/sessionInfo_figure3.txt')
cat('Figure 4 generated; visual QA required.\n')

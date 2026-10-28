args<-commandArgs(TRUE);ROOT<-normalizePath(args[1],winslash='/');suppressPackageStartupMessages(library(ape))
rd<-function(n)read.delim(file.path(ROOT,'sources',n),check.names=FALSE,stringsAsFactors=FALSE,quote='"',na.strings=c('NA',''))
wt<-function(d,n)write.table(d,file.path(ROOT,'checks',n),sep='\t',row.names=FALSE,quote=TRUE)
tr<-read.tree(file.path(ROOT,'sources/Lsakei_229.treefile'));N<-length(tr$tip.label);stopifnot(N==229)
refs<-rd('../02_CONTROLLI/Table_S10_StrainPhlAn_Lsakei/References_182.tsv');samples<-rd('../02_CONTROLLI/Table_S10_StrainPhlAn_Lsakei/Samples_100.tsv');sp<-samples[samples$Retained_in_tree=='True',]
metagenome_labels<-rd('../02_CONTROLLI/Table_S10_StrainPhlAn_Lsakei/Metagenome_labels.tsv');stopifnot(nrow(metagenome_labels)==47,setequal(sp$Sample,metagenome_labels$original_code));mlab<-setNames(metagenome_labels$leaf_label,metagenome_labels$original_code)
stopifnot(nrow(refs)==182,nrow(sp)==47,setequal(tr$tip.label,c(refs$Accession,sp$Sample)))
key<-function(z){a<-sort(z);b<-sort(setdiff(tr$tip.label,a));if(length(a)>length(b))a<-b else if(length(a)==length(b)&&paste(a,collapse=';')>paste(b,collapse=';'))a<-b;paste(a,collapse=';')}
children<-function(t)split(t$edge[,2],t$edge[,1])
descfun<-function(t){ch<-children(t);desc<-vector('list',N+t$Nnode);rec<-function(v){if(v<=N)z<-v else z<-unlist(lapply(ch[[as.character(v)]],rec),use.names=FALSE);desc[[v]]<<-z;z};rec(setdiff(t$edge[,1],t$edge[,2])[1]);desc}
desc0<-descfun(tr);bs<-rd('../02_CONTROLLI/Table_S10_StrainPhlAn_Lsakei/Branch_support.tsv');bs$key<-vapply(strsplit(bs$Original_Newick_side_leaves,'; ',fixed=TRUE),key,character(1));stopifnot(!anyDuplicated(bs$key))
# Verify that the published support table encodes the original Newick edges.
orig<-do.call(rbind,lapply(which(tr$edge[,2]>N),function(i){v<-tr$edge[i,2];lab<-tr$node.label[v-N];vals<-as.numeric(strsplit(lab,'/',fixed=TRUE)[[1]]);if(length(vals)!=2)vals<-c(NA_real_,NA_real_);data.frame(key=key(tr$tip.label[desc0[[v]]]),SH=vals[1],UF=vals[2])}))
ix<-match(orig$key,bs$key);stopifnot(!anyNA(ix),all(is.na(orig$SH)==is.na(bs$SH_aLRT_percent[ix])),all(is.na(orig$UF)==is.na(bs$UFBoot_percent[ix])),max(abs(orig$UF-bs$UFBoot_percent[ix]),na.rm=TRUE)<1e-10,max(abs(orig$SH-bs$SH_aLRT_percent[ix]),na.rm=TRUE)<1e-10)
# Midpoint root by inserting a display node on the diameter path; branch lengths are preserved.
D<-cophenetic(tr);ab<-which(D==max(D),arr.ind=TRUE)[1,];a<-match(rownames(D)[ab[1]],tr$tip.label);b<-match(colnames(D)[ab[2]],tr$tip.label)
E<-tr$edge;L<-tr$edge.length;adj<-vector('list',N+tr$Nnode+1)
for(i in seq_len(nrow(E))){adj[[E[i,1]]]<-c(adj[[E[i,1]]],E[i,2]);adj[[E[i,2]]]<-c(adj[[E[i,2]]],E[i,1])}
prev<-rep(NA_integer_,length(adj));prev[a]<-0L;q<-a
while(length(q)){v<-q[1];q<-q[-1];for(u in adj[[v]])if(is.na(prev[u])){prev[u]<-v;q<-c(q,u)}}
path<-b;while(tail(path,1)!=a)path<-c(path,prev[tail(path,1)]);path<-rev(path)
geti<-function(x,y)which((E[,1]==x&E[,2]==y)|(E[,1]==y&E[,2]==x))
pidx<-vapply(seq_len(length(path)-1),function(i)geti(path[i],path[i+1]),integer(1));lens<-L[pidx];half<-max(D)/2;k<-which(cumsum(lens)>=half)[1];before<-if(k==1)0 else sum(lens[seq_len(k-1)]);left<-half-before;right<-lens[k]-left;u<-path[k];v<-path[k+1];ri<-N+tr$Nnode+1
E<-rbind(E[-pidx[k],,drop=FALSE],c(ri,u),c(ri,v));L<-c(L[-pidx[k]],left,right)
adj<-vector('list',ri);for(i in seq_len(nrow(E))){adj[[E[i,1]]]<-c(adj[[E[i,1]]],E[i,2]);adj[[E[i,2]]]<-c(adj[[E[i,2]]],E[i,1])}
newE<-matrix(integer(),ncol=2);newL<-numeric();walk<-function(x,par=0){for(y in adj[[x]])if(y!=par){newE<<-rbind(newE,c(x,y));jj<-which((E[,1]==x&E[,2]==y)|(E[,1]==y&E[,2]==x));newL<<-c(newL,L[jj]);walk(y,x)}};walk(ri)
internal<-c(ri,setdiff(sort(unique(newE[newE>N])),ri));mp<-setNames(c(seq_len(N),N+seq_along(internal)),c(seq_len(N),internal));tt<-list(edge=matrix(unname(mp[as.character(newE)]),ncol=2),edge.length=newL,tip.label=tr$tip.label,Nnode=length(internal));class(tt)<-'phylo';tt<-reorder(tt,'cladewise');tt<-ladderize(tt,right=FALSE)
errdist<-max(abs(cophenetic(tt)[rownames(D),colnames(D)]-D));stopifnot(errdist<1e-10)
ch<-children(tt);desc<-descfun(tt);root<-N+1;depth<-rep(0,N+tt$Nnode);walkdepth<-function(x){for(y in ch[[as.character(x)]]){ii<-which(tt$edge[,2]==y);depth[y]<<-depth[x]+tt$edge.length[ii];walkdepth(y)}};walkdepth(root)
# Preserve the eight individually displayed Italian references.
italian<-refs$Accession[which(refs$Source_geography=='Italy')];stopifnot(length(italian)==8)
refids<-match(refs$Accession,tt$tip.label);prot<-match(italian,tt$tip.label)
candidates<-which(vapply(desc,function(z)length(z)>=2&&all(z%in%refids)&&!any(z%in%prot),logical(1)))
collapsed<-integer();visible<-integer();term<-integer()
visit<-function(x){visible<<-c(visible,x);if(x%in%candidates){collapsed<<-c(collapsed,x);term<<-c(term,x)}else if(x<=N){term<<-c(term,x)}else for(y in ch[[as.character(x)]])visit(y)};visit(root)
stopifnot(all(vapply(desc[collapsed],length,integer(1))>=2),all(match(sp$Sample,tt$tip.label)%in%term))
represented<-unlist(lapply(term,function(x)desc[[x]]),use.names=FALSE);stopifnot(length(represented)==229,!anyDuplicated(represented),setequal(represented,1:229))
# Radial phylogram. The hidden descendants contribute only to triangle radial extents.
nt<-length(term);reference_individual<-term[term<=N & !term%in%match(sp$Sample,tt$tip.label)];weights<-ifelse(term%in%reference_individual,2,1);step<-2*pi/sum(weights);theta<-rep(NA_real_,length(depth));theta[term]<-(cumsum(weights)-weights/2)*step+pi/2
atheta<-function(x){if(!is.na(theta[x]))return(theta[x]);ys<-vapply(ch[[as.character(x)]],atheta,numeric(1));theta[x]<<-mean(range(ys));theta[x]};atheta(root)
catcodes<-c(meat='M',vegetal='V',animal_host='A',dairy='D',starter='T',cereal='C',environment='E',seafood='F');print(unique(refs$Source_macro_category));stopifnot(all(refs$Source_macro_category%in%names(catcodes)))
refcode<-setNames(unname(catcodes[refs$Source_macro_category]),refs$Accession)
sourcecols<-c(M='#C62828',V='#2E7D32',A='#8E44AD',D='#E5A600',T='#008C95',C='#A65E2E',E='#5F6B73',F='#1976D2')
source_textcols<-c(M='#B71C1C',V='#1B5E20',A='#70318A',D='#8A6400',T='#006E75',C='#85431C',E='#465158',F='#1558A6')
# One biological colour key applies to both references and metagenomes.
stopifnot(all(sp$Raw_material%in%c('pork','beef')))
spcol<-setNames(rep(unname(sourcecols['M']),nrow(sp)),sp$Sample)
# Apply the source locality labels used in the figure.
meta<-rd('../02_CONTROLLI/Table_S3_Samples_accessions_QC/Samples_100.tsv');sp$PGI_area<-meta$PGI_area[match(sp$Sample,meta$sample_name)]
sap<-c('SAP1','SAP2');spg<-c('SPG1','SPG2');stopifnot(all(sap%in%sp$Sample),all(spg%in%sp$Sample))
stopifnot(all(grepl('^inside PGI area',sp$PGI_area[match(sap,sp$Sample)])),all(grepl('^outside PGI area',sp$PGI_area[match(spg,sp$Sample)])))
sap_label<-'Ciauscolo-type - within PGI area';spg_label<-'Soft salami (Perugia) - outside PGI area'
sp$Display_label<-sp$Product;sp$Display_label[sp$Sample%in%sap]<-sap_label;sp$Display_label[sp$Sample%in%spg]<-spg_label
stopifnot(length(unique(spcol))==1,all(spcol==unname(sourcecols['M'])))
sp$Source_category<-'meat';sp$Label_colour<-unname(source_textcols['M']);sp$Leaf_label<-unname(mlab[sp$Sample]);sp$Display_colour<-unname(spcol[sp$Sample]);wt(sp[,c('Leaf_label','Sample','Product','Country','Raw_material','PGI_area','Display_label','Source_category','Display_colour','Label_colour')],'Figure6_sample_display_coding.tsv')
lab<-setNames(character(nt),term);tipcol<-setNames(character(nt),term);maxr<-setNames(numeric(nt),term)
cltab<-data.frame()
for(x in term){tips<-tt$tip.label[desc[[x]]];maxr[as.character(x)]<-max(depth[desc[[x]]]);if(x%in%collapsed){counts<-table(factor(refcode[tips],levels=names(sourcecols)));nz<-counts[counts>0];lab[as.character(x)]<-paste0(length(tips),' | ',paste0(names(nz),as.integer(nz),collapse=' '));win<-which(counts==max(counts));tipcol[as.character(x)]<-if(length(win)==1)sourcecols[win]else '#B6B6B6';cltab<-rbind(cltab,data.frame(display_node=x,n=length(tips),label=lab[as.character(x)],references=paste(sort(tips),collapse='; ')))}else{lab[as.character(x)]<-if(tips%in%names(mlab))unname(mlab[tips])else tips;tipcol[as.character(x)]<-if(tips%in%names(spcol))spcol[tips]else sourcecols[refcode[tips]]}}
# Label and terminal-edge colours follow the existing metadata categories.
dark_by_fill<-setNames(unname(source_textcols),unname(sourcecols))
leaf_label_colour<-setNames(vapply(term,function(x){id<-as.character(x);if(x<=N&&tt$tip.label[x]%in%names(spcol))return(unname(source_textcols['M']));co<-unname(dark_by_fill[tipcol[id]]);if(is.na(co))'#666666'else co},character(1)),as.character(term))
stopifnot(!anyNA(leaf_label_colour),length(leaf_label_colour)==length(term))
# Reference labels use only fields in the provided Table S10.
labelmeta<-rd('reference_labels_verified.tsv');stopifnot(nrow(labelmeta)==182,setequal(labelmeta$Accession,refs$Accession))
for(x in term[term<=N & term%in%refids]){j<-match(tt$tip.label[x],labelmeta$Accession);lab[as.character(x)]<-paste0(labelmeta$Figure_label_header[j],'\n',labelmeta$Figure_label_name[j])}
nd<-visible[visible>N];nd<-setdiff(nd,root);nodekeys<-vapply(nd,function(x)key(tt$tip.label[desc[[x]]]),character(1));bix<-match(nodekeys,bs$key);stopifnot(!anyNA(bix));uf<-bs$UFBoot_percent[bix];sh<-bs$SH_aLRT_percent[bix]
# A midpoint root splits one original edge in two: label that bipartition only once.
dot<-!is.na(uf)&!duplicated(nodekeys);bin<-ifelse(uf<75,'<75',ifelse(uf<=95,'75-95','>95'));dotcols<-c('<75'='#C8C8C8','75-95'='#686868','>95'='#101010')
wt(data.frame(node=nd,key=nodekeys,SH_aLRT=sh,UFBoot=uf,bin=bin,dot_drawn=dot),'Figure6_display_support.tsv');wt(cltab,'Figure6_collapsed_clades.tsv')
wt(data.frame(display_node=term,label=unname(lab),original_tip=ifelse(term<=N,tt$tip.label[pmin(term,N)],NA_character_),type=ifelse(term%in%collapsed,'collapsed references',ifelse(tt$tip.label[pmin(term,N)]%in%sp$Sample,'metagenome','individual reference')),n_descendants=vapply(desc[term],length,integer(1)),marker_colour=unname(tipcol[as.character(term)]),label_colour=unname(leaf_label_colour[as.character(term)]),terminal_branch_colour=unname(leaf_label_colour[as.character(term)])),'Figure6_display_tips.tsv')
wt(data.frame(tips=229,samples=47,references=182,display_terminals=nt,collapsed_clades=length(collapsed),min_collapse=min(cltab$n),max_patristic_error=errdist,Italian_refs_visible=all(prot%in%term),support_matches=TRUE),'Figure6_checks.tsv')
write.tree(tt,file.path(ROOT,'checks/midpoint_display_only.treefile'))
# Separate scale and metadata legends below the circular tree.
plot_tree<-function(){layout(matrix(c(1,2),2,1),heights=c(1,.28));par(mar=c(.5,.5,.5,.5),family='Arial',xpd=NA)
R<-max(depth);rr<-R*1.10;lim<-R*2.0;plot.new();plot.window(xlim=c(-lim,lim),ylim=c(-lim,lim),asp=1)
xy<-function(r,a)cbind(r*cos(a),r*sin(a))
for(x in visible[visible>N&!visible%in%collapsed]){ys<-ch[[as.character(x)]];aa<-range(theta[ys]);aaseq<-seq(aa[1],aa[2],length.out=max(4,ceiling(diff(aa)/.01)));co<-xy(depth[x],aaseq);lines(co,col='#777777',lwd=.6);for(y in ys){co<-xy(c(depth[x],depth[y]),theta[y]);lines(co,col=if(y%in%term)unname(leaf_label_colour[as.character(y)])else'#777777',lwd=.6)}}
for(x in term){a<-theta[x];kk<-as.character(x);if(x%in%collapsed){aaseq<-seq(a-step*.34,a+step*.34,length.out=15);co<-rbind(xy(depth[x],a),xy(maxr[kk],aaseq),xy(depth[x],a));polygon(co,col=adjustcolor(tipcol[kk],alpha.f=.88),border='#777777',lwd=.6)}
end<-maxr[kk];co<-xy(c(end,rr*.99),a);lines(co,col=adjustcolor(unname(leaf_label_colour[kk]),alpha.f=.22),lwd=.45,lty=3)
if(x<=N){co<-xy(end,a);if(tt$tip.label[x]%in%names(spcol))points(co,pch=15,col=spcol[tt$tip.label[x]],cex=.55)else points(co,pch=21,bg=tipcol[kk],col=if(tt$tip.label[x]%in%italian)'#2565AD'else'#777777',lwd=if(tt$tip.label[x]%in%italian)1.1 else .5,cex=.55)}
deg<-(a*180/pi)%%360;left<-deg>90&&deg<270;co<-xy(rr,a);isciau<-x<=N&&tt$tip.label[x]%in%c('SAP1','SAP2','SPG1','SPG2');text(co[1],co[2],labels=lab[kk],srt=if(left)deg+180 else deg,adj=if(left)1 else 0,cex=if(x%in%collapsed).65 else if(x%in%reference_individual).60 else .60,font=if(isciau)2 else 1,col=unname(leaf_label_colour[kk]))}
for(j in which(dot)){x<-nd[j];co<-xy(depth[x],theta[x]);points(co,pch=16,col=dotcols[bin[j]],cex=.43)}
# Branch-length scale, in the open corner.
scale<-0.1;segments(-lim*.9,-lim*.91,-lim*.9+scale,-lim*.91,lwd=1);text(-lim*.9+scale/2,-lim*.97,'0.1 substitutions/site',cex=.65)
par(mar=c(0,0,0,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
text(.01,.95,'UFBoot (%)',adj=0,cex=.72,font=2);for(j in 1:3){points(.21+(j-1)*.19,.95,pch=16,col=dotcols[j],cex=.8);text(.235+(j-1)*.19,.95,c('<75','75-95','>95')[j],adj=0,cex=.72)}
# Shapes separate metagenomes and references; colour has a single source meaning.
points(.02,.68,pch=15,col=sourcecols['M'],cex=.8);text(.04,.68,'Metagenomic consensus (n = 47)',adj=0,cex=.72)
points(.52,.68,pch=21,bg='white',col='#555555',cex=.8);text(.54,.68,'Reference genome (n = 182)',adj=0,cex=.72)
points(.79,.95,pch=21,bg='white',col='#2565AD',lwd=1.1,cex=.7);text(.81,.95,'Italian reference',adj=0,cex=.63)
text(.01,.43,'Source category (all tips; collapsed references)',adj=0,cex=.66,font=2)
sl<-c('M: meat','V: vegetal','A: animal host','D: dairy','T: starter','C: cereal','E: environment','F: seafood');for(j in 1:8){x<-.02+((j-1)%%4)*.245;y<-.27-((j-1)%/%4)*.14;points(x,y,pch=21,bg=sourcecols[j],col=source_textcols[j],lwd=.55,cex=.75);text(x+.015,y,sl[j],adj=0,cex=.70,col=unname(source_textcols[j]))}
}
png(type='cairo',file.path(ROOT,'outputs/Figure6_circular.png'),width=320,height=320,units='mm',res=600,bg='white');plot_tree();dev.off()
cairo_pdf(file.path(ROOT,'outputs/Figure6_circular.pdf'),width=320/25.4,height=320/25.4,family='Arial');plot_tree();dev.off()
svg(file.path(ROOT,'outputs/Figure6_circular.svg'),width=320/25.4,height=320/25.4,family='Arial',bg='white');plot_tree();dev.off()
png(type='cairo',file.path(ROOT,'checks/Figure6_circular_preview.png'),width=320,height=320,units='mm',res=160,bg='white');plot_tree();dev.off()
cat('CIRCULAR TREE PASS\n');print(rd('../02_CONTROLLI/Table_S10_StrainPhlAn_Lsakei/Group_support.tsv')[,c(1,6,7,8)]);print(cltab[,1:3])


quit(save='no',status=0,runLast=FALSE)

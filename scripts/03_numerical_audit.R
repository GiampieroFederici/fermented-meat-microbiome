args<-commandArgs(TRUE);A<-normalizePath(args[1],winslash='/');C<-file.path(A,'02_CONTROLLI');suppressPackageStartupMessages(library(vegan));options(stringsAsFactors=FALSE)
rd<-function(n,s){p<-list.files(C,pattern=paste0('^Table_S',n,'_'),full.names=TRUE);read.delim(file.path(p,paste0(s,'.tsv')),check.names=FALSE,quote='"',na.strings=c('','NA'))}
res<-data.frame();record<-function(test,error,tol=1e-7,detail=''){res<<-rbind(res,data.frame(test=test,max_error=error,tolerance=tol,pass=is.finite(error)&&error<=tol,detail=detail))}
mat<-function(n,s){d<-rd(n,s);m<-as.matrix(d[-1]);storage.mode(m)<-'double';rownames(m)<-d[[1]];m}
K<-mat(4,'Kraken_100');M<-mat(4,'MetaPhlAn_99');ann<-rd(3,'Samples_100');rownames(ann)<-ann$sample_name;ann<-ann[rownames(K),]
clr<-function(P,eps=NULL){P<-P[,colSums(P)>0,drop=FALSE];P<-P/rowSums(P);if(is.null(eps))eps<-min(P[P>0])/2;L<-log(P+eps);L-rowMeans(L)}
X<-clr(K);Y<-clr(M);HK<-diversity(K);HM<-diversity(M);RK<-rowSums(K>0);RM<-rowSums(M>0)
# Omnibus tests and all102 existing pairwise comparisons.
om<-rd(6,'Omnibus_KW');calc<-numeric(nrow(om));pp<-calc
for(i in seq_len(nrow(om))){v<-if(om$metric[i]=='shannon')HK else RK;r<-kruskal.test(v~as.factor(ann[[om$factor[i]]]));calc[i]<-unname(r$statistic);pp[i]<-r$p.value}
record('S6 six Kruskal-Wallis H',max(abs(calc-om$H)));record('S6 six Kruskal-Wallis p',max(abs(pp-om$p)),1e-12);record('S6 omnibus BH six-test family',max(abs(p.adjust(pp,'BH')-om$q_BH)),1e-12)
for(metric in c('Shannon','richness')){tab<-rd(6,paste0('Conover_',metric));v<-if(metric=='Shannon')HK else RK;pe<-numeric(nrow(tab));te<-pe
for(g in unique(tab$factor)){ii<-which(tab$factor==g);groups<-as.factor(ann[[g]]);N<-length(v);k<-nlevels(groups);ranks<-rank(v);S2<-sum((ranks-mean(ranks))^2)/(N-1);H<-unname(kruskal.test(v~groups)$statistic);means<-tapply(ranks,groups,mean);ns<-table(groups)
for(j in ii){t<-abs(means[tab$group1[j]]-means[tab$group2[j]])/sqrt(S2*(N-1-H)/(N-k)*(1/ns[tab$group1[j]]+1/ns[tab$group2[j]]));te[j]<-t;pe[j]<-2*pt(-t,N-k)}}
record(paste('S6 Conover statistic',metric),max(abs(te-abs(tab$statistic))));record(paste('S6 Conover p',metric),max(abs(pe-tab$p)),1e-10);adj<-unsplit(lapply(split(pe,tab$factor),p.adjust,method='BH'),tab$factor);record(paste('S6 Conover BH',metric),max(abs(adj-tab$p.adj)),1e-10)}
# Depth correlations and richness agreement.
dep<-rd(6,'Depth_samples');dc<-rd(6,'Depth_correlations')
for(i in seq_len(nrow(dc))){v<-if(dc$metric[i]=='shannon')HK[dep$sample_name] else RK[dep$sample_name];cc<-cor.test(v,dep$total_reads,method='spearman',exact=FALSE);record(paste('S6 depth rho',dc$metric[i]),abs(unname(cc$estimate)-dc$rho[i]),1e-10);record(paste('S6 depth p',dc$metric[i]),abs(cc$p.value-dc$p[i]),1e-10)}
common<-intersect(rownames(K),rownames(M));rich<-rd(5,'Richness_tests');cc<-cor.test(RK[common],RM[common],method='spearman',exact=FALSE);record('S5 richness Spearman',abs(unname(cc$estimate)-rich$spearman_rho),1e-10);ww<-wilcox.test(RK[common],RM[common],paired=TRUE,exact=FALSE);record('S5 richness Wilcoxon V',abs(unname(ww$statistic)-rich$wilcoxon_V),1e-10);record('S5 richness Wilcoxon p',abs(ww$p.value-rich$wilcoxon_p),1e-10)
# Global PCA and6-sample PCA, signs may differ.
pca<-prcomp(X);coords<-rd(7,'PCA_Kraken');for(j in 1:2){v<-pca$x[coords$sample_name,j];target<-coords[[paste0('PC',j)]];sgn<-if(cor(v,target)<0)-1 else 1;record(paste('S7 PCA scores PC',j),max(abs(v*sgn-target)))}
variance<-pca$sdev^2/sum(pca$sdev^2)*100;write.table(data.frame(PC=seq_along(variance),variance_percent=variance),file.path(C,'PCA_VARIANCE_AUDIT.tsv'),sep='\t',row.names=FALSE)
cp<-rd(4,'Figure2C_coordinates');cp<-cp[cp$sample%in%rownames(X),];pc6<-prcomp(X[cp$sample,]);vv<-rd(4,'Figure2C_variance');vv<-vv[vv$axis%in%c('PC1','PC2'),];record('S4 six-sample PCA variance',max(abs((100*pc6$sdev^2/sum(pc6$sdev^2))[1:2]-vv$variance_percent)))
for(j in 1:2){v<-pc6$x[,j];target<-cp[[paste0('PC',j)]];sgn<-if(cor(v,target)<0)-1 else 1;record(paste('S4 six-sample PCA scores',j),max(abs(v*sgn-target)))}
# R2/pseudoF are determined exactly by sums of squares, independent of permutations.
ss<-function(Z,g){Z<-scale(Z,center=TRUE,scale=FALSE);N<-nrow(Z);groups<-split(seq_len(N),g);b<-sum(vapply(groups,function(i)length(i)*sum(colMeans(Z[i,,drop=FALSE])^2),numeric(1)));t<-sum(Z^2);k<-length(groups);c(R2=b/t,F=(b/(k-1))/((t-b)/(N-k)))}
st<-rd(7,'PERMANOVA_PERMDISP')
for(i in seq_len(nrow(st))){Z<-if(st$tool[i]=='kraken')X else Y;aa<-ann[rownames(Z),];v<-ss(Z,aa[[st$factor[i]]]);record(paste('S7 R2',st$tool[i],st$factor[i]),abs(v['R2']-st$PERMANOVA_R2[i]),1e-10);record(paste('S7 pseudoF',st$tool[i],st$factor[i]),abs(v['F']-st$PERMANOVA_F[i]),1e-9);h<-kruskal.test(if(st$tool[i]=='kraken')HK else HM,g=aa[[st$factor[i]]]);record(paste('S7 Shannon H',st$tool[i],st$factor[i]),abs(unname(h$statistic)-st$Shannon_H[i]))}
# Country screen from original canonical CLR.
screen<-rd(8,'Country_screening');p<-vapply(screen$species,function(s)kruskal.test(X[,s],g=ann$country)$p.value,numeric(1));record('S8 all463 country screening p',max(abs(p-screen$p)),1e-10);record('S8 BH across463 taxa',max(abs(p.adjust(p,'BH')-screen$q)),1e-10)
# Product-level country exact permutation tests:2520 assignments.
prod<-sort(unique(ann$product));Z1<-t(vapply(prod,function(p)colMeans(X[ann$product==p,,drop=FALSE]),numeric(ncol(X))));Pmean<-t(vapply(prod,function(p)colMeans(K[ann$product==p,,drop=FALSE]),numeric(ncol(K))));Z2<-clr(Pmean);country<-vapply(prod,function(p)unique(ann$country[ann$product==p]),character(1));exact<-rd(9,'Country_product_units')
for(i in 1:2){Z<-list(Z1,Z2)[[i]];v<-ss(Z,country);record(paste('S9 country product-units R2',i),abs(v['R2']-exact$R2[i]),1e-10);record(paste('S9 country product-units F',i),abs(v['F']-exact$pseudoF[i]),1e-9);G<-tcrossprod(scale(Z,center=TRUE,scale=FALSE));t<-sum(diag(G));vals<-numeric();for(italy in combn(1:10,3,simplify=FALSE)){rest<-setdiff(1:10,italy);for(spain in combn(rest,5,simplify=FALSE)){china<-setdiff(rest,spain);b<-sum(G[italy,italy])/3+sum(G[spain,spain])/5+sum(G[china,china])/2;vals<-c(vals,(b/2)/((t-b)/7))}};p<-mean(vals>=v['F']-1e-10);record(paste('S9 all2520 exact assignments p',i),abs(p-exact$p_exact[i]),1e-12)}
# Miao exclusion.
no<-rd(9,'Without_Miao');Z<-clr(K[ann$batch!='miao',]);for(i in seq_len(nrow(no))){v<-ss(Z,ann[rownames(Z),no$factor[i]]);record(paste('S9 without Miao R2',no$factor[i]),abs(v['R2']-no$R2[i]),1e-10);record(paste('S9 without Miao F',no$factor[i]),abs(v['F']-no$pseudoF[i]),1e-9)}
# Zero-replacement and filter scope math.
zr<-rd(9,'Zero_rule_and_filter');P<-K;delta<-1/ncol(P)^2;for(i in 1:nrow(P)){zero<-P[i,]==0;P[i,zero]<-delta;P[i,!zero]<-P[i,!zero]*(1-sum(zero)*delta)};L<-log(P);Z<-L-rowMeans(L);sel<-grepl('Multiplicative',zr$analysis)
for(i in which(sel)){v<-ss(Z,ann[[zr$factor[i]]]);record(paste('S9 multiplicative replacement R2',zr$factor[i]),abs(v['R2']-zr$PERMANOVA_R2[i]),1e-10);record(paste('S9 multiplicative replacement F',zr$factor[i]),abs(v['F']-zr$PERMANOVA_F[i]),1e-9)}
write.table(res,file.path(C,'NUMERICAL_AUDIT_R.tsv'),sep='\t',row.names=FALSE,quote=TRUE);cat('R CHECKS',nrow(res),'PASS',sum(res$pass),'\n');print(res[!res$pass,]);cat('PCA',variance[1:2],'\n');writeLines(capture.output(sessionInfo()),file.path(C,'sessionInfo_AUDIT_R.txt'))

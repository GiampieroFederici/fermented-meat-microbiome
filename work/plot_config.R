# Shared names, relative ordering and colour mappings for every figure.
product_order<-c('Ciauscolo','Felino salami','Piedmont salami','cecina','cured tongue','dry-aged beef','DFS chorizo','DFS salchichon','Chinese sausage (miao)','Mutton sausage')
product_labels<-setNames(c('Ciauscolo-type','Felino salami','Piedmont salami','Cecina','Cured tongue','Dry-aged beef','Chorizo','Salchichon','Miao sour meat','Mutton sausage'),product_order)
product_colours<-setNames(c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499','#661100'),product_order)
# Preserve the original paper order, collapsing the two Ciauscolo labels only
# as required by the approved ten-group R design.
product_order<-c('Piedmont salami','Ciauscolo','Felino salami','Chinese sausage (miao)','DFS salchichon','DFS chorizo','cecina','dry-aged beef','cured tongue','Mutton sausage')
product_colours<-c('Piedmont salami'='#1F9E89',Ciauscolo='#C1272D','Felino salami'='#2E7D32','Chinese sausage (miao)'='#9AA0A6','DFS salchichon'='#E8A33D','DFS chorizo'='#C97B1A',cecina='#8C564B','dry-aged beef'='#5D4037','cured tongue'='#D2A679','Mutton sausage'='#5D6D7E')
country_order<-c('Italy','Spain','China')
country_colours<-c(Italy='#C1272D',Spain='#F5A623',China='#2E6DA4')
country_shapes<-c(Italy=16,Spain=15,China=17)
raw_order<-c('pork','beef','mutton')
raw_colours<-c(pork='#E4A0B7',beef='#8C564B',mutton='#7F7F7F')
technique_colours<-c(Amplicon='#4477AA',Isolate='#228833',Shotgun='#AA3377')
continent_colours<-c(Europe='#0072B2',Asia='#E69F00','South America'='#CC79A7',Africa='#009E73','North America'='#D55E00')
review_product_colours<-c(Sausage='#4E79A7',Salami='#E15759',Ham='#F28E2B','Whole-muscle'='#59A14F',Other='#BAB0AC')
review_raw_colours<-c(Pork='#C7522B',Beef='#8C510A',Poultry='#E9C46A',Mixed='#5B8E7D',Other='#6A4C93')
theme_paper<-function(size=9)ggplot2::theme_bw(base_size=size,base_family='sans')+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major=ggplot2::element_line(colour='grey93',linewidth=.2),axis.text=ggplot2::element_text(colour='black',size=8),legend.text=ggplot2::element_text(size=8),plot.title=ggplot2::element_text(face='bold',size=size+1),plot.tag=ggplot2::element_text(face='bold',size=12),legend.key.size=grid::unit(3,'mm'),legend.title=ggplot2::element_text(face='bold'),plot.margin=ggplot2::margin(5,7,5,5))
save_figure<-function(p,name,w,h,dpi=400,max_h=235) {
  # Reflow at the actual A4 insertion width; do not shrink type after export.
  w<-min(w,180)
  captions<-read.delim('figure_captions.tsv',check.names=FALSE,stringsAsFactors=FALSE)
  row<-captions[captions$figure==name,,drop=FALSE]
  stopifnot(nrow(row)==1)
  # Captions are editable italic Word paragraphs, never rasterised into the figure.
  # Reflow (rather than scale) leaves at least 42 mm on A4 for the external caption.
  # Limit the plot height to fit one A4 page with a separate caption.
  h<-min(h,max_h)
  p$patches$annotation$caption<-NULL
  stopifnot(is.null(p$patches$annotation$caption))
  figure_dir<-Sys.getenv('FIGURE_DIR','../figure_revision_A4')
  if(grepl('^Figure_S[0-9]',name))figure_dir<-Sys.getenv('SUPPLEMENT_FIGURE_DIR','../supplementary_figures')
  dir.create(figure_dir,showWarnings=FALSE,recursive=TRUE)
  writeLines(row$caption,file.path(figure_dir,row$caption_file))
  if(!grepl('_page[12]$',name))ggplot2::ggsave(file.path(figure_dir,paste0(name,'.pdf')),p,width=w,height=h,units='mm',device=grDevices::cairo_pdf,bg='white')
  ggplot2::ggsave(file.path(figure_dir,paste0(name,'.png')),p,width=w,height=h,units='mm',dpi=dpi,bg='white')
  cat('Saved',name,'\n')
}
profiler_colours<-c('Kraken2/Bracken'='#1F78B4','MetaPhlAn 4'='#E6AB02')
# Arial throughout: panel titles 10 pt bold; other text 9 pt;
# geom_text/annotate sizes are in mm, so 9 pt = 9/.pt.
text_mm<-9/ggplot2::.pt
theme_arial<-function()ggplot2::theme(text=ggplot2::element_text(family='Arial',size=9),axis.text=ggplot2::element_text(family='Arial',size=9,colour='black'),
  axis.title=ggplot2::element_text(family='Arial',size=9),legend.text=ggplot2::element_text(family='Arial',size=9),legend.title=ggplot2::element_text(family='Arial',size=9,face='plain'),
  strip.text=ggplot2::element_text(family='Arial',size=9),plot.subtitle=ggplot2::element_text(family='Arial',size=9),plot.title=ggplot2::element_text(family='Arial',size=10,face='bold'))

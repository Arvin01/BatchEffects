# Standard packages
ll <- c('dplyr','magrittr','tidyr','ggplot2','cowplot','ggrepel','GGally','broom',
        'stringr','reshape2','gridExtra','grid','RColorBrewer','BiocInstaller','githubinstall',
        'latex2exp')
sapply(ll,function(l) require(l,character.only = T))
# Bioconductor packages
ip <- rownames(installed.packages())
lb <- c('Biobase','genefilter','annotationTools','hgfocus.db','hgu95av2.db','Homo.sapiens','sva','limma')
for (k in lb) {
  if(k %in% ip) { require(k,character.only=T)} 
  else { biocLite(k); require(k,character.only = T)  }
}
# Note! Data set requires a Github install
library(GSE5859)

options(max.print = 75)

################################################
######## ------ (1) CONFOUNDING ------ #########
################################################

# Load in the GSE5895 data
data(GSE5859)

# We need to remove observations for both the control gene index and duplicate values
control.idx <- grep('AFFX',featureNames(e))
dup.idx <- which(abs(cor(exprs(e)))>0.9999,arr.ind=TRUE) %>% tbl_df %>% 
  mutate(test=(row==col)) %>% filter(!test) %>% use_series(row) %>% min # Use the correlation matrix to quickly identify

# Assign the gene expression data
geneExp <- exprs(e)
# Assign the gene indo
geneInfo <- pData(e) %>% tbl_df %>% mutate(date2=format(date,'%Y')) %>% 
  mutate(ethnicity=recode_factor(ethnicity,ASN='CHB/S',HAN='CHB/S'))
# Drop the duplicated person data
geneExp <- geneExp[,-dup.idx]
geneInfo <- geneInfo[-dup.idx,]
# Drop the control genes
geneExp <- geneExp[-control.idx,]

# Remember the convention for expression data is the each column represents a person! So we actually have 8793 features per person!
gene.names <- rownames(geneExp)
id.names <- colnames(geneExp)
# Assign out ethnicity to a vector
eth <- geneInfo$ethnicity

# Make a plot showing how the sequencing of ethnicities is linked to the date
long.date <- with(geneInfo,table(ethnicity,date2)) %>% data.frame %>% mutate(Freq=ifelse(Freq==0,NA,Freq))
# gg it
gg.date <- 
  ggplot(long.date %>% mutate(Freq=ifelse(is.na(Freq),0,Freq)),aes(x=date2,y=Freq,fill=ethnicity)) + xlab('') + ylab('# of samples') +
  geom_bar(stat='identity',position=position_dodge(),color='black') +
  theme(legend.position=c(0.9,0.85)) + scale_fill_discrete(name='Ethnicity') +
  geom_text(aes(x=factor(2002),y=c(75),label=c('CHB/S=Asian')),color=gg_color_hue(2)[1],nudge_x = 0.4) +
  geom_text(aes(x=factor(2002),y=c(65),label=c('CEU=Caucasian')),color=gg_color_hue(2)[2],nudge_x = 0.4) +
  labs(subtitle='Sequencing date for microarray assays')


# Make sure we can replicate the paper's findings: http://www.nature.com/ng/journal/v39/n2/extref/ng1955-S1.pdf
geneExp[which(gene.names %in% '207245_at'),] %>% as.numeric %>% tbl_df %>% mutate(eth=eth) %>% group_by(eth) %>% summarise(mu=mean(value))

# --- Volcano plot: Asian and European --- #

# Get the t-tests accross the rows
v1.dat <- rowttests(geneExp,eth) %>% data.frame %>% mutate(probe=rownames(.)) %>%
  tbl_df %>% mutate(plog10=-log10(p.value),rid=1:nrow(.))

# # Get the associated 4 most expressed genes
# v1.outliers <- v1.dat %>% mutate(prank=rank(-plog10),dmrank=rank(-abs(dm))) %>% filter(prank<11 | dmrank<11) %>% arrange(prank)
# # Seartch the Affymetrix Human Genome U133 Set annotation data (chip hgu133a) for the extreme outliers
# v1.extreme <- v1.outliers %>% filter(abs(dm)>=1 & round(plog10)>=29)
# v1.extreme <- v1.extreme %>% mutate(gene=gene.names[v1.extreme$rid]) %>% 
#   cbind(.,select(hgu133a.db,keys=.$gene,columns=c('SYMBOL','GENENAME'),keytype = 'PROBEID'))

# Get the Bonferoni cutoff
bonfer <- -log10(0.01/nrow(geneExp))
# Get the Sidak cut off
sidak <- -log10(1-(1-0.01)**(1/nrow(v1.dat)))
# Get number of non-sig
nbonfer <- nrow(v1.dat) - nrow(filter(v1.dat,plog10>bonfer))

# Define dm/pval cutoff
pl=30;dm=c(-1,1)

# Volcano plot with Bonferoni cut-off
v1.count <- v1.dat %>% mutate(bf=ifelse(plog10>bonfer,1,0),pos=ifelse(dm>0.5,'pos','not'),neg=ifelse(dm< -0.5,'neg','not')) %>% 
  group_by(bf,pos,neg) %>% count %>% filter(bf==1 & (pos=='pos' | neg=='neg') ) %>% data.frame %>% mutate(dm=c(-1.5,1.5),plog10=c(40,40))

gg.v1 <-
  ggplot(v1.dat %>% mutate(bonfer=ifelse(plog10>bonfer & abs(dm)>0.5,0,1)),aes(x=dm,y=plog10)) + 
  geom_point(aes(color=factor(bonfer)),shape=1,alpha=0.5) +
  xlab('Mean diff. (CEU vs. CHB/S)') + ylab('-log10(P-value)') +
  geom_hline(yintercept = bonfer,color='orange') + 
  geom_text(aes(x=-1.5,y=8.5,label='Bonferonni cutoff'),color='orange') +
  scale_color_manual(name='',values=c('orange','grey')) + theme(legend.position='none') +
  labs(subtitle='Gene-wise t-tests between ethnicities') +
  geom_vline(xintercept=c(-0.5,0.5),linetype=2) + 
  geom_text(data=v1.count,aes(label=n))

# Calculate the empirical cumulative dist
gg.cdist <- 
ggplot(v1.dat %>% mutate(cum=percent_rank(plog10)) %>% arrange(cum),aes(x=plog10,y=cum*100)) +
  geom_line() + geom_vline(xintercept = bonfer,color='orange',linetype=2) + xlab('-plog10') +
  geom_hline(yintercept = nbonfer/nrow(v1.dat)*100,color='blue',linetype=2) +
  scale_y_continuous('Cumulative share (%)',sec.axis=sec_axis(~.*nrow(v1.dat)/100,name='# of genes')) +
  geom_segment(aes(x=max(v1.dat$plog10),xend=max(v1.dat$plog10),y=nbonfer/nrow(v1.dat)*100,yend=100),
                        arrow=arrow(length=unit(0.5,'cm')),color='blue') +
  geom_text(aes(x=max(v1.dat$plog10)-6,y=nbonfer/nrow(v1.dat)*100+5,
                label=paste(nrow(v1.dat)-nbonfer,'genes',sep=' ')),color='blue') +
  labs(subtitle='Number of t-tests which were significant')


# Find the genes which pass the bonferonni cutoff!
bidx <- v1.dat %>% filter(plog10>bonfer) %>% use_series(rid)

# --------------------------------------------------------- #
# Compare the gene extression for Caucasians 2002 versus 2003

# Get the column ids for the Caucasian people sequenced 2002/2003
cauc.info <- geneInfo %>% mutate(rid=1:nrow(.)) %>% filter(ethnicity=='CEU' & date2<=2003)
widx <- cauc.info$rid
# Run the t-tests on these genes too!
v2.dat <- rowttests(geneExp[,widx],factor(cauc.info$date2)) %>% data.frame %>% mutate(probe=rownames(.)) %>% tbl_df %>%
                    mutate(plog10=-log10(p.value))
# Make a volcanco plot for the CEU v CEU!
v2.count <- v2.dat %>% mutate(bf=ifelse(plog10>bonfer,1,0),pos=ifelse(dm>0.5,'pos','not'),neg=ifelse(dm< -0.5,'neg','not')) %>% 
  group_by(bf,pos,neg) %>% count %>% filter(bf==1 & (pos=='pos' | neg=='neg') ) %>% data.frame %>% 
  mutate(dm=c(-1.0,1.5),plog10=c(12.5,12.5))

gg.v2 <-
  ggplot(v2.dat %>% mutate(bonfer=ifelse(plog10>bonfer & abs(dm)>0.5,0,1)),aes(x=dm,y=plog10)) + 
  geom_point(aes(color=factor(bonfer)),shape=1,alpha=0.5) +
  xlab('Mean diff. (CEU 2002 vs. 2003)') + ylab('-log10(P-value)') +
  geom_hline(yintercept = bonfer,color='darkgreen') + 
  geom_text(aes(x=0.9,y=5,label='Bonferonni cutoff'),color='darkgreen') +
  scale_color_manual(name='',values=c('darkgreen','grey')) + theme(legend.position='none') +
  labs(subtitle='Gene-wise t-tests for Caucasians with year as control') + 
  geom_vline(xintercept = c(-0.5,0.5),linetype=2) +
  geom_text(data=v2.count,aes(label=n))


# Join the tests
v12.dat <- 
  v1.dat %>% transmute(pl1=plog10,dm1=dm,probe=probe) %>%
  left_join(transmute(v2.dat,pl2=plog10,dm2=dm,probe=probe),by='probe')

# Find the genes which show expression differences above the bonferroni...
bf.dat <- 
  v12.dat %>% mutate(pl1=pl1-bonfer,pl2=pl2-bonfer) %>% filter(pl1>0 & pl2>0) %>%
  mutate(dm_abs=(abs(dm1)+abs(dm2))/2)

# Get some genes
xtr.gene <-  bf.dat %>% arrange(desc(dm_abs)) %>% mutate(r=1:nrow(.)) %>% filter(r %in% c(1,14)) %>% use_series(probe)
xtr.idx <- which(gene.names %in% xtr.gene)
# Look up the gene!
xtr.symbol <- select(hgfocus.db,keys=xtr.gene,keytype = 'PROBEID',columns=c('SYMBOL'))

# Get the gene expression data for the CEU/ASN
all.xtr <- geneExp[xtr.idx,] %>% t %>% tbl_df %>% mutate(eth=eth,date=geneInfo$date2) 
race.xtr <- all.xtr %>% gather(var,val,-eth,-date) %>%
  mutate(date2=cut(as.numeric(date),breaks=c(2002,2003,2004,2007),right=F,labels=c('2002','2003','2004+')))
# Plot the differences
gg.xtr <- 
ggplot(race.xtr %>% mutate(var=factor(var,levels=xtr.symbol$PROBEID,labels=xtr.symbol$SYMBOL)),
       aes(x=date2,y=val,color=eth,fill=eth)) + 
  geom_jitter(size=3,alpha=0.75,aes(shape=date2)) + facet_wrap(~var,ncol=2) + 
  scale_shape_manual(values=c(21,22,25),guide=F) +
  xlab('') + ylab('Gene expression') + theme(legend.position = 'bottom', legend.direction = "horizontal") +
  scale_color_discrete(name='Ethnicity') + scale_fill_discrete(name='Ethnicity') + 
  labs(subtitle='Visual demonstation of batch effects')

# Show where the 4 extreme genes stand!
gg.genes <- 
ggplot(bf.dat,aes(x=pl1+bonfer,y=pl2+bonfer)) + 
  geom_point(aes(color=dm_abs)) + scale_color_continuous(low='blue',high='red',name='') +
  xlab('-plog10 (CEU vs. CHB/S)') + ylab('-plog10 (CEU vs. CEU)') +
  geom_text_repel(aes(label=probe),size=4,color='darkgreen',nudge_y =1,
    data=bf.dat %>% filter(probe %in% xtr.gene) %>% mutate(probe=factor(probe,levels=xtr.symbol$PROBEID,labels=xtr.symbol$SYMBOL))) +
  theme(legend.position = 'none') + 
  labs(subtitle='Differentially expressed genes by race/year')

# Add on the two genes to V1 and V2
probe.v1 <- v1.dat %>% filter(probe %in% xtr.gene) %>% mutate(symbol=xtr.symbol$SYMBOL)
probe.v2 <- v2.dat %>% filter(probe %in% xtr.gene) %>% mutate(symbol=xtr.symbol$SYMBOL)
gg.v1 <- gg.v1 + geom_text_repel(data=probe.v1,aes(x=dm,y=plog10,label=symbol),nudge_y = -5)
gg.v2 <- gg.v2 + geom_text_repel(data=probe.v2,aes(x=dm,y=plog10,label=symbol),nudge_y = -2)

# Now check to see whether the 2004/5 white people are differentially expressed to the asian counterparts
widx.0405 <- which(id.names %in% (geneInfo  %>% filter(ethnicity=='CEU' & date2>2003) %>% use_series(filename)))
asn.idx <- which(id.names %in% (geneInfo  %>% filter(ethnicity=='CHB/S') %>% use_series(filename)))
# Run the row-tests
ndiff.0405 <- rowttests(geneExp[,c(widx.0405,asn.idx)],eth[c(widx.0405,asn.idx)]) %>% tbl_df %>% 
  mutate(plog10=-log10(p.value)) %>% filter(plog10>bonfer)
# Now do sampling on the remaining guys
nsim <- 500
nsig <- rep(NA,nsim)
widx.0203 <- which(id.names %in% (geneInfo  %>% filter(ethnicity=='CEU' & date2<=2003) %>% use_series(filename)))
# Loop
for(k in 1:nsim){
  set.seed(k)
  rand.col <- sample(widx.0203,length(widx.0405),replace=F) # Randomly select 16 columns
  nsig[k] <- rowttests(geneExp[,c(rand.col,asn.idx)],eth[c(rand.col,asn.idx)]) %>% tbl_df %>% 
    mutate(plog10=-log10(p.value)) %>% filter(plog10>bonfer) %>% nrow
}
# Now make a histogram chart and compare to the ndiff.0405
gg.mc <-
ggplot(data.frame(nsig),aes(x=nsig,y=..count..)) + 
  geom_histogram(color='blue',fill='grey',bins=35) + 
  geom_vline(xintercept=nrow(ndiff.0405),color='orange',linetype=2) +
  geom_text(aes(x=nrow(ndiff.0405)+400,y=20,label=paste('CEU 04/05',nrow(ndiff.0405),sep='=')),color='orange') + 
  xlab('# of genes > Bonferroni') + ylab('Frequency') + 
  scale_x_continuous(breaks=seq(0,2500,500),labels=seq(0,2500,500),limits=c(0,2500)) + 
  labs(subtitle='MC simulations of 16 CEUs from 2002/03')

####################################################
######## ------ (2) CORRELATION EXP ------ #########
####################################################

# Get the gene correlation between individuals
gene.cor <- cor(geneExp)
# Normalize the rows to have zero means (i.e. normalized gene expression)
y <- geneExp - rowMeans(geneExp)
# Get the chromosome location for our genes
map2gene <- mapIds(hgfocus.db, keys=gene.names,column="ENTREZID", keytype="PROBEID",multiVals="first")
if(length(map2gene==nrow(geneExp))) {
  map2chr <- mapIds(Homo.sapiens, keys=map2gene,column="TXCHROM", keytype="ENTREZID",multiVals="first")
} else {print('HOLD THE GENETIC PHONE!')}
# Which which variables are in chromosome y
chry.idx <- str_match(map2chr,'chrY')[,1] %>% is.na %>% not %>% which
# Now that we know the gene rows associated with Chromosome Y, we can see how their expression varies between person
chry.Exp <- colMeans(y[chry.idx,])
# hist(chry.Exp)
# Clearlay anyone above 0 is 'male' and below is female
gender.idx <- ifelse(as.numeric(chry.Exp)>0,'M','F') %>% as.factor
geneInfo$gender <- gender.idx
# Now check for differential expression in white people 2002/2003 not on Y-chromosome -chry.idx
v3.dat <- rowttests(geneExp[,widx.0203],gender.idx[widx.0203]) %>% data.frame %>% 
  mutate(probe=rownames(.)) %>% tbl_df %>% mutate(plog10=-log10(p.value))
v3.labels <- v3.dat %>% mutate(rid=1:nrow(.),chr=as.character(map2chr)) %>% filter(plog10>bonfer & abs(dm)>0.5)
gg.v3 <-
ggplot(v3.dat %>% mutate(bonfer=ifelse(plog10<=bonfer,1,0)),aes(x=dm,y=plog10)) + 
  geom_point(aes(color=factor(bonfer)),shape=1,alpha=0.5) +
  xlab('Mean diff. (CEU 2002-3 Men vs. Women)') + ylab('-log10(P-value)') +
  geom_hline(yintercept = bonfer,color='red') +
  geom_text(aes(x=-3.5,y=9.5,label='Bonferonni cutoff'),color='red') +
  scale_color_manual(name='',values=c('red','grey')) + theme(legend.position='none') +
  geom_text_repel(aes(x=dm,y=plog10,label=chr),data=v3.labels,color='red') + 
  labs(subtitle='Gene-wise t-tests between genders')

# ------------ CORRELATION OF GENE EXPRESSION ---------------- #

library(pheatmap)

# Create test matrix
test = cor(y)
colnames(test) = paste('id',1:nrow(test),sep='')
rownames(test) = paste('id',1:nrow(test),sep='')

# Generate annotations for rows and columns
annotation_col = data.frame(
  EthnicType = geneInfo$ethnicity %>% factor(labels=c('CHB','CEU'))
)
rownames(annotation_col) = colnames(test)

annotation_row = data.frame(
  YearSeq = geneInfo$date2 %>% factor
)
rownames(annotation_row) = rownames(test)

library(RColorBrewer)

gr=colorRampPalette(brewer.pal(3,"Greens"))(5)
rr=colorRampPalette(brewer.pal(3,"Reds"))(6)

# Specify colors
YearSeq <- gr
names(YearSeq) <- 2002:2006
ann_colors = list(
  EthnicType = c(CEU = rr[2], CHB = rr[4]),
  YearSeq = YearSeq
)

# Create chart
gg.heat <- pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row,
         annotation_colors = ann_colors,cluster_rows = F,cluster_cols = F,
         show_rownames = F,show_colnames = F,annotation_names_col	= F,annotation_names_row = F,
         color=brewer.pal(4,"Blues"))

lblue <-  brewer.pal(4,"Blues")[1]
hblue <-  brewer.pal(4,"Blues")[4]
test.dat <- data.frame(x=seq(-0.5,1,0.01),y=seq(-0.5,1,0.01)) %>% tbl_df %>% mutate(cc=ntile(x,4),cc2=-7/20+3/20*cc)
test.gg <- 
  ggplot(data.frame(cc=test.dat$cc2 %>% unique,x=1:4)) + 
  geom_tile(aes(x=x,y=cc,fill=factor(cc))) +
  scale_fill_manual(name='',values=brewer.pal(4,'Blues'),
                    labels=c(-0.5,0,0.5,1)) +
  theme(legend.direction = 'horizontal',legend.position = 'bottom')

# And grid it
gg.heat.grid <- 
  plot_grid(gg.heat$gtable[3:5,6] %>% gtable_squash_cols(cols=c(1:5)),
          plot_grid(gg.heat$gtable[3:4,1:4] %>% gtable_squash_cols(cols=c(1,2,4:6)),
                    get_legend(test.gg),nrow=2,rel_heights = c(7,1)),
          ncol=2,rel_widths = c(1,15)) +
  labs(subtitle='            Correlation of genetic expression by person')
# dev.off()
gg.heat.grid
###################################################
######## ------ (3) EXPLORATORY EDA------ #########
###################################################

# Get the SVD for the gene-normalized data
y.svd  <- svd(y)

# Get the share of variance explained by the factors
y.share <- with(y.svd,data.frame(nshare=d^2/sum(d^2)*100,npc=1:length(d))) %>% tbl_df %>% mutate(cumshare=cumsum(nshare))

gg.svd <- 
ggplot(y.share,aes(x=npc,y=nshare)) + geom_point(shape=21) +
  geom_line(color='blue',linetype=2) + xlab('# of components') + ylab('Explained (%)') +
  geom_point(aes(x=npc,y=cumshare/(100/max(nshare))),color='red',shape=21) +
  geom_line(aes(x=npc,y=cumshare/(100/max(nshare))),color='red',linetype=2) +
  scale_y_continuous(breaks=seq(0,20,5),limits=c(0,20),
                     sec.axis = sec_axis(~.*(100/max(.)),name='Cumulative (%)')) + 
  labs(subtitle='% of variation explained by ranked PCs')

# ------ MDS PLOT ------- #

# Use first two PCs for multidemnsional scalings
y.12 <- y.svd$v[,1:4] %>% set_colnames(paste('PC',1:ncol(.),sep=''))
y.12.tbl <- data.frame(y.12,eth=eth,date=geneInfo$date2,datemon=format(geneInfo$date,'%b-%y')) %>% tbl_df

# ------ SVM ------- #
svm.fun <- function(my.data) {
  # Fit the svm model
  svm.model <- svm(type ~ ., data=my.data, type='C-classification', kernel='linear',scale=FALSE)
  # get parameters of hiperplane
  w <- t(svm.model$coefs) %*% svm.model$SV
  b <- -svm.model$rho
  # in this 2D case the hyperplane is the line w[1,1]*x1 + w[1,2]*x2 + b = 0
  int <- -b/w[1,2]
  slope <- -w[1,1]/w[1,2]
  return(c(int,slope))
  # abline(a=-b/w[1,2], b=-w[1,1]/w[1,2], col="blue", lty=3)
}

library(e1071)
library(kernlab)

# Fit a linear SVM for ethnicity
eth.data <- with(y.12.tbl,data.frame(x1=PC1,x2=PC2,type=eth))
eth.line <- svm.fun(eth.data)

# Fit a linear SVM for date
date.data <- with(y.12.tbl,data.frame(x1=PC1,x2=PC2,type=ifelse(date %in% c(2002:2004),'early','late')))
date.line <- svm.fun(date.data)



# By ethnicity/date
gg.svd.eth.date <-
  ggplot(y.12.tbl,aes(x=PC1,y=PC2,color=date,shape=eth)) + geom_point(size=3) + 
  scale_color_discrete(name='Date: ') + scale_shape_manual(name='Ethnicity',values=c(1,3)) + 
  labs(subtitle='SVM with linear separator (ethnicity/date<2005) \nMDS clustering by ranked PC1/PC2') +
  geom_abline(intercept=eth.line[1],slope=eth.line[2],linetype=2) +
  geom_abline(intercept=date.line[1],slope=date.line[2],linetype=4)



# ----- Intra-PCA distribution ------ #

dm.levels <-  unique(y.12.tbl$datemon)[y.12.tbl$datemon %>% unique %>% paste('01',sep='-') %>% as.Date('%b-%y-%d') %>% order]
# Get long format
long.pc <- y.12.tbl %>% dplyr::select(-eth,-date) %>% gather(pc,val,-datemon) %>%
                mutate(datemon=factor(datemon,levels=rev(dm.levels)))
gg.datemon.pc <-  
  ggplot(long.pc %>% filter(pc %in% c('PC1','PC2')),aes(y=val,x=datemon,color=datemon)) + 
  geom_boxplot() + facet_wrap(~pc,ncol=2) + coord_flip() + 
  theme(legend.position='none',axis.text.y = element_text(size=8)) + xlab('') + ylab('') + 
  labs(subtitle='Box plot for first 2 PCs by year-month of sequence')
 
# ----- SHARE OF EACH PC EXPLAINED BY THE MONTHS ------ #

all.pc <- y.svd$v %>% set_colnames(paste('PC',1:ncol(.),sep='')) %>% tbl_df %>% 
                mutate(datemon=format(geneInfo$date,'%b-%y')) %>% gather(pc,val,-datemon)
ols.pc <- all.pc %>% group_by(pc) %>% do(ols=lm(val~factor(datemon),data=.))
# Get the stats
reg.stats <- glance(ols.pc %>% rowwise,ols) %>%
      transmute(ar2=adj.r.squared,fs=p.value,sig=ifelse(fs<0.05,1,0),pcnum=gsub('PC','',pc) %>% as.numeric) %>% arrange(pcnum)
# plot!
gg.pc.ols <- 
ggplot(reg.stats,aes(x=pcnum,y=ar2*100)) + 
  geom_point(aes(color=factor(sig))) + 
  xlab('PC number') + ylab('Adj. R-squared (%)') + 
  scale_color_discrete(name='F-stat',labels=c('<5%','>5%')) +
  theme(legend.position = c(0.7,0.8)) + 
  labs(subtitle='Regression of PCs on year-month dummies')


############################################################
######## ------ (4) Adjusting Linear Models ------ #########
############################################################

library(xtable)

# Get the controls
yy <- geneInfo$date %>% format('%Y')
ym <- geneInfo$date %>% format('%b%y')

# Prepare an example
x <- geneExp[1,] # Some gene
batch <- as.factor(yy)
contrasts(batch) <- contr.sum(levels(batch))
batch.mat <- model.matrix(~batch)[, -1, drop = FALSE]
batch.dat <- data.frame(x=x,int=1,batch.mat,row.names = NULL)
beta <- as.matrix(coefficients(lm(x~.-1,data=batch.dat))[-1],ncol=1)
clean <- as.matrix(x) - batch.mat %*% beta

# Check
q=removeBatchEffect(geneExp[1:2,] %>% set_colnames(NULL),batch=yy)
lapply(list(q[1,],clean,x), function(ll) ll %>% as.numeric %>% head)

# Toy example
y <- matrix(rnorm(10*8),10,8)
y[,3:4] <- y[,3:4] + 5
batch <- rep(LETTERS[c(2,1,4,3)],each=2)
y.ctr <- rbe.contrs(y,batch)
y.ctr$beta %>% apply(1,scale) %>% t %>% apply(2,mean)

# Print contrast to Latex
# attributes(batch)[[3]] %>% set_colnames(paste('Contrast',1:ncol(.),sep=' ')) %>% xtable %>% print()

# Now run the remove Batch effects on them all
rbe.geneExp <- removeBatchEffect(geneExp,batch=ym)

# Get the t-tests accross the rows
v4.dat <- rowttests(rbe.geneExp,eth) %>% data.frame %>% mutate(probe=rownames(.)) %>%
  tbl_df %>% mutate(plog10=-log10(p.value),rid=1:nrow(.))

# Volcano plot with Bonferoni cut-off
gg.v4 <-
  ggplot(v4.dat %>% mutate(bonfer=ifelse(plog10<=bonfer,1,0)),aes(x=dm,y=plog10)) + 
  geom_point(aes(color=factor(bonfer)),shape=1,alpha=0.5) +
  geom_hline(yintercept = bonfer,color='purple') + 
  xlab('Mean diff. (CEU vs. CHB/S)') + ylab('-log10(P-value)') +
  geom_text(aes(x=0,y=5.0,label='Bonferonni cutoff'),color='purple') +
  scale_color_manual(name='',values=c('grey','purple')) + theme(legend.position='none') + 
  labs(subtitle='Gene t-tests after linear contrast adjust')

# Get the contrasts
rbe.ctrs <- rbe.contrs(geneExp,batch=ym)
# average batch effect
rbe.avg <- rbe.ctrs$beta %>% apply(2,mean) %>% as.numeric
rbe.df <- data.frame(ctr=rbe.avg,batch=rbe.ctrs$batch %>% colnames,date=unique(ym)[-length(unique(ym))])
# ggplot(rbe.df,aes(x=factor(date),y=ctr,color=date)) + 
#   geom_point(show.legend = F) + coord_flip()
  
rel <- rbe.ctrs$beta %>% apply(1,scale) %>% t %>% tbl_df
ymf <- levels(factor(ym))[-length(unique(ym))]
yml <- as.Date(paste(ymf,'01',sep=''),format='%b%y%d') %>% sort %>% format('%b%y')

rbe.df <- data.frame(ctr=rel %>% summarise_each(funs(mean(.))) %>% as.numeric,
                     batch=factor(ymf,levels=yml))
rbe.df <- geneInfo %>% mutate(batch=format(date,'%b%y'),is.white=ifelse(ethnicity=='CEU',1,0)) %>% 
  group_by(batch) %>% summarise(white.share=mean(is.white)) %>% 
  filter(batch %in% ymf) %>% mutate(batch=factor(batch,levels=yml)) %>% 
  left_join(rbe.df,by='batch')

gg.rbe <- ggplot(rbe.df %>% mutate(white.share=paste(round(white.share,0)*100,'%',sep='')),
       aes(y=batch,x=ctr,color=factor(white.share))) + geom_vline(xintercept = 0,linetype=2) +
  geom_point(size=3) + xlab('Contrast coefficients') +
  ylab('') + labs(subtitle='Standardized coefficients averaged over all genes') +
  scale_color_discrete(name='Share of CEU') + theme(legend.position = c(0.75,0.2))

###########################################
######## ------ (5) ComBat ------ #########
###########################################

# Get the batch vector and intercept model matrix
batch <- ym
modcombat = model.matrix(~1, data=geneInfo)
# Run ComBat
combat_edata = ComBat(dat=geneExp, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# Get the t-tests accross the rows
v5.dat <- rowttests(combat_edata,eth) %>% data.frame %>% mutate(probe=rownames(.)) %>%
  tbl_df %>% mutate(plog10=-log10(p.value),rid=1:nrow(.))
gg.v5 <- 
  ggplot(v5.dat %>% mutate(bonfer=ifelse(plog10<=bonfer,1,0)),aes(x=dm,y=plog10)) + 
  geom_point(aes(color=factor(bonfer)),shape=1,alpha=0.5) +
  geom_hline(yintercept = bonfer,color='cyan3') + 
  xlab('Mean diff. (CEU vs. CHB/S)') + ylab('-log10(P-value)') +
  geom_text(aes(x=0,y=5.0,label='Bonferonni cutoff'),color='cyan3') +
  scale_color_manual(name='',values=c('grey','cyan3')) + theme(legend.position='none') + 
  labs(subtitle='By ethnicity: gene t-tests after ComBat') +
  scale_x_continuous(limits=c(-0.5,0.5)) + geom_vline(xintercept = c(-0.5,0.5),linetype=2)
# Robustness check by gender
v6.dat <- rowttests(combat_edata,gender.idx) %>% data.frame %>% mutate(probe=rownames(.)) %>%
  tbl_df %>% mutate(plog10=-log10(p.value),rid=1:nrow(.))
probe.check <- v6.dat %>% filter(plog10>bonfer & abs(dm)>0.5) %>% use_series(probe)
m2g <- mapIds(hgfocus.db, keys=probe.check,column="ENTREZID", keytype="PROBEID",multiVals="first")
m2c <- mapIds(Homo.sapiens, keys=m2g,column="TXCHROM", keytype="ENTREZID",multiVals="first")
gp <- data.frame(probe=names(m2g),chr=m2c) %>% left_join(v6.dat %>% filter(plog10>bonfer))
# See what we've lost!
lost.sex <- v3.labels$probe[(v3.labels$probe %in% probe.check) %>% not %>% which]
lostg <- mapIds(hgfocus.db, keys=lost.sex,column="ENTREZID", keytype="PROBEID",multiVals="first")
lostc <- mapIds(Homo.sapiens, keys=lostg,column="TXCHROM", keytype="ENTREZID",multiVals="first")
lp <- data.frame(probe=lost.sex,chr=lostc) %>% left_join(v6.dat %>% filter(plog10>bonfer))

gg.v6 <- 
  ggplot(v6.dat %>% mutate(bonfer=ifelse(plog10>bonfer & abs(dm)>0.5,1,0)),aes(x=dm,y=plog10)) + 
  geom_point(aes(color=factor(bonfer)),shape=1,alpha=0.5) +
  geom_hline(yintercept = bonfer,color='indianred') + 
  xlab('Mean diff. (Men vs. Women)') + ylab('-log10(P-value)') +
  geom_text(aes(x=-3,y=12.0,label='Bonferonni cutoff'),color='indianred') +
  scale_color_manual(name='',values=c('grey','indianred')) + theme(legend.position='none') + 
  labs(subtitle='By gender: gene t-tests after ComBat \nBlue genes no longer significant after adjustment') + 
  geom_text_repel(aes(x=dm,y=plog10,label=chr),data=gp,color='indianred') + 
  geom_text_repel(aes(x=dm,y=plog10,label=chr),data=lp,color='blue',nudge_x=1) +
  geom_vline(xintercept = c(-0.5,0.5),linetype=2)

# mapIds(Homo.sapiens, keys=m2g[which(m2c=='chr5')],column=c('GENENAME'), keytype="ENTREZID",multiVals="first")

#################################################
######## ------ (X) GG PLOTTING ------ #########
################################################

dev.off()
# Aggregate plot genetic analysis
agg1 <- plot_grid(gg.v1,gg.v2,gg.date,gg.xtr + theme(legend.position = 'none'),labels=LETTERS[1:4])
# Aggregate plots for confounding
agg2 <- plot_grid(gg.mc,gg.svd.eth.date,gg.datemon.pc,gg.heat.grid,
                  ncol=2,nrow=2,labels=LETTERS[1:4])

# # Aggregate plot for Explaratory PCA
# agg3 <- plot_grid(gg.v3,gg.svd.eth.date,gg.svd,gg.pc.ols,labels=LETTERS[1:4],ncol=2,nrow=2)
# # Aggregate plot for linear constrats
# agg4 <- plot_grid(gg.v4,gg.rbe,labels=LETTERS[1:2],ncol=2)

# Aggregate plots for ComBat
agg5 <- plot_grid(gg.v5,gg.v6,labels=LETTERS[1:2],ncol=2)

# Save the list for markdown
rmd.list <- list(agg1=agg1,
                 agg2=agg2,
                 agg5=agg5)
save(rmd.list,file='rmd_data.RData')
load('rmd_data.RData')


## R code for article "Is There a Generational Shift in Preferences for Forest Carbon Sequestration 
#vs Preservation of Agricultural Landscapes?"

# load packages 
library(labelled)	 # labeling data
library(MASS)	     # ordinal  

# MCA 
library(FactoMineR)  # MCA analysis
library(factoextra)  # MCA visulization
library(RColorBrewer)
library(poLCA)   	 # latent class analysis for categorical variables
library(cowplot)  	 # plot_grid
library(stringr) 	 # str_to_title, convert first to capital
library(prettyGraphs)# add.alpha
library(ggrepel)     # 
library(plotrix) 	 # image and contour,sizeplot
library(ggtern)	 	 # kde2d.weighted 
library(ggplot2)
library(ggpubr)		 # combine figures,map_dbl


# define the root directories
root<- "C:/Users/a34543/Documents/Norce"
res.root<-paste0(root,"/wave16_manipulated photo/analysis_022021/")

# LOAD data
load(file=paste0(root,"/wave16_manipulated photo/processed data/wave16_clean0.Rdata"),ver=TRUE) # dat

## extract motivation variables
mt.start<-which(names(dat)=="aesthetic")
mt.end<-which(names(dat)=="carbon")

motiv<-dat[,c(mt.start:mt.end)]
motiv.id<-cbind(newid=dat$newid,motiv)

## extract supplementary variables
# upbringing variables
g.start<-which(names(dat)=="growup.spruce")
g.end<-which(names(dat)=="growup.none")

gro<-dat[,c(g.start:g.end)] 
# rename
names(gro)<-as.factor(paste0("g.",substr(colnames(gro),8,12)))
head(gro)

# variables about profession 
id<-which(names(dat)=="newid")

wk.start<-which(names(dat)=="work.agriculture")
wk.end<-which(names(dat)=="work.nonrelated")
work<-dat[,c(wk.start:wk.end)]
colnames(work)<-as.factor(paste0("w.",substr(colnames(work),6,14)))
work<-work[,!names(work)%in%c("vetikke")]
head(work)
# combine data 
mot.gw<-cbind(motiv.id,gro,work)

# supplementary variables
var<-c("newid","gender","region","birth.k6","norind","concernClimate",
		"image.a","image.b","image.c",
		"residenceType2","edu.k3","r14.income")
supp<-dat[,var1]

mgw.data<-merge(mot.gw,supp,by="newid",all.x=TRUE)

summary(mgw.data) 


#######################################
######### 1. RUN MCA Analysis #########
#######################################

n.q<-which(names(mgw.data[,-1])=="g.spruc") # position of the last supplementary variables used for MCA analysis
res.mca<-MCA(mgw.data[,-1],quali.sup=n.q:ncol(mgw.data[,-1]),graph=FALSE) 
summary(res.mca)

## different visulizations of MCA results 
# confidence ellipses around suppl variables
p1<-plotellipses(res.mca,keepvar=c("birth.k6"),invisible=c("ind"),type="p",ylim=c(-0.25,0.5),xlim=c(-0.25,0.5))+labs(title="Year of birth")
p1
p2<-plotellipses(res.mca,keepvar=c("concernClimate"),invisible=c("ind"),ylim=c(-0.25,0.5),xlim=c(-0.25,0.5))+labs(title="Climate change concern")
p2
p3<-plotellipses(res.mca,keepvar=c("image.a"),invisible=c("ind"),ylim=c(-0.25,0.5),xlim=c(-0.25,1))+labs(title="Image chosen")#,ylim=c(-0.2,0.3)) # this is the candidate plot for paper
p3
p4<-plotellipses(res.mca,keepvar=c("image.c"),invisible=c("ind"),ylim=c(-0.25,0.5),xlim=c(-0.25,1))+labs(title="Image chosen")#,ylim=c(-0.2,0.3)) # this is the candidate plot for paper
p4
ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)

plotellipses(res.mca,keepvar=c("residenceType2"),invisible=c("ind"),graph.type = c("classic"))# ggplot
plotellipses(res.mca,keepvar=c("edu.k3"),invisible=c("ind"),graph.type = c("classic"))# ggplot

### Table A2 #########
# correlation with Principal dimensions
quali<-function(n,dat){
	Dim<-dimdesc(dat,proba=0.05)[[n]]$quali # proba=.05, significance threshold considered to characterize the dimension
	Dim<-data.frame(Dim)
	Dim$Variables<-rownames(Dim)
	Dim<-Dim[,c("R2","Variables")]
	row.names(Dim) <- NULL
	Dim$R2<-round(Dim$R2,2)
	Dim$Variables<-gsub("norind","party",Dim$Variables)
	Dim
}
dim1<-quali(1,res.mca)
dim2<-quali(2,res.mca)
dim3<-quali(3,res.mca)

# link var.star with var1
var.star<-unique(c(rownames(dim1),rownames(dim2),rownames(dim3)))

mot<-rownames(res.mca$var$coord)
n1<-unlist(lapply(gregexpr(pattern="_",mot),min))
mot<-unique(substr(mot,1,n-1))

# Table A2
myfun<-function(Dim,n,dat){ 
	d1<-round(dimdesc(dat,proba=0.05)[[n]]$category,3)
	nx<-unlist(lapply(gregexpr(pattern="=",rownames(d1)),min))
	name<-substr(rownames(d1),nx+1,nchar(rownames(d1)))
	name.var<-substr(rownames(d1),1,nx-1)
	d1<-data.frame(dimension=n,Categories=name, Variables=name.var,d1)
	d1$Variables<-gsub("norind","party",d1$Variables)
	#d1<-d1[abs(d1$Estimate)>=cutpoint,]
	row.names(d1) <- NULL
	d1<-d1[,1:4]
	
	d1<-merge(d1,Dim,by="Variables",all.x=TRUE)
	d1<-d1[,c("Categories","Estimate","Variables","R2")]
	d1<-d1[order(d1$R2,d1$Variables,decreasing=TRUE),]
	d1$new.var<-duplicated(d1$Variables)
	d1<-d1[complete.cases(d1),]
	#d1[c(grep("TRUE",d1$new.var)),c("Variables","R2")]<-""
	d1$vartype<-with(d1,ifelse(Variables %in% mot, 1,0))
	d1$Variables<-paste(toupper(substr(d1$Variables, 1, 1)), substr(d1$Variables, 2, nchar(d1$Variables)), sep="") # first letter upper
	#d1$Estimate<- round(d1$Estimate,3)
	
	d1<-d1[order(d1$vartype,d1$R2,d1$Variables,decreasing=TRUE),]
	d1<-d1[,!names(d1) %in% c("new.var")]
}

x1<-myfun(dim1,1,res.mca) 
x2<-myfun(dim2,2,res.mca)

dim(x1);dim(x2)


##############################################
## 2. K-means clustering based on MCA  #######
##############################################

# extract mca results
d.mca<- data.frame(get_mca_ind(res.mca)$coord)

# select optimal cluster
d.mca<-d.mca[,1:5] # first 5 dimensions, explained 85% of variance

k.values<- 1:10
wss_values<-sapply(k.values,FUN=function(x)wss(x,d.mca))

# determine number of clusters
plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")

# based on the plot to inspect optimal k
# k.optimal=3 is more consistent with HCPC result
k.optimal<-4
k.optimal<-3

seed<-runif(1,max=10000)
set.seed(seed)
k_best<-kmeans(d.mca,centers=k.optimal, nstart=25) # can change number of dimensions used as inputs

# ordered centers by running two times
centers_ny<-k_best$centers[order(k_best$centers[,1],decreasing=FALSE),]
k_best<-kmeans(d.mca,centers=centers_ny, nstart=25) # can change number of dimensions used as inputs
round(proportions(table(k_best$cluster)),2)

d.mca$clust.kmeans<-as.factor(k_best$cluster)

# Inspect k-means cluster's separability 
k1<-fviz_cluster(k_best,axes=c(1,2),data=d.mca[,colnames(d.mca)!=c("clust.kmeans")],geom="point",stand=FALSE) #
k1 # inspect the dimensions
k12<-ggplot(d.mca, aes(x=Dim.1,y=Dim.2,colour=as.factor(clust.kmeans)))+geom_point()
k12
ggarrange(k1,k12,nrow=2)


###############################################
## 3. Hierachical clustering from MCA results ##
###############################################
# data: res.mca, mgw.data, motiv.id

set.seed(12345)
hcpc.mca<-HCPC(res.mca, nb.clust=-1,consol=1,min=3,max=3,iter.max=10,order=TRUE) # nb.clust=-1, automatically cut 

# make sure it is the same order (res.mca also use mgw.data)
d.hcpc<-cbind(newid=mgw.data$newid,hcpc.mca$data.clust)   
round(proportions(table(d.hcpc$clust)),2)
names(d.hcpc)[names(d.hcpc)=="clust"]<-"clust.hcpc"

## comb combines classification from hcpc, kmeans and respondent id variables
# first save result form hcpc
comb<-d.hcpc[,c(colnames(motiv.id),"clust.hcpc")]
comb<-merge(motiv.id,comb[,c("newid","clust.hcpc")],by="newid",all.x=TRUE)
comb[,colnames(motiv)]<-data.frame(sapply(comb[,colnames(motiv)],FUN=function(x) as.integer(x)+1)) 

## combine k-means clustering results
# same order: d.mca/res.mca use mgw.data  
d.mca<-cbind(mgw.data[,colnames(motiv.id)],d.mca)

comb<-merge(comb,d.mca[,c("newid","Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","clust.kmeans")]) # equal to d1.mca

# compare cluster differences between HCPC and k-means
with(comb,round(proportions(table(clust.hcpc)),2))
with(comb,round(proportions(table(clust.kmeans)),2))



####################################
##### 4. Latent Class Regression ###
####################################

# process the data
dat1<-dat
dat1[,colnames(motiv)]<-data.frame(sapply(dat1[,colnames(motiv)],FUN=function(x)as.integer(x)+1))
dat1$residenceType2<-as.character(dat1$residenceType2)
dat1$residenceType2<-with(dat1, ifelse(!residenceType2 %in% c("countryside","urban"),"not_asked",residenceType2))
table(dat1$residenceType2)
dat1$residenceType2<-as.factor(dat1$residenceType2)
dat1$birth.k6<- ordered(dat1$birth.k6,levels=levels(dat1$birth.k6))


# 
type<-data.frame(foto=c("beite","gjen","gran"),Canvas.c=c(1,2,3))

# use photo choice from image.c (or canvas.c)
dat2<-merge(dat1,type,by.x="image.c",by.y="foto")
head(dat2)

## LCR formula, 7 motivational variables and 1 photo choice variable as item response indicators
# third photo choice (canvas.c) is the 8th item response indicators
f1<-cbind(aesthetic,access,forestry,grazing,biodiversity, culture, carbon,Canvas.c)~  birth0.k7


## Generate initial values from HCPC or K-means for LCR analysis in poLCA
n1<-grep("clust.hcpc",colnames(comb))
n2<-grep("clust.kmeans",colnames(comb))

start.hcpc<-start.lca(comb,n1)    # results of HCPC as starting values for LCR 
start.kmeans<-start.lca(comb,n2)  # results of k-means as starting values for LCR

# note p.start perhaps can be used as initial value for glca
results<-NULL
n<-1  # change n can be a high number so that ; run multiple times to see the stability of the results
for (i in 1:n){
	#seed<-runif(1,max=10000) # use this to check stability of results
	seed<-7826.186  		  # for reproduction of results in the paper 
	set.seed(seed) 
	nc<-3
	lc0<-poLCA(f1,dat2,nclass=nc, maxiter=5000,probs.start=start.hcpc,nrep=20,graphs=F,na.rm=FALSE,verbose=TRUE,calc.se=TRUE)
	p.start<-poLCA.reorder(lc0$probs.start,order(lc0$P,decreasing=TRUE)) # this order needs to be change manually
	lc2<-poLCA(f1,dat2,nclass=nc, maxiter=5000,probs.start=p.start,nrep=30,graphs=TRUE,na.rm=FALSE,verbose=TRUE,calc.se=TRUE) # continue in line 
	lc<-lc2 # 
	clas<-c(apply(lc2$posterior,2,mean))
	clas<-clas[order(clas)]
	results<-rbind(results,data.frame(seeds=seed,cla1=clas[1],cla2=clas[2],cla3=clas[3])) # mean
	print(i)
}



## Dignostics 
# birth cohort effect from logistic model of the LCR
lc2$coeff     

# Proportions of class membership by LCR
apply(lc2$posterior,2,mean) # 

# Check the fitting results
fit<-lc$predcell
with(fit,plot(observed,expected));lines(x = c(0,100), y = c(0,100))


## combine results from LCA clustering#######
# lc2 from LCR analysis below, note use same dat name as in poLCA (dat2)
lca<-data.frame(newid=dat2$newid,class=lc2$predclass,lc2$posterior)

names(lca)[names(lca)%in%c("X1","X2","X3")]<-c("prob.c1","prob.c2","prob.c3")
names(lca)[names(lca)=="class"]<-"clust.lca"
# note class 1 in lca may NOT same as classes in kmeans or hcpc
proportions(table(lca$clust.lca))

head(lca)

# change names after having combined results from lca
comb1<- merge(comb,lca,by="newid",all.x=TRUE)
comb1<-merge(comb1,dat1[,c("newid","birth0.k7","birth.k6","residenceType2","image.a","image.b","image.c")],all.x=TRUE)
comb1$clust.hcpc<- as.factor(comb1$clust.hcpc)
comb1$clust.lca<- as.factor(comb1$clust.lca)

head(comb1)

## compare cluster separability between latent class regression (LCA), k-means and Hierarchical Clustering on Principal Components(hcpc)
# results indicate that LCA results in better separability

p.lca<-ggplot(comb1, aes(x=Dim.1,y=Dim.2,colour=clust.lca))+geom_point()
p.kmeans<-ggplot(comb1, aes(x=Dim.1,y=Dim.2,colour=clust.kmeans))+geom_point()
p.hcpc<-ggplot(comb1, aes(x=Dim.1,y=Dim.2,colour=clust.hcpc))+geom_point()
ggarrange(p.lca,p.kmeans,p.hcpc,nrow=3)

## LCA with weighted posterior probabilities predicted by LCA.
# LCA predict class probability, hence can account for uncertainty.
lca1<-ggplot(comb1, aes(x=Dim.1,y=Dim.2,colour=prob.c1))+geom_point()
lca2<-ggplot(comb1, aes(x=Dim.1,y=Dim.2,colour=prob.c2))+geom_point()
lca3<-ggplot(comb1, aes(x=Dim.1,y=Dim.2,colour=prob.c3))+geom_point()

ggarrange(lca1,lca2,lca3,nrow=3)


######################################
####### 5. Plot the results ##########
######################################

#save(res.mca,lc2,comb1,file=paste0(res.root,"res_mca_lca2.Rdata"))
load(file=paste0(res.root,"res_mca_lca2.Rdata"),ver=TRUE)

###### Figure 3 ##############

# Figure 3a: plot motivational variables
plot.motiv<-function(res.mca){
d.var<-as.data.frame(res.mca$var$coord)[,1:2]  # extract results from MCA analysis
d.var<-cbind(d.var,contr.d1=res.mca$var$contrib[,1],contr.d2=res.mca$var$contrib[,2])
d.var$labs<-row.names(d.var)

d.var$labs<-gsub("FALSE",0,d.var$labs)
d.var$labs<-gsub("TRUE",1,d.var$labs)
d.var$labs<-as.factor(d.var$labs)
names(d.var)[1:2]<-c("dim1","dim2")
n<-unlist(lapply(gregexpr(pattern="_",d.var$labs),min))
d.var$Variables<-substr(d.var$labs,1,n-1)
d.var$contr.d<-with(d.var,pmax(contr.d1,contr.d2))
# first letter in capital
d.var$Variables<-paste(toupper(substr(d.var$Variables, 1, 1)), substr(d.var$Variables, 2, nchar(d.var$Variables)), sep="")

p1<-ggplot(d.var, aes(dim1, dim2, label = labs,size=contr.d,colour=Variables,shape=Variables)) + #
		  geom_text_repel(max.overlaps =20,show.legend = FALSE) +
		  #geom_label_repel(box.padding = 0.5, max.overlaps =20) +
		  geom_point(size=2.5) + 
		  xlim(-.6,2.5)+ylim(-.5,.92)+
		  scale_shape_manual(values=c(0:2,15:17,9,10,13))+
		  #scale_colour_gradientn(colours = rainbow(5))+
		  xlab("")+
		  ylab("")+
		  geom_hline(yintercept=0,linetype="dashed")+
		  geom_vline(xintercept=0,linetype="dashed")+
		  theme_classic(base_size = 14)+
		  theme(legend.position="right",legend.key.size = unit(.5, 'cm'))+
		  labs(col="",shape="")+
		  ggtitle("A. Motivation variables")+
		  theme(plot.title = element_text(size = 12),
				  plot.margin = margin(b=-0.4,t=.1,unit="cm"))# title size
p1
}

fig3a<-plot.motiv(res.mca)
fig3a


## Figure 3b: plot supplementary variables
plot.supp<-function(res.mca,dat){

d2<-as.data.frame(as.data.frame(cbind(res.mca$quali.sup$coord[,1:2],res.mca$quali.sup$cos2[,1:2],res.mca$quali.sup$v.test[,1:2])))
names(d2)<-c("d1.coord","d2.coord","d1.cos2","d2.cos2","d1.v","d2.v")

d2$labs<-row.names(d2)
d2$labs<-gsub("FALSE",0,d2$labs)
d2$labs<-gsub("TRUE",1,d2$labs)
m<-grep("^g\\.+",rownames(d2)) # matching .
d2$labs[m]<-substr(d2$labs[m],3,nchar(d2$labs[m]))

# dat is from raw data
d2$group<-c(rep("Upbringing",12),
				rep("profession",12),
				#rep("Gender",2),
				#rep("Region",length(levels(dat$region))),
				rep("Birth cohort",6),
				rep("Political party",length(levels(dat$norind))),				
				rep("Climate concern",length(levels(dat$concernClimate))),
				rep("Photo choice",length(levels(dat$image.a))),
				rep("Photo choice",length(levels(dat$image.b))),
				rep("Photo choice",length(levels(dat$image.c))),
				rep("Residence",3),
				rep("Education",length(levels(dat$edu.k3))),
				rep("Income",length(levels(dat$r14.income))))
d3<-d2

# only display most important ones for visibility
threshold<-2.576# 90% inertia threshold, v.test>=1.96 is similar to P value 95%

d3<-with(d3,d3[abs(d1.v)>=threshold | abs(d2.v)>=threshold,])

# rename some categories so that it fits better in the figure
d3$labs[grep("countryside",rownames(d3))]<-"rural"

n<-grep("image",rownames(d3))
d3$labs[n]<-substr(rownames(d3[n,]),7,nchar(rownames(d3[n,])))

d3$labs<-gsub("beite","grass",d3$labs)
d3$labs<-gsub("gran","spruce",d3$labs)
d3$labs<-gsub("gjen","reforest",d3$labs)

d3$labs<-gsub("grazi","graz",d3$labs)
d3$labs<-gsub("agric","agri",d3$labs)
#d3$labs<-gsub("19","",d3$labs)
d3$labs<-gsub("b_spruce","",d3$labs)
d3$labs<-gsub("a_spruce","",d3$labs)
d3$labs<-gsub("a_grass","",d3$labs)
d3$labs<-gsub("b_grass","",d3$labs)
d3$labs<-gsub("b_reforest","",d3$labs)
d3$labs<-gsub("a_reforest","",d3$labs)

#dim(d3)
#head(d3)

fig3b<-ggplot(d3, aes(d1.coord, d2.coord, label = labs,colour=group,shape=group),size=20) +
	  geom_text_repel(max.overlaps =40,show.legend = FALSE) +
	  #geom_label_repel(box.padding = 0.5, max.overlaps =50,show.legend = FALSE) +
	  geom_point(size=2) + 
	  xlim(-.6,2.25)+ylim(-.5,.8)+
	  scale_shape_manual(values=c(9,1,25,15:17,8,13,14))+ # symbol signs
	  xlab("")+# xlab(paste0("Dim1 (",round(res.mca$eig[1,2],1),"%)"))+ expression(symbol('\256'))
	  ylab("")+#ylab(paste0("Dim2 (",round(res.mca$eig[2,2],1),"%)"))+
	  geom_hline(yintercept=0,linetype="dashed")+
	  geom_vline(xintercept=0,linetype="dashed")+
	  theme_classic(base_size = 14)+
	  labs(col="",shape="")+ theme(legend.position="right")+
	  ggtitle("B. Supplementary variables")+
	  theme(plot.title = element_text(size = 12),plot.margin = margin(b=-0.5,t=0,unit="cm")) # l=left, b=bottom
	  
fig3b 
}

fig3b<-plot.supp(res.mca,dat)
fig3b

## Figure 3c: Kernel density plot 
n.col<-6
par(mar=c(3,3,1.2,.5)) # same setting as baseplot so that the figures can be aligned

baseplot<-function(d){
	par(mar=c(3,3,1.2,.5),oma=c(0,0,0,5))
	with(d,image(kde2d.weighted(x=Dim.1,y=Dim.2,w=prob.c3,n=100),xlim=c(-.6,2.2),ylim=c(-.8,.8),col=colorRampPalette(c("white","blue"))(n.col),
						xlab="",ylab="",main="",labels=FALSE,tick=FALSE))
	axis(side=1,labels=TRUE,cex.axis=.9,col.axis="grey30",lwd.ticks=.01,mgp=c(2,.5,0),tck=-0.01)
	axis(side=2,las=2,cex.axis=.9, mgp=c(2,.5,0),col.axis="grey30",lwd.ticks=.01,tck=-0.01) # add own axis, mgp=c(3,.75,0) to adjust label position
	with(d,image(kde2d.weighted(x=Dim.1,y=Dim.2,w=prob.c1,n=100),col=add.alpha(colorRampPalette(c("white","red"))(n.col),.5),add=T))
	with(d,image(kde2d.weighted(x=Dim.1,y=Dim.2,w=prob.c2,n=100),col=add.alpha(colorRampPalette(c("white","green"))(n.col),0.5),add=T))
	with(d,contour(kde2d.weighted(x=Dim.1,y=Dim.2,w=prob.c1,n=100),col="grey",add=T,nlevels=5))
	with(d,contour(kde2d.weighted(x=Dim.1,y=Dim.2,w=prob.c2,n=100),col="green",add=T,nlevels=5))
	with(d,contour(kde2d.weighted(x=Dim.1,y=Dim.2,w=prob.c3,n=100),col="blue",add=T,nlevels=5))

	box()
	abline(h=0,lty="dashed"); abline(v=0,lty="dashed"); box()

	with(d,sizeplot(Dim.1,y=Dim.2,pch=16,add=T,scale=0.2,col="grey50")) # the size of dots show the overlapping observations
	#text(-.4,-.7,"Class 1",cex=1.0) #Agriculture \n
	#text(-.4,.7,"Class 2",cex=1.0) #Recreation \n
	#text(1,.7,"Class 3",cex=1.0) #Forestation \n
	title(main=list("C. Class density weighted by posterior probabilities", cex=1,font=1),adj = .40, line = .5)
	legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
       c("Individuals","Class 1","Class 2","Class 3"), pch=19,col=c("grey60","pink","green","#3333FF"),bg="black",
	   y.intersp=6,pt.cex=c(1.3,2,2,2),cex=.9) #y.intersp adjust legend vertical space
	
}

## combine all three subfigures in Figure 3
baseplot(comb1)
fig3c<-recordPlot()

Figure3<-plot_grid(fig3a,fig3b,fig3c,nrow=3,rel_heights=c(1.3,1.95,1.3),align=c("v"))
Figure3

# add common x, y labels and dimension arrows					
Figure3 + annotate("segment", 
           x = 0.03, xend = 0.03, y=.3,  yend = 0.7,arrow=arrow(length=unit(0.3, "cm"),ends="both"))+
		   annotate("text", x = .02, y = 0.22, label = "More grazing",angle=90)+
		   annotate("text", x = .02, y = 0.8, label = "More aesthetics",angle=90)+
		   annotate("text", x = .01, y = 0.5, label = paste0("Dim2(",round(res.mca$eig[2,2],1),"%)"),angle=90)+
		   # x axis
		   annotate("segment", x = 0.3, xend = 0.6, y=0.02,  yend = 0.02,arrow=arrow(length=unit(0.3, "cm"),ends="both"))+
		   annotate("text", x = .74, y = 0.025, label = "High envir. concern")+
		   annotate("text", x = .17, y = 0.025, label = "Low envir. concern")+	
		   annotate("text", x = .44, y = 0.035, label = paste0("Dim1(",round(res.mca$eig[1,2],1),"%)"))



#######Figure 4 ########

## Fig. 4A. class profile 
# note lcmodel now include image.a as item indicator
lcmodel <- reshape2::melt(lc$probs, level=2)

#lcmodel$Category<-with(lcmodel,ifelse(Var2=="Pr(1)","Not selected","Selected"))
perc<-paste0(round(lc$P*100),"%")
label<-c("Agriculture","Recreation","Forestation")
lcmodel$Var1<-str_to_title(lcmodel$Var1)
lcmodel$L2<-str_to_title(lcmodel$L2)
lcmodel$Var3<-with(lcmodel,paste0(Var1,ifelse(Var1=="Class 1: ",label[1],ifelse(Var1=="Class 2: ", label[2],label[3]))))


lcmodel$Indicators<-as.factor(lcmodel$Var2)
lcmodel$Indicators<-with(lcmodel,ifelse(Var2=="Pr(1)"& L2=="Canvas.c","Pastureland",
								 ifelse(Var2=="Pr(2)"& L2=="Canvas.c","Reforest",
								 ifelse(Var2=="Pr(3)"& L2=="Canvas.c","Spruce",
								 ifelse(Var2=="Pr(1)"& L2!="Canvas.c","Not selected","Selected"))))) 
# fixed the order of ggplot
lcmodel$Indicators<-factor(as.factor(lcmodel$Indicators), levels=c("Not selected","Selected","Pastureland","Reforest","Spruce"))
table(lcmodel$Indicators)
lcmodel$L2<-factor(as.factor(lcmodel$L2),levels=c("Access","Aesthetic","Biodiversity","Carbon","Culture","Forestry","Grazing","Canvas.c"))

head(lcmodel)
model<-lcmodel
fig4a<-ggplot(model, aes(fill=Indicators, y=value,x=L2))+		# Var2
		scale_fill_manual(values=c("lightgrey", "seagreen3","#88CCEE" ,"#999933" ,"#CC6677"))+
		geom_bar(position="stack", stat="identity") +
		facet_wrap(~Var3)+theme_bw()+
		theme(strip.background = element_rect(fill="white"))+
		theme(axis.text.x = element_text(angle = 90,size=11,hjust=0.95,vjust=.1),
			  strip.text = element_text(size=11),plot.margin = margin(b=-0.4,t=0.1,unit="cm"))+ # l=left, b=bottom)+ 
		ylab("Item response probabilities")+
		xlab("")+ggtitle("A. Class profiles and labels")
fig4a


## Fig. 4B: prediction the relationship between birth and class membership
# birth as the predictor of class membership

mat<-function(coefs){
	pidmat<- cbind(1,c(1:7)) #beite
	exb <- exp(pidmat %*% coefs)
	m1<-cbind(1,exb)/(1+rowSums(exb))
}

m1<-mat(lc$coeff)

birth<-unique(dat1[,c("birth0.k7","birth.k7")])
birth<-birth[order(birth$birth0.k7),]

dd<-as.data.frame(rbind(cbind(birth=1:7,coefs=m1[,1]),cbind(birth=1:7,coefs=m1[,2]),cbind(birth=1:7,coefs=m1[,3])))
dd$Class<-c(rep("Class 1",7),rep("Class 2",7),rep("Class 3",7))
dd<-merge(dd,birth,by.x="birth",by.y="birth0.k7",all.x=TRUE)

# plot figure 4b
fig4b<-ggplot(dd,aes(x=birth.k7,y=coefs,group=Class,color=Class))+
				geom_line(aes(linetype=Class),size=1)+
				geom_point(aes(shape=Class),size=2)+
				ylab("Predicted class probabilities")+ xlab("")+
				ggtitle("B. Birth as a predictor of candidate affinity class")+
				theme_bw()+ylim(0,1)+
				theme(axis.text.x = element_text(size=11,angle = 30,vjust=.6,hjust=.5),
					  plot.margin = margin(b=-.5,t=0.1,unit="cm")) # l=left, b=bottom))
fig4b


ggarrange(fig4a,fig4b,nrow=2)







###########################
## user-defined functions##
###########################


# function to assign kmean or hierclustering results to poLCA
start.lca<-function(d,coln){
	aesthetic<-proportions(table(d[,coln],d[,"aesthetic"]),1)
	access<-proportions(table(d[,coln],d[,"access"]),1)
	forestry<-proportions(table(d[,coln],d[,"forestry"]),1)
	grazing<-proportions(table(d[,coln],d[,"grazing"]),1)
	biodiversity<-proportions(table(d[,coln],d[,"biodiversity"]),1)
	culture<-proportions(table(d[,coln],d[,"culture"]),1)
	carbon<-proportions(table(d[,coln],d[,"carbon"]),1)


	# assign kmean or hierclustering results to poLCA
	prob<-list(aesthetic,access,forestry,grazing,biodiversity,culture,carbon) # hcpc
	prob
}

# kernel density 
kde2d.weighted <- function (x, y, w, h, n = n, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
      stop("data vectors must be the same length")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}

# add alpha
add.alpha <- function(COLORS, ALPHA){
 if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
 RGB <- col2rgb(COLORS, alpha=TRUE)
 RGB[4,] <- round(RGB[4,]*ALPHA)
 NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
 return(NEW.COLORS)
}

## plot k-mean cluster 
wss<- function(k,d){
  kmeans(d,centers=k, nstart=10)$tot.withinss # centers=number of class, nstart=multiple initial configurations and reports the best one
}

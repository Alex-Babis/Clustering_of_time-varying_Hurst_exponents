library(fGarch)
library(pracma)
library(data.table)
library(cluster)
library(dendextend)
library(lubridate)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(igraph)


setwd("C:/Users/alexb/Desktop/Teaching/cvicenie6")

### DATA PREPARATION ###
# https://www.ecb.europa.eu/stats/policy_and_exchange_rates/euro_reference_exchange_rates/html/index.en.html
raw_data <- read.table("eurofxref-hist.csv", 
                                 header = TRUE, 
                                 sep = ",",
                                 stringsAsFactors = FALSE,
                                 na.strings = "N/A")
head(raw_data)

# reversing order 
raw_data <- raw_data[nrow(raw_data):1, ]

# find beginning of 2018
beginning.data <- which(grepl(raw_data[,1], pattern = "2018-"))[1]

# find ending of 2020
end.data <-  tail(which(grepl(raw_data[,1], pattern = "2020-")), 1)

# filter out period of 2018-2020 (3 years)
data <- raw_data[beginning.data:end.data, ]  

str(data)

colSums(is.na(data))

# deleting columns with too much NA
data <- data[ , colSums(is.na(data)) == 0] 

# store datetime 
xr.date <- data[, 1] 

# deleting datetime from data
data <- data[, -1] 

# first difference of logarithms of data
data <- log(data[2:nrow(data), ]) - log(data[1:(nrow(data)-1), ]) 

summary(data)

# deleting first obs of data because of differencing data (first one is NA)
xr.date <- xr.date[-1]

# third data column is constant so we get rid of it
data <- data[, -3] 

# plotting data
for (ind in 1:ncol(data)) 
  plot(data[, ind], type = "l", main = paste(ind, colnames(data)[ind]))


### COMPUTING TIME-VARYING HURST EXPONENT ON STANDARDIZED RESIDUALS FROM AR+GARCH MODELS ###

# function returning p-value from all test performed by `garchFit` 
# function using summary(garch_fit)
pvalues <- function(y)
{
  m <- capture.output(summary(y))
  prvy <- which( m == "Standardised Residuals Tests:") + 2
  posledny <- which( m %like% "LM Arch Test") 
  m <- m[prvy:posledny]
  mm <- strsplit(m, split= " ")
  mmm <- lapply(mm,function(x){i <- x != ""
  a <- x[i]
  a <- a[sum(i)]
  return(a)
  })
  pvalue <- as.numeric(unlist(mmm))
  names(pvalue) <- c("normtest","normtest","R10","R15","R20","R^210",
                     "R^215","R^220","archtest")
  return(pvalue)
}

# function returning model with all tests for residuals above 5% and
# with the highest BIC
AGorder <- function(x)
{
  # AR order that will be tried
  ar <- c(0,1,2,3) 
  
  # garch model that will be tried
  gar <- rbind(c(1,0), 
               c(1,1), c(2,0), 
               c(2,1), c(1,2), c(3,0)) 
  
  pval_tab <- matrix(0, nrow = length(ar), ncol = nrow(gar))
  BIC_tab <- matrix(0, nrow = length(ar), ncol = nrow(gar))
  
  # looping through all possible combination of specified ar order and
  # garch order
  for (i in 1:length(ar))
    for (j in 1:nrow(gar)) {
      form <- as.formula(paste0("~arma(",ar[i],",0)+garch(",gar[j,1],",",gar[j,2],")"))
      fit_ret <- garchFit(formula = form,data = x,trace = FALSE)
      m <- capture.output(summary(fit_ret))
      p <- pvalues(fit_ret)
      pval_tab[i, j] <- min(p[3:9])
      BIC_tab[i, j] <- as.numeric(strsplit(m[which( m == "Information Criterion Statistics:") + 2], split = " ")[[1]][2])
    }
  
  p.val.max <- max(pval_tab)
  
  BIC_tab[pval_tab < 0.05] <- BIC_tab[1, 1] + 1
  
  ind <- arrayInd(which.min(BIC_tab), dim(BIC_tab))
  i <- ind[1,1]
  j <- ind[1,2]
  form <- as.formula(paste0("~arma(",ar[i],",0)+garch(",gar[j,1],",",gar[j,2],")"))
  fit_ret <- garchFit(formula = form,data = x,trace = FALSE)
  
  # returning form of model, standardized residuals and maximum of p-values 
  # for final model 
  return(list(form,
              as.numeric(residuals(fit_ret, standardize = TRUE)),
              p.val.max))
}

left.out <- c(4)

garch_residuals <- data

# looping through data - for each exchange rate finding best suitable model
for (i in 1:ncol(data)) {
  if (!(i %in% left.out)) {
    model <- AGorder(data[, i])
    print(i)
    print(colnames(data)[i])
    print(model[[3]])
    print(model[[1]])
    garch_residuals[, i] <- model[[2]]
  }
}

# left out exchange rate with small values, so its treated separately
i <- 4
model <- AGorder(10 * data[, i])
print(model[[1]])

garch_residuals[, i] <- model[[2]]/10

# exchange rates of IDR and ILS left out due to not single model was
# suitable for either one of them
garch_residuals <- garch_residuals[,-c(20,21)]

# function calculating time-varying Hurst exponent from data 
# by rolling window approach with parameter rw corresponding to that
timehe <- function(data ,rw)
{
  n <- nrow(data) # lenght of time series of exchange rates
  m <- ncol(data) # number of exchange rate
  list_he <- list() 
  
  pb <- txtProgressBar(min = 0, max = m, style=3)
  
  for ( i in 1:m)
  { 
    resmodel <- data[, i]  
    he <- c() 
    for (j in 1:(n-rw))
    {
      # calculation of hurst exponent (all 4 types)for specific time window from 
      # standardized residuals
      h <- hurstexp(resmodel[j:(j+rw)], d = 50, display = FALSE)
      h <- unlist(h)
      he <- cbind(he,h)
    }
    list_he[[i]] <- he
    
    setTxtProgressBar(pb, i)
  }
  return(list_he)
}

# calculating time-varying hurst exponents for windonw leght 252 days
HE_all <- timehe(garch_residuals, 252)


select_HE <- function(ind){
  HE <- c()
  for( i in 1:length(HE_all))
  {
    # creating matrix with 1 type of hurst exponent for all exchange rates
    HE <- rbind(HE, HE_all[[i]][ind,]) 
  }
  
  rownames(HE) <- colnames(garch_residuals)
  return(HE)
}


### Plots ###
xr.date <-  as.Date(xr.date, "%Y-%m-%d")

# preparing data for ggplots
tables.ts <- list()
for(i in 1:ncol(data))
{
  hexp <- rbind(matrix(NA,126,5),
                t(HE_all[[i]]),
                matrix(NA,126,5)
  )
  
  hexp <- hexp[, -ncol(hexp)]
  df <- data.frame(Time = year(xr.date)+ (yday(xr.date)-1)/365,
                   data = as.numeric(data[,i]),
                   Residuals = as.numeric(garch_residuals[,i]),
                   Hurst_exponent = as.numeric(hexp[,1]),
                   Corrected = as.numeric(hexp[,2]),
                   Empirical = as.numeric(hexp[,3]),
                   `Corrected empirical` = as.numeric(hexp[,4]))
  colnames(df)[c(2,4)] <- c(colnames(data)[i],"Hurst exponents")
  tables.ts[[i]] <- df
}

# function ploting 2 exchange rates, their standardized residuals and
# time-varying Hurst exponents to one plot
graph.he.ts <- function(x)
{
  plot.label <- "Hurst_exp"
  dir.create(file.path("grafy/", plot.label), showWarnings = FALSE) 
  
  tab.x <- tables.ts[[x[1]]]
  
  tab.x.long <- melt(tab.x[,1:4], id.vars = "Time", variable.name = "ts")
  
  colnames(tab.x)[4] <- "Simplified"
  
  facet_labeller_left <- function(variable, value) {
    c(colnames(tab.x)[2],"","")
  }
  
  p1 <- ggplot(data = tab.x.long, aes(x = Time, y = value, group = ts))+
    geom_line()+
    facet_wrap(~ts, nrow = 3, ncol = 1,scales = "free",
               labeller = as_labeller(facet_labeller_left))+
    geom_hline(data = data.frame(yint=0.5,ts=as.factor("Hurst exponents")), 
               aes(yintercept = yint), linetype = "dashed",alpha=0.4)+
    geom_line(data = data.frame(xint =tab.x$Time ,yint=tab.x$Simplified,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Simplified"))+
    geom_line(data = data.frame(xint =tab.x$Time ,yint=tab.x$Corrected,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Corrected"))+
    geom_line(data = data.frame(xint =tab.x$Time ,yint=tab.x$Empirical,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Empirical"))+
    geom_line(data = data.frame(xint =tab.x$Time ,yint=tab.x$Corrected.empirical,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Corrected empirical"))+
    theme(legend.position = c(0.103, 0.21),
          legend.background=element_blank(),
          strip.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15)) +
    scale_color_manual(values = c("Simplified" = "black",
                                  "Corrected" = "#FF6B6B",
                                  "Empirical" ="#4ECDC4",
                                  "Corrected empirical" = "grey"))
  
  p1$labels$colour <- ""
  p1$labels$y <- ""
  
  tab.y <- tables.ts[[x[2]]]
  
  tab.y.long <- melt(tab.y[,1:4], id.vars = "Time", variable.name = "ts")
  
  colnames(tab.y)[4] <- "Simplified"
  
  facet_labeller_right <- function(variable, value) {
    c(colnames(tab.y)[2],"","")
  }
  
  p2 <-ggplot(data = tab.y.long, aes(x = Time, y = value, group = ts))+
    geom_line()+
    facet_wrap(~ts, nrow = 3, ncol = 1,scales = "free",
               labeller = as_labeller(facet_labeller_right))+
    geom_hline(data = data.frame(yint=0.5,ts=as.factor("Hurst exponents")), 
               aes(yintercept = yint), linetype = "dashed",alpha=0.4)+
    geom_line(data = data.frame(xint =tab.y$Time ,yint=tab.y$Simplified,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Simplified"))+
    geom_line(data = data.frame(xint =tab.y$Time ,yint=tab.y$Corrected,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Corrected"))+
    geom_line(data = data.frame(xint =tab.y$Time ,yint=tab.y$Empirical,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Empirical"))+
    geom_line(data = data.frame(xint =tab.y$Time ,yint=tab.y$Corrected.empirical,
                                ts=as.factor("Hurst exponents")), 
              aes(x = xint,y = yint, color = "Corrected empirical"))+
    theme(legend.position = c(0.103, 0.21),
          legend.background=element_blank(),
          strip.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15))  +
    scale_color_manual(values = c("Simplified" = "black",
                                  "Corrected" = "#FF6B6B",
                                  "Empirical" ="#4ECDC4",
                                  "Corrected empirical" = "grey"))
  
  p2$labels$colour <- ""
  p2$labels$y <- ""
  
  ggplot2::ggsave(paste0("grafy/",plot.label,"/",plot.label,"priebeh1.png"),
                  grid.arrange(p1, p2, ncol=2),width = 2400,
                  height = 1200,units = "px",
                  scale = 2,
                  dpi = 300)
  
  
  ggplot2::ggsave(paste0("grafy/",plot.label,"/",plot.label,"priebeh2.tiff"),
                  grid.arrange(p1, p2, ncol=2),width = 2400,
                  height = 1200,units = "px",
                  scale = 2,
                  dpi = 300)
}
# PNL and CZK exchange rate plot
graph.he.ts(c(7,3))

# ploting 4 time series of exchgane rate to one plot 
graph.ts.exchangerate <- function(x)
{
  plot.label <- "Grafy_priebeh"
  
  df <- data.frame(Time = year(xr.date)+ (yday(xr.date)-1)/365,
                   data1 = as.numeric(data[,x[1]]),
                   data2 = as.numeric(data[,x[2]]),
                   data3 = as.numeric(data[,x[3]]),
                   data4 = as.numeric(data[,x[4]])
  )
  colnames(df) <- c("Time",colnames(data)[x])
  
  df.long <- melt(df, id.vars = "Time", variable.name = "ts")
  
  p <- ggplot(data = df.long, aes(x = Time, y = value, group = ts))+
    geom_line()+
    facet_wrap(~ts, nrow = 2, ncol = 2,scales = "free")+
    theme(strip.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15)
    )
  
  p$labels$y <- ""
  p$labels$x <- ""

  ggplot2::ggsave(paste0("grafy/",plot.label,"/",plot.label,"priebeh",paste0(x,collapse = ""),".png"),
                  p,width = 2400,
                  height = 1200,units = "px",
                  scale = 2,
                  dpi = 300)
  
  ggplot2::ggsave(paste0("grafy/",plot.label,"/",plot.label,"priebeh",paste0(x,collapse = ""),"_2.tiff"),
                  p,width = 2400,
                  height = 1200,units = "px",
                  scale = 2,
                  dpi = 300)
}

graph.ts.exchangerate(c(1,3,18,28))


### CLUSTERING ###

# setting possible number of cluster and determining optimal via sulhouette criterium
number.of.clusters <- function(hc, plot = FALSE){
  pocet <- 2:5
  sil_mean <- rep(0, length(pocet))
  for (i in 1:length(pocet)){
    cut_hc <- cutree(hc, k = pocet[i]) 
    clust1 <- cut_hc
    sil <- silhouette(cut_hc, d)
    sil_mean[i] <- mean(sil[,3])}
  if(plot) plot(pocet, sil_mean, type = "b")
  return(pocet[which.max(sil_mean)])
}

clust <- list()
hurst.names <- c("Simplified","Corrected",
                 "Empirical","Corrected empirical")

# clustering time-varying Hurst exponents + plotting results
for(i in 1:4){
  print(paste("* * *", i, "* * *"))
  HE <- select_HE(i)
  matplot(t(HE), type="l", lty = 1)
  d <- dist(as.data.frame(HE))
  hc <- hclust(d^2 , method = "ward.D")
  plot(hc, hang = -1, cex = 0.6)
  opt <- number.of.clusters(hc, plot = TRUE)
  
  pal <- c("#FF6B6B", "#4ECDC4", "#556270","orange","red","blue","grey")
  plot.label <- "HC"
  
  dend <- as.dendrogram(hc)
  dend <- dend %>%
    color_labels(dend, k = opt, col = pal[1:opt]) %>%
    color_branches( col = pal[1:opt], k = opt) 
  
  
  tiff(filename=paste0("grafy/",plot.label,"/",plot.label,"priebeh",paste0(i,collapse = ""),"d2.tiff"),
       pointsize=10, width=2600, height=1600, res = 300)
  dend%>%
    plot(yaxt = "n")
  dend %>% rect.dendrogram(k=opt, 
                           border = 8, lty = 5, lwd = 2)
  dev.off()
  
  cut_hc <- cutree(hc, k = opt)  
  clust[[i]] <- cut_hc
}


# determining number of times each pair ended up in the same cluster
matches <- matrix(0, nrow = ncol(garch_residuals), ncol = ncol(garch_residuals))
for (i in 1:length(clust)) matches <- matches + outer(clust[[i]], clust[[i]], "==") 

# creating graph from matches matrix
g <- graph_from_adjacency_matrix(matches, weighted = TRUE, diag = FALSE, mode = "undirected")

# table of distribution of weights
table(E(g)$weight)

# plotting results via igraph
set.seed(135)
lay <- layout_nicely(g)
V(g)$label.cex <- 2.5
E(g)$lty <- c(2,1,2,1)[E(g)$weight]
width.e <- c(1,2,5,5)[E(g)$weight]
plot.label <- "siet"

# outputting plots
png(filename=paste0("grafy/",plot.label,"/",plot.label,"priebeh_vsetky.png"),
    pointsize=10, width=3200, height=3200,res = 300)
plot(g, 
     edge.width = width.e, 
     edge.color =  c("lightgrey", "grey", "#4ECDC4", "#FF6B6B")[E(g)$weight],
     layout = lay,
     vertex.size = 0,
     vertex.label.color = "black")
dev.off()


# creating graph from matches matrix only for weights equalt to 4
g4 <- graph_from_adjacency_matrix(matches == 4, diag = FALSE, mode = "undirected")

# plotting results
V(g4)$label.cex <- 2.5

plot.label <- "siet"

png(filename=paste0("grafy/",plot.label,"/",plot.label,"priebeh_iba4.png"),
    pointsize=10, width=3200, height=3200,res = 300)
plot(g4, 
     edge.width = 5, 
     edge.color = "#FF6B6B",
     layout = lay,
     vertex.size = 0,
     vertex.label.color = "black")
dev.off()


# creating graph from matches matrix only for weights higher or equalt to 3
matches3 <- matches
matches3[matches < 3] <- 0 

g3 <- graph_from_adjacency_matrix(matches3, weighted = TRUE,
                                  diag = FALSE, mode = "undirected")

# plotting results
V(g3)$label.cex <- 2.5
E(g3)$lty <- c(0,0,2,1)[E(g3)$weight]
width.e <- c(0,0,5,5)[E(g3)$weight]
plot.label <- "siet"

png(filename=paste0("grafy/",plot.label,"/",plot.label,"priebeh_viac_3.png"),
    pointsize=10, width=3200, height=3200,res = 300)
plot(g3, 
     edge.width = width.e, 
     edge.color =  c("lightgrey", "grey", "#4ECDC4", "#FF6B6B")[E(g3)$weight],
     layout = lay,
     vertex.size = 0,
     vertex.label.color = "black")
dev.off()

# computing components of G4 graph 
com <- components(g4)
komponenty <- which(com$csize > 1)
for(i in 1:length(komponenty)){
  print(paste("* * * ", i, "* * * "))
  print(names(which(com$membership==komponenty[i])))
}

# computing components of G3 graph
com <- components(g3)
komponenty <- which(com$csize > 1)
for(i in 1:length(komponenty)){
  print(paste("* * * ", i, "* * * "))
  print(names(which(com$membership==komponenty[i])))
}

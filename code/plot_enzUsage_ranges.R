# devtools::install_github("hrbrmstr/ggalt")
library(ggplot2)
library(ggalt)
#theme_set(theme_classic())
current        <- '/Users/ivand/Documents/GitHub/heme_production_yeast/results'
fileName       <- paste(current,'/enzUsageRanges_hemeGenes.txt',sep='')
enzUsageRanges <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
#Transform computationally-introduced negative values
enzUsageRanges$maxUsages[which(enzUsageRanges$maxUsages<0)] <- 0
enzUsageRanges$minUsages[which(enzUsageRanges$minUsages<0)] <- 0
#
inconsistent <- which(enzUsageRanges$maxUsages< enzUsageRanges$minUsages)
enzUsageRanges$maxUsages[inconsistent] <- enzUsageRanges$minUsages[inconsistent]
#Recalculate range
enzUsageRanges$range   <- enzUsageRanges$maxUsages-enzUsageRanges$minUsages
#Get relative FV range
enzUsageRanges$percVar <- enzUsageRanges$range/enzUsageRanges$maxUsages
#Order dataset acording to descendent range
df                     <- enzUsageRanges[order(-enzUsageRanges$range,-enzUsageRanges$maxUsages),]
#df                     <- df[x>0,]
#df                     <- df[df$percVar<1,]
df$ordered             <- seq(1,nrow(df))
#df$enzNames <- factor(enzUsageRanges$enzNames, levels=as.character(enzUsageRanges$ordered)) 
# health$Area <- factor(health$Area)
gg <- ggplot(df) + geom_dumbbell(aes(x=df$minUsage, xend=df$maxUsage, y=df$ordered),
                 colour="grey80",size=0.5,colour_x = 'blue',size_x= 0.7,colour_xend= 'red',size_xend=0.5) + 
      labs(x = "Usage ranges [mmol/gDw]", y = "Enzymes") +
  theme(title = element_text(size = rel(1.4))) + scale_x_log10() + #scale_x_log10(limits = c(1E-8,1E-4))
  geom_text(aes(x = 6E-13, y = ordered[nrow(df)-30]),label = "minUsage",color = "black", vjust = -1,size = 7) + 
  geom_text(aes(x = 1E-5, y = ordered[nrow(df)-30]),label = "maxUsage",color = "black", vjust = -1,size=7) + 
  theme(axis.text.x = element_text(color="black", size=14),axis.text.y = element_text(color="black", size=14),
  plot.background = element_rect(fill = "white"),panel.background = element_blank(),axis.line = element_line(colour = "black"),
  panel.grid.major = element_line('gray'))
plot(gg)

#Get a histogram 
#Remove NaN and Inf values
dataset <- df
to.keep <- !(is.na(dataset$percVar))
dataset <- dataset[to.keep,]
#plot data
bars    <- ggplot (dataset, aes(x=dataset$percVar)) + 
           geom_histogram(binwidth = 0.05,color="black", fill="light blue") + 
           labs(y = 'Frequency', x='Relative enzUSage variability range') + theme_classic(base_size = 14)
plot(bars)
#Get a CDF plot
to.keep     <- !(is.na(dataset$percVar))
dataset     <- dataset[to.keep,]
cdfP        <- ggplot (data=dataset, aes(x= dataset$percVar)) + 
  stat_ecdf(geom = "step",color = 'blue') + 
  labs(y = 'Relative frequency', x='Relative enzUSage variability range') +xlim(0,1) +
  theme_classic(base_size = 14)

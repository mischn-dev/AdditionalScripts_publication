library(corrplot)
library(lubridate)
library(Stack)
we = read.csv("Kon_Evapotranspiration96-19_180nFk.csv")
# the radiation is not correct for the first years until 2009 - might be the column mismatched - but the calculation of ET0 is fine!
dic = read.csv("evapotranspiration_calc.csv")
dic$Date = as.Date(dic$Date)

we$Kanal = as.Date(we$Kanal)
we$year = year(we$Kanal)

wed = we[we$year<2010,]
wek = we[we$year>2010,]

weed = merge(dic, wed, by.x = "Date", by.y = "Kanal")
wed$Radiation = weed$Radiation.x

we = rbind(wed, wek)

## corelation plot
png("correlation_plotClimateData.png", res = 100, width = 1250, height = 1250)
par(mfrow=c(1,2))
corrplot(cor(na.omit(wek)[,c(2:11)]), method="circle", type="lower", tl.pos = "lt", tl.col = "black", tl.srt=45)
corrplot(cor(na.omit(wek)[,c(2:11)]), method="number", type="upper", add=T, tl.pos = "n",  tl.col = "black", tl.srt=45, col="black")
text(0,12,"A",cex=2,font=2)

corrplot(cor(na.omit(wed)[,c(2:11)]), method="circle", type="lower", tl.pos = "lt", tl.col = "black", tl.srt=45)
corrplot(cor(na.omit(wed)[,c(2:11)]), method="number", type="upper", add=T, tl.pos = "n",  tl.col = "black", tl.srt=45, col="black")
text(0,12,"B",cex=2,font=2)
dev.off()


we$year = as.factor(year(we$Kanal))
aov(we$ETa~we$year+we$Humidity_min+we$Temp_Max+we$Windspeed+we$SWC+we$precipitation+we$Temp_min+we$TSum+we$GrowthStage+we$Radiation)
tuk = TukeyHSD(aov(we$ETa~we$year))
ty = as.data.frame(tuk$`we$year`)
ty[ty$`p adj`<0.05,]

# running a linear mixed model instead of a linear model
library(lmerTest)
we2 = na.omit(we)
summary(lmer(we2$ETa ~ (1|we2$Humidity_min) + we2$GrowthStage + (1|we2$Temp_Max) + (1|we2$Radiation) + we2$year + (1|we2$SWC) + we2$Drought))
# create groups for the years 
library(agricolae)
library(ggplot2)
year = HSD.test(lm(we2$ETa ~ we2$Humidity_min + we2$GrowthStage + we2$Temp_Max + we2$Radiation + we2$year + we2$SWC + we2$Drought, data = we2), "we2$year", group=T)
qa = year$means
colnames(qa)[1] = "ETa"
qa$year = row.names(qa)
qa = merge(qa, year$groups, by = "row.names")
ggplot(qa, aes(x =year, y = ETa, color = groups)) + geom_point() + geom_errorbar(aes(ymin = Q25, ymax = Q75)) + stat_summary(geom = 'text', label = qa$groups, fun = max, vjust = 3) + xlab("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.6, size = 12), axis.text.y = element_text(size =12),strip.text.x = element_text(size = 12, face="bold.italic"), legend.position = "none")
ggsave("ETa_year_grouping_ggplot.png", width = 15, height = 10, dpi = 200, units = "in", scale = 0.6)

png("ETa_year_grouping.png", width = 15, height = 5, res=100, units = "in")
par(bg = "white")
plot(year, main="", ylab = "ETa")
dev.off()

# deviation of ETa to ETc
we2$ETdev = we2$ETa-we2$Etc
year = HSD.test(lm(we2$ETdev ~ we2$Humidity_min + we2$GrowthStage + we2$Temp_Max + we2$Radiation + we2$year + we2$SWC + we2$Drought, data = we2), "we2$year", group=T)
qa = year$means
colnames(qa)[1] = "ETa"
qa$year = row.names(qa)
qa = merge(qa, year$groups, by = "row.names")
ggplot(qa, aes(x =year, y = ETa, color = groups)) + geom_point() + geom_errorbar(aes(ymin = -std, ymax = Q75)) + stat_summary(geom = 'text', label = qa$groups, fun = max, vjust = 3) + xlab("") + ylab(expression(delta~ET[textstyle("A")]~~ET[textstyle("C")])) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.6, size = 12), axis.text.y = element_text(size =12),strip.text.x = element_text(size = 12, face="bold.italic"), legend.position = "none")
ggsave("ETa_year_grouping_deltaET.png", width = 15, height = 10, dpi = 200, units = "in", scale = 0.6)



## create plots for the deviation of the mean for all major climatic relevant measurements
# eta
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
                              qa[qa$V1==i,2] = sum(we[we$year==i, 22])
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = sum(wed$ETa, na.rm = T) / length(unique(wed$year))
  qa[qa$V1>2010,3] = sum(wek$ETa, na.rm = T) / length(unique(wek$year))
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "Evapotranspiration sum of growth period [mm]"
  
  set = qa

# tempmax
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = mean(we[we$year==i, 2], na.rm = T)
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = mean(wed$Temp_Max, na.rm = T)
  qa[qa$V1>2010,3] = mean(wek$Temp_Max, na.rm = T)
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "maximum daily Temperature [?C]"
  
  set = rbind(set,qa)

  # tempmin
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = mean(we[we$year==i, 3], na.rm = T)
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = mean(wed$Temp_min, na.rm = T)
  qa[qa$V1>2010,3] = mean(wek$Temp_min, na.rm = T)
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "minimum daily Temperature [?C]"
  
  set = rbind(set,qa)
  
  # tempmean
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = mean(we[we$year==i, 8], na.rm = T)
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = mean(wed$Temp_Mean, na.rm = T)
  qa[qa$V1>2010,3] = mean(wek$Temp_Mean, na.rm = T)
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "average daily Temperature [?C]"
  
  set = rbind(set,qa)
    
  # radiation - the radiation mena has to be calculated seperated for the two weather stations and periods, the variation ist just too high
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = sum(we[we$year==i, 5], na.rm = T)
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = sum(wed$Radiation, na.rm = T) / length(unique(wed$year))
  qa[qa$V1>2010,3] = sum(wek$Radiation, na.rm = T) / length(unique(wek$year))
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "Radiation sum in growth period [W/m?]"
  
  set = rbind(set,qa)
  
  # humidity # the humidity mean has to be calucalted seperate for the two weather stations - the variation is just too high
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = mean(we[we$year==i, 6])
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = mean(wed$Humidity_min, na.rm = T)
  qa[qa$V1>2010,3] = mean(wek$Humidity_min, na.rm = T)
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "minimum Humidity [%]"
  
  set = rbind(set,qa)
  
# precipitation
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = sum(we[we$year==i, 7])
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = sum(wed$precipitation, na.rm = T) / length(unique(wed$year))
  qa[qa$V1>2010,3] = sum(wek$precipitation, na.rm = T) / length(unique(wek$year))
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "Precipitation sum in growth period [mm]"
  
  set = rbind(set,qa)
  
# days of drought
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = sum(we[we$year==i, 21]=="stressed")
  }
  qa$V3= sum(we$Drought=="stressed") / length(unique(we$year))
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "days of Drought [d]"
  
  set = rbind(set,qa)
  
# diverence of eta and etc 
  we$etd = we$Etc-we$ETa
  qa = as.data.frame(matrix(ncol=2, nrow=length(unique(we$year))))
  qa$V1 = unique(we$year)
  for (i in unique(we$year)){
    qa[qa$V1==i,2] = sum(we[we$year==i, 24])
  }
  qa$V3 = 0
  qa$V1 = as.integer(as.character(qa$V1))
  qa[qa$V1<2010,3] = sum(we[we$year<2010,]$etd, na.rm = T) / length(unique(we[we$year<2010,]$year))
  qa[qa$V1>2010,3] = sum(we[we$year<2010,]$etd, na.rm = T) / length(unique(we[we$year>2010,]$year))
  qa$V1 = as.factor(qa$V1)
  qa$V4 = "Sum of ETa to ETc variation of growth period [mm]"
  
  set = rbind(set,qa)
  
  

library(ggplot2)
p1 = ggplot(set, aes(x = V1, y = V2)) +
  #geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
  geom_segment(aes(xend = V1, yend = V3), alpha = .8) +  # alpha to fade lines
  geom_point()  + xlab("Year") + ylab("") + 
  geom_point(aes(y = V3), shape = 1) + facet_wrap(~ V4, scales = "free_y") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.6), strip.text.x = element_text(size = 12, face="bold.italic"))# Add theme for cleaner look

ggsave("Ubersicht_Wetterdaten_Jahrweise.png", p1, dpi=300, width = 16, height = 8, units = "in")



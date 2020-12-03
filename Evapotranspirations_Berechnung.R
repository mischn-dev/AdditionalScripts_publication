#### calculation of the evapotranspiration with penman-monteith equation
library(data.table)

for ( i in list.files("../Rohdaten/")){
  
    filename1 = paste0("../Rohdaten/", i)
    we = read.csv2(filename1, header=T)
    print(paste("Processing file", i))
    
    # remove the empty columns and replace the headers
    we = we[,-grep("X",colnames(we))]
    we = apply(we,2,as.character)
    
    colnames(we) = c("Kanal", "Luftdruck", "Temp_DruckSensor", "Solarspannung", "Batteriespannung", "Wechselspannung", "relFeuchte", "Luftemperatur", "Windgeschw", 
                     "Windrichtung", "StdWindgeschw", "StdWindrichtung", "ErdoberflTemp", "Temp3", "Temp4", "Temp5", "Globalstrahlung", "PAR_Strahlung", "UVB_Strahlung", 
                     "NiederschlagMenge", "NiederschlagArt", "Niederschlagstatus", "Niederschlagsdauer", "Taupunkt")
    we = we[-c(1,2),] # remove 1st + 2nd row
    
    # there might be some values with "[]" surrounding them, these have to be removed
    we[,20] = gsub("\\[|\\]", "", we[,20])

# save and reload the file
    outfile1 = paste0("../cleaned_files/", i)
    write.csv2(we, outfile1, quote = F, row.names = F)
    we = fread(outfile1, sep=";", dec=",")

    # convert the date entries to date and remove the time values
    we[,1] = as.Date(as.data.frame(we[,1])[,1], "%d.%m.%Y %H:%M")
    # remove the last row - 24:00 of the last day will be allocated to next month, not wanted
    we = we[-nrow(we),]
    # group by the date with data.table function to new dataframe ## base values
    we2 = we[,.(Temp_Max=max(Luftemperatur), Temp_min = min(Luftemperatur), Windspeed = mean(Windgeschw), Radiation = mean(Globalstrahlung)*0.0864, Humidity_max = max(relFeuchte), Humidity_min = min(relFeuchte), precipitation=sum(NiederschlagMenge)), by=Kanal]
    we2$Temp_Mean = (we2$Temp_Max + we2$Temp_min)/2
    # slope of saturation vapor pressure curve (delta)
    we2$delta = (4098+(0.6108*exp((17.27*we2$Temp_Mean)/(we2$Temp_Mean+237.3))))/(we2$Temp_Mean +237.3)^2
    # atmospheric pressure 
    z = 150 # height above sea level
    # psychometric constant y
    we2$psych_constant = 0.000665 * 101.3 * ((293 - 0.0065 * z) / 293)^5.26 
    # delta term radiation
    we2$radiationTerm = we2$delta / (we2$delta + we2$psych_constant * (1 + 0.34 * we2$Windspeed))
    # delta term wind
    we2$windTerm =  we2$psych_constant / (we2$delta + we2$psych_constant * (1 + 0.34 * we2$Windspeed))
    # delta term temperature
    we2$TempTerm = (900 / (we2$Temp_Mean + 273)) * we2$Windspeed
    ## mean saturation vapor pressure
    # e(Tmax)
    we2$eTmax = 0.6108 * exp((17.27 * we2$Temp_Max) / (we2$Temp_Max + 237.3))
    # e(Tmin)
    we2$eTmin = 0.6108 * exp((17.27 * we2$Temp_min) / (we2$Temp_min * 237.3))
    # es
    we2$es = (we2$eTmax + we2$eTmin) / 2
    # ea
    we2$ea = (we2$eTmin * (we2$Humidity_max/100) + we2$eTmax * (we2$Humidity_min / 100)) / 2
    ## radiation calculation 
    # inverse rekultive distance earth to sun - dr
    we2$doy = as.integer(strftime(we2$Kanal, format = "%j")) # calculate the day of the year
    we2$dr = 1+ 0.033 * cos((2* pi) / 365 * we2$doy)
    # solar declination d
    we2$SolarD = 0.409 * sin(((2*pi) / 365) * we2$doy - 1.39)
    # latitude to radius global position of klein altendorf
    ld = (pi / 180) * 50
    # sunset hour angle
    we2$ws = acos(-tan(ld) * tan(we2$SolarD))
    # extra trerrestrial radiation (Ra)
    we2$Ra = ((24*60) / pi) * 0.0820 * we2$dr * ((we2$ws * sin(ld) * sin(we2$SolarD)) + (cos(ld) * cos(we2$SolarD) * sin(we2$ws)))
    # clear sky solar radiation Rso
    we2$Rso = (0.75 + 2* exp(1) * 10^-5 * z) * we2$Ra
    # net solar short wave radiation
    we2$Rns = (1-0.23) * we2$Radiation # 0.23 is the albeldo
    # net solar long wave radiation
    we2$Rnl = (4.903 * 10^-9) * (((we2$Temp_Max + 273.16)^4 + (we2$Temp_min + 273.16)^4) / 2) * (0.34 - 0.14 * sqrt(we2$ea)) * (1.35 * (we2$Radiation / we2$Rso) - 0.35)
    # net radiation Rn as equivalent to evaporation 
    we2$Rng = 0.408 *  (we2$Rns - we2$Rnl)
    
    # Radiation term
    we2$ETrad = we2$radiationTerm * we2$Rng
    # Wind term
    we2$ETwind = we2$windTerm * we2$TempTerm * (we2$es - we2$ea)
    # Evapotranspiration value ETo
    we2$ETo = we2$ETwind + we2$ETrad
    
    ## save the create file in the Wetterdaten folder
    
    fwrite(we2, "Klein_Altendorf_2011-2019.csv", row.names = F, append = T, quote = F)
    
}

# file cleaning - remove duplicate dates from the file
wetter = read.csv("Klein_Altendorf_2011-2019.csv")
dup = wetter$Kanal[duplicated(wetter$Kanal)]
for (i in dup){
               qw = which(wetter$Kanal==i)[1]
               wetter = wetter[-qw,]
}

fwrite(wetter, "Klein_Altendorf_2011-2019.csv", row.names = F, append = T, quote = F)

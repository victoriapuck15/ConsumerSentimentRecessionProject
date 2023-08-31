library(reshape2)
library(fredr)
library(tidyverse)
library(zoo)
library(xgboost)
library(Matrix)
library(caret)
library(readxl)
library(httr)
library(lubridate)


####################################
# Start

url <- "https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/coincident/coincident-revised.xls"
download.file(url = url, mode = "wb", destfile = "statecoincidentdata.xls")

statecoin<- read_excel("statecoincidentdata.xls")
statecoin<- statecoin[c(4:length(statecoin$Date)),]
Date<- as.character(statecoin$Date[-1])
statecoin<- as.data.frame(statecoin[,-1])
statecoin<- apply(statecoin, 2, diff)
statecoin<- as.data.frame(cbind(Date, statecoin)) %>% 
  melt(id="Date")
statecoin$Date<- as.yearmon(statecoin$Date, format="%Y-%m-%d")

ab    <- (c("AL",
            "AK", "AZ", "KS", "UT", "CO", "CT",
            "DE", "FL", "GA", "HI", "ID", "IL",
            "IN", "IA", "AR", "KY", "LA", "ME",
            "MD", "MA", "MI", "MN", "MS", "MO",
            "MT", "NE", "NV", "NH", "NJ", "NM",
            "NY", "NC", "ND", "OH", "OK", "OR",
            "PA", "RI", "SC", "SD", "TN", "TX",
            "CA", "VT", "VA", "WA", "WV", "WI",
            "WY", "US", "DC"))
st    <- c("Alabama",
           "Alaska", "Arizona", "Kansas",
           "Utah", "Colorado", "Connecticut",
           "Delaware", "Florida", "Georgia",
           "Hawaii", "Idaho", "Illinois",
           "Indiana", "Iowa", "Arkansas",
           "Kentucky", "Louisiana", "Maine",
           "Maryland", "Massachusetts", "Michigan",
           "Minnesota", "Mississippi", "Missouri",
           "Montana", "Nebraska", "Nevada",
           "New Hampshire", "New Jersey", "New Mexico",
           "New York", "North Carolina", "North Dakota",
           "Ohio", "Oklahoma", "Oregon",
           "Pennsylvania", "Rhode Island", "South Carolina",
           "South Dakota", "Tennessee", "Texas",
           "California", "Vermont", "Virginia",
           "Washington", "West Virginia", "Wisconsin",
           "Wyoming", "United States", "District of Columbia")
colp<- as.data.frame(cbind(st, ab))
names(colp)<- c("State", "variable")

statecoin <- left_join(statecoin, colp, by="variable")
statecoin<- statecoin[c(1,4,3)]
names(statecoin)<- c("Date", "State", "Coincident")


recessioning<- function(x) {
  
  y<- statecoin %>% 
    filter(State==x)
  x<- y %>% 
    select(Coincident)
  x<- as.numeric(x$Coincident)
  
  
  peaks<- list()
  troughs<- list()
  
  for(i in 3:(length(y$Date)-6)) {
    peaks[i]<- if(x[i]<0) {
      0} else if(x[i+1]>0) {
        0} else if(sum(x[i],x[i-1],x[i-2])<0) {
          0} else if(sum(x[i+1],x[i+2],x[i+3])>0) {
            0} else if(sum(x[i+4],x[i+5],x[i+6])>0) {
              0} else {1}
    
    troughs[i]<- if(x[i]>0) {
      0} else if(x[i+1]<0) {
        0} else if(sum(x[i],x[i-1],x[i-2])>0) {
          0} else if(sum(x[i+1],x[i+2],x[i+3])<0) {
            0} else if(sum(x[i+4],x[i+5],x[i+6])<0) {
              0} else {1}
  }
  
  peaks<- as.data.frame(t(as.data.frame(peaks[3:(length(y$Date)-6)])))
  names(peaks)<- "Peaks"
  rownames(peaks)<- NULL
  date<- y$Date[3:(length(y$Date)-6)]
  
  
  troughs<- as.data.frame(t(as.data.frame(troughs[3:(length(y$Date)-6)])))
  names(troughs)<- "Troughs"
  rownames(troughs)<- NULL
  
  numbers<- cbind(date, peaks, troughs)
  
  tabling2<- numbers %>% 
    mutate(Peaks=as.numeric(Peaks)) %>% 
    mutate(Troughs=as.numeric(Troughs)) %>% 
    mutate(Index= 1:length(numbers$date))
  
  
  # Condition 4, There Cannot be two troughs in a quarter or two peaks in a quarter
  for(i in 3:length(tabling2$date)) {
    if(tabling2$Troughs[i]==1 & tabling2$Troughs[i-2]==1) {
      tabling2$Troughs[i]<- 0
    }
  }
  for(i in 3:length(tabling2$date)) {
    if(tabling2$Peaks[i]==1 & tabling2$Peaks[i-2]==1) {
      tabling2$Peaks[i]<- 0
    }
  }
  # Condition 5, Choosing between two troughs in between peaks
  peaktablecheck<- tabling2 %>% 
    filter(Peaks==1)
  troughtablecheck<- tabling2 %>% 
    filter(Troughs==1)
  
  # Finding Competing Peaks
  if(peaktablecheck$Index[1]<troughtablecheck$Index[1]) {
    peaktablecheck$ind2[1]<- 0
  }
  for(i in 1:length(troughtablecheck$Trough)) {
    for(j in 1:length(peaktablecheck$Peak)) {
      if(peaktablecheck$Index[j]>troughtablecheck$Index[length(troughtablecheck$Trough)]) {
        peaktablecheck$ind2[j]<- length(troughtablecheck$Trough)
      } else if(peaktablecheck$Index[j]>troughtablecheck$Index[i] & peaktablecheck$Index[j]<troughtablecheck$Index[i+1]) {
        peaktablecheck$ind2[j]<- i
      }}
  }
  
  # Removing competing Peak losers
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  # This is working to determine competing troughs
  if(troughtablecheck$Index[1]<peaktablecheck$Index[1]) {
    troughtablecheck$ind2[1]<- 0
  }
  for(i in 1:length(peaktablecheck$Peak)) {
    for(j in 1:length(troughtablecheck$Trough)) {
      if(troughtablecheck$Index[j]>peaktablecheck$Index[length(peaktablecheck$Peak)]) {
        troughtablecheck$ind2[j]<- length(peaktablecheck$Peak)
      } else if(troughtablecheck$Index[j]>peaktablecheck$Index[i] & troughtablecheck$Index[j]<peaktablecheck$Index[i+1]) {
        troughtablecheck$ind2[j]<- i
      }}
  }
  
  # Removing competing trough losers, Peak Losers, multiple iterations
  for(i in 1:length(unique(troughtablecheck$ind2))) {
    if(is.na(troughtablecheck$ind2[i+1])) {
      print("")
    } else if(troughtablecheck$ind2[i] == troughtablecheck$ind2[i+1]) {
      if(sum(x[troughtablecheck$Index[i]:troughtablecheck$Index[i+1]])>0) {
        troughtablecheck<- troughtablecheck[-(i+1),]
      } else {
        troughtablecheck<- troughtablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(troughtablecheck$ind2))) {
    if(is.na(troughtablecheck$ind2[i+1])) {
      print("")
    } else if(troughtablecheck$ind2[i] == troughtablecheck$ind2[i+1]) {
      if(sum(x[troughtablecheck$Index[i]:troughtablecheck$Index[i+1]])>0) {
        troughtablecheck<- troughtablecheck[-(i+1),]
      } else {
        troughtablecheck<- troughtablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(peaktablecheck$ind2))) {
    if(is.na(peaktablecheck$ind2[i+1])) {
      print("")
    } else if(peaktablecheck$ind2[i] == peaktablecheck$ind2[i+1]) {
      if(sum(x[peaktablecheck$Index[i]:peaktablecheck$Index[i+1]])<0) {
        peaktablecheck<- peaktablecheck[-(i+1),]
      } else {
        peaktablecheck<- peaktablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(troughtablecheck$ind2))) {
    if(is.na(troughtablecheck$ind2[i+1])) {
      print("")
    } else if(troughtablecheck$ind2[i] == troughtablecheck$ind2[i+1]) {
      if(sum(x[troughtablecheck$Index[i]:troughtablecheck$Index[i+1]])>0) {
        troughtablecheck<- troughtablecheck[-(i+1),]
      } else {
        troughtablecheck<- troughtablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(troughtablecheck$ind2))) {
    if(is.na(troughtablecheck$ind2[i+1])) {
      print("")
    } else if(troughtablecheck$ind2[i] == troughtablecheck$ind2[i+1]) {
      if(sum(x[troughtablecheck$Index[i]:troughtablecheck$Index[i+1]])>0) {
        troughtablecheck<- troughtablecheck[-(i+1),]
      } else {
        troughtablecheck<- troughtablecheck[-(i),]
      }} 
  }
  for(i in 1:length(unique(troughtablecheck$ind2))) {
    if(is.na(troughtablecheck$ind2[i+1])) {
      print("")
    } else if(troughtablecheck$ind2[i] == troughtablecheck$ind2[i+1]) {
      if(sum(x[troughtablecheck$Index[i]:troughtablecheck$Index[i+1]])>0) {
        troughtablecheck<- troughtablecheck[-(i+1),]
      } else {
        troughtablecheck<- troughtablecheck[-(i),]
      }} 
  }
  # Above iterations should get rid of all competing points
  tabling3<- peaktablecheck %>% 
    filter(Peaks==1) %>% 
    select(date) %>% 
    mutate(date=as.character(date))
  names(tabling3)<- "Peak"
  tabling4<- troughtablecheck %>% 
    filter(Troughs==1) %>% 
    select(date) %>% 
    mutate(date=as.character(date))
  names(tabling4)<- "Trough"
  
  startdate<- as.data.frame(as.yearmon("May 1979", format="%b %Y"))
  startdate2<- as.data.frame("May 1979")
  nas<- as.data.frame(NA)
  names(nas)<- "Trough"
  names(startdate)<- "Peak"
  names(startdate2)<- "Peak"
  if(nrow(tabling4)>nrow(tabling3)) {
    tabling3<- as.data.frame(rbind(startdate2, tabling3))
  }
  if(nrow(tabling3)>nrow(tabling4)) {
    tabling4<- as.data.frame(rbind(tabling4, nas))
  }
  
  tabling3<- tabling3 %>% 
    mutate(Peak=as.yearmon(Peak, "%b %Y"))
  tabling4<- tabling4 %>% 
    mutate(Trough=as.yearmon(Trough, "%b %Y"))
  
  if(tabling3[1,]>tabling4[1,]) {
    tabling3<- as.data.frame(rbind(startdate, tabling3))
    tabling4<- as.data.frame(rbind(tabling4, nas))
  }
  
  staterecessiontable<- as.data.frame(cbind(tabling3, tabling4))
  names(staterecessiontable)<- c("Peaks", "Troughs")
  staterecessiontable<- staterecessiontable %>% 
    mutate(Troughs= as.yearmon(Troughs, format="%b %Y"))
  currentdate<- as.yearmon(Sys.Date(), format="%Y-%m-%d")
  if(is.na(staterecessiontable$Troughs[length(staterecessiontable$Troughs)])) {
    staterecessiontable$Troughs[length(staterecessiontable$Troughs)]<- currentdate
  }
  staterecessiontable<<- staterecessiontable %>% 
    mutate(Months= abs((Peaks-Troughs)*12)) %>% 
    filter(Months>=5)
  
}





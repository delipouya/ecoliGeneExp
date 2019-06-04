##  in this script we will preprocess and clean the data

##  input: Gene-expression matrix + array-information matrix
##  output: list of splitted gene-expression matrix based on 
##      Aerobic/anaerobic/Anaerobic/Nitrate conditions


## loading data
rawData = read.csv("Data/GeneExp_data.csv", header = T)
arrayInfo = read.csv("Data/array_info.csv", header = T)


## data preprocessing

# row: gene , columns: conditions 
row.names(rawData) = rawData[,'X']
rawData = subset(rawData, select=-X)


## finding indexs based on Aerobic status 

aerobic_index = arrayInfo$Aerobic.Anaerobic %in% c("Aerobic", "aerobic")
anaerobic_index = arrayInfo$Aerobic.Anaerobic %in% c("anaerobic" ,"Anaerobic")
Anaerobic_Nitrate_index = arrayInfo$Aerobic.Anaerobic == "Anaerobic/Nitrate"


## splitting rawData based on aerobic/anaerobic/Anaerobic-Nitrate
rawData_aerobic = rawData[,aerobic_index]
rawData_Anaerobic = rawData[,anaerobic_index]
rawData_Anaerobic_Nitrate = rawData[,Anaerobic_Nitrate_index]

listOfGeneExp <- list(rawData_aerobic, rawData_Anaerobic, rawData_Anaerobic_Nitrate)
names(listOfGeneExp) <- c('rawData_aerobic', 'rawData_Anaerobic', 'rawData_Anaerobic_Nitrate')

lapply(listOfGeneExp, dim)
saveRDS(listOfGeneExp, 'Data/listOfGeneExp.rds')




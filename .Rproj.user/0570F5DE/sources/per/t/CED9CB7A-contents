# Read Lung exp data (Gene ID) excel file as matrix
LungExpData_ID <- as.matrix(read.table(file="GSE10072.txt", sep="", header=T))
# Read Lung exp data (Gene Symbol) excel file as matrix
LungExpData_GS <- as.matrix(read.table(file="PID10072.txt", sep="", header=T))

# Takes the data file object from the R environment and saves the R data object into data/
usethis::use_data(LungExpData_ID, overwrite=TRUE)
usethis::use_data(LungExpData_GS, overwrite=TRUE)


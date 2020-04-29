load('Data/Data.Rdata') # this line is to load the workspace
source('functions_demo.R') #  this line is to load the functions

#Model 1: cox model with random effects for patients
#######patient
#the data should be ordered based on the ID of patients before getting the C-index
###not run
#train_data = train_data[with(train_data, order(train_data$VisitLink)), ]
#test_data = test_data[with(test_data, order(test_data$VisitLink)), ]
#get the  C-index, running time can be long
###not run
#cindex_coxph_gamma_patient = gammaFrailty.up(train_data,test_data)
#Because the running time can be over 30 minutes, just skip running this line

#Show saved results
#C-index for train
cindex_coxph_gamma_patient[[1]]
#C-index for test
cindex_coxph_gamma_patient[[2]]
#Estimation of random effect
cindex_coxph_gamma_patient[[3]]
# It can be seen that there are over 4000 patients and each of them has a estimated frailty value

#Model 2:Parametric regression model with exponential baseline function, no random effect
#get the C-index
cindex_expreg_tt= expreg.u(train_data,test_data)
#Just skip running this line

#Show saved results
#C-index for train
cindex_expreg_tt[[1]]
#C-index for test
cindex_expreg_tt[[2]]

# For many different models compared in the project, each of them can get the C-index
# for train and test data, which indicate the prediction performance of the models.


################################################################
#The estimated frailty values can actually indicate the readmission risk of each patients
#As a example, the frailty values and their corresponding time to readmission can be plotted.

#Plot two in one column
par(mfrow = c(2,1))

#Plot for selected patients of their frailty values, indicating the random effects
p_patient = log(as.matrix(frailty_patient[c(22,237,536,1113,1266,1509,1524,2279,3803,3985)]))
barplot(p_patient[,1] ,xlab = "Patient indicators",ylab = "Frailties in patients"
                      ,cex.axis =1,cex=2,cex.names = 1,cex.lab = 1.2
                      ,ylim = c(-0.07,0.06)
                      ,col = c("red","white","white","green","red"
                      ,"blue","green","white","red","green"))
box(bty="l")


###Demo finished, Thank you!


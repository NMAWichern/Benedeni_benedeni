#R script for error propagation of clumped isotopes based temperature calibration

#Uses Monte Carlo simulations to propagate the non-linear errors associated with
#the nonlinear temperature calibration for clumped isotopes.

#Part of the supplementary data to Wichern et al., 
#"The potential of high-resolution stable isotope records in the bivalve 
#Angulus benedeni shells to investigate Pliocene seasonality".

#For further explanation, see sect. 3.6 the main text.

#Written by Niels J. de Winter



setwd("[INSERT DIRECTORY")                                          #set the directory from which the data is to be
                                                                    #loaded into R

#Define parameters
N = 10^5                                                            #Number of simulations
a_mean = 0.0397                                                     #Slope of T equation (Meinecke et al. 2021)
a_sd = 0.0004                                                       #SD on slope of T equation (Meinecke et al. 2021)
b_mean = 0.152                                                      #intercept of T equation (Meinecke et al. 2021)
b_sd = 0.004                                                        #SD on intercept of T equation (Meinecke et al. 2021)

#Read data
D_mean = read.table("D47_cold.txt", header=F)                       #read in D47 data
D_sd = read.table("SD_cold.txt", header=F)                          #read in D47 sd data

D_mean <- D_mean[[1]]                                               #formatting; may not be necessary depending on data format
D_sd <- D_sd[[1]]                                                   #formatting; may not be necessary depending on data format

#Monte Carlo simulation
a = rnorm(N, a_mean, a_sd)                                          #Simulate a values, type '?rnorm for more information
b = rnorm(N, b_mean, b_sd)                                          #Simulate b values


D = matrix(ncol = length(D_mean), nrow = N)                         #Create empty D47 matrix
Temp = matrix(ncol = length(D_mean), nrow = N)                      #Create empty temperature matrix
for(i in 1:length(D_mean)) {
    D[, i] <- rnorm(N, D_mean[i], D_sd[i])                          #simulate D47 values
    Temp[, i] = sqrt(a * 10 ^6 / (D[, i] - b))                      #simulate temperature values
}

#Calculate required parameters
T_mean = mean(Temp, na.rm = TRUE)                                   #mean temperature
T_sd = sd(Temp, na.rm = TRUE)                                       #sd on temperature
T_se = T_sd / sqrt(length(D_mean))                                  #standard error on temperature
T_95CL <- qt(1 - 0.05 / 2, length(D_mean)) * T_se                   #95% confidence level on temperature

T_mean_C <- T_mean-273.15                                           #transfer temperature to celcius


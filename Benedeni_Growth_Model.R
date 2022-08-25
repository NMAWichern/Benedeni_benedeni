#R script to estimate maximum shell height of a specimen of the bivalve Angulus
#benedeni benedeni using Von Bertalanffy and Gompertz growth functions.

#uses R base function nls(), which calculates the best fit of user-supplied data
#to an user-defined function based on an nonlinear least squares algorithm.

#Part of the supplementary data to Wichern et al. (202x), 
#"The potential of high-resolution stable isotope records in the bivalve 
#Angulus benedeni shells to investigate Pliocene seasonality".

#For further explanation, see sect. 3.7, 4.4.2 of the main text, 
#or use ?nls for more information on the nls() function.

#Written by Nina M.A. Wichern


#define data                               
age <- 1:9                                                  #define estimated age [yr]
ht <- c(4.02,9.01,11.81,14.9,18.23,21.89,23.9,26.8,28.16)   #define estimated shell height at each year [mm]
dat <- data.frame(age, ht)                                  #combine age and height into a dataframe

#Von Bertalanffy Growth Function (VBFG)
h.asymp_VBGF<-45                                            #estimate initial asymptotic height
k_VBGF<-0.09                                                #estimate initial growth coefficient
c_VBGF<-0                                                   #estimate initial t0  

nlc=nls.control(maxiter=1000)                               #set maximum number of iterations. 
VBGF_nls <-nls(ht~a*(1-exp(-b*(age-c))),                    #set VBGF as the function to fit the data to
               start=list(a=h.asymp_VBGF,                   #and run nonlinear least-squares algorithm
                          b=k_VBGF,
                          c=c_VBGF),control=nlc)

summary(VBGF_nls)                                           #summary of estimated parameters


#Gompertz Growth Function
h.asymp_G<-40                                               #estimate initial asymptotic height
b_G<-0.9                                                    #estimate initial Gompertz parameter b
c_G<-0.22                                                   #estimate initial Gompertz parameter c

nlc=nls.control(maxiter=1000)                               #set maximum number of iterations
G_nls <-nls(ht~a*exp(-b*exp(-c*age)),                       #set G as the function to fit the data to
            start=list(a=h.asymp_G,                         #and run nonlinear least-squares algorithm
                       b=b_G,
                       c=c_G),control=nlc)

summary(G_nls)                                              #summary of estimated parameters


#Plotting the growth functions
plot(age, ht, xlab= "Age [yr]", ylab= "Shell height [mm]",  #plot age-shell height datapoints
     xlim = c(0, 10), ylim = c(0, 30), pch=19)

lines(age,predict(VBGF_nls), col = "chocolate1", lwd = 2.5) #add VBGF growth curve to plot

lines(age,predict(G_nls), col = "royalblue", lwd = 2.5)     #add Gompertz growth curve to plot


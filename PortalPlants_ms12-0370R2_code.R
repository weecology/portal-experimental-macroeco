######## This is where the code for ms 12-0370R2 lives
######## "An experimental test of the response of macroecological patterns…”
######## Supp, Xiao, Ernest and White, Ecology, 2013, doi: 10.1890/12-0370.1
######## Code written by Sarah R. Supp and Xiao Xiao

# set working directory
wd = "pathname"
setwd(wd)

#import source code that contains the functions necessary to clean the code, generate the macroecological
#patterns, grab the parameters of interest, run the statistical analyses, and to generate the figures used
#in the manuscript.
source("PortalPlants_fxns.R")

# load winter and summer annual plant data and subset out bad data
winter = read.csv("PortalWinterAnnuals_1995_2009.csv")
winter = subset(winter, plot!=1 & plot!=9 & plot!=24)   #omit spectabs removal (1,9) and misshapen plot (24)

summer = read.csv("PortalSummerAnnuals_1995_2009.csv")
summer = subset(summer, plot!=1 & plot!=9 & plot!=24)   #omit spectabs removal (1,9) and misshapen plot (24)

# add spatial areas needed to get averages for SARs
winter_new = add_areas(winter)
summer_new = add_areas(summer)
  summer_new = subset(summer_new, year != 1997 & year != 1998) # yrs where Unid Spp > 10% community

########## MACROECOLOGICAL PATTERNS ##########
##### The below functions summarize the plant community data, construct the three macroecological patterns,
##### and get the parameters used to characterize each pattern as output for analysis

##### RADs
# reshape data for analysis
winter_wide = reshape_data(winter_new)
summer_wide = reshape_data(summer_new)

#get descriptive parameters for SADs for each plot-year combination
years = c(1995:2009)
winter_rad = RAD_data(winter_wide, years) #returns all the parameters
summer_rad = RAD_data(summer_wide, years) 


##### SARs
# aggregate data at spatial scales below the plot
areas = c(0.25, 0.5, 1, 2, 4)
years = c(1995 : 2009)
winter_richness = Create_richness_dataframe(winter_new, areas, years)
summer_richness = Create_richness_dataframe(summer_new, areas, years)

SAR_winter_means = Create_means(winter_richness)
SAR_summer_means = Create_means(summer_richness)

# get SAR parameters for both seasonal communities, which will be used in the statistical analyses below
# get slope and intercept for winter annuals, using means at each spatial scale
plots = unique(sort(SAR_winter_means[,3]))
years = unique(sort(SAR_winter_means[,2]))
winter_sar = SAR_slope(SAR_winter_means, plots, years)

# get slope and intercept for summer annuals, using means at each spatial scale
plots = unique(sort(SAR_summer_means[,3]))
years = unique(sort(SAR_summer_means[,2]))
summer_sar = SAR_slope(SAR_summer_means, plots, years)


##### STRs
years = c(1995 : 2009)   # Winter Annuals
timespans = c(1 : 15)
winter_time = Create_time_dataframe(winter_new, timespans, years)
years = c(1999 : 2009)   #summer annuals -- shorter time series because 1997-98 data has >10% unidentified individuals
timespans = c(1 : 11)
summer_time = Create_time_dataframe(summer_new, timespans, years)

STR_winter_means = Create_time_means(winter_time)
STR_summer_means = Create_time_means(summer_time)

# get slopes and intercepts for STRs using means at each timeframe
plots = unique(STR_winter_means[,2])
winter_str = time_slope(STR_winter_means, plots)

plots = unique(STR_summer_means[,2])
summer_str = time_slope(STR_summer_means, plots)

########## STATISTICAL ANALYSIS ##########
#For this paper, we did both standard statistical analyses for the parameters (species composition,
# S, N, SAD mu, SAD sigma, SAR intercept, SAR slope, STR slope) to look for statistically significant
# differences and equivalence testing (SAD mu, SAD sigma, SAR intercept, SAR slope, STR slope) to look
# for statistically signficant similarities between the macroecological patterns in the experimental 
#treatments. 

##### Species composition - pCCA analysis
#Compositional differences among rodent treatments were characterized with partially constrained 
#correspondence analysis (pCCA; Oksanen et al. 2010) and permutational significance tests were used to 
#determine significance of the pCCA axes. We square root transformed the abundance data and controlled 
#for the effect of year.

#preliminary setup of variables    
# Winter
names(winter_wide) <- c(names(winter_wide)[1:3],sapply(names(winter_wide)[-(1:3)],function(x) strsplit(x, "abun.")[[1]][2]))
win.sp <- sqrt(as.matrix(winter_wide[,-(1:3)])) #sqrt transform to smush super-abun (Erod texa)
win.trmt <- as.factor(winter_wide$trmt)
win.yr <- as.factor(winter_wide$year)
win.plot <- as.factor(winter_wide$plot)  
win.abu <- colSums(win.sp)

# Summer annuals
names(summer_wide) <- c(names(summer_wide)[1:3],sapply(names(summer_wide)[-(1:3)],function(x) strsplit(x, "abun.")[[1]][2]))
sum.sp <- sqrt(as.matrix(summer_wide[,-(1:3)])) #sqrt transform to smush super-abun (Erod texa)
sum.trmt <- as.factor(summer_wide$trmt)
sum.yr <- as.factor(summer_wide$year)
sum.plot <- as.factor(summer_wide$plot)   
sum.abu <- colSums(sum.sp)

#partially constrained correspondence analysis
win.pcca <- cca(win.sp ~ win.trmt + Condition(win.yr))   #winter annuals
sum.pcca <- cca(sum.sp ~ sum.trmt + Condition(sum.yr))   #summer annuals

### SIGNIFICANCE TEST ON pCCA AXES
# winter
anova(win.pcca)
permutest(win.pcca,permutations=500) # should be similar to anova on pcca
anova(win.pcca,strata=win.yr) # more conservative test           

# summer           
anova(sum.pcca)
permutest(sum.pcca,permutations=500)
anova(sum.pcca,strata=sum.yr) # more conservative test

### ADONIS - another test to look for compositional differences, with similar results to above.
# winter
win.spp.canb = vegdist(win.sp, method = "canb")
win.canb = adonis(win.spp.canb ~ win.trmt, permutation=1000)

# summer
sum.spp.canb = vegdist(sum.sp, method = "canb")
sum.canb = adonis(sum.spp.canb ~ sum.trmt, permutations = 1000)



##### STATE VARIABLES: species richness (S) and total abundance (N)
# get just the S and N columns from the whole-plot scale (all 16 quadrats)
winter_sn = winter_richness[winter_richness$area == 4, c(1, 2, 4, 5)]
summer_sn = summer_richness[summer_richness$area == 4, c(1, 2, 4, 5)]

# testing for significant differences
# linear mixed effects models - generates values in Appendix D, Tables S2 and S3
summer_s=lmer_analysis(summer_sn$S, summer_sn)
## E-C: .1653   R-C: .7532   E-R: .2675
summer_n=lmer_analysis((summer_sn$N)^(1/3), summer_sn)
## E-C: .3095  R-C: .734   E-R: .4745
winter_s=lmer_analysis((winter_sn$S)^0.75, winter_sn)
## E-C: .0012  R-C: .981   E-R: .0013
winter_n=lmer_analysis((winter_sn$N)^(1/3), winter_sn)
## E-C: .0138  R-C: .2541    E-R: .152
FDR_sig(as.numeric(c(summer_s, summer_n))) #0
FDR_sig(as.numeric(c(winter_s, winter_n))) #3

# testing for significant equivalence
# linear mixed effects models - values not included in the paper. 
# Equivalence testing explained in Appendix B -- some values may slightly differ from published because it is an approximation
summer_s=lmer_equi(summer_sn$S, summer_sn, bound_level = 0.1)
## E-C: 0.5636   R-C: 0.9221   E-R: 0.6898
summer_n=lmer_equi(summer_sn$N, summer_sn, pw = 1/3, bound_level = 0.1)
## E-C: 0.3742  R-C: 0.6484   E-R: 0.5031
winter_s=lmer_equi(winter_sn$S, winter_sn, pw = 0.75, bound_level = 0.1)
## E-C: 0.0681  R-C: 0.9699   E-R: 0.0724
winter_n=lmer_equi(winter_sn$N, winter_sn, pw = 1/3, bound_level = 0.1)
## E-C: 0.0182  R-C: 0.2811    E-R: 0.1821
FDR_equi(as.numeric(c(summer_s, summer_n))) #0
FDR_equi(as.numeric(c(winter_s, winter_n))) #0



##### SADs
#The maximum likelihood (MLE) of the Poisson log-normal parameters, µ (mean) and s (standard deviation), 
#were estimated with R function “poilogMLE” from package “poilog” (Grøtan & Engen 2008). 
#Since µ took both positive and negative values, we used its exponentiated form, exp (µ).
# linear mixed effects models

# mu and sigma
# testing for significant differences
# linear mixed effects models - generates values in Appendix D, Tables S2 and S3
summer_rad_mu = lmer_analysis((exp(summer_rad$mu))^0.5, summer_rad)
## E-C: .2006  R-C: .7243  E-R: .1091
summer_rad_sig = lmer_analysis(summer_rad$sigma^0.5, summer_rad)
## E-C: .4862  R-C: .6853  E-R: .7298
winter_rad_mu =  lmer_analysis((exp(winter_rad$mu))^0.3, winter_rad)
## E-C: .1645  R-C: .7826  E-R: .2475
winter_rad_sig =  lmer_analysis(winter_rad$sigma^0.1, winter_rad)
## E-C: 2e-4  R-C: .1804  E-R: .0094

#testing for significant equivalence
# linear mixed effects models - generates values in Appendix D, Tables S4 and S5
#  -- some values may slightly differ from published because it is an approximation
summer_rad_mu_equi = lmer_equi(exp(summer_rad$mu), summer_rad, pw = 0.5)
## E-C: 0.1356  R-C: 0.2751  E-R: 0.0943
summer_rad_sig_equi = lmer_equi(summer_rad$sigma, summer_rad, pw = 0.5)
## E-C: 0.7149  R-C: 0.8536  E-R: 0.8068
winter_rad_mu_equi =  lmer_equi(exp(winter_rad$mu), winter_rad, pw = 0.3)
## E-C: 0.0661  R-C: 0.1771  E-R: 0.0861
winter_rad_sig_equi =  lmer_equi(winter_rad$sigma, winter_rad, pw = 0.1)
## E-C: 0.021  R-C: 0.7587  E-R: 0.2483



##### SARs
# intercept and slope
# testing for significant differences
# linear mixed effects models - generates values in Appendix D, Tables S2 and S3
summer_sar_slope = lmer_analysis(summer_sar$slope, summer_sar)
## E-C: .4385  R-C: .9324  E-R: .3949
summer_sar_intercept = lmer_analysis(summer_sar$intercept^2.5, summer_sar)
## E-C: .3374  R-C: .4102  E-R: .811
winter_sar_slope =  lmer_analysis(winter_sar$slope, winter_sar)
## E-C: .9083  R-C: .437  E-R: .3999
winter_sar_intercept =  lmer_analysis(winter_sar$intercept^2, winter_sar) # not included in the tables
## E-C: 2e-4  R-C: .4094  E-R: .0024

#testing for significant equivalence
# linear mixed effects models - generates values in Appendix D, Tables S4 and S5
#  -- some values may slightly differ from published because it is an approximation
summer_sar_slope_equi = lmer_equi(summer_sar$slope, summer_sar)
## E-C: .4641  R-C: .6905  E-R: .4396
summer_sar_intercept_equi = lmer_equi(summer_sar$intercept, summer_sar, pw = 2.5)
## E-C: 0.7822  R-C: 0.879  E-R: 0.9358
winter_sar_slope_equi =  lmer_equi(winter_sar$slope, winter_sar)
## E-C: 0.7607  R-C: 0.6176  E-R: 0.546
winter_sar_intercept_equi =  lmer_equi(winter_sar$intercept, winter_sar, pw = 2) #not included in the tables
## E-C: 0.0142  R-C: 0.8953  E-R: 0.088


##### STRs
# slope
# ANOVA
# testing for significant differences
summer_str_slope = str_anova(summer_str$slope, summer_str)
## E-C: 0.5929169  R-C: 0.600602  E-R: 0.3250864
winter_str_slope = str_anova(winter_str$slope, winter_str)
## E-C: 0.02944050  R-C:0.686883  E-R: 0.06011208


# testing for significant equivalence
# ANOVA -- generates values in Appendix D, tables S
# -- some values may slightly differ from published because it is an approximation
summer_str_slope_equi = str_equi(summer_str$slope, summer_str)
## E-C: 0.6065354  R-C: 0.672158  E-R: 0.4690185
winter_str_slope_equi = str_equi(winter_str$slope, winter_str)
## E-C: 0.1302223  R-C: 0.7819378  E-R: 0.2154747


# We used FALSE DISCOVERY RATE CONTROL (Benjamini & Hochberg 1995; Garcia 2004) to correct for 
# multiple statistical tests within each seasonal community. 
# False Discovery Rate Control on the conventional and equivalence statistical tests
FDR_sig(as.numeric(c(summer_sar_slope, summer_sar_intercept, summer_rad_mu, summer_rad_sig,
                     summer_str_slope)))  #0
FDR_sig(as.numeric(c(winter_sar_slope, winter_sar_intercept, winter_rad_mu, winter_rad_sig,
                     winter_str_slope)))  #4

FDR_equi(as.numeric(c(summer_sar_slope_equi, summer_sar_intercept_equi, summer_rad_mu_equi,
                      summer_rad_sig_equi, summer_str_slope_equi)))  #0
FDR_equi(as.numeric(c(winter_sar_slope_equi, winter_sar_intercept_equi, winter_rad_mu_equi,
                      winter_rad_sig_equi, winter_str_slope_equi)))  #0


########## FIGURE 1: Possible responses of three macroecological patterns to manipulated seed predation ##########
#assuming that the manipulation has no effect on species richness (S) and total abundance (N). Please 
#note that each macroecological pattern varies with manipulations that impact species composition 
#(blue dotted line) despite fixed S and N.

# set up plotting window for FIGURE 1
pdf("Fig1_conceptual.pdf", 7, 3, pointsize = 10)
par(mfcol=c(1,3), oma = c(2,3.75,2,1), mar=c(1,1,1,1))

# plot SAD
abun1 = c(5000, 500, 50, 10, 7, 6, 6, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1)
par1 = as.list(poilogMLE(abun1, startVals = c(mean(log(abun1)), sd(log(abun1))))$par)
abun2 = c(1000, 848, 750, 650, 500, 300, 275, 250, 250, 250, 150, 100, 90, 80, 40, 30, 20, 10, 10, 5, 1, 1 )
par2 = as.list(poilogMLE(abun2, startVals = c(mean(log(abun2)), sd(log(abun2))))$par)
#plot relationships
plot(c(1:length(abun1)), abun1, log = 'y', type = 'l', lwd = 4, xlim = c(1,22), lty = 3, ylim = c(1,5000),
     col = 'dodgerblue2', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
lines(c(1:length(abun2)), abun2,  type = 'l', lwd = 2, col = 'black')
mtext('log(N)', side = 2, line = .5)
mtext('Rank', side = 1, line = .25)
mtext('SAD', side = 3, line = .5)

# plot SAR 
areas = c(1, 2, 4, 8, 16, 32)
S1 = c(8, 10, 20, 35, 40, 50)
S2 = c(2, 5, 10, 15, 30, 50)
reg1 = lm(log10(S1)~log10(areas))
reg2 = lm(log10(S2)~log10(areas))
plot(NA, NA, xlab = '', ylab = '', log = 'xy', xlim = c(1,32), ylim = c(1,55), xaxt = 'n', yaxt = 'n')
par(new=TRUE)
abline(reg1, lwd = 4, col = 'dodgerblue2', lty = 3)
abline(reg2, lwd = 2, col = 'black')
mtext('log(S)', side = 2, line = .5)
mtext('log(Area)', side = 1, line = .25)
mtext('SAR', side = 3, line = .5)

# plot STR 
times = c(1:15)
S4 = c(7, 10, 16, 21, 25, 27, 30, 35, 37, 42, 49, 53, 59, 66, 71)
S5 = c(7, 9, 11, 15, 17, 20, 21, 22, 23, 24, 26, 27, 30, 31, 33)
reg4 = lm(log10(S4)~log10(times))
reg5 = lm(log10(S5)~log10(times))
plot(NA, NA, xlab = '', ylab = '', log = 'xy', xlim = c(1,15), ylim = c(1,75), xaxt = 'n', yaxt = 'n')
par(new=TRUE)
abline(reg4, lwd = 4, col = 'dodgerblue2', lty = 3)
abline(reg5, lwd = 2, col = 'black')
mtext('log(S)', side = 2, line = .5)
mtext('log(Time)', side = 1, line = .25)
mtext('STR', side = 3, line = .5)
legend('bottomright', c('experiment','     ' ,'control'), bty = 'n', lty = c(1,1, 3), lwd =c(2,1,4), seg.len = 3, col = c('black', 'white', 'dodgerblue2'))

dev.off()



########## FIGURE 2: Statistical differences among the parameters were only detected in the winter annual ##########
# community when experimental manipulation (C = control, K = kangaroo rat removal, R = total rodent removal)
# also impacted species richness and total abundance. Top panels display results from standard statistical
# tests (linear mixed effects models - SAD, SAR; ANOVA – STR) for significant differences and lower panels
# display results from equivalence tests. Points represent the mean difference in parameter estimation 
# between two treatments, and whiskers indicate 95% confidence intervals (CIs; top) and 90% CI (bottom) of
# the difference in parameter estimates. Because parameter estimates differ in magnitude for different 
# patterns, all values and their CIs are standardized with respect to their designated range of 
# equivalence in both the upper and lower panels for better visualization.

summer_point = as.data.frame(matrix(nrow = 7, ncol = 9))
summer_equi = as.data.frame(matrix(nrow = 7, ncol = 9))
winter_point = as.data.frame(matrix(nrow = 7, ncol = 9))
winter_equi = as.data.frame(matrix(nrow = 7, ncol = 9))

################ SAR ######################
summer_sar_slope = lmer_quan(summer_sar$slope, summer_sar)
summer_point[1, ] = as.numeric(unlist(summer_sar_slope[c(1, 3, 5)]))
summer_equi[1, ] = as.numeric(unlist(summer_sar_slope[c(2, 4, 6)]))

summer_sar_intercept = lmer_quan(summer_sar$intercept, summer_sar, pw = 2.5)
summer_point[2, ] = as.numeric(unlist(summer_sar_intercept[c(1, 3, 5)]))
summer_equi[2, ] = as.numeric(unlist(summer_sar_intercept[c(2, 4, 6)]))

winter_sar_slope =  lmer_quan(winter_sar$slope, winter_sar)
winter_point[1, ] = as.numeric(unlist(winter_sar_slope[c(1, 3, 5)]))
winter_equi[1, ] = as.numeric(unlist(winter_sar_slope[c(2, 4, 6)]))

winter_sar_intercept =  lmer_quan(winter_sar$intercept, winter_sar, pw = 2)
winter_point[2, ] = as.numeric(unlist(winter_sar_intercept[c(1, 3, 5)]))
winter_equi[2, ] = as.numeric(unlist(winter_sar_intercept[c(2, 4, 6)]))

################## RAD #####################
summer_rad_mu = lmer_quan(exp(summer_rad$mu), summer_rad, pw = 0.5)
summer_point[3, ] = as.numeric(unlist(summer_rad_mu[c(1, 3, 5)]))
summer_equi[3, ] = as.numeric(unlist(summer_rad_mu[c(2, 4, 6)]))

summer_rad_sig = lmer_quan(summer_rad$sig, summer_rad, pw = 0.5)
summer_point[4, ] = as.numeric(unlist(summer_rad_sig[c(1, 3, 5)]))
summer_equi[4, ] = as.numeric(unlist(summer_rad_sig[c(2, 4, 6)]))

winter_rad_mu = lmer_quan(exp(winter_rad$mu), winter_rad, pw = 0.3)
winter_point[3, ] = as.numeric(unlist(winter_rad_mu[c(1, 3, 5)]))
winter_equi[3, ] = as.numeric(unlist(winter_rad_mu[c(2, 4, 6)]))

winter_rad_sig = lmer_quan(winter_rad$sig, winter_rad, pw = 0.1)
winter_point[4, ] = as.numeric(unlist(winter_rad_sig[c(1, 3, 5)]))
winter_equi[4, ] = as.numeric(unlist(winter_rad_sig[c(2, 4, 6)]))

################### STR ######################
summer_str_slope = str_quan(summer_str$slope, summer_str)
summer_point[5, ] = as.numeric(unlist(summer_str_slope[c(1, 3, 5)]))
summer_equi[5, ] = as.numeric(unlist(summer_str_slope[c(2, 4, 6)]))

winter_str_slope = str_quan(winter_str$slope, winter_str)
winter_point[5, ] = as.numeric(unlist(winter_str_slope[c(1, 3, 5)]))
winter_equi[5, ] = as.numeric(unlist(winter_str_slope[c(2, 4, 6)]))

#################### SN ##########################
summer_s = lmer_quan(summer_sn$S, summer_sn, bound_level = 0.1)
summer_point[6, ] = as.numeric(unlist(summer_s[c(1, 3, 5)]))
summer_equi[6, ] = as.numeric(unlist(summer_s[c(2, 4, 6)]))

summer_n = lmer_quan(summer_sn$N, summer_sn, pw = 1/3, bound_level = 0.1)
summer_point[7, ] = as.numeric(unlist(summer_n[c(1, 3, 5)]))
summer_equi[7, ] = as.numeric(unlist(summer_n[c(2, 4, 6)]))

winter_s = lmer_quan(winter_sn$S, winter_sn, pw = 0.75, bound_level = 0.1)
winter_point[6, ] = as.numeric(unlist(winter_s[c(1, 3, 5)]))
winter_equi[6, ] = as.numeric(unlist(winter_s[c(2, 4, 6)]))

winter_n = lmer_quan(winter_sn$N, winter_sn, pw = 1/3, bound_level = 0.1)
winter_point[7, ] = as.numeric(unlist(winter_n[c(1, 3, 5)]))
winter_equi[7, ] = as.numeric(unlist(winter_n[c(2, 4, 6)]))

# make summary figure of results
pdf("Fig2_summary.pdf", 7, 5, paper = "letter", pointsize = 10)
XMIN = floor(min(summer_point, summer_equi, winter_point, winter_equi))
XMAX = ceiling(max(summer_point, summer_equi, winter_point, winter_equi))
par(mfrow = c(2, 2), oma = c(4, 3.2, 2, 1), mar = c(0, 0, 0, 0))
sub_plot(summer_point, "point", yax = T)
symbol_list = c(15, 16, 21:25)
legend("bottomleft", legend = c("S", "N", c(expression(paste("SAD exp(",mu,")"))), 
                                c(expression(paste("SAD ", sigma))), "SAR slope", "SAR intercept", "STR slope"), 
       pch = symbol_list, pt.bg = "white", lty = "solid", lwd = c(3, 3, rep(1, 5)), cex = 0.9)
sub_plot(winter_point, "point")
sub_plot(summer_equi, "equi", yax = T, xax = T)
sub_plot(winter_equi, "equi", xax = T)

mtext("Summer                               Winter", outer = TRUE, cex = 1.5, line = 0.4)
dev.off()



########## APPENDIX B: PLOT ALL THE PATTERNS AND THE FUNCTIONS USED TO CHARACTERIZE THEM ##########
# The code below should print 6 separate pdf files for each of the patterns in each of the seasonal communities.
#Each pattern is labeled with the experimental plot identification number and year combination. 
#Black points represent the plotted data and the red lines represent the function used to fit the data. 
#SADs were characterized using the Poisson log-normal distribution and we plot the data as rank abundance 
#distributions (RADs) for visual ease. The x-axis is rank and the y-axis is abundance. SADs and STRs were 
#characterized using power-laws. For SARs, the x-axis is the area sampled (0.25, 0.5, 1, 2, 4) in square 
#meters, and the y-axis is mean abundance at each spatial scale. For STRs, the x-axis is the timespan 
#sampled in years (winter, 1-15; summer, 1-11), and the y-axis is mean abundance for each timespan. 
#Experimental plot identification numbers refer to experimental treatment as follows: 
#Controls (2, 4, 8, 11, 12, 14, 17, 22), Kangaroo rat removals (3, 6, 13, 15, 18, 19, 20, 21), 
#and total rodent removals (5, 7, 10, 16, 23). 

### Winter RADs
pdf("WinterRADs.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(7,5), mar=c(1,2,3,1), oma=c(1,0,0,0))

years = c(1995:2009)
RAD_plot(winter_wide, years) #move the function b.c it takes a long time to run

dev.off()

### Summer RADs
pdf("SummerRADs.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(7,5), mar=c(1,2,3,1), oma=c(1,0,0,0))
years = c(1995:2009)
RAD_plot(summer_wide, years)

dev.off()

### Winter SARs
pdf("WinterSARs.pdf",7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(7,5), mar=c(1,2,3,1), oma=c(1,0,0,0))

years = unique(sort(SAR_winter_means$year))
plots = unique(sort(SAR_winter_means$plot))
plot_SARs(SAR_winter_means, years, plots)

dev.off()


### summer SARs
pdf("SummerSARs.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(7,5), mar=c(1,2,3,1), oma=c(1,0,0,0))
years = unique(sort(SAR_summer_means$year))
plots = unique(sort(SAR_summer_means$plot))
plot_SARs(SAR_summer_means, years, plots)

dev.off()


### Winter STRs
pdf("WinterSTRs.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(7,5), mar=c(1,2,3,1), oma=c(1,0,0,0))
plots = unique(STR_winter_means[,2])
plot_STRs(STR_winter_means, plots)

dev.off()

### Summer STRs
pdf("SummerSTRs.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(7,5), mar=c(1,2,3,1), oma=c(1,0,0,0))
plots = unique(STR_summer_means[,2])
plot_STRs(STR_summer_means, plots)

dev.off()



########## APPENDIX C: DETAILS ON THE METHODS AND RESULTS OF EQUIVALENCE TESTING ##########
#Figure C1. Visual depiction of equivalence test ranges. We deemed ranges within patterns equivalent 
# +/- 5% for all response variables. These ranges translate to roughly 25% deviation in the abundance 
# of the most abundant species for SADs (left) and a 20% deviation in species richness for SARs and STRs 
# at all scales (middle, right), which we felt represented reasonable fluctuation for claiming equivalence.
pdf("FigC1_equivalence.pdf", 7, 3, pointsize = 10)
par(mfrow=c(1,3), mar=c(1,2,3,1), oma=c(1,0,0,0))

YEAR = 1999

## RAD
summer_rad_year = summer_rad[summer_rad$year == YEAR, ]
summer_rad_trmt = add_trmt(summer_rad_year)
mu_mean = mean(c(mean(summer_rad_trmt$mu[summer_rad_trmt$trmt == "C"]), 
                 mean(summer_rad_trmt$mu[summer_rad_trmt$trmt == "E"]), 
                 mean(summer_rad_trmt$mu[summer_rad_trmt$trmt == "R"])))
sig_mean = mean(c(mean(summer_rad_trmt$sig[summer_rad_trmt$trmt == "C"]), 
                  mean(summer_rad_trmt$sig[summer_rad_trmt$trmt == "E"]), 
                  mean(summer_rad_trmt$sig[summer_rad_trmt$trmt == "R"])))  
S = round(mean(summer_rad_trmt$S))
N = round(mean(summer_rad_trmt$N))
rad_plot_eq(N, S, mu_mean, sig_mean)

## SAR
summer_sar_year = summer_sar[summer_sar$year == YEAR, ]
summer_sar_trmt = add_trmt(summer_sar_year)
sar_slope = mean(c(mean(summer_sar_trmt$slope[summer_sar_trmt$trmt == "C"]), 
                   mean(summer_sar_trmt$slope[summer_sar_trmt$trmt == "E"]), 
                   mean(summer_sar_trmt$slope[summer_sar_trmt$trmt == "R"]))) 
sar_intercept = mean(c(mean(summer_sar_trmt$intercept[summer_sar_trmt$trmt == "C"]), 
                       mean(summer_sar_trmt$intercept[summer_sar_trmt$trmt == "E"]), 
                       mean(summer_sar_trmt$intercept[summer_sar_trmt$trmt == "R"]))) 
sar_str_plot(sar_slope, sar_intercept, "SAR")

## STR
summer_str_trmt = add_trmt(summer_str)
str_slope = mean(c(mean(summer_str_trmt$slope[summer_str_trmt$trmt == "C"]), 
                   mean(summer_str_trmt$slope[summer_str_trmt$trmt == "E"]), 
                   mean(summer_str_trmt$slope[summer_str_trmt$trmt == "R"]))) 
str_intercept = mean(c(mean(summer_str_trmt$intercept[summer_str_trmt$trmt == "C"]), 
                       mean(summer_str_trmt$intercept[summer_str_trmt$trmt == "E"]), 
                       mean(summer_str_trmt$intercept[summer_str_trmt$trmt == "R"]))) 
sar_str_plot(str_slope, str_intercept, "STR")

dev.off()

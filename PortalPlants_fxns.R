# Portal data fxns for use with ms 12-0370R2 in Ecology; Supp, Xiao, Ernest, and White
# Code was written collaboratively by Sarah Supp and Xiao Xiao

#import libraries needed for functions
library(BiodiversityR)
library(car)
library(CCA)
library(equivalence)
library(gplots)
library(graphics)
library(languageR)
library(lme4)
library(nlme)
library(plotrix)
library(poilog)
library(vegan)
library(VGAM)

########## FUNCTIONS FOR CLEANING UP THE DATA ######################################################################

add_areas=function(dat){
  # This code adds areas to the Portal dataset to enable later calculation of community structure at different spatial scales
  # add "two", neighboring pairs of quadrats
  for (i in 1:nrow(dat)){
    if(dat[i,3]==11 | dat[i,3]==13) {
      dat$two[i]=1}
    else if(dat[i,3]==31 | dat[i,3]==33){
      dat$two[i]=2 }
    else if(dat[i,3]==51 | dat[i,3]==53){
      dat$two[i]=3 }
    else if(dat[i,3]==71 | dat[i,3]==73){
      dat$two[i]=4 }
    else if(dat[i,3]==15 | dat[i,3]==17){
      dat$two[i]=5 }
    else if(dat[i,3]==35 | dat[i,3]==37){
      dat$two[i]=6 }
    else if(dat[i,3]==55 | dat[i,3]==57){
      dat$two[i]=7 }
    else {dat$two[i]=8} }
  # add corner
  for (i in 1:nrow(dat)) {
    if (dat[i,3]==11 | dat[i,3]==13 | dat[i,3]==31 | dat[i,3]==33 ){
      dat$corner[i] = 1}
    else if (dat[i,3]==15 | dat[i,3]==17 | dat[i,3]==35 | dat[i,3]==37 ){
      dat$corner[i] = 2}
    else if (dat[i,3]==51 | dat[i,3]==53 | dat[i,3]==71 | dat[i,3]==73 ){
      dat$corner[i] = 3}
    else { dat$corner[i] = 4}
  }
  # add bisection
  for (i in 1:nrow(dat)){
    if(dat[i,3]==11 | dat[i,3]==13 | dat[i,3]==15 | dat[i,3]==17 |
      dat[i,3]==31 | dat[i,3]==33 | dat[i,3]==35 | dat[i,3]==37) {
      dat$bisect[i]=1}
    else {dat$bisect[i]=2} }
  # add treatment
  for (i in 1:nrow(dat)){
    if(dat[i,2]==2 | dat[i,2]==4 | dat[i,2]==8 | dat[i,2]==11 |
      dat[i,2]==12 | dat[i,2]==14 | dat[i,2]==17 | dat[i,2]==22) {
      dat$trmt[i] = 1}
    else if (dat[i,2]==3 | dat[i,2]==6 | dat[i,2]==13 | dat[i,2]==15 |
      dat[i,2]==18 | dat[i,2]==19 | dat[i,2]==20 | dat[i,2]==21) {
      dat$trmt[i] = 2}
    else {dat$trmt[i] = 3} }
  return(dat)
}


reshape_data = function(dat){
  # puts data in "wide" format aggregated BY PLOT for use with certain functions, esp. related to compostional analyses
  # and the rank-abundance distribution
  dat2 = aggregate(dat$abundance,by=list(year=dat$year, plot=dat$plot, trmt=dat$trmt,species=dat$species),FUN=sum)
  #names the new column
  names(dat2)[5]="abun"
  dat3 = reshape(dat2,idvar=c("year","plot","trmt"),timevar="species",direction="wide")
  dat3[is.na(dat3)] = 0
  return(dat3)
}


add_trmt=function(dat_old){
  # Adds treatment codes to the data (C = control, R = kangaroo rat removal, E = total rodent exclosures)
  dat=dat_old
  for (i in 1:dim(dat)[1]){
    if (dat$plot[i] %in% c(2,4,8,11,12,14,17,22)){
      dat$trmt[i]="C"}
    else if (dat$plot[i] %in% c(5,7,10,16,23)){
      dat$trmt[i]="E"}
    else {dat$trmt[i]="R"}
  }
  dat$trmt=as.factor(dat$trmt)
  dat$plot=as.factor(dat$plot)
  return(dat)
}



########## FUNCTIONS FOR MACROECOLOGICAL PATTERNS AND PARAMETERS #########################################################################
##### RADs #####

poilog_cdf = function(mu, sigma, x){
  # Function to compute the cdf of poilog by adding up the pmfs. XX.
  x_list = seq(x)
  cdf = sum(dpoilog(x_list, mu, sigma) / (1 - dpoilog(0, mu, sigma))) # Left-truncated
  return(cdf)
}


rad_poilog = function(x_vec, plot, year){
  #to get back parameter values mu and sigma for RADs
  par_mle = poilogMLE(x_vec, startVals = c(mu = mean(log(x_vec)),
                                           sig = sd(log(x_vec))))$par
  mu = as.numeric(par_mle)[1]
  sig = as.numeric(par_mle)[2]
  S = length(x_vec)
  N = sum(x_vec)
  
  parameters = c(year, plot, S, N, mu, sig)
  return(parameters)
}


RAD_data=function(dat,year) {
  ## Records data on the mu and sigma values for each plot's species abundance distribution in each year
  #set up dataframe to store results
  abundance = data.frame("year"=1, "plot"=1, "S"=1, "N"=1, "mu"=1, "sigma"=1)
  outcount = 0
  plot=unique(sort(dat[,2]))   #generate a list of all the plots
  for (iYr in 1:length(year)) {   #loop thru years
    for (iPlot in 1:length(plot)){   #loop thru plots
      tmp <- which(dat[,2]==plot[iPlot] & dat[,1]==year[iYr])
      if (nrow(dat[tmp,]) > 0){   #make sure data exists for this plot-year combo
        abuns = as.numeric(dat[tmp,c(4:ncol(dat))])
        abuns = sort(abuns[abuns > 0], decreasing = T)
        if (length(abuns) >= 5) {  #Don't plot years when S < 5
          # Get Rank Abundance parameters using species data and poilog prediction
          param = rad_poilog(abuns, plot[iPlot], year[iYr])
          #record parameter values
          outcount = outcount + 1
          abundance[outcount,] <- param
        }}}}
  return (abundance)
}


##### SARs #####

Create_richness_dataframe = function(dat, areas, years){
  # Generates a dataframe with richness at 5 scales below the scale of a whole Portal plot.
  wa_richness = data.frame("year"=1, "plot"=1, "area"=1, "S"=1, "N"=1)
  for (A in 1:length(areas)){
    if (A == 1){
      dat$myLabels = paste(dat[,9], dat[,2], dat[,8], dat[,7], dat[,6], dat[,3], sep=",")}
    else if (A == 2){
      dat$myLabels = paste(dat[,9], dat[,2], dat[,8], dat[,7], dat[,6], sep=",")}
    else if (A == 3){
      dat$myLabels = paste(dat[,9], dat[,2], dat[,8], dat[,7], sep=",") }
    else if (A == 4){
      dat$myLabels = paste(dat[,9], dat[,2], dat[,8], sep=",")}
    else{
      dat$myLabels = paste(dat[,9], dat[,2], sep=",")}
    
    plots = unique(dat$myLabels)
    dat_r = richness_at_areas(dat, areas[A], plots, years)
    wa_richness = rbind(wa_richness, dat_r)
  }
  return(wa_richness[-1,])
}


richness_at_areas = function(dat, area, plots, years){
  # sums species richness (S) and total abundance (N) at a specific plot, year, area within that plot
  richness = data.frame("year"=1, "plot"=1, "area"=1, "S"=1, "N"=1)
  outCount = 0
  for (iPlot in 1:length(plots)){
    for(iYr in 1:length(years)){
      outCount<-outCount+1
      tmp<-which(dat$myLabels==plots[iPlot] & dat[,1]==years[iYr])
      if (nrow(dat[tmp,]) == 0) {
        S = 0
        N = 0
        plot = unique(dat[which(dat$myLabels==plots[iPlot]),2])}
      else {
        S<-length(unique(dat$species[tmp]))   # count species and individuals
        N<-sum(dat$abundance[tmp],na.rm=T)
        plot = unique(dat$plot[tmp])}
      richness[outCount,]<-c(years[iYr], plot, area, S, N)
    }}
  return(richness)}


Create_means=function(dat){
  # Get mean richness at each sub-plot spatial scale for the species-area relationship (SAR)
  #aggregate data by each area within a given plot-year combination
  means=aggregate(dat$S,by=list(area=dat$area,year=dat$year,plot=dat$plot),FUN=mean)
  names(means)[4]="S"
  N=aggregate(dat$N,by=list(area=dat$area,year=dat$year,plot=dat$plot),FUN=mean)
  means[5]=N[,4]
  names(means)[5]="N"
  return(means)
}


SAR_slope=function(dat, plot, year){
  # Get SAR slope, z, for each plot-year combination using power law
  # dataframe to store results
  zValue = data.frame("plot" = 1, "year" = 1, "intercept" = 1, "slope" = 1, "R2adj" = 1)
  outcount = 0
  area = unique(sort(dat[,1]))
  # generate slopes for each plot-year combo, record in dataframe
  for(i in 1:length(year)){ # loop thru years
    for (j in 1:length(plot)){  # loop thru plots
      tmp<-which(dat[,3]==plot[j] & dat[,2]==year[i])
      if (max(dat[tmp,4]) > 5){
        A = log10(dat[tmp,1])
        S = log10(dat[tmp,4])
        Z = lm(S ~ A)  #get slope and fit of points to the line (adj.r.squared)
        outcount = outcount + 1
        zValue[outcount,]<-c(plot[j], year[i], Z$coefficients[1], Z$coefficients[2], summary(Z)$adj.r.squared)
      }}}
  return(zValue)
}



##### STRs #####

Create_time_dataframe=function(dat, timespans, years){
  # Generates a dataframe for STRs with richness at all timescales for each plot, records data in a dataframe
  wa_time_richness = data.frame("plot"=1,"timespan"=1,"S"=1, "N"=1)
  lastYr=max(years)
  firstYr=min(years)
  plots = unique(sort(dat$plot))
  r = richness_in_time(dat, plots, lastYr, firstYr, timespans)
  wa_time_richness = rbind(wa_time_richness, r)
  
  return(wa_time_richness[-1,])
}


richness_in_time = function(dat, plots, lastYr, firstYr, timespans){
  # Calculates richness within a given plot for all possible timespans using a moving-window approach
  t_richness = data.frame("plot"=1, "timespan"=1, "S"=1, "N"=1)
  outCount = 0
  for (iPlot in 1:length(plots)){
    for(tSpan in 1:length(timespans)){
      for(startYr in firstYr:(lastYr-timespans[tSpan]+1)){
        outCount<-outCount+1
        tmp<-which(dat$plot==plots[iPlot] & dat$year>=startYr & dat$year<(startYr+timespans[tSpan]))
        if (nrow(dat[tmp,]) == 0) {
          S = 0
          N = 0
          plot = unique(dat[which(dat$plot==plots[iPlot]),2])}
        else {
          S<-length(unique(dat$species[tmp]))     # count species and individuals
          N<-sum(dat$abundance[tmp],na.rm=T)
          plot = unique(dat$plot[tmp])}
        t_richness[outCount,]<-c(plot, timespans[tSpan], S, N)
      }}}
  return(t_richness)}


Create_time_means=function(dat){
  # Gets mean richness at each timespan in each plot
  means = aggregate(dat$S,by=list(timespan=dat$timespan,plot=dat$plot),FUN=mean)
  names(means)[3]="S"
  N = aggregate(dat$N,by=list(timespan=dat$timespan,plot=dat$plot),FUN=mean)
  means[4] = N[,3]
  names(means)[4] = "N"
  return(means)
}


time_slope=function(dat, plot){
  # Gets STR slopes, w, for each plot using power law
  #set up dataframe to store results
  wValue = data.frame("plot" = 1, "intercept" = 1, "slope" = 1, "adjR2" = 1)
  outcount = 0
  for(iPlot in 1:length(plot)){   #loop thru plots
    tmp<-which(dat[,2]==plot[iPlot])
    Tspan = log10(dat[tmp,1])
    S = log10(dat[tmp,3])
    W = lm(S ~ Tspan)  #get slope and fit of points to the line (adj.r.squared)
    outcount=outcount + 1
    wValue[outcount,] <- c(plot[iPlot], W$coefficients[1], W$coefficients[2], summary(W)$adj.r.squared)
  }
  return(wValue)
}


########## FUNCTIONS FOR PLOTTING THE DATA #########################################################################

rad_poilog_plot = function(x_vec, plot, year){
  # Function to plot the RAD with the fitted poilog curve overlaid. XX
  # use with RAD_plot.r
  par_mle = poilogMLE(x_vec, startVals = c(mu = mean(log(x_vec)),
                                           sig = sd(log(x_vec))))$par
  mu = as.numeric(par_mle)[1]
  sig = as.numeric(par_mle)[2]
  S = length(x_vec)
  cdf_list = NULL
  for (i in 1:max(x_vec)){
    p = poilog_cdf(mu, sig, i)
    cdf_list = c(cdf_list, (1 - p) * S + 0.5)
  }
  # Plot black curve for fitted poilog. XX
  plot(cdf_list, seq(max(x_vec)), type = "l", lwd = 2, lty = 2, col = "red", ylim = c(1, max(x_vec)),
       xlab = "", ylab = "")
  # Plot empirical data as red points
  points(seq(S), sort(x_vec, decreasing = T), pch = 19, col = "black")
  mtext(side = 3, paste("plot", plot, "in", year, sep = " "), cex = .75)
}


RAD_plot=function(dat,year) {
  # plot RADs
  plot=unique(sort(dat[,2]))   #generate a list of the plots
  for (iYr in 1:length(year)) {   #loop thru years
    for (iPlot in 1:length(plot)){   #loop thru plots
      tmp <- which(dat[,2]==plot[iPlot] & dat[,1]==year[iYr])
      if (nrow(dat[tmp,]) > 0){  #make sure data exists for this plot-year combination
        abuns = as.numeric(dat[tmp,c(4:ncol(dat))])
        abuns = sort(abuns[abuns > 0], decreasing = T)
        if (length(abuns) >= 5) {  #Don't plot years when S < 5
          #Make Rank Abundance plot using species data and poilog prediction
          rad_poilog_plot(abuns, plot[iPlot], year[iYr])
        }}}}}


plot_SARs=function(dat, year, plot){
  # plot SARs
  for(iYr in 1:length(year)) {   #loop thru years
    for(iPlot in 1:length(plot)) {   #loop thru plots
      tmp<-which(dat[,3]==plot[iPlot] & dat[,2]==year[iYr])
      if (max(dat[tmp,4])>=5) {  #Don't plot years when S<5
        S = dat[tmp,4]
        A = dat[tmp,1]
        Z = lm(log10(S) ~ log10(A))
        plot(NA,NA,log="xy",xlab="Area",ylab="Species Richness",ylim=c(.1*min(dat[tmp,4]),max(dat[tmp,4])),
             xlim=c(0.25,4),main=paste("plot ",plot[iPlot], " in ", year[iYr],sep=""))
        par(new=TRUE)
        plot(dat[tmp,1],dat[tmp,4],log="xy",xlab="",ylab="", pch = 19, cex = 1.5,
             ylim=c(.1*min(dat[tmp,4]),max(dat[tmp,4])),xlim=c(0.25,4),type="p")
        abline(Z, col = "red", lty = 2, cex = 2, lwd = 2)
      }}}
}


plot_STRs=function(dat, plot){
  ## STR plotting function
  for(iPlot in 1:length(plot)) {
    tmp<-which(dat[,2]==plot[iPlot])
    Tspan = log10(dat[tmp,1])
    S = log10(dat[tmp,3])
    W = lm(S ~ Tspan)
    plot(NA,NA,log="xy",xlab="Timespan",ylab="Species Richness",ylim=c(min(dat[,3]),max(dat[,3])),
         xlim=c(1,15))
    par(new=TRUE)
    plot(dat[tmp,1],dat[tmp,3],log="xy",xlab="",ylab="", pch = 19,
         ylim=c(min(dat[,3]),max(dat[,3])),xlim=c(1,15),type="b")
    abline(W, col = "red", lty = 2, lwd = 2)
    mtext(side = 3, paste("plot", plot[iPlot], sep = " "), line = 0, cex = 0.75)
  }}



########## FUNCTIONS FOR STATISTICAL ANALYSIS #########################################################################


lmer_analysis = function(var_name, dat){
  #linear mixed effects model code used with RADs and SARs to look for statistical differences
  dat = add_trmt(dat)
  dat$year = as.factor(dat$year)
  dat$trmt = relevel(dat$trmt, ref = "C")
  ## Change: no transformation. Instead, check normality first.
  ans1 = lmer(var_name~trmt+(1|plot)+(1|year)+(1|trmt:year),data=dat)
  qqnorm(resid(ans1))
  qqline(resid(ans1))
  p1 = pvals.fnc(ans1, nsm = 10000, addPlot = F, withMCMC = F)$fixed
  EC = as.numeric(p1[2, 6])
  RC = as.numeric(p1[3, 6])
  dat$trmt = relevel(dat$trmt, ref = "E")
  ans2 = lmer(var_name~trmt+(1|plot)+(1|year)+(1|trmt:year),data=dat)
  p2 = pvals.fnc(ans2, nsm = 10000, addPlot = F, withMCMC = F)$fixed
  ER = as.numeric(p2[3, 6])
  return (list(EC = EC, RC = RC, ER = ER))
}

lmer_equi = function(var_name, dat, add = 0, pw = 1, bound_level = 0.05){
  #linear mixed effects model code used with RADs and SARs to look for statistical equivalence
  ## Define bound before transformation.
  dat = add_trmt(dat)
  dat$year = as.factor(dat$year)
  dat$trmt = relevel(dat$trmt, ref = "C")
  avg = mean(c(mean(var_name[dat$trmt=="E"]), mean(var_name[dat$trmt=="C"]),
               mean(var_name[dat$trmt=="R"])))
  y = (var_name + add)^pw
  bound1 = (avg * (1 + bound_level) + add)^pw - (avg + add)^pw  ## Compute post-transformation bounds
  bound2 = (avg * (1 - bound_level) + add)^pw - (avg + add)^pw
  bound_hi = max(bound1, bound2)
  bound_lo = min(bound1, bound2)
  ans1 = lmer(y~trmt+(1|plot)+(1|year)+(1|trmt:year),data=dat)
  mc1 = pvals.fnc(ans1, nsim = 10000, addPlot = F, withMCMC = T)
  EC = length(which(mc1$mcmc$trmtE <= bound_hi & mc1$mcmc$trmtE >= bound_lo)) / 10000
  RC = length(which(mc1$mcmc$trmtR <= bound_hi & mc1$mcmc$trmtR >= bound_lo)) / 10000
  dat$trmt = relevel(dat$trmt, ref = "E")
  ans2 = lmer(y~trmt+(1|plot)+(1|year)+(1|trmt:year),data=dat)
  mc2 = pvals.fnc(ans2, nsim = 10000, addPlot = F, withMCMC = T)
  ER = length(which(mc2$mcmc$trmtR <= bound_hi & mc2$mcmc$trmtR >= bound_lo)) / 10000
  return (list(EC = EC, RC = RC, ER = ER))
}

str_anova = function(var_name, dat){
  #ANOVA model code to be used with STRs to look for statistical differences
  dat = add_trmt(dat)
  dat$trmt = relevel(dat$trmt, ref = "C")
  ans1 = lm(var_name~trmt, data = dat)
  qqnorm(resid(ans1))
  qqline(resid(ans1))
  EC = coef(summary(ans1))[2, 4]
  RC = coef(summary(ans1))[3, 4]
  dat$trmt = relevel(dat$trmt, ref = "E")
  ans2 = lm(var_name~trmt, data = dat)
  ER = coef(summary(ans2))[3, 4]
  return (list(EC = EC, RC = RC, ER = ER))
}

bound_prob = function(row_in_sum, dof, bound){
  ## Returns the probability that the parameter estimate is within (-bound, bound)
  var_mean = as.numeric(row_in_sum[1])
  var_sd = as.numeric(row_in_sum[2])
  var_prob = pt((bound - var_mean) / var_sd, dof) - pt((-bound - var_mean) / var_sd, dof)
  return (var_prob)
}

str_equi = function(var_name, dat, bound_level = 0.05){
  #to be used with STRs to look for statistical equivalence
  require(equivalence)
  dat = add_trmt(dat)
  dat$trmt = relevel(dat$trmt, ref = "C")
  ans1 = lm(var_name~trmt, data = dat)
  bound = bound_level * (coef(summary(ans1))[1, 1] * 3 + coef(summary(ans1))[2, 1]
                         + coef(summary(ans1))[3, 1]) / 3
  dof = summary(ans1)$df[2]
  EC = bound_prob(coef(summary(ans1))[2, ], dof, bound)
  RC = bound_prob(coef(summary(ans1))[3, ], dof, bound)
  dat$trmt = relevel(dat$trmt, ref = "E")
  ans2 = lm(var_name~trmt, data = dat)
  ER = bound_prob(coef(summary(ans2))[3, ], dof, bound)
  return (list(EC = EC, RC = RC, ER = ER))
}

FDR_sig = function(p_vec){
  #False discovery rate correction for statistical differences
  p = sort(p_vec)
  i = length(p)
  alpha = 0.05
  while(p[i] > i * alpha / length(p) && i >= 1){
    i = i-1
  }
  return (i)
}

FDR_equi = function(p_vec){
  #False discovery rate correction for statistical equivalence
  p = sort(1 - p_vec)
  i = length(p)
  alpha_equi = 0.1
  while(p[i] > i * alpha_equi / length(p) && i >= 1){
    i = i - 1
  }
  return (i)
}


########## FUNCTIONS FOR GENERATING FIGURES IN THE PAPER ##############################################################

pred_sad=function(S, mu, sig){
  #get SAD predictions using poilog
  library(poilog)
  abd=vector(mode="numeric",length=S)
  rank=rev(1:S)
  for (i in 1:S){
    y.cdf=function(x) ppoilog(x, mu, sig)-(rank[i]-0.5)/S
    if (y.cdf(1)>0){abd[i]=1}
    else {
      abd[i]=round(uniroot(y.cdf,lower=1,upper=N)$root)
    }
  }
  abd=rev(abd)
  return (abd)
}


ppoilog = function(x, mu, sig){
  #cdf poilog
  x_list = seq(x)
  return (sum(dpoilog(x_list, mu, sig)) / (1 - dpoilog(0, mu, sig)))
}


get_abuns = function(abun_list){
  #To get mu and sigma from a list of empirical abundances, run:
  dat_par = as.list(poilogMLE(abun_list, startVals = c(mean(log(abun_list)), sd(log(abun_list))))$par)
  S = length(abun_list)
  mu = dat_par$mu
  sig = dat_par$sig
  plog_abunds = pred_sad(S, mu, sig)
  return(plog_abunds)
}


abd_solver=function(N, S, mu, sig){
  library(poilog)
  abd=vector(mode="numeric",length=S)
  rank=rev(1:S)
  for (i in 1:S){
    y.cdf=function(x) ppoilog(x, mu, sig)-(rank[i]-0.5)/S
    if(y.cdf(1)>0) {abd[i]=1}
    else {
      abd[i]=round(uniroot(y.cdf,lower=1,upper=N)$root)
    }
  }
  return(rev(abd))
}


rad_plot_eq = function(N, S, mu, sigma){
  #generates RAD plot for figure C1 -- equivalence testing -- in Appendix C
  rank = rev(1:S)
  abd = abd_solver(N, S, mu, sigma)
  abd_upper = abd_solver(N, S, mu * 1.05, sigma * 1.05)
  abd_lower = abd_solver(N, S, mu * 0.95, sigma * 0.95)
  plot(rank, abd, type = "b", pch = 19, lwd = 2, main = "RAD", xlab = "Rank", ylim = c(0, max(abd_upper)), 
       ylab = "Abundance", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex = 0.5)
  lines(rank, abd_lower, col = "#778899", lwd = 2, lty = "dashed", type = "b", pch = 19, cex = 0.5)
  lines(rank, abd_upper,  col = "#778899", lwd = 2, lty = "dashed", type = "b", pch = 19, cex = 0.5)
}


sar_str_plot = function(slope, intercept, plot_type){
  #generates SAR and STR plots for figure C1 -- equivalence testing -- in Appendix C
  if (plot_type == "SAR"){
    xlabp = "Area"
    xlower = 0.25
    xupper = 4}
  else {
    xlabp = "Timespan"
    xlower = 1
    xupper = 15}
  curve(exp(intercept)* x ^ slope, from = xlower, to = xupper, main = plot_type, lwd = 2,
        xlab = xlabp, ylab = "Number of species", log = "xy", cex.lab = 1.5, cex.main = 1.5, 
        cex.axis = 1.5)
  curve(exp(intercept * 0.95) * x ^ (slope * 0.95), add = TRUE, lwd = 2, lty = "dashed", col = "#778899")
  curve(exp(intercept * 1.05) * x ^ (slope * 1.05), add = TRUE, lwd = 2, lty = "dashed", col = "#778899")
}

rescale = function(value_list, bound_hi, bound_lo){
  ## Rescale values (i.e., mean, lower CI and upper CI) with respect to bounds
  ## Bounds might not be symmetrical
  value_new = NULL
  for (i in 1:length(value_list)){
    item = value_list[i]
    if (item < 0){
      value_new = c(value_new, item / abs(bound_lo))
    }
    else {
      value_new = c(value_new, item / abs(bound_hi))
    }
  }
  return (value_new)
}

lmer_quan = function(var_name, dat, add = 0, pw = 1, bound_level = 0.05){
  ## Get the two CIs rescaled with respect to bounds
  dat = add_trmt(dat)
  dat$year = as.factor(dat$year)
  dat$trmt = relevel(dat$trmt, ref = "C")
  y = (var_name + add)^pw
  avg = mean(c(mean(var_name[dat$trmt=="E"]), mean(var_name[dat$trmt=="C"]),
               mean(var_name[dat$trmt=="R"])))
  bound1 = (avg * (1 + bound_level) + add)^pw - (avg + add)^pw  ## Compute post-transformation bounds
  bound2 = (avg * (1 - bound_level) + add)^pw - (avg + add)^pw
  bound_hi = max(bound1, bound2)
  bound_lo = min(bound1, bound2)
  ans1 = lmer(y~trmt+(1|plot)+(1|year)+(1|trmt:year),data=dat)
  p1 = pvals.fnc(ans1, nsm = 10000, addPlot = F, withMCMC = T)
  EC_1 = rescale(as.numeric(c(p1$fixed[2, 1], quantile(p1$mcmc$trmtE, c(0.025, 0.975)))), bound_hi, bound_lo)
  EC_2 = rescale(as.numeric(c(p1$fixed[2, 1], quantile(p1$mcmc$trmtE, c(0.05, 0.95)))), bound_hi, bound_lo)
  RC_1 = rescale(as.numeric(c(p1$fixed[3, 1], quantile(p1$mcmc$trmtR, c(0.025, 0.975)))), bound_hi, bound_lo)
  RC_2 = rescale(as.numeric(c(p1$fixed[3, 1], quantile(p1$mcmc$trmtR, c(0.05, 0.95)))), bound_hi, bound_lo)
  dat$trmt = relevel(dat$trmt, ref = "R")
  ans2 = lmer(y~trmt+(1|plot)+(1|year)+(1|trmt:year),data=dat)
  p2 = pvals.fnc(ans2, nsm = 10000, addPlot = F, withMCMC = T)
  ER_1 = rescale(as.numeric(c(p2$fixed[3, 1], quantile(p2$mcmc$trmtE, c(0.025, 0.975)))), bound_hi, bound_lo)
  ER_2 = rescale(as.numeric(c(p2$fixed[3, 1], quantile(p2$mcmc$trmtE, c(0.05, 0.95)))), bound_hi, bound_lo)
  return (list(EC1 = EC_1, EC2 = EC_2, RC1 = RC_1, RC2 = RC_2,
               ER1 = ER_1, ER2 = ER_2))
}

CI = function(row_in_sum, dof, level = 0.05){
  ## Returns the CI for a variable based on a row taken from coef(summary(lm.model))
  var_mean = as.numeric(row_in_sum[1])
  var_sd = as.numeric(row_in_sum[2])
  t_critical = abs(qt(level / 2, dof))
  return (c(var_mean, var_mean - var_sd*t_critical, var_mean + var_sd*t_critical))
}

str_quan = function(var_name, dat, bound_level = 0.05){
  #STR quantiles
  dat = add_trmt(dat)
  dat$trmt = relevel(dat$trmt, ref = "C")
  ans1 = lm(var_name~trmt, data = dat)
  bound = bound_level * (coef(summary(ans1))[1, 1] * 3 + coef(summary(ans1))[2, 1]
                         + coef(summary(ans1))[3, 1]) / 3
  dof = summary(ans1)$df[2]
  EC_1 = rescale(CI(coef(summary(ans1))[2, ], dof), bound, -bound)
  EC_2 = rescale(CI(coef(summary(ans1))[2, ], dof, level = 0.1), bound, -bound)
  RC_1 = rescale(CI(coef(summary(ans1))[3, ], dof), bound, -bound)
  RC_2 = rescale(CI(coef(summary(ans1))[3, ], dof, level = 0.1), bound, -bound)
  dat$trmt = relevel(dat$trmt, ref = "R")
  ans2 = lm(var_name~trmt, data = dat)
  ER_1 = rescale(CI(coef(summary(ans2))[3, ], dof), bound, -bound)
  ER_2 = rescale(CI(coef(summary(ans2))[3, ], dof, level = 0.1), bound, -bound)
  return (list(EC1 = EC_1, EC2 = EC_2, RC1 = RC_1, RC2 = RC_2,
               ER1 = ER_1, ER2 = ER_2))
}


sub_plot = function(dat, typ, xax = F, yax = F){
  #Generates panels for Figure 2 in the paper -- summary of results
  require(gplots)
  ## Function to plot a sub-plot
  plot(0, 0, ylim = c(0, 15), type = "n", xlim = c(XMIN, XMAX),
       xaxt = "n", xlab = "", yaxt = "n", ylab = "", lwd = 1.2, cex.lab = 1, cex.axis = 1)
  if (xax){
    axis(side = 1, at = seq(-15, 5, 5), labels = T, cex.axis = 1.2)
    mtext("Rescaled difference", side = 1, line = 2.7)
  }
  if (yax){
    axis(side = 2, at = c(12, 3), labels = c("R-C", "K-C"), las = 1,
         tick = FALSE, cex.axis = 1.2, hadj = 1)
  }
  col_list_sn = c("#043A6B", "#A40004")
  col_list_par = c("#408DD2", "#FE3F44")
  symbol_list = c(15, 16, 21:25)
  y = 15
  for (i in 1:2){
    for (j in 1:7){
      wid = 1
      if (j == 1|j == 2){
        wid = 3
        col_plot = col_list_sn[i]
      }
      else {col_plot = col_list_par[i]}
      plotCI(x = dat[j, 3*i-2], y = y, ui = dat[j, 3*i], li = dat[j, 3*i-1],
             err = "x", col = col_plot, type = "b", gap = 0, add = TRUE, pch = symbol_list[j],
             cex = 1.4, lwd = wid, pt.bg = "white")
      y = y - 1
    }
    y = y - 2
  }
  if (typ == "point") {
    abline(v = 0, lty = "dashed")
  }
  else {
    abline(v = 1, lty = "dashed")
    abline(v = -1, lty = "dashed")
  }
}






# ############################################ EXTRANEOUS CODE?? #################################
# 
# # Get R2 from comparing observed abundances with predicted abundances
# obs_pred_rsquare = function(obs, pred){
#   return(1 - sum((obs - pred) ^ 2) / sum((obs - mean(obs)) ^ 2))
# }
# 
# ## Function: Calculate lambda_1 using (3) from Harte et al. 2009
# get_lambda_sad=function(S,N){
#   if (N<=0){
#     print("Error: N must be greater than 0.")
#     return(NA)}
#   else if (S<=0){
#     print("Error: S must be greater than 0.")
#     return(NA)}
#   else if (S/N>=1){
#     print("Error: N must be greater than S.")
#     return(NA)}
#   else {
#     ## Solve for lambda 1 in Harte et al. 2009 based on (3)
#     m=seq(1:N)
#     y=function(x) S/N*sum(x^m)-sum(x^m/m)
#     bound=10^-15  ## Define lower bound for x
#     ## Upper limit for p is 1.1. Haven't found any p value above 1.1 yet. 
#     p=uniroot(y,lower=bound,upper=1.005,tol=.Machine$double.eps^0.5)$root   ## Parameter for log-series
#     lambda_sad=-log(p)
#     return(list(p=p,lambda=lambda_sad))
#   }
# }
# # 
# # 
# # ######## Add Trmt to SAR_richness 
# # # add treatment
# # add_trmt = function(dat){
# #   for (i in 1:nrow(dat)){
# #     if(dat[i,3]==2 | dat[i,3]==4 | dat[i,3]==8 | dat[i,3]==11 |
# #       dat[i,3]==12 | dat[i,3]==14 | dat[i,3]==17 | dat[i,3]==22) {
# #       dat$trmt[i] = 1}
# #     else if (dat[i,3]==3 | dat[i,3]==6 | dat[i,3]==13 | dat[i,3]==15 |
# #       dat[i,3]==18 | dat[i,3]==19 | dat[i,3]==20 | dat[i,3]==21) {
# #       dat$trmt[i] = 2}
# #     else {dat$trmt[i] = 3} }
# #   return(dat)
# # }
# # 
# # # Compares models of line fit for SARs. Does power-law function work the best?
# # model_comp = function(dat){
# #   
# #   years = unique(sort(dat$year))
# #   plots = unique(sort(dat$plot))
# #   models = data.frame("year"=1, "plot"=1, "power"=1, "exponential"=1,
# #                       "negative_exp"=1, "monod"=1)
# #   outcount = 0
# #   for (iYr in 1:length(years)){
# #     for (iPlot in 1:length(plots)){
# #       tmp<-which(dat$plot == plots[iPlot] & dat$year == years[iYr])
# #       if (max(dat[tmp,4])>4){
# #         A_S_dat = cbind(dat[tmp,1],dat[tmp,4])
# #         source("C://Documents and Settings//sarah//My Documents//Active Research Projects//Rodent-Plant-MaxEnt//MacroEco_stuff_2011_2_23//SAR_comp.R")
# #         comp = SAR_comp(A_S_dat)
# #         outcount = outcount + 1
# #         models[outcount,] <- c(years[iYr], plots[iPlot], comp)
# #       }}}
# #   return(models)
# # }
# # 
# # # puts data in "wide" format CATEGORIZED BY TREATMENT for use with certain functions
# # reshape_data_by_trmt = function(dat){
# #   
# #   dat2 = aggregate(dat$abundance,by=list(year=dat$year, trmt=dat$trmt,
# #                                          species=dat$species),FUN=sum)
# #   names(dat2)[4]="abun"
# #   
# #   dat3 = reshape(dat2,idvar=c("year","trmt"),timevar="species",direction="wide")
# #   dat3[is.na(dat3)] = 0
# #   
# #   dummy = c(rep(NA,(nrow(dat3))))
# #   dat3[,ncol(dat3)+1] = dummy
# #   names(dat3)[ncol(dat3)] = "dummy"
# #   dat4 = dat3[,c(1,2,ncol(dat),3:ncol(dat)-1)]
# #   return(dat4)
# # }
# # 
# # converts data from raw abundance to relative abundance
# convert_relAbund = function(dat){
#   
#   columns = 4:ncol(dat)  # where species data begins at col 4. We don't want to change site info
#   
#   for (irow in 1:nrow(dat)) {
#     tot = sum(dat[irow,columns])   # total abundance in a row
#     for (icol in 1:length(columns)){
#       raw = dat[irow, columns[icol]]
#       rel = raw/tot                 # relative abundance of a cell in a row
#       dat[irow,columns[icol]] = rel   # fill new value into existing dataframe
#     }}
#   return (dat)
# }
# 
# #### Rank-abundance distribution relAbund, rank outputs
# RAD_relAbund=function(dat,year) {
#   require(vegan)
#   require(BiodiversityR)
#   
#   plot=unique(sort(dat[,2]))
#   #trmtcolor=c("darkgreen","blue","red")[as.numeric(dat$trmt)]
#   abundance = data.frame("year"=1, "plot"=1, "S"=1, "N"=1, "p"=1, "lambda"=1, "R2" = 1)
#   outcount = 0
#   
#   for (iYr in 1:length(year)) {
#     for (iPlot in 1:length(plot)){
#       tmp <- which(dat[,2]==plot[iPlot] & dat[,1]==year[iYr])
#       if (nrow(dat[tmp,]) > 0){
#         s = as.numeric(dat[tmp,c(4:ncol(dat))])
#         s = s[s>0]
#         if (length(s) >= 5) {  #Don't plot years when S < 5
#           #Make Rank Abundance plot using species data
#           #RankAbun <- rankabundance(dat[tmp,4:ncol(dat)])
#           #rankabunplot(RankAbun,scale='logabun',scaledx=F, specnames=c(1:3), ylim = c(1,max(s)),
#           #             col=trmtcolor[tmp],xlim=c(1,length(s)),main=paste("RAD",plot[iPlot],year[iYr],sep=","))
#           #save observed abundance and rank values
#           S = length(s)
#           N = sum(s)
#           #p = get_lambda_sad(S,N) #from get lambda_sad_new.r
#           # get R2 of obs v. pred # FIXME FIXME FIXME FIXME FIXME!!!!!!!!
#           obs = sort(s)
#           logn_vals = as.list(poilogMLE(s, startVals = c(mean(log(s)), sd(log(s))))$par)
#           pred = sum(log(dpoilog(s, logn_vals$mu, logn_vals$sig)))
#           R2 = obs_pred_rsquare(obs, pred)
#           # Record data in data table
#           outcount = outcount + 1
#           abundance[outcount,] <- c(year[iYr], plot[iPlot], S, N, p$p, p$lambda, R2)
#         }}}}
#   return(abundance)
# }
# 
# #### creates a two column matrix with column 1 = ID and col 2 = abun of a species
# generate_abun_list=function(dat,year) {
#   require(vegan)
#   require(BiodiversityR)
#   
#   plot=unique(sort(dat[,2]))
#   abundance_list = data.frame("ID_list" = 0, "s" = 0)
#   
#   for (iYr in 1:length(year)) {
#     for (iPlot in 1:length(plot)){
#       tmp <- which(dat[,2]==plot[iPlot] & dat[,1]==year[iYr])
#       if (nrow(dat[tmp,]) > 0){
#         s = as.numeric(dat[tmp,c(4:ncol(dat))])
#         s = s[s>0]
#         if (length(s) >= 5) {  #Don't plot years when S < 5
#           s = sort(s)
#           ID = as.numeric(paste(plot[iPlot],year[iYr],sep="."))
#           ID_list = rep(ID,length(s))
#           #copy to matrix
#           merge_dat = cbind(ID_list, s)
#           abundance_list = rbind(abundance_list, merge_dat)
#         }}}}
#   abundance_list = abundance_list[-1,] # delete meaningless first row
#   return(abundance_list)
# }


# # ########################################## Distance Decay across years ##########
# # distance_by_year = function(dat, plot_list, year_list){
# #   
# #   require(vegan)
# #   data_matrix = data.frame("plot" = 1, "bc" = 1, "jaccard" = 1, "canb" = 1,
# #                            "horn" = 1, "yr1" = 1, "yr2" = 1)
# #   outcount = 0
# #   
# #   for (iPlot in 1:length(plot_list)){
# #     for (iYr in 1:length(year_list)){
# #       tmp1 <- which(dat[,2]==plot_list[iPlot] & dat$year==year_list[iYr])
# #       tmp2 <- which(dat[,2]==plot_list[iPlot] & dat$year==year_list[iYr+1])
# #       # don't compare years that have no species
# #       if (nrow(dat[tmp1,]) == 0 | nrow(dat[tmp2,]) == 0) {
# #         bc_dist = NA
# #         j_dist = NA
# #         canb_dist = NA
# #         horn_dist = NA }
# #       # Calculate similarity metrics when both quads have > 0 species
# #       else {
# #         yr_to_yr = rbind(dat[tmp1,],dat[tmp2,])
# #         bc_dist = vegdist(yr_to_yr[,4:ncol(yr_to_yr)], method = "bray")    # all range 0 (same) to 1 (different)
# #         j_dist = vegdist(yr_to_yr[,4:ncol(yr_to_yr)], method = "jaccard")       #CHOOSE BEST OPTION(S)!!
# #         canb_dist = vegdist(yr_to_yr[,4:ncol(yr_to_yr)],method = "canberra")
# #         horn_dist = vegdist(yr_to_yr[,4:ncol(yr_to_yr)],method = "horn")
# #       }
# #       outcount <- outcount + 1
# #       data_matrix[outcount,] = c(plot_list[iPlot], bc_dist, j_dist,
# #                                  canb_dist, horn_dist, year[iYr], year[iYr+1])
# #     }}
# #   return(data_matrix)
# # }
# # 
# # ################ adonis by year
# # adonis_by_year = function(dat, year_list){
# #   
# #   require(vegan)
# #   adonis_table = data.frame("year" = 1, "bray" = 1, "horn" = 1, "canb" = 1, "jaccard" = 1)
# #   outcount = 0
# #   
# #   for (y in 1:length(year_list)){
# #     tmp <- which(dat$year==year_list[y])  #make temporal window
# #     spp_data = subset(dat[tmp,], select=c(4:ncol(dat))) #This used to create the distance matrix
# #     env_data = subset(dat[tmp,], select=c(1:3))  #This is for the adonis function
# #     if (nrow(spp_data) > 10) { # make sure over half the plots have something growing!
# #       # make distance matrices
# #       spp_bc = vegdist(spp_data, method = "bray")
# #       spp_horn = vegdist(spp_data, method = "horn")
# #       spp_canb = vegdist(spp_data, method = "canb")
# #       spp_jac = vegdist(spp_data, method = "jaccard")
# #       # analyze distance matrices
# #       bc = adonis(spp_bc ~ trmt, data = env_data, permutation=1000)
# #       horn = adonis(spp_horn ~ trmt, data = env_data, permutation=1000)
# #       canb = adonis(spp_canb ~ trmt, data = env_data, permutation=1000)
# #       jac = adonis(spp_canb ~ trmt, data = env_data, permutation=1000)
# #       # record output in table
# #       outcount <- outcount + 1
# #       adonis_table[outcount,] = c(year_list[y], bc$aov.tab[1,6], horn$aov.tab[1,6], 
# #                                   canb$aov.tab[1,6], jac$aov.tab[1,6])
# #       outcount <- outcount + 1
# #       adonis_table[outcount,] = c("r2", bc$aov.tab[1,5], horn$aov.tab[1,5], 
# #                                   canb$aov.tab[1,5], jac$aov.tab[1,5])
# #     }}
# #   return(adonis_table)
# # }
# # 
# # ######### adonis of each treatment pair
# # adonis_by_pairs = function(dat){
# #   
# #   require(vegan)
# #   adonis_pair_table = data.frame("method" = 1, "CR" = 1, "CE" = 1, "RE" = 1)
# #   outcount = 0
# #   
# #   methods = c("bray", "horn", "canb", "jaccard")
# #   
# #   for (m in 1:length(methods)){
# #     # subset out pairs of treatments
#     dat_CoRe = subset(dat, trmt != 3)
#     dat_CoEx = subset(dat, trmt != 2)
#     dat_ReEx = subset(dat, trmt != 1)
#     # make species and env tables
#     spp_CoRe = vegdist(dat_CoRe[4:ncol(dat_CoRe)], method = methods[m]) 
#     env_CoRe = subset(dat_CoRe, select=c(1:3)) 
#     spp_CoEx = vegdist(dat_CoEx[4:ncol(dat_CoEx)], method = methods[m]) 
#     env_CoEx = subset(dat_CoEx, select=c(1:3)) 
#     spp_ReEx = vegdist(dat_ReEx[4:ncol(dat_ReEx)], method = methods[m]) 
#     env_ReEx = subset(dat_ReEx, select=c(1:3)) 
#     # do the adonis
#     CR = adonis(spp_CoRe ~ trmt, data = env_CoRe, permutations = 1000, strata=env_CoRe$year)
#     CE = adonis(spp_CoEx ~ trmt, data = env_CoEx, permutations = 1000, strata=env_CoEx$year)
#     RE = adonis(spp_ReEx ~ trmt, data = env_ReEx, permutations = 1000, strata=env_ReEx$year)
#     # fill the table
#     outcount = outcount + 1
#     adonis_pair_table[outcount,] = c(methods[m], CR$aov.tab[1,6], CE$aov.tab[1,6], RE$aov.tab[1,6])
#     outcount = outcount + 1
#     adonis_pair_table[outcount,] = c("r2", CR$aov.tab[1,5], CE$aov.tab[1,5], RE$aov.tab[1,5])
#   }
#   return(adonis_pair_table)
# }
# 
# ## power transform function
# transform_to_normal = function(vector) {
#   require(car)
#   yj1 = powerTransform(vector, family = "yjPower")
#   yj = as.numeric(yj1[7])   
#   transformed_vector = as.numeric()
#   for (i in 1:length(vector)){
#     trans = yj * vector[i]
#     transformed_vector <- append(trans, transformed_vector)
#     tranformed_vector = yj*vector }
#   return(transformed_vector)
# } 
# 
# ## Comparing the fit of four SAR functional forms with 2 parameters
# ## See Table S2 in Guilhaumon et al. 2008
# ## Simple function to calculate residual sum of squares
# res2 = function(model, dat){
#   S = dat[, 2]
#   pred = predict(model)
#   res = pred - S
#   return (sum(res^2))
# }
# ## input: a dataframe with two columns (A and S)
# ## ouput: residual sum of squares for the four models
# SAR_comp = function(dat){
#   A = dat[, 1]
#   S = dat[, 2]
#   ## Power
#   mod_power = lm(log(S) ~ log(A))
#   ## Cannot use the above function b/c log scale
#   c = exp(mod_power$coefficients[1])
#   z = mod_power$coefficients[2]
#   pred = c * A^z
#   res2_power = sum((pred - S)^2)
#   ## Exponential
#   mod_exp = lm(S ~ log(A))
#   res2_exp = res2(mod_exp, dat)
#   ## Negative exponential
#   mod_neg = nls(S ~ c * (1 - exp(- z * A)), start = list(c = log(max(S)), z = 0.1),
#                 control = nls.control(maxiter = 2000, warnOnly = TRUE))
#   res2_neg = res2(mod_neg, dat)
#   ## Monod
#   mod_monod = nls(S ~ (c * A) / (z + A), start = list(c = max(S), z = 0.1),
#                   control = nls.control(maxiter = 2000, warnOnly = TRUE))
#   res2_monod = res2(mod_monod, dat)
#   return (c(res2_power, res2_exp, res2_neg, res2_monod))
# }

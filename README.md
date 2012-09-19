portal-experimental-macroeco
============================

Code for reproducing the results from "An experimental test of the response of macroecological patterns to altered species interactions" by Sarah R. Supp, Xiao Xiao, S. K. Morgan Ernest, and Ethan P. White 2013 published in Ecology, doi: 10.1890/12-0370.1 Code was written by Sarah R. Supp and Xiao Xiao.

The code and data in this repository allow for the analyses and figures in the paper to be fully replicated using a subset of the published Portal dataset which includes annual plant data from 1995-2009.

Requirements: R 2.x and the following packages: Biodiversity R, car, CCA, gplots, graphics, languageR, lme4, nlme, plotrix, poilog, vegan, VGAM and the file containing functions specific to this code, PortalPlants_fxns.R. 

The analyses can then be replicated by changing the working directory at the top of the file PortalPlants_ms12-0370R2.R to the location on your computer where you have stored the .R and .csv files.

Please note that the pvalues generated for Appendix D in the published paper were done using R 2.12.2. Because of approximations, the values for SAD sigma and mu may differ slightly (around the 10th decimal place) from Appendix D, Tables S2 and S3. Because the equivalence testing also uses approximations, there may be very small differences in the exact values generated compared to Appendix D, tables S4 and S5.

It should take approximately 30 minutes to run all the code from start to finish. Figures should output as pdfs in your working directory.

Data use: Data is provided in this supplement for the purposes of replication. If you wish to use the data for additional research, they should be obtained from the published source (Ecological Archives E090-118-D1; S. K. Morgan Ernest, Thomas J. Valone, and James H. Brown. 2009. Long-term monitoring and experimental manipulation of a Chihuahuan Desert ecosystem near Portal, Arizona, USA. Ecology 90:1708. doi:10.1890/08-1222.1)

Included Files
============================

* PortalPlants_ms12-0370R2.R script -- cleans up the data, constructs the macroecological patterns, pulls out descriptive parameters of these patterns, runs the statistical analyses, and outputs figures.  
* PortalPlants_fxns.R script -- holds the relevant functions for executing the PortalPlants_ms12-0370R2.R script.  
* PortalSummerAnnuals_1995_2009.csv -- summer annual plant community abundance data
* PortalWinterAnnuals_1995_2009.csv -- winter annual plant community abundance data

License
============================

This code is available under a BSD 2-Clause License.

Copyright (c) 2012 Sarah R. Supp. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact information
============================
Sarah Supp's email: sarah@weecology.org and sarah.supp@usu.edu

Sarah's website: http://weecology.org/people/sarahsupp/Sarah_Supp/

author: Corey James MYERS (corey.james.myers@cern.ch)

# Postprocessing macros for Xi(1820)- in pp 13TeV #

- cinvmass4.c
- overhead.C
- pinvmass7.c

## Description ##

The following ReadMe is an explanation of how to use the Xi(1820) macros for pp 13 TeV data. This ReadMe may be updated at a later data. The following codes are still under construction for optimization, but are still viable.

- overhead.C: This macro is a template for many of the conventions used in the next 2 codes, such as multiplicity percentile binning and pt binning. It also is important for definitions of normalization ranges and naming conventions.

- cinvmass4.c: This macro should be the first macro run for analysis of the Xi(1820) data. It’s primary purpose is to take the 3-d histograms of the invariant mass plots of Lambda+K0s and Lambda+K(-+) and use them to construct background subtracted histograms of mixed-event or like-charge data subtracted from raw data. This data is organized by its normalization range, pt binning, and multiplicity percentile binning. Also, several histograms are constructed for analysis to Monte Carlo data. These histograms are useful for efficiency corrections and resolution calculations.

- pinvmass7.c: This should be the second macro run for the analysis of the Xi(1820) data. Its primary purpose is to calculate and fit the invariant mass plots constructed from cinvmass4.c. It does this though multiple functions and procedures.
The simplified procedure is as follows: (this is not how the code is structured exactly, only a logical understanding)



1.) Collect background-subtracted histogram from cinvmass4.c that corresponds to specific decay channel, normalization range, background distribution, pt bin, and multiplicity percentile bin.
<br> 2.) Due to the fact that a residual background may be present, a polynomial function is used to fit the residual background while ignoring the signal.
3.) The polynomial used to fit the residual background is looped though several degrees of polynomial (1st to 4th degree polynomial). This is called the “version” loop.
4.) The polynomial used to fit the residual background is looped though several fit ranges (1.70-2.0, 1.76-1.9, 1.76-2.0 GeV/c^2). This is called the “range” loop.
5.) After the residual background is fit with the polynomial (technically also fit with Breit-Wigner + polynomial), the difference between the polynomial fit the background-subtracted histogram is calculated and used to fill another histogram. This histogram would contain only the signal we are looking for. Note the nomenclature for these plots is that the plots before the polynomial subtraction are called “sub” plots, while the plots after the polynomial subtraction are called “linear” plots. 
6.) The “linear” plots are fitted with a Voigtian function (a convolution between a Breit-Wigner function for the resonance and a gaussian for the defector resolution). The sigma of the Voigtian is either fixed to 2 MeV/c^2, or free to moved between 1-3 MeV/c^2 (more details in the MC section later). This is called the “type” loop.
7.) The yield is calculated 2 times as either the integral of the Voigtian function from section 6, or the bin counting of the histogram. These 2 yields are stored in arrays for further use in later sections. (This section is still under construction to fully account for statistical errors and systematic procedures)
8.) The mean and width for the voigtian fit (with associated errors) are stored for further use in later sections.
9.) Both “sub” and “linear” plots are printed to pdfs for storage.
10.) Repeat steps 1-9 for every background, normalization range, “version”, “range”, and “type”.
11.)  For the specific decay channel, multiplicity, and pt used, calculate the average yield, mean, and width (with errors) for each background, normalization range, “version”, “range” and “type” that satisfies the following conditions:
A.) The selection of background, normalization range, “version”, “range” and “type” has been selected to offer either the lowest chi-squared fit for the data, or simply fit the data the best, and is selected as the “default” conditions. For Lambda+K(-+), the “default” selection is like-charge (mixed-event is kept as separate measurement, but similar “default” conditions), normalization range 2, 2nd-degree polynomial “version 2” or “v2”, 1.70 - 2.0 GeV/c^2 fit range “range 0”, Voigtian with fixed sigma of 2 MeV/c^2 “type 0”. For Lambda+K0s, the “default” selection is mixed-event, normalization range 0, “v2”, “range 0”, “type 0”.
B.) The selection of normalization range, “version”, “range”, and “type” is only one variation from the “default” selection. For example: Lambda+K(-+), normalization range 2, !“v3”!, “range 0”, “type 0” is included in the calculations. Lambda+K(-+), normalization range !3!, !”v4”!, “range !2!”, “type !1!” is not included in the calculations.
12.) The final yield is calculated as the average between the integral and bin counting. The systematic error is calculated as the difference between the average and integral or bin counting yields. (This section is still under construction)
12. A.) The average of the final yields that satisfy conditions 11A and B and used to calculate the “total” average yield. The standard deviation of the final yields that satisfy conditions 11A and B is set as the systematic uncertainty, with the systematic uncertainty from 12 added in quadrature, to set the “total” systematic uncertainty. (This section is still under construction)
13.) The average of the means that satisfy conditions 11A and B are used to calculate “total” mean. The standard deviation of the means that satisfy conditions 11A and B is set as the “total” systematic uncertainty. The errors of the means from the fits that satiety conditions 11A and B are added in quadrature to calculate the “total” statistical error. (This section is still under construction)
14.) The average of the widths that satisfy conditions 11A and B are used to calculate “total” width. The standard deviation of the widths that satisfy conditions 11A and B is set as the “total” systematic uncertainty. The errors of the widths from the fits that satiety conditions 11A and B are added in quadrature to calculate the “total” statistical error. (This section is still under construction)
15.) The “total” averaged yields, means, and widths are outputted to a text file (StatisticalDataValues.txt) for storage.
16.) The contributions of each variation of normalization range, “version”, “range”, and “type” to each systematic uncertainty is outputted to a text file (Percentage.txt) for storage.
17.) Repeat steps 1-16 for each decay channel, background, multiplicity, and pt bin.
18.) If high-multiplicity triggered data is available, repeat steps 1-17, but only one multiplicity (0-100% is used as the stand in for the correct 0-0.1% multiplicity percentile) is used. (This section is still under construction)

In the event that Monte Carlo data is available, then a few changes to the previous procedures are implemented. (This section is still under construction)
Around section 1.) Monte Carlo resolution data is fitted with a gaussian. The sigma is stored for later use. (This section is still under construction)
Around section 1.)  Monte Carlo acceptance * efficiency data is fitted with a linear function to estimate the acceptance * efficiency factor (constant) and any mass effects (linear). This information is stored for later use. (This section is still under construction)
Around section 6.) The Voigtian fit function is modified to account for any correction effect found in the slope of the acceptance * efficiency calculations. See funVoightMC. Correction is approximated by including (Voigtian * (1+0.15*(slope)*(Erf((x-mean)/(0.15*1.1283416)))) ).
Around section 6.) The sigma of the voigtian is fixed to the sigma of the gaussian of the resolution histogram (not 2 Mev/c^2). This is set as the “default” fit. 
Around section 15.) The “total” averaged yields and errors are corrected for Monte Carlo efficiency corrections.

In the event that different PID cuts are selected, then a few changes too the previous procedures are implemented. (This section is still under construction) 
Around section 1.) A “cutloop” is implemented that corresponds different PID cuts to different histograms created from cinvmass4.c. Currently, only 2 cuts for Lambda+K(-+) are used, but more are expected.
Around sections 6-8.) The arrays that store the yields, means and widths are modified to store them as different “type” variations. For example: type0, normal PID cut, Voigtian with fixed sigma. type1, normal PID cut, Voigtian with free sigma. type2, cut1 PID, Voigtian with fixed sigma. type3, cut1 PID, Voigtian with free sigma. type4, cut2 PID, Voigtian with fixed sigma. type5, cut2 PID, Voigtian with free sigma.
Around sections 11-13.) Because the “default” section is the same, only one variation form each PID change is accounted for. So, only cut1 PID (type2) and cut2 PID (type4) that match all other selection criteria are accounted for.

There are many more details of the codes that will be explained at a later date.

To use the codes, use the following procedure:
1.) Make sure your data is in a format similar to the ones used in AddTaskRare_pp13.C
2.) Make sure the data and codes overhead.C, cinvmass4.c, and pinvmass7.c are in the directory. 
3.) Make sure you have the correct directory format to store the output and input from cinvmass4.c and pinvmass7.c.
4.) Make sure you have a correct assortment of switches in cinvmass4.c and pinvmass7.c to allow the code to work. (This section is still under construction). For example, for civmass4.c
int hmtflag=1;//if hmtflag=1, look for hmt events.
int cutflag=0;//if cutflag=1, begin alternate cut process.
//int cutmax=0;
int cutloop=0;
int cutstart=0;
int MCflag=0;//if MCflag=1, MC information is run
For pinvmass7.c:
only1flag=1;//using only1flag will remove the average and lambda+KX mixed event calcualtions
MCflag=0;
cutflag=0;//1,2
hmtflag=1;//if hmtflag==1, activate loop
writeplots=1;//if writeplots=1, writes mass and sub plots
//int writeall=0;//if writeall=1, writes allmass plots
fitflag=1;//if fitflag=1, fits are preformed
fitsubtract=1;//if fitsubtract=1, fits are subtracted
typemax=0;
newflag=1;
Use1NR=0;//1
UseLinear=1;
Nsigmaflag=0; //if Nsigmaflag=1, use Nsigma calculations
const int versionmax=4;//6+1
writetxt=1;//if writetxt=1, txt file is written
Use1PT=0;
UseNoPoly=1;
Use1Type=0;//if 1, Only one type is used
UseSigmaRange=1;//if 1, use sigma as range
UseSigmaRange2=0;
UseVoigt=1;//if 1, use voigt
Fixmean=0;//if 1, fix mean for fit functions
onechange=1;//if onechange=1, use only the systematic calcualtions that are "one" change from default (free sigma + 3rd, free sigma + 4th etc)
5.) If all conditions are meet, then in aliroot:
6.) [type] .x $HOME/Desktop/xiTest/cinvmass4.c(0,”Lambdak0") [enter]
7.) [type] .q [enter] //some of the canvases are resused between codes, so exiting and returning prevents some errors.
8.) Check if root files were created correctly and filled. Return to aliroot
9.) [type] .x $HOME/Desktop/xiTest/pinvmass7.c(0,”Lambdak0") [enter]

The ReadMe will be expanded at a later date to include more detailed explanations of the code and procedures.

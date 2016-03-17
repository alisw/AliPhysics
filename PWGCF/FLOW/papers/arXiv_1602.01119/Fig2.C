// For the lists of output files see the bottom of this macro.

static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

void Fig2(Int_t rWrite = 0, Int_t rPerformance = 0)
{
 myOptions();
 gROOT->ForceStyle();   
 TDatime now;
 int iDate = now.GetDate();
 int iYear=iDate/10000;
 int iMonth=(iDate%10000)/100;
 int iDay=iDate%100;
 char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
 char cStamp1[25],cStamp2[25];
 sprintf(cStamp1,"%i %s %i",iDay,cMonth[iMonth-1],iYear);
 sprintf(cStamp2,"%i/%.2d/%i",iDay,iMonth,iYear);
    
    
    // Durham format
    /*
     Int_t nPoints = fv22Gap10_cent05->GetN();
     Double_t x = 0.;
     Double_t y = 0.;
     Double_t yErr = 0.;
     printf("\p_{t}-binCenter  v_{2}{2}     +-stat.error \n");
     for(Int_t p=0;p<nPoints;p++)
     {
     fv22Gap10_cent05->GetPoint(p,x,y);
     yErr = fv22Gap10_cent05->GetErrorY(p);
     printf("%f  %f  +-%f\n",x,y,yErr);
     } // end of for(Int_t p=0;p<nPoints;p++)
     cout<<endl;
     */
    
    
    Double_t PtBin[] = {0.3, 0.5, 0.7, 0.9, 1.125, 1.375, 1.625, 1.875, 2.25, 2.75, 3.25, 3.75, 4.5, 5.5};
    Double_t PtBinErr[20] = {0.};
    
    Double_t v22Gap10Run2_cent05[]={0.015355, 0.0249461, 0.0282724, 0.0347966, 0.0405938, 0.0494957, 0.0572995, 0.0586614, 0.0630034, 0.0666916, 0.0812318, 0.0489493, 0.0795362};
    Double_t v22Gap10Run2_cent05_Err[]={0.00103295, 0.00119824, 0.000904021, 0.000761121, 0.00122489, 0.00197264, 0.0030834, 0.00287649, 0.00275123, 0.00366634, 0.0102748, 0.0116112, 0.0115641};
    Double_t v22Gap10Run2_cent05_Sys[]={0.000568137, 0.000923005, 0.00104608, 0.00128747, 0.00101484, 0.00123739, 0.00143249, 0.00146653, 0.00157508, 0.00166729, 0.00203079, 0.00122373, 0.0019884};
    Double_t v22Gap10Run2_cent05CombErr[]={0.00117889, 0.00151252, 0.00138258, 0.00149563, 0.00159068, 0.00232862, 0.00339991, 0.00322876, 0.0031702, 0.00402764, 0.0104735, 0.0116755, 0.0117338};
 
    TGraphErrors *fv22Gap10_cent05 = new TGraphErrors(13,PtBin,v22Gap10Run2_cent05,PtBinErr,v22Gap10Run2_cent05CombErr);
    fv22Gap10_cent05->SetMarkerStyle(20);
    fv22Gap10_cent05->SetMarkerSize(1.00);
    fv22Gap10_cent05->SetMarkerColor(kRed+1);
    fv22Gap10_cent05->SetLineColor(kRed+1);
    fv22Gap10_cent05->SetLineWidth(2);
    
     Double_t v22Gap10Run2_cent3040[20]={0.0476818, 0.0770824, 0.102426, 0.124908, 0.143984, 0.16946, 0.191963, 0.206824, 0.222508, 0.239456, 0.249393, 0.244855, 0.198};
     Double_t v22Gap10Run2_cent3040_Err[20]={0.000747631, 0.000668875, 0.000674352, 0.001619, 0.00205859, 0.00190391, 0.00348401, 0.00223386, 0.00221814, 0.00415282, 0.00267322, 0.00964771, 0.00823685};
     Double_t v22Gap10Run2_cent3040_Sys[20]={0.00143045, 0.00231247, 0.00307277, 0.00374723, 0.00172781, 0.00203352, 0.00230356, 0.00248188, 0.0026701, 0.00287348, 0.00299271, 0.00293827, 0.002376};
     Double_t v22Gap10Run2_cent3040CombErr[20]={0.00161405, 0.00240726, 0.0031459, 0.00408202, 0.00268758, 0.00278569, 0.00417668, 0.00333914, 0.00347125, 0.00505003, 0.00401278, 0.0100852, 0.0085727};
    
    TGraphErrors *fv22Gap10_cent3040 = new TGraphErrors(13,PtBin,v22Gap10Run2_cent3040,PtBinErr,v22Gap10Run2_cent3040CombErr);
    fv22Gap10_cent3040->SetMarkerStyle(20);
    fv22Gap10_cent3040->SetMarkerSize(1.00);
    fv22Gap10_cent3040->SetMarkerColor(kRed+1);
    fv22Gap10_cent3040->SetLineColor(kRed+1);
    fv22Gap10_cent3040->SetLineWidth(2);

    Double_t v32Gap10Run2_cent05[20]={0.0061003, 0.013519, 0.0191132, 0.0273897, 0.0347368, 0.0428637, 0.0511093, 0.0476368, 0.0677752, 0.0830096, 0.110166, 0.110179, 0.0765003};
    Double_t v32Gap10Run2_cent05_Err[20]={0.0015357, 0.00187404, 0.00127069, 0.00161112, 0.00206893, 0.00100997, 0.00240842, 0.00372999, 0.00517277, 0.00732543, 0.00460911, 0.00992513, 0.017338};
    Double_t v32Gap10Run2_cent05_Sys[20]={0.000176909, 0.000392051, 0.000554282, 0.000794301, 0.00100737, 0.00124305, 0.00148217, 0.00138147, 0.00196548, 0.00240728, 0.00319482, 0.00319519, 0.00221851};
    Double_t v32Gap10Run2_cent05CombErr[20]={0.00154586, 0.00191461, 0.00138632, 0.00179628, 0.00230114, 0.00160162, 0.00282796, 0.00397759, 0.0055336, 0.00771083, 0.0056081, 0.0104268, 0.0174793};
    
    TGraphErrors *fv32Gap10_cent05 = new TGraphErrors(13,PtBin,v32Gap10Run2_cent05,PtBinErr,v32Gap10Run2_cent05CombErr);
    fv32Gap10_cent05->SetMarkerStyle(21);
    fv32Gap10_cent05->SetMarkerSize(1.000);
    fv32Gap10_cent05->SetMarkerColor(kBlue+1);
    fv32Gap10_cent05->SetLineColor(kBlue+1);
    fv32Gap10_cent05->SetLineWidth(2);
    
    Double_t v32Gap10Run2_cent3040[20]={0.0110293, 0.0241038, 0.0346258, 0.0428405, 0.0518479, 0.063612, 0.0746615, 0.0855818, 0.101463, 0.105487, 0.135066, 0.119134, 0.102021};
    Double_t v32Gap10Run2_cent3040_Err[20]={0.00138643, 0.00164461, 0.00228094, 0.00166515, 0.0024811, 0.00276232, 0.00411546, 0.00546156, 0.00624752, 0.00711386, 0.00358371, 0.0153258, 0.0196387};
    Double_t v32Gap10Run2_cent3040_Sys[20]={0.000209557, 0.000457972, 0.00065789, 0.00081397, 0.000985111, 0.00120863, 0.00141857, 0.00162605, 0.00192779, 0.00200425, 0.00256626, 0.00226354, 0.0019384};
    Double_t v32Gap10Run2_cent3040CombErr[20]={0.00140218, 0.00170719, 0.00237392, 0.00185345, 0.00266951, 0.00301516, 0.00435308, 0.00569848, 0.00653819, 0.00739081, 0.0044078, 0.015492, 0.0197342};
    
    TGraphErrors *fv32Gap10_cent3040 = new TGraphErrors(13,PtBin,v32Gap10Run2_cent3040,PtBinErr,v32Gap10Run2_cent3040CombErr);
    fv32Gap10_cent3040->SetMarkerStyle(21);
    fv32Gap10_cent3040->SetMarkerSize(1.000);
    fv32Gap10_cent3040->SetMarkerColor(kBlue+1);
    fv32Gap10_cent3040->SetLineColor(kBlue+1);
    fv32Gap10_cent3040->SetLineWidth(2);
    
    Double_t v42Gap10Run2_cent05[20]={0.00298618, 0.0079179, 0.00924839, 0.0107313, 0.0195288, 0.0127875, 0.0327479, 0.0424276, 0.0524524, 0.0488274, 0.0466606, 0.0552678, 0.0812232};
    Double_t v42Gap10Run2_cent05_Err[20]={0.00125167, 0.00208424, 0.0025715, 0.00397315, 0.00381933, 0.00437264, 0.00492297, 0.00800186, 0.00769408, 0.0171918, 0.0200979, 0.016969, 0.0363349};
    Double_t v42Gap10Run2_cent05_Sys[20]={5.37513e-05, 0.000142522, 0.000166471, 0.000193163, 0.000351519, 0.000230175, 0.000589463, 0.000763697, 0.000944143, 0.000878893, 0.00083989, 0.00099482, 0.00146202};
    Double_t v42Gap10Run2_cent05CombErr[20]={0.00125282, 0.00208911, 0.00257688, 0.00397784, 0.00383547, 0.0043787, 0.00495814, 0.00803822, 0.00775179, 0.0172142, 0.0201154, 0.0169981, 0.0363643};
    
    TGraphErrors *fv42Gap10_cent05 = new TGraphErrors(13,PtBin,v42Gap10Run2_cent05,PtBinErr,v42Gap10Run2_cent05CombErr);
    fv42Gap10_cent05->SetMarkerStyle(33);
    fv42Gap10_cent05->SetMarkerSize(1.500);
    fv42Gap10_cent05->SetMarkerColor(kGreen+3);
    fv42Gap10_cent05->SetLineColor(kGreen+3);
    fv42Gap10_cent05->SetLineWidth(2);
    
    Double_t v42Gap10Run2_cent3040[20]={0.00462444, 0.00629752, 0.0188221, 0.0148352, 0.0301152, 0.0511986, 0.0590605, 0.0648285, 0.0888981, 0.103704, 0.0654188, 0.0886088, 0.0813909};
    Double_t v42Gap10Run2_cent3040_Err[20]={0.00383527, 0.00319158, 0.00563652, 0.00621685, 0.00778161, 0.00876169, 0.0142432, 0.0119483, 0.0122702, 0.0143192, 0.0229749, 0.043931, 0.0554728};
    Double_t v42Gap10Run2_cent3040_Sys[20]={8.78644e-05, 0.000119653, 0.000357619, 0.000281868, 0.000572189, 0.000972774, 0.00112215, 0.00123174, 0.00168906, 0.00197037, 0.00124296, 0.00168357, 0.00154643};
    Double_t v42Gap10Run2_cent3040CombErr[20]={0.00383627, 0.00319383, 0.00564785, 0.00622324, 0.00780262, 0.00881552, 0.0142874, 0.0120116, 0.0123859, 0.0144542, 0.0230085, 0.0439633, 0.0554943};
    
    TGraphErrors *fv42Gap10_cent3040 = new TGraphErrors(13,PtBin,v42Gap10Run2_cent3040,PtBinErr,v42Gap10Run2_cent3040CombErr);
    fv42Gap10_cent3040->SetMarkerStyle(33);
    fv42Gap10_cent3040->SetMarkerSize(1.500);
    fv42Gap10_cent3040->SetMarkerColor(kGreen+3);
    fv42Gap10_cent3040->SetLineColor(kGreen+3);
    fv42Gap10_cent3040->SetLineWidth(2);
    

    //10-20
    Double_t v24Run2_cent1020[20]={0.0288237, 0.0436239, 0.0586879, 0.0697475, 0.0845251, 0.0962882, 0.107704, 0.118858, 0.133664, 0.147125, 0.152484, 0.153556, 0.136461};
    Double_t v24Run2_cent1020_Err[20]={0.000470525, 0.000454916, 0.000620154, 0.00100315, 0.000535354, 0.00110078, 0.00135911, 0.00150855, 0.00259672, 0.00241045, 0.00318281, 0.00611818, 0.0071198};
    Double_t v24Run2_cent1020_Sys[20]={0.000547651, 0.000828855, 0.00111507, 0.0013252, 0.00160598, 0.00182948, 0.00204638, 0.00225829, 0.00253962, 0.00279537, 0.00289719, 0.00291756, 0.00259277};
    Double_t v24Run2_cent1020CombErr[20]={0.000722022, 0.000945489, 0.00127592, 0.00166207, 0.00169286, 0.00213511, 0.00245659, 0.00271581, 0.00363216, 0.00369112, 0.00430395, 0.00677822, 0.0075772};
    
    TGraphErrors *fv24_cent1020 = new TGraphErrors(13,PtBin,v24Run2_cent1020,PtBinErr,v24Run2_cent1020CombErr);
    fv24_cent1020->SetMarkerStyle(24);
    fv24_cent1020->SetMarkerColor(kCyan+2);
    fv24_cent1020->SetMarkerSize(1.);
    fv24_cent1020->SetLineColor(kCyan+2);
    fv24_cent1020->SetLineWidth(2);

    
    //20-30
    Double_t v24Run2_cent2030[20]={0.0373426, 0.0600879, 0.0777636, 0.0951139, 0.112949, 0.130113, 0.139457, 0.16252, 0.173051, 0.187569, 0.191713, 0.192826, 0.17515};
    Double_t v24Run2_cent2030_Err[20]={0.000526843, 0.00071135, 0.00102737, 0.000999068, 0.000947932, 0.00113189, 0.00171729, 0.0018113, 0.0024593, 0.00414407, 0.00501406, 0.00361673, 0.0118352};
    Double_t v24Run2_cent2030_Sys[20]={0.00115762, 0.00186273, 0.00241067, 0.00294853, 0.00112949, 0.00130113, 0.00139457, 0.0016252, 0.00173051, 0.00187569, 0.00191713, 0.00192826, 0.0017515};
    Double_t v24Run2_cent2030CombErr[20]={0.00127187, 0.00199393, 0.00262046, 0.00311319, 0.00147456, 0.00172456, 0.00221222, 0.00243353, 0.00300713, 0.00454879, 0.00536807, 0.00409865, 0.0119641};
    
    TGraphErrors *fv24_cent2030 = new TGraphErrors(13,PtBin,v24Run2_cent2030,PtBinErr,v24Run2_cent2030CombErr);
    fv24_cent2030->SetMarkerStyle(25);
    fv24_cent2030->SetMarkerColor(kMagenta+1);
    fv24_cent2030->SetMarkerSize(1.);
    fv24_cent2030->SetLineColor(kMagenta+1);
    fv24_cent2030->SetLineWidth(2);
    
    //30-40
    Double_t v24Run2_cent3040[20]={0.0402175, 0.0661117, 0.0895367, 0.110332, 0.128524, 0.145721, 0.166351, 0.179628, 0.193559, 0.191064, 0.20225, 0.18186, 0.165222};
    Double_t v24Run2_cent3040_Err[20]={0.000734471, 0.00105856, 0.00150956, 0.00114635, 0.00204589, 0.00262661, 0.0023523, 0.00263422, 0.00276848, 0.00598472, 0.00846769, 0.0123023, 0.0125843};
    Double_t v24Run2_cent3040_Sys[20]={0.00124674, 0.00204946, 0.00277564, 0.00342031, 0.00115672, 0.00131149, 0.00149716, 0.00161665, 0.00174203, 0.00171958, 0.00182025, 0.00163674, 0.001487};
    Double_t v24Run2_cent3040CombErr[20]={0.001447, 0.0023067, 0.00315958, 0.0036073, 0.00235025, 0.00293583, 0.00278833, 0.00309074, 0.00327096, 0.00622686, 0.00866112, 0.0124107, 0.0126719};
    
    TGraphErrors *fv24_cent3040 = new TGraphErrors(13,PtBin,v24Run2_cent3040,PtBinErr,v24Run2_cent3040CombErr);
    fv24_cent3040->SetMarkerStyle(27);
    fv24_cent3040->SetMarkerColor(1);
    fv24_cent3040->SetMarkerSize(1.5);
    fv24_cent3040->SetLineColor(1);
    fv24_cent3040->SetLineWidth(2);
    
  
    //===================================================================================================================
    // Figure a:
    //===================================================================================================================
    //  v2{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v2_etaGap10 = 17;
    Double_t xSP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.900000,1.100000,
        1.300000,1.500000,1.700000,1.900000,2.200000,2.600000,3.000000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.013672,0.015745,0.019944,0.024169,0.026921,0.030296,0.034916,0.038472,
        0.043103,0.046626,0.052386,0.057132,0.060070,0.068419,0.061459,0.056816,0.050311};
    Double_t xErrSP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v2_etaGap10[nPointsSP_0005ALICE_v2_etaGap10] = {0.000813,0.000804,0.000852,0.000930,0.001029,0.001140,0.000939,0.001155,
        0.001410,0.001736,0.002131,0.002620,0.002502,0.003759,0.005615,0.006617,0.014242};
    TGraphErrors *GrSP_0005ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v2_etaGap10,xSP_0005ALICE_v2_etaGap10,ySP_0005ALICE_v2_etaGap10,xErrSP_0005ALICE_v2_etaGap10,yErrSP_0005ALICE_v2_etaGap10);
    GrSP_0005ALICE_v2_etaGap10->SetLineColor(kRed-6);
    GrSP_0005ALICE_v2_etaGap10->SetLineWidth(1);
    GrSP_0005ALICE_v2_etaGap10->SetFillColor(kRed-6);
    GrSP_0005ALICE_v2_etaGap10->SetFillStyle(3002);
    
    //  v3{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v3_etaGap10 = 16;
    Double_t xSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.900000,1.100000,
        1.300000,1.500000,1.700000,1.900000,2.300000,2.900000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.006303,0.009800,0.011143,0.014246,0.017628,0.019437,0.028412,0.030580,
        0.038730,0.045653,0.052469,0.062303,0.071522,0.082223,0.083373,0.076951};
    Double_t xErrSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.001012,0.000996,0.001062,0.001158,0.001277,0.001415,0.001158,0.001422,
        0.001734,0.002124,0.002610,0.003206,0.002724,0.004977,0.008123,0.017180};
    TGraphErrors *GrSP_0005ALICE_v3_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v3_etaGap10,xSP_0005ALICE_v3_etaGap10,ySP_0005ALICE_v3_etaGap10,xErrSP_0005ALICE_v3_etaGap10,yErrSP_0005ALICE_v3_etaGap10);
    GrSP_0005ALICE_v3_etaGap10->SetLineColor(kBlue-4);
    GrSP_0005ALICE_v3_etaGap10->SetLineWidth(1);
    GrSP_0005ALICE_v3_etaGap10->SetFillColor(kBlue-4);
    GrSP_0005ALICE_v3_etaGap10->SetFillStyle(3002);
   
    
    //  v4{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v4_etaGap10 = 11;
    Double_t xSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.300000,0.500000,0.700000,0.950000,1.250000,1.550000,1.850000,2.300000,2.900000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.002042,0.002556,0.009693,0.013286,0.016780,0.027865,0.031797,0.051101,0.060164,
        0.095985,0.094607};
    Double_t xErrSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.001460,0.001624,0.001930,0.002021,0.002737,0.003717,0.005042,0.005564,0.010160,
        0.016472,0.035083};
    TGraphErrors *GrSP_0005ALICE_v4_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v4_etaGap10,xSP_0005ALICE_v4_etaGap10,ySP_0005ALICE_v4_etaGap10, xErrSP_0005ALICE_v4_etaGap10,yErrSP_0005ALICE_v4_etaGap10);
    GrSP_0005ALICE_v4_etaGap10->SetLineColor(kGreen-3);
    GrSP_0005ALICE_v4_etaGap10->SetLineWidth(1);
    GrSP_0005ALICE_v4_etaGap10->SetFillColor(kGreen-3);
    GrSP_0005ALICE_v4_etaGap10->SetFillStyle(3002);

   
    
    //  v2{SP}(pt) for 30-40%, rapidity gap = 1.0:
    const Int_t nPointsSP_3040ALICE_v2_etaGap10 = 22;
    Double_t xSP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,
        0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.300000,2.500000,2.800000,3.200000,
        3.600000,4.000000,4.400000,4.800000};
    Double_t ySP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.039107,0.053955,0.068788,0.083037,0.095749,
        0.107080,0.119547,0.127823,0.143002,0.160203,0.180128,0.191935,0.201014,0.210955,0.220726,0.227580,0.244640,0.230650,
        0.234024,0.220566,0.217360,0.208241};
    Double_t xErrSP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.000746,0.000744,0.000790,0.000867,0.000960,
        0.001069,0.001187,0.001321,0.001091,0.001339,0.001642,0.002016,0.002483,0.003050,0.003736,0.004558,0.004267,0.006108,
        0.008517,0.011491,0.017248,0.017496};
    TGraphErrors *GrSP_3040ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_3040ALICE_v2_etaGap10,xSP_3040ALICE_v2_etaGap10,ySP_3040ALICE_v2_etaGap10,xErrSP_3040ALICE_v2_etaGap10, yErrSP_3040ALICE_v2_etaGap10);
    GrSP_3040ALICE_v2_etaGap10->SetLineColor(kRed-6);
    GrSP_3040ALICE_v2_etaGap10->SetLineWidth(1);
    GrSP_3040ALICE_v2_etaGap10->SetFillColor(kRed-6);
    GrSP_3040ALICE_v2_etaGap10->SetFillStyle(3013);
    
    //  v3{SP}(pt) for 30-40%, rapidity gap = 1.0:
    const Int_t nPointsSP_3040ALICE_v3_etaGap10 = 18;
    Double_t xSP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,
        0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.009546,0.013907,0.013673,0.021384,0.029748,
        0.035769,0.037629,0.044245,0.053656,0.061235,0.063584,0.080152,0.089683,0.089203,0.094393,0.115172,0.107330,0.062049};
    Double_t xErrSP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.};
    Double_t yErrSP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.001753,0.001759,0.001898,0.002071,0.002291,
        0.002545,0.002838,0.003136,0.002585,0.003170,0.003875,0.004736,0.005828,0.007141,0.006755,0.009900,0.010855,0.023521};
    TGraphErrors *GrSP_3040ALICE_v3_etaGap10 = new TGraphErrors(nPointsSP_3040ALICE_v3_etaGap10,xSP_3040ALICE_v3_etaGap10,ySP_3040ALICE_v3_etaGap10,xErrSP_3040ALICE_v3_etaGap10,yErrSP_3040ALICE_v3_etaGap10);
    GrSP_3040ALICE_v3_etaGap10->SetLineColor(kBlue-4);
    GrSP_3040ALICE_v3_etaGap10->SetLineWidth(1);
    GrSP_3040ALICE_v3_etaGap10->SetFillColor(kBlue-4);
    GrSP_3040ALICE_v3_etaGap10->SetFillStyle(3013);
    
    
    //  v4{SP}(pt) for 30-40%, rapidity gap = 1.0:
    const Int_t nPointsSP_3040ALICE_v4_etaGap10 = 12;
    Double_t xSP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.150000,1.450000,1.750000,2.050000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.002125,0.005223,0.017491,0.020758,0.026932,0.043126,0.055596,0.055244,0.081049,0.097779,0.074055,0.158711};
    Double_t xErrSP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.};
    Double_t yErrSP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.002391,0.002677,0.003227,0.003998,0.004216,0.005692,0.007729,0.010472,0.012836,0.018727,0.020749,0.044065};
    TGraphErrors *GrSP_3040ALICE_v4_etaGap10 = new TGraphErrors(nPointsSP_3040ALICE_v4_etaGap10,xSP_3040ALICE_v4_etaGap10,ySP_3040ALICE_v4_etaGap10,xErrSP_3040ALICE_v4_etaGap10,yErrSP_3040ALICE_v4_etaGap10);
    GrSP_3040ALICE_v4_etaGap10->SetLineColor(kGreen-3);
    GrSP_3040ALICE_v4_etaGap10->SetLineWidth(1);
    GrSP_3040ALICE_v4_etaGap10->SetFillColor(kGreen-3);
    GrSP_3040ALICE_v4_etaGap10->SetFillStyle(3013);

    

    // Systematical errors reported separatly for pt < 2 GeV and pt > 2 GeV
    // (see also bottom of this macro for the details!)
    // pt < 2 GeV:
    Double_t a1 = 0.01; // centrality estimators, pt < 2 GeV
    Double_t a2 = 0.01; // polarities, pt < 2 GeV
    Double_t a3 = 0.01; // # of TPC clusters, pt < 2 GeV
    Double_t a4 = 0.01; // DCAxy for POIs, pt < 2 GeV
    Double_t a5 = 0.01; // DCAz for POIs, pt < 2 GeV
    Double_t a6 = 0.02; // ESD global, pt < 2 GeV
    Double_t a7 = 0.04; // TPC tracks constrained to SPD vertex, pt < 2 GeV
    Double_t sysErrorRelativeCombinedPtSmallerThen2GeV = pow(pow(a1,2.)+pow(a2,2.)+pow(a3,2.)+pow(a4,2.)
                                                             +pow(a5,2.)+pow(a6,2.)+pow(a7,2.),0.5); // = 5 %
    // pt > 2 GeV:
    Double_t b1 = 0.01; // centrality estimators, pt > 2 GeV
    Double_t b2 = 0.00; // polarities, pt > 2 GeV
    Double_t b3 = 0.01; // # of TPC clusters, pt > 2 GeV
    Double_t b4 = 0.00; // DCAxy for POIs, pt > 2 GeV
    Double_t b5 = 0.00; // DCAz for POIs, pt > 2 GeV
    Double_t b6 = 0.00; // ESD global, pt > 2 GeV
    Double_t b7 = 0.03; // TPC tracks constrained to SPD vertex, pt > 2 GeV
    Double_t sysErrorRelativeCombinedPtLargerThen2GeV = pow(pow(b1,2.)+pow(b2,2.)+pow(b3,2.)+pow(b4,2.) +pow(b5,2.)+pow(b6,2.)+pow(b7,2.),0.5); // = 3.3 %
    

    // v2{4}(pt) for 10-20%:
    Double_t xQC4_1020[] = {0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
        1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
        2.100000,2.300000,2.500000,2.700000,2.900000,3.100000,3.300000,3.500000,3.700000,3.900000,
        4.150000,4.450000,4.750000,5.100000,5.500000,6.000000,6.600000,7.300000,8.200000,9.200000,
        10.200000,11.200000,12.350000,14.000000,17.500000,25.000000,35.000000,45.000000,55.000000,65.000000,
        75.000000,85.000000,0.000000};
    Double_t yQC4_1020[] = {0.017774,0.023225,0.031787,0.040604,0.048850,0.056175,0.063071,0.069683,0.075780,
        0.081264,0.086234,0.091964,0.096743,0.102048,0.107092,0.111773,0.116437,0.121063,0.124777,
        0.130815,0.138226,0.143429,0.148270,0.151187,0.152727,0.152423,0.154532,0.150250,0.145989,
        0.143389,0.140081,0.125633,0.117222,0.106690,0.098739,0.091938,0.077646,0.069070,0.066492,
        0.059181,0.064943,0.060974,0.041624,0.009953,0.031629,-0.352207,-0.011664,-0.024524,0.412278,
        1.105233,0.121755,0.000000};
    Double_t xErrQC4_1020[53] = {0.};
    Double_t yErrQC4_1020[] = {0.000223,0.000165,0.000186,0.000216,0.000245,0.000277,0.000306,0.000338,0.000367,
        0.000397,0.000429,0.000462,0.000496,0.000538,0.000580,0.000625,0.000679,0.000742,0.000811,
        0.000662,0.000814,0.000999,0.001204,0.001419,0.001641,0.001876,0.002129,0.002412,0.002731,
        0.002594,0.003127,0.003728,0.003946,0.004898,0.005079,0.006499,0.007411,0.008809,0.011717,
        0.014599,0.017963,0.020450,0.021298,0.021012,0.027724,0.031261,0.068300,0.131392,0.201637,
        0.293752,1.056439,0.000000}; // Statistical errors only
    // Combining in quadrature statistical and systematical errors:
    Double_t yErrQC4_1020_Combined[53] = {0.};
    for(Int_t p=0;p<53;p++)
    {
        Double_t relativeSysError_1020 = 0.;
        if(xQC4_1020[p]<2.0)
        {
            relativeSysError_1020 = sysErrorRelativeCombinedPtSmallerThen2GeV;
        } else
        {
            relativeSysError_1020 = sysErrorRelativeCombinedPtLargerThen2GeV;
        }
        yErrQC4_1020_Combined[p] = pow(pow(relativeSysError_1020*yQC4_1020[p],2.)+pow(yErrQC4_1020[p],2.),0.5);
    } // end of for(Int_t p=0;p<53;p++)
    Int_t nPointsQC4_1020 = sizeof(xQC4_1020)/sizeof(Double_t);
    TGraphErrors *QC4_1020 = new TGraphErrors(31,xQC4_1020,yQC4_1020,xErrQC4_1020,yErrQC4_1020_Combined);
    QC4_1020->SetMarkerStyle(kFullSquare);
    QC4_1020->SetMarkerColor(kCyan+2);
    QC4_1020->SetLineColor(kCyan+2);
    QC4_1020->SetLineWidth(1);
    QC4_1020->SetFillColor(kCyan+2);
    QC4_1020->SetFillStyle(3005);
    
    //===================================================================================================================
    
    //===================================================================================================================
    // v2{4}(pt) for 20-30%:
    Double_t xQC4_2030[] = {0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
        1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
        2.100000,2.300000,2.500000,2.700000,2.900000,3.150000,3.450000,3.750000,4.050000,4.350000,
        4.650000,5.050000,5.550000,6.300000,7.400000,9.000000,11.500000,14.750000,18.250000,25.000000,
        35.000000,45.000000,55.000000,65.000000,75.000000,0.000000,0.000000};
    Double_t yQC4_2030[] = {0.023211,0.030399,0.041714,0.053495,0.064452,0.074603,0.083592,0.092134,0.100406,
        0.107896,0.115173,0.122144,0.128718,0.135064,0.141218,0.148384,0.153673,0.158816,0.164991,
        0.170478,0.179390,0.184598,0.188347,0.193158,0.191820,0.188337,0.185207,0.175270,0.165661,
        0.160623,0.146838,0.126817,0.113601,0.109363,0.095271,0.070494,0.079178,0.103382,0.053497,
        0.014589,0.604099,0.653579,-0.075916,-1.194077,0.000000,0.000000};
    Double_t xErrQC4_2030[47] = {0.};
    Double_t yErrQC4_2030[] = {0.000224,0.000165,0.000182,0.000207,0.000235,0.000265,0.000293,0.000324,0.000354,
        0.000385,0.000415,0.000451,0.000490,0.000531,0.000577,0.000633,0.000699,0.000774,0.000859,
        0.000711,0.000877,0.001065,0.001255,0.001455,0.001394,0.001682,0.002010,0.002434,0.002915,
        0.003497,0.003359,0.004392,0.004284,0.005854,0.007387,0.010793,0.017269,0.028365,0.028359,
        0.044697,0.071455,0.228457,0.245604,2.378602,0.000000,0.000000}; // Statistical errors only
    // Combining in quadrature statistical and systematical errors:
    Double_t yErrQC4_2030_Combined[47] = {0.};
    for(Int_t p=0;p<47;p++)
    {
        Double_t relativeSysError_2030 = 0.;
        if(xQC4_2030[p]<2.0)
        {
            relativeSysError_2030 = sysErrorRelativeCombinedPtSmallerThen2GeV;
        } else
        {
            relativeSysError_2030 = sysErrorRelativeCombinedPtLargerThen2GeV;
        }
        yErrQC4_2030_Combined[p] = pow(pow(relativeSysError_2030*yQC4_2030[p],2.)+pow(yErrQC4_2030[p],2.),0.5);
    } // end of for(Int_t p=0;p<47;p++)
    Int_t nPointsQC4_2030 = sizeof(xQC4_2030)/sizeof(Double_t);
    TGraphErrors *QC4_2030 = new TGraphErrors(30,xQC4_2030,yQC4_2030,xErrQC4_2030,yErrQC4_2030_Combined);
    QC4_2030->SetMarkerStyle(kFullSquare);
    QC4_2030->SetMarkerColor(kMagenta-9);
    QC4_2030->SetLineColor(kMagenta-9);
    QC4_2030->SetLineWidth(1);
    QC4_2030->SetFillColor(kMagenta-4);
    QC4_2030->SetFillStyle(3012);

    //===================================================================================================================
    
    //===================================================================================================================
    // v2{4}(pt) for 30-40%:
    Double_t xQC4_3040[] = {0.150000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
        1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
        2.050000,2.150000,2.250000,2.350000,2.450000,2.550000,2.700000,2.900000,3.100000,3.300000,
        3.500000,3.750000,4.050000,4.350000,4.700000,5.150000,5.650000,6.300000,7.100000,7.900000,
        8.800000,9.800000,11.150000,13.500000,17.500000,25.000000,35.000000,45.000000,55.000000,65.000000,
        75.000000,0.000000,0.000000};
    Double_t yQC4_3040[] = {0.025999,0.034005,0.047331,0.060976,0.073396,0.084921,0.095273,0.105345,0.115212,
        0.123443,0.131304,0.139384,0.146887,0.153697,0.161578,0.166773,0.173150,0.177996,0.183631,
        0.187558,0.191383,0.195856,0.199685,0.201283,0.202841,0.204668,0.205053,0.201223,0.199740,
        0.196353,0.192957,0.183932,0.172645,0.157642,0.146704,0.129871,0.126335,0.119367,0.103949,
        0.097185,0.094387,0.100591,0.116564,0.043284,0.280880,0.117779,0.223997,0.749208,-0.209985,
        0.919698,0.000000,0.000000};
    Double_t xErrQC4_3040[53] = {0.};
    Double_t yErrQC4_3040[] = {0.000296,0.000213,0.000230,0.000260,0.000291,0.000327,0.000362,0.000402,0.000442,
        0.000479,0.000527,0.000576,0.000631,0.000696,0.000767,0.000850,0.000951,0.001066,0.001188,
        0.001327,0.001470,0.001619,0.001773,0.001920,0.002072,0.001626,0.001861,0.002125,0.002396,
        0.002722,0.002603,0.003156,0.003788,0.004041,0.004586,0.005884,0.006230,0.008552,0.011334,
        0.013120,0.017022,0.017409,0.021085,0.026571,0.017655,0.036215,0.070872,0.415189,0.925482,
        0.780913,0.000000,0.000000}; // Statistical errors only
    // Combining in quadrature statistical and systematical errors:
    Double_t yErrQC4_3040_Combined[53] = {0.};
    for(Int_t p=0;p<53;p++)
    {
        Double_t relativeSysError_3040 = 0.;
        if(xQC4_3040[p]<2.0)
        {
            relativeSysError_3040 = sysErrorRelativeCombinedPtSmallerThen2GeV;
        } else
        {
            relativeSysError_3040 = sysErrorRelativeCombinedPtLargerThen2GeV;
        }
        yErrQC4_3040_Combined[p] = pow(pow(relativeSysError_3040*yQC4_3040[p],2.)+pow(yErrQC4_3040[p],2.),0.5); 
    } // end of for(Int_t p=0;p<53;p++) 
    Int_t nPointsQC4_3040 = sizeof(xQC4_3040)/sizeof(Double_t);                                      
    TGraphErrors *QC4_3040 = new TGraphErrors(34,xQC4_3040,yQC4_3040,xErrQC4_3040,yErrQC4_3040_Combined);
    QC4_3040->SetMarkerStyle(kFullSquare);
    QC4_3040->SetMarkerColor(kGray+1);
    QC4_3040->SetLineColor(kGray+1);
    QC4_3040->SetLineWidth(1);
    QC4_3040->SetFillColor(kGray+1);
    QC4_3040->SetFillStyle(3008);

    
    Double_t PtBin[] = {0.3, 0.5, 0.7, 0.9, 1.125, 1.375, 1.625, 1.875, 2.25, 2.75, 3.25, 3.75, 4.5, 5.5};
    Double_t PtBinErr[] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05,0.05, 0.05, 0.05,0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    Double_t PtBinSys[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    
    
    // /Users/youzhou/Research/ALICE/Flow_ch/Code_ReSampling_RUN2/Run_sys_Jan12/Output/244918/Fig_sys
    Double_t RatioV24_cc2030[]= {1.03666, 1.02539, 0.987332, 0.991868, 1.00163, 1.00241, 0.95637, 1.0162, 0.98567, 0.991899, 1.00134, 1.04638, 1.06515};
    Double_t RatioV24Err_cc2030[]={0.0146759, 0.0122016, 0.0131083, 0.0104593, 0.00847948, 0.00879668, 0.0118535, 0.0115006, 0.0141548, 0.022022, 0.026447, 0.020522, 0.0731219};
    Double_t RatioV24Sys_cc2030[]={0.0103666, 0.0102539, 0.00987332, 0.00991868, 0.0100163, 0.0100241, 0.0095637, 0.010162, 0.0098567, 0.00991899, 0.0100134, 0.0104638, 0.0106515};
    Double_t RatioV24CombErr_cc2030[]={0.017968, 0.0159381, 0.0164107, 0.0144145, 0.0131236, 0.0133366, 0.0152306, 0.015347, 0.0172486, 0.0241527, 0.0282792, 0.0230357, 0.0738936};
    
    TGraphErrors *fRatioV24_cc2030 = new TGraphErrors(13,PtBin,RatioV24_cc2030,PtBinErr,RatioV24Err_cc2030);
    fRatioV24_cc2030->SetMarkerStyle(25);
    fRatioV24_cc2030->SetMarkerColor(1);
    fRatioV24_cc2030->SetLineColor(1);
    fRatioV24_cc2030->SetLineWidth(2);
    fRatioV24_cc2030->SetMarkerSize(0.8);
    
    TGraphErrors *fRatioV24sys_cc2030 = new TGraphErrors(13,PtBin,RatioV24_cc2030,PtBinSys,RatioV24Sys_cc2030);
    fRatioV24sys_cc2030->SetMarkerStyle(25);
    fRatioV24sys_cc2030->SetMarkerColor(1);
    fRatioV24sys_cc2030->SetLineColor(1);
    fRatioV24sys_cc2030->SetLineWidth(2);
    fRatioV24sys_cc2030->SetMarkerSize(0.8);
    fRatioV24sys_cc2030->SetFillStyle(1001);
    fRatioV24sys_cc2030->SetFillColor(kGray);
    
    TGraphErrors *fRatioV24Comb_cc2030 = new TGraphErrors(13,PtBin,RatioV24_cc2030,PtBinErr,RatioV24CombErr_cc2030);
    fRatioV24Comb_cc2030->SetMarkerStyle(25);
    fRatioV24Comb_cc2030->SetMarkerColor(kMagenta+1);
    fRatioV24Comb_cc2030->SetLineColor(kMagenta+1);
    fRatioV24Comb_cc2030->SetLineWidth(2);
    fRatioV24Comb_cc2030->SetMarkerSize(0.8);
    fRatioV24Comb_cc2030->SetFillStyle(1001);
    fRatioV24Comb_cc2030->SetFillColor(kGray);
    
    
    
    
    // Durham data from here %%%%%%%%%%%%%%%%%%%%%%%%
    cout << "Centrality: 0-5% (RUN2)"<<endl;
    printf("\p_{t}  v_{2}{2,|#Delta#eta|>1}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v22Gap10Run2_cent05[p],v22Gap10Run2_cent05_Err[p], v22Gap10Run2_cent05_Sys[p]);
    }
    cout<<endl;
    
    cout << "Centrality: 30-40% (RUN2)"<<endl;
    printf("\p_{t}  v_{2}{2,|#Delta#eta|>1}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v22Gap10Run2_cent3040[p],v22Gap10Run2_cent3040_Err[p], v22Gap10Run2_cent3040_Sys[p]);
    }
    cout<<endl;
    
    cout << "Centrality: 0-5% (RUN2)"<<endl;
    printf("\p_{t}  v_{3}{2,|#Delta#eta|>1}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v32Gap10Run2_cent05[p],v32Gap10Run2_cent05_Err[p], v32Gap10Run2_cent05_Sys[p]);
    }
    cout<<endl;
    
    cout << "Centrality: 30-40% (RUN2)"<<endl;
    printf("\p_{t}  v_{3}{2,|#Delta#eta|>1}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v32Gap10Run2_cent3040[p],v32Gap10Run2_cent3040_Err[p], v32Gap10Run2_cent3040_Sys[p]);
    }
    cout<<endl;

    cout << "Centrality: 0-5% (RUN2)"<<endl;
    printf("\p_{t}  v_{4}{2,|#Delta#eta|>1}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v42Gap10Run2_cent05[p],v42Gap10Run2_cent05_Err[p], v42Gap10Run2_cent05_Sys[p]);
    }
    cout<<endl;
    
    cout << "Centrality: 30-40% (RUN2)"<<endl;
    printf("\p_{t}  v_{4}{2,|#Delta#eta|>1}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v42Gap10Run2_cent3040[p],v42Gap10Run2_cent3040_Err[p], v42Gap10Run2_cent3040_Sys[p]);
    }
    cout<<endl;
    
    cout << "Centrality: 10-20% (RUN2)"<<endl;
    printf("\p_{t}  v_{2}{4}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v24Run2_cent1020[p],v24Run2_cent1020_Err[p], v24Run2_cent1020_Sys[p]);
    }
    cout<<endl;
    
    cout << "Centrality: 20-30% (RUN2)"<<endl;
    printf("\p_{t}  v_{2}{4}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v24Run2_cent2030[p],v24Run2_cent2030_Err[p], v24Run2_cent2030_Sys[p]);
    }
    cout<<endl;

    cout << "Centrality: 30-40% (RUN2)"<<endl;
    printf("\p_{t}  v_{2}{4}     +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],v24Run2_cent3040[p],v24Run2_cent3040_Err[p], v24Run2_cent3040_Sys[p]);
    }
    cout<<endl;
    
    
    cout << "Centrality: 20-30%"<<endl;
    printf("\p_{t}  ratio:v_{2}{4}(RUN2/RUN1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<13;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",PtBin[p],RatioV24_cc2030[p],RatioV24Err_cc2030[p], RatioV24Sys_cc2030[p]);
    }
    cout<<endl;
    
    // Style histogramm:
    TH1F *myBlankHisto = new TH1F("myBlankHisto","Blank Histogram",100,0.,10.);
    //myBlankHisto->SetXTitle("p_{t} (GeV/#font[72]{c})");
    myBlankHisto->SetYTitle("#it{v}_{n}{2, |#Delta #eta| > 1.0}  ");
    myBlankHisto->SetNdivisions(505,"y");
    myBlankHisto->GetYaxis()->SetTitleOffset(0.85);
    myBlankHisto->GetYaxis()->SetLabelSize(0.08);
    myBlankHisto->GetYaxis()->SetTitleSize(0.08);
    //myBlankHisto->GetYaxis()->SetRangeUser(-0.02101,0.399);
    myBlankHisto->GetXaxis()->SetRangeUser(0.,5.1);
    
    // Main canvas:
    TCanvas *myCan = new TCanvas("myCan",cStamp1,400,500);
    myCan->Divide(1,2,0.,0.);
    myCan->SetTopMargin(0.0002);
    myCan->Draw();
    
    //=================================================================================================
    
    // Figure 3a:
    myCan->cd(1);
    gPad->SetTopMargin(0.001);
    gPad->SetBottomMargin(0.);
    gPad->SetRightMargin(0.001);
    gPad->SetLeftMargin(0.17);
    TH1D* styleHistPad1 = myBlankHisto->Clone("styleHistPad1");
    styleHistPad1->GetYaxis()->SetRangeUser(-0.01,0.27);
    styleHistPad1->GetYaxis()->SetLabelSize(0.085);
    styleHistPad1->GetYaxis()->SetTitleSize(0.09);
    styleHistPad1->GetYaxis()->SetTitleOffset(0.9);
    styleHistPad1->GetXaxis()->SetNdivisions(10);
    styleHistPad1->Draw();
    
    styleHistPad1->Draw();
    
    GrSP_0005ALICE_v2_etaGap10->Draw("3same");
    GrSP_0005ALICE_v3_etaGap10->Draw("3same");
    GrSP_0005ALICE_v4_etaGap10->Draw("3same");
    
    fv22Gap10_cent05->Draw("pzsame");
    fv32Gap10_cent05->Draw("pzsame");
    fv42Gap10_cent05->Draw("pzsame");
    
    
    // Common legend:
    TLegend *myLegendPad1 = new TLegend(0.22,0.55,0.52,0.9,"5.02 TeV");
    //myLegendSetUp(myLegendPad1,0.07);
    myLegendPad1->SetTextSize(0.07);
    myLegendPad1->SetLineColor(0);
    myLegendPad1->AddEntry(fv22Gap10_cent05,"v_{2}{2, |#Delta#eta|>1}","p");
    myLegendPad1->AddEntry(fv32Gap10_cent05,"v_{3}{2, |#Delta#eta|>1}","p");
    myLegendPad1->AddEntry(fv42Gap10_cent05,"v_{4}{2, |#Delta#eta|>1}","p");
    myLegendPad1->SetMargin(0.2);
    myLegendPad1->SetHeader("5.02 TeV");
    myLegendPad1->Draw("same");
    
    
    TLegend *myLegendPad11 = new TLegend(0.58,0.55,0.94,0.9,"2.76 TeV");
    myLegendPad11->SetTextSize(0.07);
    myLegendPad11->SetLineColor(0);
    myLegendPad11->AddEntry(GrSP_0005ALICE_v2_etaGap10,"v_{2}{2, |#Delta#eta|>1}","f");
    myLegendPad11->AddEntry(GrSP_0005ALICE_v3_etaGap10,"v_{3}{2, |#Delta#eta|>1}","f");
    myLegendPad11->AddEntry(GrSP_0005ALICE_v4_etaGap10,"v_{4}{2, |#Delta#eta|>1}","f");
    myLegendPad11->SetMargin(0.2);
    myLegendPad11->SetHeader("2.76 TeV");
    myLegendPad11->Draw("same");
    
    TLatex *aa = new TLatex(0.22,0.92,"ALICE Pb-Pb");
    aa->SetNDC();
    aa->SetTextSize(0.07);
    aa->Draw();
    
    TLatex *a = new TLatex(0.81,0.92," 0-5% (a)");
    a->SetNDC();
    a->SetTextSize(0.08);
    a->Draw();
    
    //    TLatex *a11 = new TLatex(0.84,0.07,"0-5%");
    //    a11->SetNDC();
    //    a11->SetTextSize(0.07);
    //    a11->Draw();
    
    
    
    //=================================================================================================
    
    // Figure 3b:
    myCan->cd(2);
    gPad->SetTopMargin(0.);
    gPad->SetBottomMargin(0.2);
    gPad->SetRightMargin(0.001);
    gPad->SetLeftMargin(0.17);
    TH1D* styleHistPad2 = myBlankHisto->Clone("styleHistPad2");
    styleHistPad2->GetYaxis()->SetRangeUser(-0.01,0.29);
    styleHistPad2->GetYaxis()->SetLabelSize(0.08);
    styleHistPad2->GetYaxis()->SetTitleSize(0.08);
    styleHistPad2->GetYaxis()->SetTitleOffset(1.);
    styleHistPad2->GetXaxis()->SetLabelSize(0.08);
    styleHistPad2->GetXaxis()->SetTitleSize(0.08);
    styleHistPad2->GetXaxis()->SetNdivisions(10);
    styleHistPad2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    styleHistPad2->Draw();
    TLatex *b = new TLatex(0.77,0.93,"30-40% (b)");
    b->SetNDC();
    b->SetTextSize(0.07);
    b->Draw();
    
    //    TLatex *a12 = new TLatex(0.84,0.05,"30-40%");
    //    a12->SetNDC();
    //    a12->SetTextSize(0.07);
    //    a12->Draw();
    
    TLegend *myLegendPad2 = new TLegend(0.12,0.8,0.45,0.98,"Centrality 20-30%");
    myLegendSetUp(myLegendPad2,0.07);
    
    
    GrSP_3040ALICE_v2_etaGap10->Draw("3same");
    GrSP_3040ALICE_v3_etaGap10->Draw("3same");
    GrSP_3040ALICE_v4_etaGap10->Draw("3same");
    
    
    fv22Gap10_cent3040->Draw("pzsame");
    fv32Gap10_cent3040->Draw("pzsame");
    fv42Gap10_cent3040->Draw("pzsame");
    
    //myCan->SaveAs("fig2_20160129_1.eps");
    
    
    //=================================================================================================
    
    TCanvas *myCan2 = new TCanvas("myCan2",cStamp2,400,500);
    myCan2->Divide(1,2,0.,0.);
    myCan2->SetTopMargin(0.0002);
    myCan2->Draw();
    
    myCan2->cd(1);
    gPad->SetTopMargin(0.001);
    gPad->SetBottomMargin(0.);
    gPad->SetRightMargin(0.001);
    gPad->SetLeftMargin(0.17);
    TH1D* styleHistPad3 = myBlankHisto->Clone("styleHistPad3");
    styleHistPad3->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    styleHistPad3->SetYTitle("#it{v}_{2}{4} ");
    styleHistPad3->GetYaxis()->SetRangeUser(0.001,0.42);
    styleHistPad3->GetYaxis()->SetLabelSize(0.085);
    styleHistPad3->GetYaxis()->SetTitleSize(0.09);
    styleHistPad3->GetYaxis()->SetTitleOffset(0.9);
    styleHistPad3->GetXaxis()->SetNdivisions(10);
    styleHistPad3->Draw();
    TLatex *c = new TLatex(0.93,0.92,"(c)");
    c->SetNDC();
    c->SetTextSize(0.08);
    c->Draw("same");
    TLegend *legend = new TLegend(0.5,0.54,0.8,0.9, "2.76 TeV");
    legend->SetFillStyle(0); // white legend background
    legend->SetLineColor(0);
    legend->SetTextSize(gStyle->GetTextSize()*1.4);
    legend->AddEntry(QC4_1020,"10-20 %","f");
    legend->AddEntry(QC4_2030,"20-30 %","f");
    legend->AddEntry(QC4_3040,"30-40 %","f");
    legend->Draw("same");
    
    TLegend *legend1 = new TLegend(0.22,0.54,0.45,0.9, "5.02 TeV");
    legend1->SetFillStyle(0); // white legend background
    legend1->SetLineColor(0);
    legend1->SetTextSize(gStyle->GetTextSize()*1.4);
    legend1->AddEntry(fv24_cent1020,"10-20 %","p");
    legend1->AddEntry(fv24_cent2030,"20-30 %","p");
    legend1->AddEntry(fv24_cent3040,"30-40 %","p");
    legend1->Draw("same");
    
    
    QC4_1020->Draw("3same");
    QC4_2030->Draw("3same");
    QC4_3040->Draw("3same");
    
    
    fv24_cent1020->Draw("pzsame");
    fv24_cent2030->Draw("pzsame");
    fv24_cent3040->Draw("pzsame");
    
    aa->Draw();
    
    
    
    // Figure 3c:
    myCan2->cd(2);
    gPad->SetTopMargin(0.);
    gPad->SetBottomMargin(0.2);
    gPad->SetRightMargin(0.001);
    gPad->SetLeftMargin(0.17);
    TH1D* styleHistPad4 = myBlankHisto->Clone("styleHistPad3");
    styleHistPad4->SetXTitle("#it{p}_{T} (GeV/#it{c}) ");
    styleHistPad4->SetYTitle("Ratio  ");
    styleHistPad4->GetYaxis()->SetRangeUser(0.86,1.14);
    styleHistPad4->GetYaxis()->SetLabelSize(0.08);
    styleHistPad4->GetYaxis()->SetTitleSize(0.08);
    styleHistPad4->GetYaxis()->SetTitleOffset(1.);
    styleHistPad4->GetXaxis()->SetLabelSize(0.08);
    styleHistPad4->GetXaxis()->SetTitleSize(0.08);
    styleHistPad4->GetXaxis()->SetNdivisions(10);
    styleHistPad4->Draw();
    TLatex *d = new TLatex(0.93,0.93,"(d)");
    d->SetNDC();
    d->SetTextSize(0.07);
    d->Draw("same");
    //  TPaveText *d = new TPaveText(0.9,0.9,0.99,0.99,"brNDC");
    //  d->SetBorderSize(0.1);                                     //no borders
    //  d->SetFillColor(0);
    //  d->AddText("(d)");
    //  d->Draw("same");
    
    
    TLegend *legendv2 = new TLegend(0.22,0.8,0.39,0.96, "#it{v}_{2 }{4} (5.02 TeV)  / #it{v}_{2}{4} (2.76 TeV)");
    legendv2->SetFillStyle(0); // white legend background
    legendv2->SetLineColor(0);
    legendv2->SetTextSize(gStyle->GetTextSize()*1.3);
    legendv2->AddEntry(fRatioV24Comb_cc2030,"20-30 %","p");
    legendv2->Draw("same");
    
    TLine *line1 = new TLine(0,1.,5.1,1.);
    line1->SetLineStyle(2);
    line1->SetLineColor(1);
    line1->Draw("lsame");
    
    //fRatioV24sys_cc2030->Draw("2same");
    //fRatioV24_cc2030->Draw("pzsame");
    fRatioV24Comb_cc2030->Draw("pzsame");
    
    //myCan2->SaveAs("fig2_20160129_2.eps");
    
    
    return;
 
    
    if (rPerformance){
        TLatex *alice = new TLatex(0.75,0.34,"Performance");
        alice->SetNDC();
        alice->SetTextColor(myDarkRed);
        //    alice->SetTextFont(42);
        alice->SetTextSize(0.05);
        alice->SetLineWidth(2);
        alice->Draw();
        
        TText *date = new TText(0.78,0.28,cStamp2);
        date->SetNDC();
        date->SetTextFont(42);
        date->SetTextSize(0.04);
        date->Draw();
        
        TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.77,0.37,0.92,0.60);
        myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
        myPadSetUp(myPadLogo,0,0,0,0);
        myPadLogo->Draw();
        myPadLogo->cd();
        TASImage *myAliceLogo = new TASImage("alice_logo_transparent.png");
        //.TASImage *myAliceLogo = new TASImage("alice_logo.pdf");
        
        myAliceLogo->Draw();
    }
    if (rWrite == 1)  myCan->SaveAs("fig_template.pdf");
    if (rWrite == 2)  myCan->SaveAs("fig_template.png");
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
    currentLegend->SetTextFont(42);
    currentLegend->SetBorderSize(0);
    currentLegend->SetFillStyle(0);
    currentLegend->SetFillColor(0);
    currentLegend->SetMargin(0.25);
    currentLegend->SetTextSize(currentTextSize);
    currentLegend->SetEntrySeparation(0.5);
    return;
}

void myPadSetUp(TPad *currentPad, float currentLeft=0.14, float currentTop=0.01, float currentRight=0.01, float currentBottom=0.16){
    currentPad->SetLeftMargin(currentLeft);
    currentPad->SetTopMargin(currentTop);
    currentPad->SetRightMargin(currentRight);
    currentPad->SetBottomMargin(currentBottom);
    return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
                  int currentMarkerStyle=21, int currentMarkerColor=0,
                  int currentLineStyle=1, int currentLineColor=0, int currentFillColor=0, int currentFillStyle=0)
{
    currentGraph->SetMarkerSize(currentMarkerSize);
    currentGraph->SetMarkerStyle(currentMarkerStyle);
    currentGraph->SetMarkerColor(currentMarkerColor);
    currentGraph->SetLineStyle(currentLineStyle);
    currentGraph->SetLineColor(currentLineColor);
    currentGraph->SetLineWidth(2.);
    currentGraph->SetFillColor(currentFillColor);
    currentGraph->SetFillStyle(currentFillStyle);
    return;
}

void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
 // Shift original TGraphErrors along x-axis by Double_t shift.
 
 if(!ge){cout<<"!ge"<<endl;exit(0);}
 
 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 Double_t y = 0.;
 for(Int_t p=0;p<nPoints;p++)
 { 
  ge->GetPoint(p,x,y);
  x+=shift;
  ge->SetPoint(p,x,y);
 } // end of for(Int_t p=0;p<nPoints;p++)
 
} // end of void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)

void myOptions(Int_t lStat=0){
    // Set gStyle
    int font = 42;
    // From plain
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    //
    gStyle->SetDrawBorder(0);
    gStyle->SetTextFont(font);
    gStyle->SetStatFont(font);
    gStyle->SetStatFontSize(0.06);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetTickLength(0.02,"y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.07,"xyz");
    gStyle->SetLabelFont(font,"xyz"); 
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetTitleFont(font,"xyz");  
    gStyle->SetTitleOffset(1.1,"xyz");
    gStyle->SetTitleSize(0.07,"xyz");
    gStyle->SetMarkerSize(1); 
    gStyle->SetPalette(1,0); 
    if (lStat){
        gStyle->SetOptTitle(1);
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(1111);
    }
    else {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
    }
}

TH1D* StyleHistogram()
{
    // Style histogram for top panel:
    
    Int_t nPtBins = 100;
    Double_t ptMin = 0.;
    Double_t ptMax = 10.0;
    TString xTitle  = "p_{t} (GeV/#font[72]{c})";
    TString yTitle  = "#it{v}_{n}  ";
    
    TH1D* hist = new TH1D("","",nPtBins,ptMin,ptMax);
    hist->SetXTitle(xTitle.Data());
    hist->SetYTitle(yTitle.Data());
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetTitleOffset(1.10);
    hist->SetMinimum(-0.02101); // minimum y-value
    hist->SetMaximum(0.3799); // maximum y-value 
    hist->GetXaxis()->SetRangeUser(0.,plotUpTo);  
    
    return hist;
    
} // end of TH1D* StyleHistogram()



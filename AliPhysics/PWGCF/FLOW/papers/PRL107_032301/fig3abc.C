// For the lists of output files see the bottom of this macro.

static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

void fig3abc(Int_t rWrite = 0, Int_t rPerformance = 0)
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
    
    //===================================================================================================================
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
    TGraphErrors *GrSP_3040ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_3040ALICE_v2_etaGap10,xSP_3040ALICE_v2_etaGap10,
                                                              ySP_3040ALICE_v2_etaGap10,xErrSP_3040ALICE_v2_etaGap10,
                                                              yErrSP_3040ALICE_v2_etaGap10);
    myGraphSetUp(GrSP_3040ALICE_v2_etaGap10,1.0,kOpenCircle,kBlue,1,kBlack);
    
    //  v3{SP}(pt) for 30-40%, rapidity gap = 1.0:
    const Int_t nPointsSP_3040ALICE_v3_etaGap10 = 18;
    Double_t xSP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,
        0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.009546,0.013907,0.013673,0.021384,0.029748,
        0.035769,0.037629,0.044245,0.053656,0.061235,0.063584,0.080152,0.089683,0.089203,0.094393,0.115172,0.107330,0.062049};
    Double_t xErrSP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.};
    Double_t yErrSP_3040ALICE_v3_etaGap10[nPointsSP_3040ALICE_v3_etaGap10] = {0.001753,0.001759,0.001898,0.002071,0.002291,
        0.002545,0.002838,0.003136,0.002585,0.003170,0.003875,0.004736,0.005828,0.007141,0.006755,0.009900,0.010855,0.023521};
    TGraphErrors *GrSP_3040ALICE_v3_etaGap10 = new TGraphErrors(nPointsSP_3040ALICE_v3_etaGap10,xSP_3040ALICE_v3_etaGap10,ySP_3040ALICE_v3_etaGap10,
                                                              xErrSP_3040ALICE_v3_etaGap10,yErrSP_3040ALICE_v3_etaGap10);
    myGraphSetUp(GrSP_3040ALICE_v3_etaGap10,1.2,kOpenTriangleUp,kGreen+2,1,kBlack);

    //  v4{SP}(pt) for 30-40%, rapidity gap = 1.0:
    const Int_t nPointsSP_3040ALICE_v4_etaGap10 = 12;
    Double_t xSP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.150000,1.450000,1.750000,2.050000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.002125,0.005223,0.017491,0.020758,0.026932,0.043126,0.055596,0.055244,0.081049,0.097779,0.074055,0.158711};
    Double_t xErrSP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.};
    Double_t yErrSP_3040ALICE_v4_etaGap10[nPointsSP_3040ALICE_v4_etaGap10] = {0.002391,0.002677,0.003227,0.003998,0.004216,0.005692,0.007729,0.010472,0.012836,0.018727,0.020749,0.044065};
    TGraphErrors *GrSP_3040ALICE_v4_etaGap10 = new TGraphErrors(nPointsSP_3040ALICE_v4_etaGap10,xSP_3040ALICE_v4_etaGap10,ySP_3040ALICE_v4_etaGap10,
                                                              xErrSP_3040ALICE_v4_etaGap10,yErrSP_3040ALICE_v4_etaGap10);
    myGraphSetUp(GrSP_3040ALICE_v4_etaGap10,1.,kOpenSquare,kRed,1,kBlack);
        
    //  v5{SP}(pt) for 30-40%, rapidity gap = 1.0:
    const Int_t nPointsSP_3040ALICE_v5_etaGap10 = 12;
    Double_t xSP_3040ALICE_v5_etaGap10[nPointsSP_3040ALICE_v5_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.150000,1.450000,1.750000,2.050000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v5_etaGap10[nPointsSP_3040ALICE_v5_etaGap10] = {-0.001660,-0.000953,0.009396,0.020289,-0.003699,0.016010,0.045845,0.044840,0.055518,0.125920,0.014951,0.017336};
    Double_t xErrSP_3040ALICE_v5_etaGap10[nPointsSP_3040ALICE_v5_etaGap10] = {0.};
    Double_t yErrSP_3040ALICE_v5_etaGap10[nPointsSP_3040ALICE_v5_etaGap10] = {0.006239,0.007340,0.008126,0.010441,0.010926,0.014386,0.020467,0.027239,0.033738,0.051262,0.053235,0.113824};
    TGraphErrors *GrSP_3040ALICE_v5_etaGap10 = new TGraphErrors(10,xSP_3040ALICE_v5_etaGap10,ySP_3040ALICE_v5_etaGap10,
                                                              xErrSP_3040ALICE_v5_etaGap10,yErrSP_3040ALICE_v5_etaGap10);
    myGraphSetUp(GrSP_3040ALICE_v5_etaGap10,1.4,kOpenDiamond,kBlack,1,kBlack);
    //===================================================================================================================
    
    //===================================================================================================================
    //  v2{SP}(pt) for 30-40%, rapidity gap = 0.2:
    const Int_t nPointsSP_3040ALICE_v2_etaGap02 = 23;
    Double_t xSP_3040ALICE_v2_etaGap02[nPointsSP_3040ALICE_v2_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,
        0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.300000,2.500000,2.700000,2.900000,
        3.200000,3.600000,4.000000,4.400000,4.800000};
    Double_t ySP_3040ALICE_v2_etaGap02[nPointsSP_3040ALICE_v2_etaGap02] = {0.039506,0.054606,0.068971,0.084048,0.096809,
        0.109131,0.119861,0.130735,0.145594,0.161632,0.180863,0.194553,0.208946,0.218347,0.223919,0.230006,0.237168,0.243976,
        0.239478,0.237796,0.229257,0.211401,0.199559};
    Double_t xErrSP_3040ALICE_v2_etaGap02[nPointsSP_3040ALICE_v2_etaGap02] = {0.};
    Double_t yErrSP_3040ALICE_v2_etaGap02[nPointsSP_3040ALICE_v2_etaGap02] = {0.000453,0.000449,0.000471,0.000511,0.000564,
        0.000626,0.000697,0.000777,0.000643,0.000799,0.000986,0.001220,0.001512,0.001867,0.002282,0.002804,0.003396,0.004077,
        0.003726,0.005269,0.007165,0.010680,0.010646};
    TGraphErrors *GrSP_3040ALICE_v2_etaGap02 = new TGraphErrors(nPointsSP_3040ALICE_v2_etaGap02,xSP_3040ALICE_v2_etaGap02,
                                                              ySP_3040ALICE_v2_etaGap02,xErrSP_3040ALICE_v2_etaGap02,yErrSP_3040ALICE_v2_etaGap02);
    myGraphSetUp(GrSP_3040ALICE_v2_etaGap02,1.,kFullCircle,kBlue,1,kBlack);
    
    //  v3{SP}(pt) for 30-40%, rapidity gap = 0.2:
    const Int_t nPointsSP_3040ALICE_v3_etaGap02 = 18;
    Double_t xSP_3040ALICE_v3_etaGap02[nPointsSP_3040ALICE_v3_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,
        0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v3_etaGap02[nPointsSP_3040ALICE_v3_etaGap02] = {0.008459,0.013237,0.017909,0.024466,0.029466,
        0.036073,0.040105,0.046120,0.051994,0.064426,0.070241,0.081317,0.092806,0.091925,0.108180,0.114513,0.116033,0.113098};
    Double_t xErrSP_3040ALICE_v3_etaGap02[nPointsSP_3040ALICE_v3_etaGap02] = {0.};
    Double_t yErrSP_3040ALICE_v3_etaGap02[nPointsSP_3040ALICE_v3_etaGap02] = {0.000868,0.000874,0.000933,0.001020,0.001132,
        0.001255,0.001395,0.001548,0.001273,0.001560,0.001911,0.002346,0.002890,0.003561,0.003350,0.004920,0.005420,0.011451};
    TGraphErrors *GrSP_3040ALICE_v3_etaGap02 = new TGraphErrors(nPointsSP_3040ALICE_v3_etaGap02,xSP_3040ALICE_v3_etaGap02,ySP_3040ALICE_v3_etaGap02,
                                                              xErrSP_3040ALICE_v3_etaGap02,yErrSP_3040ALICE_v3_etaGap02);
    myGraphSetUp(GrSP_3040ALICE_v3_etaGap02,1.2,kFullTriangleUp,kGreen+2,1,kBlack);
    
    //  v4{SP}(pt) for 30-40%, rapidity gap = 0.2:
    const Int_t nPointsSP_3040ALICE_v4_etaGap02 = 12;
    Double_t xSP_3040ALICE_v4_etaGap02[nPointsSP_3040ALICE_v4_etaGap02] = {0.300000,0.500000,0.700000,0.900000,1.150000,1.450000,1.750000,2.050000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v4_etaGap02[nPointsSP_3040ALICE_v4_etaGap02] = {0.002607,0.007077,0.014585,0.020707,0.032150,0.044317,0.052834,0.062102,0.085004,0.087401,0.082642,0.094469};
    Double_t xErrSP_3040ALICE_v4_etaGap02[nPointsSP_3040ALICE_v4_etaGap02] = {0.};
    Double_t yErrSP_3040ALICE_v4_etaGap02[nPointsSP_3040ALICE_v4_etaGap02] = {0.001134,0.001267,0.001541,0.001899,0.001991,0.002707,0.003668,0.005000,0.006097,0.008925,0.009853,0.020906};
    TGraphErrors *GrSP_3040ALICE_v4_etaGap02 = new TGraphErrors(nPointsSP_3040ALICE_v4_etaGap02,xSP_3040ALICE_v4_etaGap02,ySP_3040ALICE_v4_etaGap02,
                                                              xErrSP_3040ALICE_v4_etaGap02,yErrSP_3040ALICE_v4_etaGap02);
    myGraphSetUp(GrSP_3040ALICE_v4_etaGap02,1.,kFullSquare,kRed,1,kBlack);

    //  v5{SP}(pt) for 30-40%, rapidity gap = 0.2:
    const Int_t nPointsSP_3040ALICE_v5_etaGap02 = 12;
    Double_t xSP_3040ALICE_v5_etaGap02[nPointsSP_3040ALICE_v5_etaGap02] = {0.300000,0.500000,0.700000,0.900000,1.150000,1.450000,1.750000,2.050000,2.400000,2.800000,3.500000,4.500000};
    Double_t ySP_3040ALICE_v5_etaGap02[nPointsSP_3040ALICE_v5_etaGap02] = {0.005428,0.002807,0.003034,0.008406,0.011691,0.016794,0.033101,0.025018,0.030943,0.045353,0.071972,0.064883};
    Double_t xErrSP_3040ALICE_v5_etaGap02[nPointsSP_3040ALICE_v5_etaGap02] = {0.};
    Double_t yErrSP_3040ALICE_v5_etaGap02[nPointsSP_3040ALICE_v5_etaGap02] = {0.002227,0.002588,0.003172,0.003872,0.004061,0.005520,0.007585,0.010182,0.012424,0.018294,0.020202,0.042575};
    TGraphErrors *GrSP_3040ALICE_v5_etaGap02 = new TGraphErrors(nPointsSP_3040ALICE_v5_etaGap02,xSP_3040ALICE_v5_etaGap02,ySP_3040ALICE_v5_etaGap02,
                                                              xErrSP_3040ALICE_v5_etaGap02,yErrSP_3040ALICE_v5_etaGap02);
    myGraphSetUp(GrSP_3040ALICE_v5_etaGap02,1.4,33,kBlack,1,kBlack);
    //===================================================================================================================  
    
    //===================================================================================================================
    // Schenke: 
    // ideal v3(pt) for 30-40%:
    Double_t xIdealSchenke3040_v3[] = {0.007776,0.041058,0.101300,0.189163,0.305635,0.452111,0.630493,0.843352,
        1.094190,1.387867,1.731373,2.135325,2.617293,3.210890,4.002174};
    Double_t yIdealSchenke3040_v3[] = {1.479520e-04,1.124939e-04,6.584102e-04,3.220716e-03,8.219176e-03,1.547538e-02,2.635454e-02,4.111273e-02,
        6.115409e-02,8.484866e-02,1.116099e-01,1.404768e-01,1.692029e-01,1.963983e-01,2.213624e-01};
    Double_t xErrIdealSchenke3040_v3[15] = {0.};
    Double_t yErrIdealSchenke3040_v3[] = {4.047338e-05,1.796386e-05,5.504791e-05,1.993759e-04,4.913364e-04,8.945152e-04,1.473105e-03,2.234133e-03,
        3.251702e-03,4.456732e-03,5.843234e-03,7.403563e-03,9.146037e-03,1.122705e-02,1.408298e-02};
    Int_t nPointsIdealSchenke3040_v3 = sizeof(xIdealSchenke3040_v3)/sizeof(Double_t);                                      
    TGraphErrors *GrIdealSchenke3040_v3 = new TGraphErrors(nPointsIdealSchenke3040_v3,xIdealSchenke3040_v3,yIdealSchenke3040_v3,
                                                         xErrIdealSchenke3040_v3,yErrIdealSchenke3040_v3);
    myGraphSetUp(GrIdealSchenke3040_v3,1.,33,kBlack,kDashed,kGreen+2,kGreen-10,1001);

    // viscous v3(pt) for 30-40%:
    Double_t xViscousSchenke3040_v3[] = {0.007776,0.041058,0.101300,0.189163,0.305635,0.452111,0.630493,0.843352,
        1.094190,1.387867,1.731373,2.135325,2.617293,3.210890,4.002174};
    Double_t yViscousSchenke3040_v3[] = {1.810421e-04,1.088683e-04,4.874871e-04,2.624849e-03,6.650429e-03,1.214362e-02,2.019854e-02,3.120146e-02,
        4.615661e-02,6.383687e-02,8.328145e-02,1.029784e-01,1.200315e-01,1.310700e-01,1.294474e-01};
    Double_t xErrViscousSchenke3040_v3[15] = {0.};
    Double_t yErrViscousSchenke3040_v3[] = {2.895564e-05,9.863071e-06,4.089573e-05,1.529722e-04,3.555346e-04,6.109573e-04,9.704488e-04,1.454140e-03,
        2.115684e-03,2.906832e-03,3.816101e-03,4.841144e-03,6.011145e-03,7.466822e-03,9.715151e-03};
    Int_t nPointsViscousSchenke3040_v3 = sizeof(xViscousSchenke3040_v3)/sizeof(Double_t);                                      
    TGraphErrors *GrViscousSchenke3040_v3 = new TGraphErrors(nPointsViscousSchenke3040_v3,xViscousSchenke3040_v3,yViscousSchenke3040_v3,
                                                           xErrViscousSchenke3040_v3,yErrViscousSchenke3040_v3);
    myGraphSetUp(GrViscousSchenke3040_v3,1.,33,kBlack,1,kGreen+2,kGreen-10,1001);
    //===================================================================================================================	
    
    //===================================================================================================================
    // Schenke: 
    // ideal v2(pt) for 30-40%:
    Double_t xIdealSchenke3040_v2[] = {0.007776,0.041058,0.101300,0.189163,0.305635,0.452111,0.630493,0.843352,
        1.094190,1.387867,1.731373,2.135325,2.617293,3.210890,4.002174};
    Double_t yIdealSchenke3040_v2[] = {2.558414e-05,9.867626e-04,5.661819e-03,1.649046e-02,3.377584e-02,5.560800e-02,8.121079e-02,1.086796e-01,
        1.375457e-01,1.667749e-01,1.964467e-01,2.260053e-01,2.560743e-01,2.863471e-01,3.177713e-01};
    Double_t xErrIdealSchenke3040_v2[15] = {0.};
    Double_t yErrIdealSchenke3040_v2[] = {2.497460e-05,4.656790e-05,2.650216e-04,7.813507e-04,1.647879e-03,2.800632e-03,4.085668e-03,5.417398e-03,
        6.708719e-03,7.947069e-03,9.168744e-03,1.043873e-02,1.204135e-02,1.416794e-02,1.715206e-02};
    Int_t nPointsIdealSchenke3040_v2 = sizeof(xIdealSchenke3040_v2)/sizeof(Double_t);                                      
    TGraphErrors *GrIdealSchenke3040_v2 = new TGraphErrors(nPointsIdealSchenke3040_v2,xIdealSchenke3040_v2,yIdealSchenke3040_v2,
                                                         xErrIdealSchenke3040_v2,yErrIdealSchenke3040_v2);
    myGraphSetUp(GrIdealSchenke3040_v2,1.,33,kBlack,kDashed,kBlue,kBlue-10,1001);
    
    // viscous v2(pt) for 30-40%:
    Double_t xViscousSchenke3040_v2[] = {0.007776,0.041058,0.101300,0.189163,0.305635,0.452111,0.630493,0.843352,
        1.094190,1.387867,1.731373,2.135325,2.617293,3.210890,4.002174};
    Double_t yViscousSchenke3040_v2[] = {-2.625537e-05,9.441659e-04,5.499670e-03,1.580113e-02,3.170383e-02,5.125177e-02,7.392717e-02,9.807519e-02,
        1.229014e-01,1.475803e-01,1.719447e-01,1.952542e-01,2.174093e-01,2.368600e-01,2.511733e-01};
    Double_t xErrViscousSchenke3040_v2[15] = {0.};
    Double_t yErrViscousSchenke3040_v2[] = {2.033230e-05,3.826003e-05,2.169308e-04,6.176627e-04,1.242795e-03,2.027983e-03,2.914887e-03,3.833102e-03,
        4.703648e-03,5.516849e-03,6.311937e-03,7.179270e-03,8.381866e-03,1.020568e-02,1.325437e-02};
    Int_t nPointsViscousSchenke3040_v2 = sizeof(xViscousSchenke3040_v2)/sizeof(Double_t);                                      
    TGraphErrors *GrViscousSchenke3040_v2 = new TGraphErrors(nPointsViscousSchenke3040_v2,xViscousSchenke3040_v2,yViscousSchenke3040_v2,
                                                           xErrViscousSchenke3040_v2,yErrViscousSchenke3040_v2);

    myGraphSetUp(GrViscousSchenke3040_v2,1.,33,kBlack,1,kBlue,kBlue-10,1001);
    //===================================================================================================================
    
    //===================================================================================================================
    // Figure 3b:   
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
    TGraphErrors *GrSP_0005ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v2_etaGap10,xSP_0005ALICE_v2_etaGap10,ySP_0005ALICE_v2_etaGap10,
                                                              xErrSP_0005ALICE_v2_etaGap10,yErrSP_0005ALICE_v2_etaGap10);
    myGraphSetUp(GrSP_0005ALICE_v2_etaGap10,1.,kOpenCircle,kBlue,1,kBlack,1.,kBlack);

    //  v3{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v3_etaGap10 = 16;
    Double_t xSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.900000,1.100000,
        1.300000,1.500000,1.700000,1.900000,2.300000,2.900000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.006303,0.009800,0.011143,0.014246,0.017628,0.019437,0.028412,0.030580,
        0.038730,0.045653,0.052469,0.062303,0.071522,0.082223,0.083373,0.076951};
    Double_t xErrSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v3_etaGap10[nPointsSP_0005ALICE_v3_etaGap10] = {0.001012,0.000996,0.001062,0.001158,0.001277,0.001415,0.001158,0.001422,
        0.001734,0.002124,0.002610,0.003206,0.002724,0.004977,0.008123,0.017180};
    TGraphErrors *GrSP_0005ALICE_v3_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v3_etaGap10,xSP_0005ALICE_v3_etaGap10,ySP_0005ALICE_v3_etaGap10,
                                                              xErrSP_0005ALICE_v3_etaGap10,yErrSP_0005ALICE_v3_etaGap10);
    //ShiftAlongXaxis(SP_0005ALICE_v3_etaGap10,-0.025);    
    myGraphSetUp(GrSP_0005ALICE_v3_etaGap10,1.2,kOpenTriangleUp,kGreen+2,1.,kBlack);    

    //  v4{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v4_etaGap10 = 11;
    Double_t xSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.300000,0.500000,0.700000,0.950000,1.250000,1.550000,1.850000,2.300000,2.900000,3.650000,4.550000};
    Double_t ySP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.002042,0.002556,0.009693,0.013286,0.016780,0.027865,0.031797,0.051101,0.060164,
        0.095985,0.094607};
    Double_t xErrSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v4_etaGap10[nPointsSP_0005ALICE_v4_etaGap10] = {0.001460,0.001624,0.001930,0.002021,0.002737,0.003717,0.005042,0.005564,0.010160,
        0.016472,0.035083};
    TGraphErrors *GrSP_0005ALICE_v4_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v4_etaGap10,xSP_0005ALICE_v4_etaGap10,ySP_0005ALICE_v4_etaGap10,
                                                              xErrSP_0005ALICE_v4_etaGap10,yErrSP_0005ALICE_v4_etaGap10);
    //ShiftAlongXaxis(SP_0005ALICE_v4_etaGap10,0.025); 
    myGraphSetUp(GrSP_0005ALICE_v4_etaGap10,1.,kOpenSquare,kRed,1.,kBlack);    
    
    //  v5{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v5_etaGap10 = 12;
    Double_t xSP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.100000,1.300000,1.600000,2.000000,2.400000,
        2.800000,3.500000,4.500000};
    Double_t ySP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.002016,0.003409,0.004029,0.002665,0.002765,0.003042,0.013241,0.015430,0.031845,
        0.031373,0.068504,0.017964};
    Double_t xErrSP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.};
    Double_t yErrSP_0005ALICE_v5_etaGap10[nPointsSP_0005ALICE_v5_etaGap10] = {0.001260,0.001386,0.001696,0.002101,0.002560,0.003119,0.002970,0.004472,0.006802,
        0.010073,0.011899,0.027756};
    TGraphErrors *GrSP_0005ALICE_v5_etaGap10 = new TGraphErrors(nPointsSP_0005ALICE_v5_etaGap10,xSP_0005ALICE_v5_etaGap10,ySP_0005ALICE_v5_etaGap10,
                                                              xErrSP_0005ALICE_v5_etaGap10,yErrSP_0005ALICE_v5_etaGap10);
    myGraphSetUp(GrSP_0005ALICE_v5_etaGap10,1.4,kOpenDiamond,kBlack,1.,kBlack);    
    //===================================================================================================================
    
    //===================================================================================================================
    //  v2{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v2_etaGap02 = 26;
    Double_t xSP_0005ALICE_v2_etaGap02[nPointsSP_0005ALICE_v2_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
        1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.850000,1.950000,
        2.100000,2.300000,2.500000,2.700000,2.900000,3.250000,3.750000,4.500000};
    Double_t ySP_0005ALICE_v2_etaGap02[nPointsSP_0005ALICE_v2_etaGap02] = {0.012030,0.016538,0.020312,0.024888,0.028227,0.032597,0.035418,0.039010,
        0.040737,0.044450,0.045448,0.048067,0.050659,0.055360,0.056007,0.057680,0.059254,0.062111,
        0.062703,0.064158,0.066930,0.074260,0.073086,0.073927,0.080837,0.064785};
    Double_t xErrSP_0005ALICE_v2_etaGap02[nPointsSP_0005ALICE_v2_etaGap02] = {0.};
    Double_t yErrSP_0005ALICE_v2_etaGap02[nPointsSP_0005ALICE_v2_etaGap02] = {0.000436,0.000430,0.000453,0.000492,0.000543,0.000596,0.000668,0.000737,
        0.000823,0.000912,0.001015,0.001121,0.001251,0.001386,0.001545,0.001714,0.001907,0.002121,
        0.001744,0.002149,0.002632,0.003226,0.003930,0.003436,0.005465,0.006891};
    TGraphErrors *GrSP_0005ALICE_v2_etaGap02 = new TGraphErrors(nPointsSP_0005ALICE_v2_etaGap02,xSP_0005ALICE_v2_etaGap02,ySP_0005ALICE_v2_etaGap02,
                                                              xErrSP_0005ALICE_v2_etaGap02,yErrSP_0005ALICE_v2_etaGap02);
    myGraphSetUp(GrSP_0005ALICE_v2_etaGap02,1.,kFullCircle,kBlue,1.,kBlack);    
    
    //  v3{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v3_etaGap02 = 20;
    Double_t xSP_0005ALICE_v3_etaGap02[nPointsSP_0005ALICE_v3_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
        1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.300000,2.600000,3.000000,3.500000,
        4.100000,4.700000};
    Double_t ySP_0005ALICE_v3_etaGap02[nPointsSP_0005ALICE_v3_etaGap02] = {0.005264,0.008365,0.011069,0.014565,0.015858,0.021911,0.026581,0.029673,
        0.034343,0.039329,0.044124,0.054249,0.059637,0.066966,0.076106,0.080674,0.091843,0.089706,
        0.092376,0.054630};
    Double_t xErrSP_0005ALICE_v3_etaGap02[nPointsSP_0005ALICE_v3_etaGap02] = {0.};
    Double_t yErrSP_0005ALICE_v3_etaGap02[nPointsSP_0005ALICE_v3_etaGap02] = {0.000596,0.000590,0.000626,0.000678,0.000751,0.000823,0.000909,0.001005,
        0.000831,0.001014,0.001251,0.001539,0.001895,0.002331,0.002882,0.002720,0.004047,0.005268,
        0.008874,0.014420};
    TGraphErrors *GrSP_0005ALICE_v3_etaGap02 = new TGraphErrors(nPointsSP_0005ALICE_v3_etaGap02,xSP_0005ALICE_v3_etaGap02,ySP_0005ALICE_v3_etaGap02,
                                                              xErrSP_0005ALICE_v3_etaGap02,yErrSP_0005ALICE_v3_etaGap02);
    myGraphSetUp(GrSP_0005ALICE_v3_etaGap02,1.2,kFullTriangleUp,kGreen+2,1.,kBlack);    
    
    //  v4{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v4_etaGap02 = 20;
    Double_t xSP_0005ALICE_v4_etaGap02[nPointsSP_0005ALICE_v4_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
        1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.300000,2.600000,3.000000,3.500000,
        4.100000,4.700000};
    Double_t ySP_0005ALICE_v4_etaGap02[nPointsSP_0005ALICE_v4_etaGap02] = {0.001858,0.002950,0.003658,0.005777,0.008362,0.012289,0.012218,0.012171,
        0.017361,0.020719,0.024185,0.030187,0.037104,0.042695,0.050235,0.056076,0.075308,0.080679,
        0.090205,0.077760};
    Double_t xErrSP_0005ALICE_v4_etaGap02[nPointsSP_0005ALICE_v4_etaGap02] = {0.};
    Double_t yErrSP_0005ALICE_v4_etaGap02[nPointsSP_0005ALICE_v4_etaGap02] = {0.000955,0.000951,0.001013,0.001098,0.001204,0.001326,0.001472,0.001634,
        0.001341,0.001641,0.002016,0.002476,0.003059,0.003750,0.004616,0.004365,0.006467,0.008390,
        0.014154,0.022657};
    TGraphErrors *GrSP_0005ALICE_v4_etaGap02 = new TGraphErrors(nPointsSP_0005ALICE_v4_etaGap02,xSP_0005ALICE_v4_etaGap02,ySP_0005ALICE_v4_etaGap02,
                                                              xErrSP_0005ALICE_v4_etaGap02,yErrSP_0005ALICE_v4_etaGap02);
    myGraphSetUp(GrSP_0005ALICE_v4_etaGap02,1.,kFullSquare,kRed,1.,kBlack);    
    
    //  v5{SP}(pt) for 00-05%:
    const Int_t nPointsSP_0005ALICE_v5_etaGap02 = 13;
    Double_t xSP_0005ALICE_v5_etaGap02[nPointsSP_0005ALICE_v5_etaGap02] = {0.300000,0.500000,0.700000,0.900000,1.200000,1.600000,2.000000,2.400000,2.800000,
        3.300000,3.900000,4.500000,4.900000};
    Double_t ySP_0005ALICE_v5_etaGap02[nPointsSP_0005ALICE_v5_etaGap02] = {0.000227,-0.000991,0.006407,0.010032,0.002788,0.008782,0.015454,0.012575,0.044690,
        0.061348,0.089657,-0.023651,-0.184366};
    Double_t xErrSP_0005ALICE_v5_etaGap02[nPointsSP_0005ALICE_v5_etaGap02] = {0.};
    Double_t yErrSP_0005ALICE_v5_etaGap02[nPointsSP_0005ALICE_v5_etaGap02] = {0.001875,0.002131,0.002423,0.003018,0.002881,0.004339,0.006567,0.009819,0.014850,
        0.019363,0.034072,0.053514,0.131967};
    TGraphErrors *GrSP_0005ALICE_v5_etaGap02 = new TGraphErrors(nPointsSP_0005ALICE_v5_etaGap02,xSP_0005ALICE_v5_etaGap02,ySP_0005ALICE_v5_etaGap02,
                                                              xErrSP_0005ALICE_v5_etaGap02,yErrSP_0005ALICE_v5_etaGap02);
    myGraphSetUp(GrSP_0005ALICE_v5_etaGap02,1.4,33,kBlack,1.,kBlack);    
    //=================================================================================================================== 

    
    //===================================================================================================================
    // Figure 3c:   
 //===================================================================================================================
 //  v2{SP}(pt) for 00-02%, eta gap = 0.2:
 const Int_t nPointsSP_0002_v2_etaGap02 = 27;
 Double_t xSP_0002_v2_etaGap02[nPointsSP_0002_v2_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
0.850000,0.950000,1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.900000,2.100000,2.300000,
2.500000,2.700000,2.900000,3.200000,3.600000,4.000000,4.400000,4.800000};
 Double_t ySP_0002_v2_etaGap02[nPointsSP_0002_v2_etaGap02] = {0.010458,0.013621,0.017307,0.020609,0.024900,0.027688,
0.029414,0.031976,0.034382,0.035864,0.038812,0.040049,0.043550,0.044058,0.045320,0.046360,0.047605,0.053355,0.056907,
0.057136,0.060536,0.058016,0.060921,0.054174,0.042666,0.046099,0.038495};
 Double_t xErrSP_0002_v2_etaGap02[nPointsSP_0002_v2_etaGap02] = {0.};
 Double_t yErrSP_0002_v2_etaGap02[nPointsSP_0002_v2_etaGap02] = {0.000310,0.000306,0.000323,0.000352,0.000386,0.000428,
0.000476,0.000529,0.000588,0.000653,0.000722,0.000802,0.000892,0.000991,0.001103,0.001224,0.001009,0.001243,0.001524,
0.001876,0.002292,0.002802,0.002624,0.003827,0.005442,0.008369,0.008633};
 TGraphErrors *GrSP_0002_v2_etaGap02 = new TGraphErrors(nPointsSP_0002_v2_etaGap02,xSP_0002_v2_etaGap02,ySP_0002_v2_etaGap02,
                                                      xErrSP_0002_v2_etaGap02,yErrSP_0002_v2_etaGap02);
 myGraphSetUp(GrSP_0002_v2_etaGap02,1.,kFullCircle,kBlue,1.,kBlack);    

 //  v2{SP}(pt) for 00-02%, eta gap = 1.0:
 const Int_t nPointsSP_0002_v2_etaGap10 = 15;
 Double_t xSP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
0.900000,1.100000,1.350000,1.650000,1.950000,2.250000,2.700000,3.500000,4.500000};
 Double_t ySP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.010171,0.013190,0.017342,0.020629,0.022617,0.026549,
0.027423,0.032261,0.037467,0.041001,0.045763,0.049327,0.049688,0.051480,0.038527};
 Double_t xErrSP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.};
 Double_t yErrSP_0002_v2_etaGap10[nPointsSP_0002_v2_etaGap10] = {0.000600,0.000590,0.000625,0.000683,0.000757,0.000839,
0.000692,0.000848,0.000888,0.001209,0.001653,0.002252,0.002465,0.003968,0.009391};
 TGraphErrors *GrSP_0002_v2_etaGap10 = new TGraphErrors(nPointsSP_0002_v2_etaGap10,xSP_0002_v2_etaGap10,ySP_0002_v2_etaGap10,
                                                  xErrSP_0002_v2_etaGap10,yErrSP_0002_v2_etaGap10);
 myGraphSetUp(GrSP_0002_v2_etaGap10,1.,kOpenCircle,kBlue,1.,kBlack);    

 //  v3{SP}(pt) for 00-02%, eta gap = 0.2:
 const Int_t nPointsSP_0002_v3_etaGap02 = 25;
 Double_t xSP_0002_v3_etaGap02[nPointsSP_0002_v3_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
0.850000,0.950000,1.050000,1.150000,1.250000,1.350000,1.450000,1.550000,1.650000,1.750000,1.900000,2.100000,2.300000,
2.500000,2.800000,3.200000,3.600000,4.000000,4.600000};
 Double_t ySP_0002_v3_etaGap02[nPointsSP_0002_v3_etaGap02] = {0.005351,0.007675,0.010356,0.013931,0.017761,0.019780,
0.023915,0.026834,0.030480,0.032369,0.037061,0.037480,0.041826,0.044985,0.048539,0.053263,0.058862,0.063724,0.070478,
0.079161,0.080104,0.088543,0.091927,0.081413,0.075981};
 Double_t xErrSP_0002_v3_etaGap02[nPointsSP_0002_v3_etaGap02] = {0.};
 Double_t yErrSP_0002_v3_etaGap02[nPointsSP_0002_v3_etaGap02] = {0.000353,0.000350,0.000370,0.000402,0.000441,0.000487,
0.000540,0.000597,0.000663,0.000734,0.000810,0.000899,0.001002,0.001112,0.001234,0.001368,0.001130,0.001391,0.001705,
0.002089,0.001980,0.002936,0.004254,0.006087,0.006718};
 TGraphErrors *GrSP_0002_v3_etaGap02 = new TGraphErrors(nPointsSP_0002_v3_etaGap02,xSP_0002_v3_etaGap02,ySP_0002_v3_etaGap02,
                                                      xErrSP_0002_v3_etaGap02,yErrSP_0002_v3_etaGap02);
 myGraphSetUp(GrSP_0002_v3_etaGap02,1.2,kFullTriangleUp,kGreen+2,1.,kBlack);    
   
 //  v3{SP}(pt) for 00-02%, eta gap = 1.0:
 const Int_t nPointsSP_0002_v3_etaGap10 = 15;
 Double_t xSP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
0.900000,1.100000,1.350000,1.650000,1.950000,2.250000,2.700000,3.500000,4.500000};
 Double_t ySP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.006592,0.007286,0.012180,0.012242,0.017416,0.018393,
0.024716,0.030980,0.037703,0.046558,0.051285,0.064613,0.074831,0.077093,0.082442};
 Double_t xErrSP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.};
 Double_t yErrSP_0002_v3_etaGap10[nPointsSP_0002_v3_etaGap10] = {0.000682,0.000676,0.000713,0.000782,0.000860,0.000953,
0.000782,0.000957,0.001002,0.001361,0.001862,0.002541,0.002767,0.004466,0.010586};
 TGraphErrors *GrSP_0002_v3_etaGap10 = new TGraphErrors(nPointsSP_0002_v3_etaGap10,xSP_0002_v3_etaGap10,ySP_0002_v3_etaGap10,
                                                      xErrSP_0002_v3_etaGap10,yErrSP_0002_v3_etaGap10);
 myGraphSetUp(GrSP_0002_v3_etaGap10,1.,kOpenTriangleUp,kGreen+2,1.,kBlack);    
 
 //  v4{SP}(pt) for 00-02%, eta gap = 0.2:
 const Int_t nPointsSP_0002_v4_etaGap02 = 20;
 Double_t xSP_0002_v4_etaGap02[nPointsSP_0002_v4_etaGap02] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.200000,2.600000,3.000000,3.400000,3.800000,4.250000,4.750000};
 Double_t ySP_0002_v4_etaGap02[nPointsSP_0002_v4_etaGap02] = {0.001803,0.003590,0.003214,0.005683,0.007166,0.007113,
0.010206,0.012110,0.016560,0.019941,0.026784,0.034169,0.036075,0.046846,0.058077,0.060063,0.072744,0.076006,
0.067621,0.044664};
 Double_t xErrSP_0002_v4_etaGap02[nPointsSP_0002_v4_etaGap02] = {0.};
 Double_t yErrSP_0002_v4_etaGap02[nPointsSP_0002_v4_etaGap02] = {0.000630,0.000623,0.000666,0.000721,0.000794,0.000881,
0.000973,0.001076,0.000882,0.001080,0.001331,0.001640,0.002016,0.001922,0.002879,0.004289,0.006267,0.009064,
0.011705,0.017076};
 TGraphErrors *GrSP_0002_v4_etaGap02 = new TGraphErrors(nPointsSP_0002_v4_etaGap02,xSP_0002_v4_etaGap02,ySP_0002_v4_etaGap02,
                                                      xErrSP_0002_v4_etaGap02,yErrSP_0002_v4_etaGap02);
 myGraphSetUp(GrSP_0002_v4_etaGap02,1.,kFullSquare,kRed,1.,kBlack);    

 //  v4{SP}(pt) for 00-02%, eta gap = 1.0:
 const Int_t nPointsSP_0002_v4_etaGap10 = 15;
 Double_t xSP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,
0.900000,1.100000,1.350000,1.650000,1.950000,2.250000,2.700000,3.500000,4.500000};
 Double_t ySP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {-0.000533,0.001167,0.002081,0.005218,0.006826,0.008440,
0.013009,0.014812,0.017125,0.030106,0.038279,0.050488,0.067640,0.071637,0.084239};
 Double_t xErrSP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {0.};
 Double_t yErrSP_0002_v4_etaGap10[nPointsSP_0002_v4_etaGap10] = {0.001427,0.001398,0.001482,0.001594,0.001758,0.001945,
0.001593,0.001951,0.002046,0.002787,0.003802,0.005182,0.005663,0.009064,0.021449};
 TGraphErrors *GrSP_0002_v4_etaGap10 = new TGraphErrors(nPointsSP_0002_v4_etaGap10,xSP_0002_v4_etaGap10,ySP_0002_v4_etaGap10,
                                                      xErrSP_0002_v4_etaGap10,yErrSP_0002_v4_etaGap10);
 myGraphSetUp(GrSP_0002_v4_etaGap10,1.,kOpenSquare,kRed,1.,kBlack);    
    
 //  v5{SP}(pt) for 00-02%, eta gap = 0.2:
 const Int_t nPointsSP_0002_v5_etaGap02 = 13;
 Double_t xSP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.300000,0.500000,0.700000,0.900000,1.100000,1.300000,1.500000,
1.700000,2.000000,2.550000,3.250000,3.950000,4.650000};
 Double_t ySP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.000570,0.002922,0.002151,0.005256,0.006287,0.005849,0.009399,
0.011420,0.012455,0.032134,0.057009,0.020607,0.013551};
 Double_t xErrSP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.};
 Double_t yErrSP_0002_v5_etaGap02[nPointsSP_0002_v5_etaGap02] = {0.001074,0.001155,0.001433,0.001725,0.002123,0.002608,0.003196,
0.003930,0.003755,0.004869,0.009719,0.018353,0.031814};
 TGraphErrors *GrSP_0002_v5_etaGap02 = new TGraphErrors(11,xSP_0002_v5_etaGap02,ySP_0002_v5_etaGap02,
                                                      xErrSP_0002_v5_etaGap02,yErrSP_0002_v5_etaGap02);
 myGraphSetUp(GrSP_0002_v5_etaGap02,1.4,33,kBlack,1.,kBlack);    

 //  v5{SP}(pt) for 00-02%, eta gap = 1.0:
 const Int_t nPointsSP_0002_v5_etaGap10 = 13;
 Double_t xSP_0002_v5_etaGap10[nPointsSP_0002_v5_etaGap10] = {0.300000,0.500000,0.700000,0.900000,1.100000,1.300000,1.500000,
1.700000,2.000000,2.550000,3.250000,3.950000,4.650000};
 Double_t ySP_0002_v5_etaGap10[nPointsSP_0002_v5_etaGap10] = {0.004709,0.002816,0.004686,0.003099,-0.000160,0.006887,0.007711,
0.014473,0.009113,0.033936,0.056363,-0.013242,0.023228};
 Double_t xErrSP_0002_v5_etaGap10[nPointsSP_0002_v5_etaGap10] = {0.};
 Double_t yErrSP_0002_v5_etaGap10[nPointsSP_0002_v5_etaGap10] = {0.001725,0.001974,0.002383,0.002966,0.003646,0.004383,0.005409,
0.006631,0.006321,0.008218,0.016372,0.030968,0.053822};
 TGraphErrors *GrSP_0002_v5_etaGap10 = new TGraphErrors(11,xSP_0002_v5_etaGap10,ySP_0002_v5_etaGap10,
                                                      xErrSP_0002_v5_etaGap10,yErrSP_0002_v5_etaGap10);
 myGraphSetUp(GrSP_0002_v5_etaGap10,1.4,kOpenDiamond,kBlack,1.,kBlack);    
 //===================================================================================================================   
    
    
 // Style histogramm:
 TH1F *myBlankHisto = new TH1F("myBlankHisto","Blank Histogram",100,0.,10.);
 //myBlankHisto->SetXTitle("p_{t} (GeV/#font[72]{c})");
 myBlankHisto->SetYTitle("v_{n}");
 myBlankHisto->SetNdivisions(505,"y");
 myBlankHisto->GetYaxis()->SetTitleOffset(0.6);
 //myBlankHisto->GetYaxis()->SetLabelSize(0.05);
 //myBlankHisto->GetYaxis()->SetTitleSize(0.07);
 //myBlankHisto->GetYaxis()->SetRangeUser(-0.02101,0.399);
 myBlankHisto->GetXaxis()->SetRangeUser(0.,5.);    
      
 // Main canvas:   
 TCanvas *myCan = new TCanvas("myCan",cStamp1,500,1000);
 myCan->Divide(1,3,0.,0.);  
 myCan->SetTopMargin(0.0002);
 myCan->Draw();

 //=================================================================================================
 
 // Figure 3a:
 myCan->cd(1);   
 gPad->SetTopMargin(0.001); 
 gPad->SetBottomMargin(0.); 
 gPad->SetRightMargin(0.001); 
 TH1D* styleHistPad1 = myBlankHisto->Clone("styleHistPad1");
 styleHistPad1->GetYaxis()->SetRangeUser(-0.02101,0.3799);  
 styleHistPad1->Draw();
 //TLatex *system = new TLatex(0.2,0.93,"ALICE Preliminary, Pb-Pb events at  #sqrt{s_{NN}} = 2.76 TeV");
 //system->SetNDC();
 //  system->SetTextFont(42);
 //system->SetTextSize(0.0475);
 //system->SetTextColor(myDarkRed);
 //system->Draw();
 TLatex *a = new TLatex(0.95,0.93,"(a)");
 a->SetNDC();
 a->SetTextSize(0.05);
 a->Draw();
 //GrIdealSchenke3040_v3->Draw("same3");
 //GrIdealSchenke3040_v3->Draw("lsameX");
 GrViscousSchenke3040_v3->Draw("same3");
 GrViscousSchenke3040_v3->Draw("lsameX");
 GrIdealSchenke3040_v3->Draw("same3");
 GrIdealSchenke3040_v3->Draw("lsameX");
 //GrIdealSchenke3040_v2->Draw("same3");
 //GrIdealSchenke3040_v2->Draw("lsameX");
 GrViscousSchenke3040_v2->Draw("same3");
 GrViscousSchenke3040_v2->Draw("lsameX");
 GrIdealSchenke3040_v2->Draw("same3");
 GrIdealSchenke3040_v2->Draw("lsameX");
 GrSP_3040ALICE_v2_etaGap02->Draw("psame");
 GrSP_3040ALICE_v2_etaGap10->Draw("psame");
 GrSP_3040ALICE_v3_etaGap02->Draw("psame");
 GrSP_3040ALICE_v3_etaGap10->Draw("psame");
 GrSP_3040ALICE_v4_etaGap02->Draw("psame");
 GrSP_3040ALICE_v4_etaGap10->Draw("psame");
 GrSP_3040ALICE_v5_etaGap02->Draw("psame");
 GrSP_3040ALICE_v5_etaGap10->Draw("psame");
 // Common legend:    
 TLegend *myLegendPad1 = new TLegend(0.12,0.50,0.48,0.98,"Centrality 30-40%");
 myLegendSetUp(myLegendPad1,0.04);
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v2_etaGap10,"v_{2}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v2_etaGap02,"v_{2}{SP} (#Delta#eta=0.2)","p");
 myLegendPad1->AddEntry(GrSP_3040ALICE_v2_etaGap02,"v_{2}{2}","p");
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v3_etaGap10,"v_{3}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v3_etaGap02,"v_{3}{SP} (#Delta#eta=0.2)","p");
 myLegendPad1->AddEntry(GrSP_3040ALICE_v3_etaGap02,"v_{3}{2}","p");
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v4_etaGap10,"v_{4}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v4_etaGap02,"v_{4}{SP} (#Delta#eta=0.2)","p");
 myLegendPad1->AddEntry(GrSP_3040ALICE_v4_etaGap02,"v_{4}{2}","p");
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v5_etaGap10,"v_{5}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad1->AddEntry(GrSP_3040ALICE_v5_etaGap02,"v_{5}{SP} (#Delta#eta=0.2)","p");
 myLegendPad1->AddEntry(GrSP_3040ALICE_v5_etaGap02,"v_{5}{2}","p");
 myLegendPad1->AddEntry(GrIdealSchenke3040_v2,"v_{2} (#eta/s = 0.0)","fl");
 myLegendPad1->AddEntry(GrViscousSchenke3040_v2,"v_{2} (#eta/s = 0.08)","fl");
 myLegendPad1->AddEntry(GrIdealSchenke3040_v3,"v_{3} (#eta/s = 0.0)","fl");
 myLegendPad1->AddEntry(GrViscousSchenke3040_v3,"v_{3} (#eta/s = 0.08)","fl");
 //myLegendPad1->SetTextSize(0.024);
 myLegendPad1->SetMargin(0.2); 
 myLegendPad1->SetHeader("Centrality 30-40%");
 //TLegendEntry *header = (TLegendEntry*)myLegendPad1->GetListOfPrimitives()->First(); 
 //header->SetTextAlign(22); // center header
 //header->SetTextSize(0.027); 
 //myLegendPad1->SetTextFont(22); // 22 = Times New Roman (bold)
 myLegendPad1->Draw("same"); 

 //=================================================================================================

 // Figure 3b:
 myCan->cd(2);   
 gPad->SetTopMargin(0.); 
 gPad->SetBottomMargin(0.); 
 gPad->SetRightMargin(0.001); 
 TH1D* styleHistPad2 = myBlankHisto->Clone("styleHistPad2");
 styleHistPad2->GetYaxis()->SetRangeUser(-0.01,0.1299);  
 styleHistPad2->Draw();
 TLatex *b = new TLatex(0.95,0.93,"(b)");
 b->SetNDC();
 b->SetTextSize(0.05);
 b->Draw();
 TLegend *myLegendPad2 = new TLegend(0.12,0.6,0.45,0.98,"Centrality 0-5%");
 myLegendSetUp(myLegendPad2,0.04);
 //myLegendPad2->AddEntry(SP_0005ALICE_v2_etaGap10,"v_{2}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad2->AddEntry(SP_0005ALICE_v2_etaGap02,"v_{2}{SP} (#Delta#eta=0.2)","p");
 //myLegendPad2->AddEntry(SP_0005ALICE_v3_etaGap10,"v_{3}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad2->AddEntry(SP_0005ALICE_v3_etaGap02,"v_{3}{SP} (#Delta#eta=0.2)","p");
 //myLegendPad2->AddEntry(SP_0005ALICE_v4_etaGap10,"v_{4}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad2->AddEntry(SP_0005ALICE_v4_etaGap02,"v_{4}{SP} (#Delta#eta=0.2)","p");
 //myLegendPad2->AddEntry(SP_0005ALICE_v5_etaGap10,"v_{5}{SP} (#Delta#eta=1.0)","p");
 //myLegendPad2->AddEntry(SP_0005ALICE_v5_etaGap02,"v_{5}{SP} (#Delta#eta=0.2)","p");
 myLegendPad2->AddEntry(GrSP_0005ALICE_v2_etaGap02,"v_{2}{2}","p");
 myLegendPad2->AddEntry(GrSP_0005ALICE_v3_etaGap02,"v_{3}{2}","p");
 myLegendPad2->AddEntry(GrSP_0005ALICE_v4_etaGap02,"v_{4}{2}","p");
 myLegendPad2->AddEntry(GrSP_0005ALICE_v5_etaGap02,"v_{5}{2}","p");
 myLegendPad2->Draw("same");
 GrSP_0005ALICE_v2_etaGap02->Draw("psame");
 GrSP_0005ALICE_v3_etaGap02->Draw("psame");
 GrSP_0005ALICE_v4_etaGap02->Draw("psame");
 GrSP_0005ALICE_v5_etaGap02->Draw("psame");
 GrSP_0005ALICE_v2_etaGap10->Draw("psame");
 GrSP_0005ALICE_v3_etaGap10->Draw("psame");
 GrSP_0005ALICE_v4_etaGap10->Draw("psame");
 GrSP_0005ALICE_v5_etaGap10->Draw("psame");
    
 //=================================================================================================

 // Figure 3c:
 myCan->cd(3);   
 gPad->SetTopMargin(0.); 
 gPad->SetBottomMargin(0.15); 
 gPad->SetRightMargin(0.001); 
 TH1D* styleHistPad3 = myBlankHisto->Clone("styleHistPad3");
 styleHistPad3->SetXTitle("p_{t} (GeV/#font[72]{c})");
 styleHistPad3->GetYaxis()->SetRangeUser(-0.01,0.1099);  
 styleHistPad3->Draw();
 TLatex *c = new TLatex(0.95,0.93,"(c)");
 c->SetNDC();
 c->SetTextSize(0.05);
 c->Draw();
 TLegend *myLegendPad3 = new TLegend(0.12,0.65,0.45,0.98,"Centrality 0-2%");
 myLegendSetUp(myLegendPad3,0.04);
 // myLegendPad3->AddEntry(SP_0002_v2_etaGap02,"v_{2}{SP} (#Delta#eta = 0.2)","p");
 myLegendPad3->AddEntry(GrSP_0002_v2_etaGap02,"v_{2}{2}","p");
 // myLegendPad3->AddEntry(SP_0002_v2_etaGap10,"v_{2}{SP} (#Delta#eta = 1.0)","p");
 // myLegendPad3->AddEntry(SP_0002_v3_etaGap02,"v_{3}{SP} (#Delta#eta = 0.2)","p");
 myLegendPad3->AddEntry(GrSP_0002_v3_etaGap02,"v_{3}{2}","p");
 // myLegendPad3->AddEntry(SP_0002_v3_etaGap10,"v_{3}{SP} (#Delta#eta = 1.0)","p");
 myLegendPad3->AddEntry(GrSP_0002_v4_etaGap02,"v_{4}{2}","p");
 // myLegendPad3->AddEntry(SP_0002_v4_etaGap02,"v_{4}{SP} (#Delta#eta = 0.2)","p");
 // myLegendPad3->AddEntry(SP_0002_v4_etaGap10,"v_{4}{SP} (#Delta#eta = 1.0)","p");
 // myLegendPad3->AddEntry(SP_0002_v5_etaGap02,"v_{5}{SP} (#Delta#eta = 0.2)","p");
 myLegendPad3->AddEntry(GrSP_0002_v5_etaGap02,"v_{5}{2}","p");
 myLegendPad3->Draw("same");
 // myLegendPad3->AddEntry(SP_0002_v5_etaGap10,"v_{5}{SP} (#Delta#eta = 1.0)","p");
 GrSP_0002_v2_etaGap02->Draw("psame");
 GrSP_0002_v2_etaGap10->Draw("psame");
 GrSP_0002_v3_etaGap02->Draw("psame");
 GrSP_0002_v3_etaGap10->Draw("psame");
 GrSP_0002_v4_etaGap02->Draw("psame");
 GrSP_0002_v4_etaGap10->Draw("psame");
 GrSP_0002_v5_etaGap02->Draw("psame");
 GrSP_0002_v5_etaGap10->Draw("psame");
 
 
 return;
 

    /*
    TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.47,0.37,0.62,0.60);
    //myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
    myPadSetUp(myPadLogo,0,0,0,0);
    myPadLogo->Draw();
    myPadLogo->cd();
    TASImage *myAliceLogo = new TASImage("alice_logo_transparent.png");
    //.TASImage *myAliceLogo = new TASImage("alice_logo.pdf");     
    myAliceLogo->Draw();
    */
    
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

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
    currentPad->SetLeftMargin(currentLeft);
    currentPad->SetTopMargin(currentTop);
    currentPad->SetRightMargin(currentRight);
    currentPad->SetBottomMargin(currentBottom);
    return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
                  int currentMarkerStyle=20, int currentMarkerColor=0,
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
    gStyle->SetStatFontSize(0.05);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetTickLength(0.02,"y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelFont(font,"xyz"); 
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetTitleFont(font,"xyz");  
    gStyle->SetTitleOffset(1.0,"xyz");  
    gStyle->SetTitleSize(0.06,"xyz");  
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
    TString yTitle  = "v_{n}";
    
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

//============================================================================================================

// OUTPUT FILES:

// Figure 3a in "Other harmonics" paper:
// (pt dependence of harmonics v2, v3, v4 and v5 in centrality 30-40%).
// Output files (all harmonics):
//  dEta = 1.0: /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/QisQb/etaGap10/merged/outputCentrality.root
//  dEta = 0.2: /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/QisQb/etaGap02/merged/outputCentrality.root
// Objects holding measurements are following TGraphErrors:
//  1.) SP_3040ALICE_v2_etaGap10: v2{SP}(pt) in centrality 30-40%, rapidity gap = 1.0;
//  2.) SP_3040ALICE_v3_etaGap10: v3{SP}(pt) in centrality 30-40%, rapidity gap = 1.0;
//  3.) SP_3040ALICE_v4_etaGap10: v4{SP}(pt) in centrality 30-40%, rapidity gap = 1.0;
//  4.) SP_3040ALICE_v5_etaGap10: v5{SP}(pt) in centrality 30-40%, rapidity gap = 1.0;
//  5.) SP_3040ALICE_v2_etaGap02: v2{SP}(pt) in centrality 30-40%, rapidity gap = 0.2;
//  6.) SP_3040ALICE_v3_etaGap02: v3{SP}(pt) in centrality 30-40%, rapidity gap = 0.2;
//  7.) SP_3040ALICE_v4_etaGap02: v4{SP}(pt) in centrality 30-40%, rapidity gap = 0.2;
//  8.) SP_3040ALICE_v5_etaGap02: v5{SP}(pt) in centrality 30-40%, rapidity gap = 0.2.

// Figure 3b in "Other harmonics" paper:
// (pt dependence of harmonics v2, v3, v4 and v5 in centrality 00-05%).
// Output files (all harmonics):
//  dEta = 1.0: /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/QisQb/etaGap10/merged/outputCentrality.root
//  dEta = 0.2: /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/QisQb/etaGap02/merged/outputCentrality.root
// Objects holding measurements are following TGraphErrors:
//  1.) SP_0005ALICE_v2_etaGap10: v2{SP}(pt) in centrality 00-05%, rapidity gap = 1.0;
//  2.) SP_0005ALICE_v3_etaGap10: v3{SP}(pt) in centrality 00-05%, rapidity gap = 1.0;
//  3.) SP_0005ALICE_v4_etaGap10: v4{SP}(pt) in centrality 00-05%, rapidity gap = 1.0;
//  4.) SP_0005ALICE_v5_etaGap10: v5{SP}(pt) in centrality 00-05%, rapidity gap = 1.0;
//  5.) SP_0005ALICE_v2_etaGap02: v2{SP}(pt) in centrality 00-05%, rapidity gap = 0.2;
//  6.) SP_0005ALICE_v3_etaGap02: v3{SP}(pt) in centrality 00-05%, rapidity gap = 0.2;
//  7.) SP_0005ALICE_v4_etaGap02: v4{SP}(pt) in centrality 00-05%, rapidity gap = 0.2;
//  8.) SP_0005ALICE_v5_etaGap02: v5{SP}(pt) in centrality 00-05%, rapidity gap = 0.2.

// Figure 3c in "Other harmonics" paper:
// (pt dependence of harmonics v2, v3, v4, v5 in centrality 00-02%).
// (Output files - see at the end of the macro).
// Objects holding measurements are following TGraphErrors:
//  1.) SP_0002_v2_etaGap02: v2{SP}(pt) in 0-2%, eta gap = 0.2;
//  2.) SP_0002_v3_etaGap02: v3{SP}(pt) in 0-2%, eta gap = 0.2;
//  3.) SP_0002_v4_etaGap02: v4{SP}(pt) in 0-2%, eta gap = 0.2;
//  4.) SP_0002_v5_etaGap02: v5{SP}(pt) in 0-2%, eta gap = 0.2;
//  5.) SP_0002_v2_etaGap10: v2{SP}(pt) in 0-2%, eta gap = 1.0;
//  6.) SP_0002_v3_etaGap10: v3{SP}(pt) in 0-2%, eta gap = 1.0;
//  7.) SP_0002_v4_etaGap10: v4{SP}(pt) in 0-2%, eta gap = 1.0;
//  8.) SP_0002_v5_etaGap10: v5{SP}(pt) in 0-2%, eta gap = 1.0.
//  eta gap = 0.2, i.e. Qa = [-0.8,-0.1], Qb = [0.1,0.8], POIs in [-0.8,-0.1], Q = Qb     
//                  OR  Qa = [-0.8,-0.1], Qb = [0.1,0.8], POIs in [0.1,0.8], Q = Qa       
//            /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/cent02/QisQb/etaGap02/merged/outputCentrality.root
//            /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/cent02/QisQa/etaGap02/merged/outputCentrality.root
//            /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/cent02/mergedQisQa_QisQb/etaGap02/outputCentrality.root
//  eta gap = 1.0, i.e. Qa = [-0.8,-0.5], Qb = [0.5,0.8], POIs in [-0.8,-0.5], Q = Qb
//                  OR  Qa = [-0.8,-0.5], Qb = [0.5,0.8], POIs in [0.5,0.8], Q = Qa
//            /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/cent02/QisQb/etaGap10/merged/outputCentrality.root
//            /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/cent02/QisQa/etaGap10/merged/outputCentrality.root
//            /data/alice2/ab/grid/21-19/trainV0/TPConly/all/sp/cent02/mergedQisQa_QisQb/etaGap10/outputCentrality.root


static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

void fig1ab(Int_t rWrite = 0, Int_t rPerformance = 0)
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
    
    Double_t xCent[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t exCent[9] = {0.};
    
    //===================================================================================================
    // Figure 1a:
    
    //=========================================================================================================================
    // <<cos(2phi1+2phi2+2phi3-3phi4-3phi5)>>:
    Double_t xf5pCorrelator[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t yf5pCorrelator[] = {-1.885920e-10,-6.829616e-10,7.614580e-10,-1.138544e-09,1.112910e-09,
        7.147351e-08,7.625261e-08,-4.045715e-09,-8.973545e-07};
    Double_t xErrf5pCorrelator[9] = {0.};
    Double_t yErrf5pCorrelator[] = {3.193590e-10,9.788524e-10,2.210425e-09,6.072454e-09,1.346291e-08,
        2.891663e-08,6.706111e-08,2.032517e-07,1.058125e-06};
    Int_t nPointsf5pCorrelator = sizeof(xf5pCorrelator)/sizeof(Double_t);         
    TGraphErrors *f5pCorrelator = new TGraphErrors(nPointsf5pCorrelator,xf5pCorrelator,yf5pCorrelator,
                                                   xErrf5pCorrelator,yErrf5pCorrelator);
    myGraphSetUp(f5pCorrelator,1.0,kFullSquare,kBlack,1,kBlack);
    //ShiftAlongXaxis(f5pCorrelator,0.0);
    //=========================================================================================================================
    
    //=========================================================================================================================
    //  v2:
    // QC2_v2 = v2{2}:
    Double_t xQC2_v2[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t yQC2_v2[] = {0.03057544, 0.04782869, 0.06730953, 0.08693288, 0.09905272, 
        0.10507712, 0.10546323, 0.10293218, 0.10292639};
    Double_t xErrQC2_v2[9] = {0.};
    Double_t yErrQC2_v2[] = {0.00006919, 0.00008033, 0.00006915, 0.00008314, 0.00010003,
        0.00012251, 0.00015668, 0.00022385, 0.00039746};
    Int_t nPointsQC2_v2 = sizeof(xQC2_v2)/sizeof(Double_t);         
    TGraphErrors *QC2_v2 = new TGraphErrors(nPointsQC2_v2,xQC2_v2,yQC2_v2,xErrQC2_v2,yErrQC2_v2);
    myGraphSetUp(QC2_v2,1.0,kFullCircle,kBlue,1,kBlack);
    //ShiftAlongXaxis(QC2_v2,1.0);
    
    // QC4_v2 = v2{4}:
    Double_t xQC4_v2[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t yQC4_v2[] = {0.01524963,0.03617246,0.05778658,0.07579155,0.08485658,0.08614990,
        0.07996274,0.06850833,0.04132845};
    Double_t xErrQC4_v2[9] = {0.};
    Double_t yErrQC4_v2[] = {0.00090087,0.00020984,0.00014058,0.00014179,0.00017103,0.00024813,
        0.00049434,0.00157658,0.02242034};
    Int_t nPointsQC4_v2 = sizeof(xQC4_v2)/sizeof(Double_t);         
    TGraphErrors *QC4_v2 = new TGraphErrors(nPointsQC4_v2,xQC4_v2,yQC4_v2,xErrQC4_v2,yErrQC4_v2);
    myGraphSetUp(QC4_v2,1.0,kOpenSquare,kBlue,1,kBlack);
    //ShiftAlongXaxis(QC4_v2,1.0);
    //=========================================================================================================================
    
    // v2{SP}, TPC ref mult, TPConly tracks, all charges:
    Double_t xSP_TPCrefMultTPConlyAll[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t ySP_TPCrefMultTPConlyAllOrg[] = {0.027582,0.044635,0.064334,0.083715,0.095463,0.100227,0.098364,0.090479,0.078384}; 
    Double_t ySP_TPCrefMultTPConlyAll[9] = {0.};
    Double_t xErrSP_TPCrefMultTPConlyAll[9] = {0.};
    
    Double_t yErrSP_TPCrefMultTPConlyAll[] = {0.000133,0.000123,0.000097,0.000115,0.000144,0.000198,0.000317,0.000755,0.009428}; // stat errors
    //Double_t yErrSP_TPCrefMultTPConlyAll[] = {0.0014,0.0011,0.0017,0.0021,0.0023,0.0025,0.0025,0.0023,0.010}; // stat + syst nonflow
    //Double_t ySystErrSP_TPCrefMultTPConlyAll[] = {0.0014,0.0011,0.0017,0.0021,0.0023,0.0025,0.0025,0.0023,0.010}; // stat + syst nonflow
    Double_t ySystErrSP_TPCrefMultTPConlyAll[9] = {0.};
    Double_t yCorrSP[9] = {0.};
    //Double_t ySP_PP[] = {0.008, 0.0086, 0.01, 0.012, 0.015, 0.02, 0.026, 0.037, 0.056}; // Vera \delta eta=0.4 
    Double_t ySP_PP[] = {0.005321,0.005928,0.007075,0.010423,0.012069,0.013821,0.020232,0.025020,0.035993}; //Ante Hijing    
    for (Int_t p=0;p<9;p++) {
        yCorrSP[p] = TMath::Sqrt(ySP_TPCrefMultTPConlyAllOrg[p]*ySP_TPCrefMultTPConlyAllOrg[p] + ySP_PP[p]*ySP_PP[p]) -
        ySP_TPCrefMultTPConlyAllOrg[p];
        cout << "yCorrSP:" << yCorrSP[p] << endl;
        ySP_TPCrefMultTPConlyAll[p] = ySP_TPCrefMultTPConlyAllOrg[p] - yCorrSP[p];
        Double_t estsys = 0.04*ySP_TPCrefMultTPConlyAll[p];
        ySystErrSP_TPCrefMultTPConlyAll[p] = TMath::Sqrt(1.*yCorrSP[p]*1.*yCorrSP[p]+estsys*estsys);
    }
    Int_t nPointsSP_TPCrefMultTPConlyAll = sizeof(xSP_TPCrefMultTPConlyAll)/sizeof(Double_t);             
    TGraphErrors *SP_TPCrefMultTPConlyAllOrg = new TGraphErrors(nPointsSP_TPCrefMultTPConlyAll,xSP_TPCrefMultTPConlyAll,
                                                                ySP_TPCrefMultTPConlyAllOrg,xErrSP_TPCrefMultTPConlyAll,
                                                                yErrSP_TPCrefMultTPConlyAll);
    myGraphSetUp(SP_TPCrefMultTPConlyAllOrg,1.0,kOpenCircle,kRed,1,kBlack);
    
    TGraphErrors *SP_TPCrefMultTPConlyAll = new TGraphErrors(nPointsSP_TPCrefMultTPConlyAll,xSP_TPCrefMultTPConlyAll,
                                                             ySP_TPCrefMultTPConlyAll,xErrSP_TPCrefMultTPConlyAll,
                                                             yErrSP_TPCrefMultTPConlyAll);
    myGraphSetUp(SP_TPCrefMultTPConlyAll,1.0,kFullCircle,kRed,1,kBlack);                                                          
    
    TGraphErrors *SP_TPCrefMultTPConlyAll2 = new TGraphErrors(nPointsSP_TPCrefMultTPConlyAll,xSP_TPCrefMultTPConlyAll,
                                                              ySP_TPCrefMultTPConlyAll,xErrSP_TPCrefMultTPConlyAll,
                                                              ySystErrSP_TPCrefMultTPConlyAll );
    myGraphSetUp(SP_TPCrefMultTPConlyAll2,1.0,kFullCircle,kRed,1,kBlack,kRed-9,1001);                                                         
    //=========================================================================================================================
    
    //=========================================================================================================================
    //  v3:
    // QC2_v3 = v3{2}:
    Double_t xQC2_v3[] = {2.5,7.5,15.,25.,35.,45.};
    Double_t yQC2_v3[] = {0.02183955,0.02539081,0.02875255,0.03241459,0.03507416,0.03730817,
        0.03889757,0.04285879,0.05030896};
    Double_t xErrQC2_v3[9] = {0.};
    Double_t yErrQC2_v3[] = {0.00005352,0.00005932,0.00004912,0.00006050,0.00007675,0.00010535,
        0.00016134,0.00027456,0.00052779};
    Int_t nPointsQC2_v3 = sizeof(xQC2_v3)/sizeof(Double_t);         
    TGraphErrors *QC2_v3 = new TGraphErrors(nPointsQC2_v3,xQC2_v3,yQC2_v3,xErrQC2_v3,yErrQC2_v3);
    myGraphSetUp(QC2_v3,1.0,kFullTriangleUp,kBlue,1,kBlack);                                                          
    //ShiftAlongXaxis(QC2_v3,0.0);
    
    // QC2_v3_same = v3{2}, same charges:
    Double_t xQC2_v3_same[] = {2.5,7.5,15.,25.,35.,45.};
    Double_t yQC2_v3_same[] = {0.02047422,0.02403312,0.02751779,0.03108505,0.03372405,0.03572696,0.03741600};
    Double_t xErrQC2_v3_same[7] = {0.};
    Double_t yErrQC2_v3_same[] = {0.00006541,0.00007080,0.00005897,0.00007455,0.00009984,0.00014546,0.00023698};
    Int_t nPointsQC2_v3_same = sizeof(xQC2_v3_same)/sizeof(Double_t);         
    TGraphErrors *QC2_v3_same = new TGraphErrors(nPointsQC2_v3_same,xQC2_v3_same,yQC2_v3_same,xErrQC2_v3_same,yErrQC2_v3_same);
    myGraphSetUp(QC2_v3_same,1.0,kStar,kBlue,1,kBlack);                                                           
    //ShiftAlongXaxis(QC2_v3_same,1.);
    
    // QC4_v3 = v3{4}:
    Double_t xQC4_v3[] = {2.5,7.5,15.,25.,35.,45.};
    Double_t yQC4_v3[] = {0.01053669,0.01257961,0.01552626,0.01506058,0.01839253,0.01770513,
        0.02810784,0.03316069};
    Double_t xErrQC4_v3[8] = {0.};
    Double_t yErrQC4_v3[] = {0.00101506,0.00115529,0.00069718,0.00148116,0.00148018,0.00352904,
        0.00763725,0.02345570};
    Int_t nPointsQC4_v3 = sizeof(xQC4_v3)/sizeof(Double_t);         
    TGraphErrors *QC4_v3 = new TGraphErrors(6,xQC4_v3,yQC4_v3,xErrQC4_v3,yErrQC4_v3);
    myGraphSetUp(QC4_v3,1.0,kOpenSquare,kBlue,1,kBlack);                                                            
    //ShiftAlongXaxis(QC4_v3,0.75);
    
    // QC4_v3_same = v3{4}, same charges:
    Double_t xQC4_v3_same[] = {2.5,7.5,15.,25.,35.,45.};
    Double_t yQC4_v3_same[] = {0.00760270,0.01289241,0.01396290,0.01593748,0.01841474,0.01355171};
    Double_t xErrQC4_v3_same[6] = {0.};
    Double_t yErrQC4_v3_same[] = {0.00433070,0.00180779,0.00151872,0.00209834,0.00279355,0.01648787};
    Int_t nPointsQC4_v3_same = sizeof(xQC4_v3_same)/sizeof(Double_t);         
    TGraphErrors *QC4_v3_same = new TGraphErrors(nPointsQC4_v3_same,xQC4_v3_same,yQC4_v3_same,xErrQC4_v3_same,yErrQC4_v3_same);
    myGraphSetUp(QC4_v3_same,1.0,kOpenSquare,kRed,1,kBlack);                                                            
    
    // ZDC_v3 = v3{ZDC} (Ilya):
    Double_t xZDC_v3[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t yZDC_v3[] = {0.2005414,-0.0008304778,0.0003468132,0.001181955,-0.0001876906,0.0005668259,-2.986464e-05,0.0007950359,0.006499858};
    Double_t xErrZDC_v3[9] = {0.};
    Double_t yErrZDC_v3[] = {0.6825955,0.002633532,0.0006703818,0.0005074313,0.0006300845,0.0009810903,0.001860773,0.004609001,0.01633596};
    Int_t nPointsZDC_v3 = sizeof(xZDC_v3)/sizeof(Double_t);         
    TGraphErrors *ZDC_v3 = new TGraphErrors(6,xZDC_v3,yZDC_v3,xErrZDC_v3,yErrZDC_v3);
    myGraphSetUp(ZDC_v3,1.0,kFullCircle,kGreen+2,1,kBlack);                                                             
    
    /*
     // ZDC_v3 = v3{ZDC} (Ilya):
     Double_t xZDC_v3[] = {2.5,7.5,15.,25.,35.,45.};
     Double_t yZDC_v3[] = {0.004703597,2.977193e-05,7.509025e-05,8.895956e-05,0.0001218988,
     5.045662e-06,0.0003617772,0.001024175,-0.0009996569};
     Double_t xErrZDC_v3[9] = {0.};
     Double_t yErrZDC_v3[] = {0.002843094,0.0003961917,9.486118e-05,7.382399e-05,9.159939e-05,
     0.0001466929,0.0002849215,0.0007678118,0.003465045};
     Int_t nPointsZDC_v3 = sizeof(xZDC_v3)/sizeof(Double_t);         
     TGraphErrors *ZDC_v3 = new TGraphErrors(nPointsZDC_v3,xZDC_v3,yZDC_v3,xErrZDC_v3,yErrZDC_v3);
     myGraphSetUp(ZDC_v3,1.0,kOpenTriangleUp,kGreen+2,1,kBlack);                                                               
     */
    
    // SP_v3_etaGap02 = v3{SP}:
    Double_t xSP_v3_etaGap02[] = {2.5,7.5,15.,25.,35.,45.};
    Double_t ySP_v3_etaGap02[] = {0.02100282,0.02441799,0.02742384,0.03061855,0.03244124,0.03295708,0.03197919};
    Double_t xErrSP_v3_etaGap02[9] = {0.};
    Double_t yErrSP_v3_etaGap02[] = {0.00007620,0.00008169,0.00006829,0.00009012,0.00011868,0.00018015,0.00031468};
    Int_t nPointsSP_v3_etaGap02 = sizeof(xSP_v3_etaGap02)/sizeof(Double_t);         
    TGraphErrors *SP_v3_etaGap02 = new TGraphErrors(nPointsSP_v3_etaGap02,xSP_v3_etaGap02,ySP_v3_etaGap02,
                                                    xErrSP_v3_etaGap02,yErrSP_v3_etaGap02);
    myGraphSetUp(SP_v3_etaGap02,1.0,kOpenTriangleUp,kBlue,1,kBlack);                                                               
    //ShiftAlongXaxis(SP_v3_etaGap02,-0.5);  
    
    // SP_v3_etaGap10 = v3{SP}:
    Double_t xSP_v3_etaGap10[] = {2.5,7.5,15.,25.,35.,45.};
    Double_t ySP_v3_etaGap10[] = {0.02070702,0.02382848,0.02672520,0.02970215,0.03109479,0.03105698,0.02908181};
    Double_t xErrSP_v3_etaGap10[7] = {0.};
    Double_t yErrSP_v3_etaGap10[] = {0.00013628,0.00014867,0.00012609,0.00016659,0.00024023,0.00038862,0.00073425};
    Int_t nPointsSP_v3_etaGap10 = sizeof(xSP_v3_etaGap10)/sizeof(Double_t);  
    TGraphErrors *SP_v3_etaGap10 = new TGraphErrors(nPointsSP_v3_etaGap10,xSP_v3_etaGap10,ySP_v3_etaGap10,
                                                    xErrSP_v3_etaGap10,yErrSP_v3_etaGap10);
    myGraphSetUp(SP_v3_etaGap10,1.2,kOpenSquare,kBlue,1,kBlack);                                                         
    //ShiftAlongXaxis(SP_v3_etaGap10,0.5);  
    
    Double_t ySystErrSPv3[7] = {0.};
    Double_t yCorrSPv3[7] = {0.};
    Double_t ySPv3_PP[] = {0.003,0.003,0.004,0.005,0.006,0.008,0.01}; //Vera pp
    
    for (Int_t p=0;p<6;p++) {
        Double_t corr = TMath::Sqrt(ySP_v3_etaGap10[p]*ySP_v3_etaGap10[p] + ySPv3_PP[p]*ySPv3_PP[p]) - ySP_v3_etaGap10[p];
        cout << "Corr v3:" << corr << endl;
        yCorrSPv3[p] = ySP_v3_etaGap10[p] - corr;
        Double_t estsys = 0.05*yCorrSPv3[p];
        ySystErrSPv3[p] = TMath::Sqrt(1.*corr*1.*corr+estsys*estsys);
    }
    TGraphErrors *SP_v3_corr = new TGraphErrors(nPointsSP_v3_etaGap10,xSP_v3_etaGap10,yCorrSPv3,
                                                xErrSP_v3_etaGap10,yErrSP_v3_etaGap10);
    myGraphSetUp(SP_v3_corr,1.,kFullSquare,kBlue,1,kBlack);                                                         
    
    TGraphErrors *SP_v3_corr2 = new TGraphErrors(nPointsSP_v3_etaGap10,xSP_v3_etaGap10,yCorrSPv3,
                                                 xErrSP_v3_etaGap10,ySystErrSPv3);
    myGraphSetUp(SP_v3_corr2,1.,kFullSquare,kBlue-9,1,kBlue,kBlue-9,1001);                                                         
    
    // SP_v4_etaGap10 = v4{SP}:
    Double_t xSP_v4_etaGap10[] = {2.5,7.5,15.,25.,35.,45.};
    Double_t ySP_v4_etaGap10[] = {0.00981464,0.01161192,0.01288244,0.01435392,0.01599797,0.01553851};
    Double_t xErrSP_v4_etaGap10[6] = {0.};
    Double_t yErrSP_v4_etaGap10[] = {0.00025098,0.00026355,0.00022703,0.00030744,0.00042727,0.00073090};
    Int_t nPointsSP_v4_etaGap10 = sizeof(xSP_v4_etaGap10)/sizeof(Double_t);  
    TGraphErrors *SP_v4_etaGap10 = new TGraphErrors(nPointsSP_v4_etaGap10,xSP_v4_etaGap10,ySP_v4_etaGap10,
                                                    xErrSP_v4_etaGap10,yErrSP_v4_etaGap10);
    myGraphSetUp(SP_v4_etaGap10,1.2,kStar,kMagenta,1,kBlack);                                                         
    ShiftAlongXaxis(SP_v4_etaGap10,1.5);  
    
    
    
    Double_t ySystErrSPv4[7] = {0.};
    Double_t yCorrSPv4[7] = {0.};
    Double_t ySPv4_PP[] = {0.0025,0.0025,0.0025,0.0025,0.0025,0.0025, 0.0025};
    
    for (Int_t p=0;p<6;p++) {
        Double_t corr = TMath::Sqrt(ySP_v4_etaGap10[p]*ySP_v4_etaGap10[p] + ySPv4_PP[p]*ySPv4_PP[p]) - ySP_v4_etaGap10[p];
        cout << "Corr v4:" << corr << endl;
        yCorrSPv4[p] = ySP_v4_etaGap10[p] - corr;
        Double_t estsys = 0.05*yCorrSPv4[p];
        ySystErrSPv4[p] = TMath::Sqrt(1.*corr*1.*corr+estsys*estsys);
    }
    TGraphErrors *SP_v4_corr = new TGraphErrors(nPointsSP_v4_etaGap10,xSP_v4_etaGap10,yCorrSPv4,
                                                xErrSP_v4_etaGap10,yErrSP_v4_etaGap10);
    myGraphSetUp(SP_v4_corr,1.2,kStar,kMagenta,1,kBlack);                                                         
    
    TGraphErrors *SP_v4_corr2 = new TGraphErrors(nPointsSP_v4_etaGap10,xSP_v4_etaGap10,yCorrSPv4,
                                                 xErrSP_v4_etaGap10,ySystErrSPv4);
    myGraphSetUp(SP_v4_corr2,1.2,kStar,kMagenta-9,1,kMagenta,kMagenta-9,1001);        
    //=========================================================================================================================   
    
    //=========================================================================================================================   
    // v3 centrality dependence prediction by Luzum and Ollitrault:
    Double_t xLuzOll[] = {2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5};
    Double_t yLuzOll[] = {0.023931343405028907,0.025578676086755334,0.02692673087997464,0.0277496106070913,
        0.02834431352143814,0.030099333306694308,0.029712119761886445,0.02910588380400307,
        0.029029451386236353,0.027650687387294634,0.02479761160079346,0.020694594742745584,
        0.014559252721547645,0.00687958280855717,-0.0011508459086818318};
    Double_t xErrLuzOll[15] = {0.};
    Double_t yErrLuzOll[15] = {0.};
    Int_t nPointsLuzOll = sizeof(xLuzOll)/sizeof(Double_t);         
    TGraphErrors *LuzOll = new TGraphErrors(nPointsLuzOll,xLuzOll,yLuzOll,xErrLuzOll,yErrLuzOll);
    myGraphSetUp(LuzOll,1.,kFullSquare,kBlack,3,kBlack);                                                         
    //=========================================================================================================================      
    
    //=========================================================================================================================
    // v3{5}:
    Double_t factor = 100.; 
    Double_t xQC5_v3_squared[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
    Double_t yQC5_v3_squared[9] = {0.};
    Double_t xErrQC5_v3_squared[9] = {0.};
    Double_t yErrQC5_v3_squared[9] = {0.};
    for(Int_t p=0;p<6;p++) 
    { 
        Double_t value = 0.; // v3{5}
        Double_t error = 0.; // v3{5} error
        Double_t mh = yf5pCorrelator[p]; 
        Double_t mhError = yErrf5pCorrelator[p];
        Double_t v24 = yQC4_v2[p];
        Double_t v24Error = yErrQC4_v2[p];
        Double_t v22 = yQC2_v2[p];
        Double_t v22Error = yErrQC2_v2[p];
        if(TMath::Abs(v24)> 0.) 
        {
            value = mh/pow(v22,3.); // v3{5}^2
            error = pow((1./pow(v22,6.))*pow(mhError,2.)+9.*(pow(mh,2.)/pow(v22,8.))*pow(v22Error,2.),0.5);
            yQC5_v3_squared[p] = value*factor;
            yErrQC5_v3_squared[p] = error*factor;
        }
    } // end of for(Int_t p=0;p<9;p++)
    Int_t nPointsQC5_v3_squared = sizeof(xQC5_v3_squared)/sizeof(Double_t);         
    TGraphErrors *QC5_v3_squared = new TGraphErrors(6,xQC5_v3_squared,yQC5_v3_squared,
                                                    xErrQC5_v3_squared,yErrQC5_v3_squared);
    myGraphSetUp(QC5_v3_squared,1.4,33,kBlack,1.,kBlack);                                                         
    ShiftAlongXaxis(QC5_v3_squared,-0.5);
    //=========================================================================================================================
    // v2 epsilons  
    Double_t r2E2Cum2WNpart[]= {0.0813921,0.122575,0.192946,0.276331,0.351085,0.421993,0.49293,0.573383,0.661991,0.709309,0.926234};
    Double_t r2E2Cum2CGCNpart[]= {0.1037700, 0.176570, 0.270957, 0.365345, 0.437808, 0.496254, 0.546311, 0.596211, 0.645011};
    Double_t xRatio_v22_e2CGC[9] = {0.};
    Double_t yRatio_v22_e2CGC[9] = {0.};
    Double_t xErrRatio_v22_e2CGC[9] = {0.};
    Double_t yErrRatio_v22_e2CGC[9] = {0.};
    Double_t ySystErrRatio_v22_e2CGC[9] = {0.};
    Double_t xRatio_v22_e2w[9] = {0.};
    Double_t yRatio_v22_e2w[9] = {0.};
    Double_t xErrRatio_v22_e2w[9] = {0.};
    Double_t yErrRatio_v22_e2w[9] = {0.};
    Double_t ySystErrRatio_v22_e2w[9] = {0.};
    
    for(Int_t p=0;p<9;p++) {
        xRatio_v22_e2CGC[p] = xSP_TPCrefMultTPConlyAll[p];
        yRatio_v22_e2CGC[p] = ySP_TPCrefMultTPConlyAll[p]/r2E2Cum2CGCNpart[p];
        cout << "v2/epslon CGC: " << yRatio_v22_e2CGC[p] << endl;
        yErrRatio_v22_e2CGC[p] = yErrSP_TPCrefMultTPConlyAll[p]/r2E2Cum2CGCNpart[p];
        ySystErrRatio_v22_e2CGC[p] = ySystErrSP_TPCrefMultTPConlyAll[p]/r2E2Cum2CGCNpart[p];
        
        
        xRatio_v22_e2w[p] = xSP_TPCrefMultTPConlyAll[p];
        yRatio_v22_e2w[p] = ySP_TPCrefMultTPConlyAll[p]/r2E2Cum2WNpart[p];
        yErrRatio_v22_e2w[p] = yErrSP_TPCrefMultTPConlyAll[p]/r2E2Cum2WNpart[p];
        ySystErrRatio_v22_e2w[p] = ySystErrSP_TPCrefMultTPConlyAll[p]/r2E2Cum2WNpart[p];
    }    
    TGraphErrors *Ratio_v22_e2CGC = new TGraphErrors(9,xRatio_v22_e2CGC,yRatio_v22_e2CGC,
                                                     xErrRatio_v22_e2CGC,yErrRatio_v22_e2CGC);
    myGraphSetUp(Ratio_v22_e2CGC,1.,kFullCircle,kRed,1.,kRed);                                                           
    
    TGraphErrors *Ratio_v22_e2w = new TGraphErrors(9,xRatio_v22_e2w,yRatio_v22_e2w,
                                                   xErrRatio_v22_e2w,yErrRatio_v22_e2w);
    myGraphSetUp(Ratio_v22_e2w,1.,kOpenCircle,kMagenta,1.,kMagenta);                                                           
    
    TGraphErrors *Ratio_v22_e2CGC2 = new TGraphErrors(9,xRatio_v22_e2CGC,yRatio_v22_e2CGC,
                                                      xErrRatio_v22_e2CGC,ySystErrRatio_v22_e2CGC);
    myGraphSetUp(Ratio_v22_e2CGC2,1.,kFullCircle,kRed-9,1.,kRed,kRed-9,1001);                                                             
    
    TGraphErrors *Ratio_v22_e2w2 = new TGraphErrors(9,xRatio_v22_e2w,yRatio_v22_e2w,
                                                    xErrRatio_v22_e2w,ySystErrRatio_v22_e2w);
    myGraphSetUp(Ratio_v22_e2w2,1.,kOpenCircle,kMagenta,1.,kMagenta,kMagenta-9,1001);                                                     
    
    //=========================================================================================================================
    
    // v3 epsilons
    Double_t r2E3Cum2CGCNpart[]= {0.0794007, 0.0954744, 0.115036, 0.141516, 0.173228, 0.212523, 0.259254, 0.302282, 0.309824};
    Double_t r2E3Cum2WNpart[]= {0.0798219,0.104613,0.135106,0.173153,0.214027,0.259469,0.308709,0.354488,0.38321,0.39706,0.377804};   
    Double_t xRatio_v32_e2CGC[7] = {0.};
    Double_t yRatio_v32_e2CGC[7] = {0.};
    Double_t xErrRatio_v32_e2CGC[7] = {0.};
    Double_t yErrRatio_v32_e2CGC[7] = {0.};
    Double_t ySystErrRatio_v32_e2CGC[7] = {0.};     
    Double_t xRatio_v32_e2w[7] = {0.};
    Double_t yRatio_v32_e2w[7] = {0.};
    Double_t xErrRatio_v32_e2w[7] = {0.};
    Double_t yErrRatio_v32_e2w[7] = {0.};
    Double_t ySystErrRatio_v32_e2w[7] = {0.};
    for(Int_t p=0;p<6;p++) {
        xRatio_v32_e2CGC[p] = xSP_v3_etaGap10[p];
        yRatio_v32_e2CGC[p] = ySP_v3_etaGap10[p]/r2E3Cum2CGCNpart[p];
        yErrRatio_v32_e2CGC[p] = yErrSP_v3_etaGap10[p]/r2E3Cum2CGCNpart[p];
        ySystErrRatio_v32_e2CGC[p] = ySystErrSPv3[p]/r2E3Cum2CGCNpart[p];
        
        xRatio_v32_e2w[p] = xSP_v3_etaGap10[p];
        yRatio_v32_e2w[p] = ySP_v3_etaGap10[p]/r2E3Cum2WNpart[p];
        yErrRatio_v32_e2w[p] = yErrSP_v3_etaGap10[p]/r2E3Cum2WNpart[p];
        ySystErrRatio_v32_e2w[p] = ySystErrSPv3[p]/r2E3Cum2WNpart[p];
        
    }
    
    TGraphErrors *Ratio_v32_e2CGC = new TGraphErrors(6,xRatio_v32_e2CGC,yRatio_v32_e2CGC,
                                                     xErrRatio_v32_e2CGC,yErrRatio_v32_e2CGC);
    myGraphSetUp(Ratio_v32_e2CGC,1.,kFullSquare,kBlue,1.,kBlue);                                                    
    
    TGraphErrors *Ratio_v32_e2w = new TGraphErrors(6,xRatio_v32_e2w,yRatio_v32_e2w,
                                                   xErrRatio_v32_e2w,yErrRatio_v32_e2w);
    myGraphSetUp(Ratio_v32_e2w,1.,kOpenSquare,kMagenta,1.,kMagenta);                                                    
    
    
    TGraphErrors *Ratio_v32_e2CGC2 = new TGraphErrors(6,xRatio_v32_e2CGC,yRatio_v32_e2CGC,
                                                      xErrRatio_v32_e2CGC,ySystErrRatio_v32_e2CGC);
    myGraphSetUp(Ratio_v32_e2CGC2,1.,kFullSquare,kBlue-9,1,kBlue,kBlue-9,1001);                                                      
    
    TGraphErrors *Ratio_v32_e2w2 = new TGraphErrors(6,xRatio_v32_e2w,yRatio_v32_e2w,
                                                    xErrRatio_v32_e2w,ySystErrRatio_v32_e2w);
    myGraphSetUp(Ratio_v32_e2w2,1.,kOpenSquare,kMagenta,1,kMagenta,kMagenta-9,1001);                                                      
    //=========================================================================================================================
    
    // Style histogramm:
    TH1F *myBlankHisto = new TH1F("myBlankHisto","Blank Histogram",100,0.,100.);
    //myBlankHisto->SetXTitle("p_{t} (GeV/#font[72]{c})");
    myBlankHisto->SetNdivisions(505,"y");
    myBlankHisto->GetYaxis()->SetTitleOffset(0.6);
    //myBlankHisto->GetYaxis()->SetLabelSize(0.05);
    //myBlankHisto->GetYaxis()->SetTitleSize(0.07);
    //myBlankHisto->GetYaxis()->SetRangeUser(-0.02101,0.399);
    myBlankHisto->GetXaxis()->SetRangeUser(0.,81.);    
    // Main canvas:   
    TCanvas *myCan = new TCanvas("myCan",cStamp1,500,750);
    myCan->Divide(1,2,0.0,0.0);   
    myCan->Draw();
    
    // Figure 1a: 
    myCan->cd(1);   
    gPad->SetBottomMargin(0.); 
    gPad->SetRightMargin(0.001); 
    TH1D* styleHistPad1 = myBlankHisto->Clone("styleHistPad1");
    styleHistPad1->GetYaxis()->SetRangeUser(-0.005,0.1099);  
    styleHistPad1->GetYaxis()->SetTitle("v_{n}");
    styleHistPad1->Draw();
    TLatex *a = new TLatex(0.95,0.93,"(a)");
    a->SetNDC();
    a->SetTextSize(0.05);
    a->Draw();
    // Legend:    
    TLegend *myLegendPad1 = new TLegend(0.45,0.44,0.67,0.74);
    myLegendSetUp(myLegendPad1,0.04); 
    //myLegendPad1->AddEntry(f5pCorrelator,"#propto v_{2}^{3} v_{3}^{2} ??","p"); 
    //myLegendPad1->AddEntry(QC2_v2,"v_{2}{2}","p");
    //myLegendPad1->AddEntry(QC4_v2,"v_{2}{4}","p");
    //myLegendPad1->AddEntry(QC2_v3,"v_{3}{2}","p");
    //myLegendPad1->AddEntry(QC2_v3_same,"v_{3}{2} (same charges)","p");
    //myLegendPad1->AddEntry(SP_v3_etaGap02,"v_{3}{SP} (#Delta#eta = 0.2)","p");
    myLegendPad1->AddEntry(SP_TPCrefMultTPConlyAll,"v_{2}{2, |#Delta#eta| > 1}","p");
    myLegendPad1->AddEntry(SP_v3_corr,"v_{3}{2, |#Delta#eta| > 1}","p");
    myLegendPad1->AddEntry(SP_v4_etaGap10,"v_{4}{2, |#Delta#eta| > 1}","p");
    myLegendPad1->AddEntry(QC4_v3,"v_{3}{4}","p");
    //myLegendPad1->AddEntry(QC4_v3_same,"v_{3}{4} (same charges)","p");
    myLegendPad1->AddEntry(ZDC_v3,"v_{3/#Psi_{RP}}","p");
    myLegendPad1->AddEntry(QC5_v3_squared,"100 #times v_{3/#Psi_{2}}^{2}","p");  
    myLegendPad1->Draw("same"); 
    // Final drawing (order is important!):
    //f5pCorrelator->Draw("psame");
    //QC2_v2->Draw("psame");
    //QC4_v2->Draw("psame");
    //QC2_v3->Draw("psame"); 
    //QC2_v3_same->Draw("psame"); 
    //SP_v3_etaGap02->Draw("psame");    
    SP_TPCrefMultTPConlyAll2->Draw("same3");
    SP_TPCrefMultTPConlyAll->Draw("lsame");
    SP_TPCrefMultTPConlyAll->Draw("psame");   
    //SP_v3_etaGap10->Draw("psame"); 
    //SP_v4_etaGap10->Draw("psame"); 
    SP_v4_corr2->Draw("same3");
    SP_v4_corr->Draw("lsame");
    SP_v4_corr->Draw("psame");
    SP_v3_corr2->Draw("same3");
    SP_v3_corr->Draw("lsame");
    SP_v3_corr->Draw("psame");
    QC4_v3->Draw("psame");
    //QC4_v3_same->Draw("psame");
    QC5_v3_squared->Draw("psame");
    ZDC_v3->Draw("psame");
    //LuzOll->Draw("lsame"); 
    
    // Figure 1b: 
    myCan->cd(2);   
    gPad->SetTopMargin(0.); 
    gPad->SetBottomMargin(0.15); 
    gPad->SetRightMargin(0.001); 
    TH1D* styleHistPad2 = myBlankHisto->Clone("styleHistPad2");
    styleHistPad2->GetYaxis()->SetRangeUser(0.,0.4799);  
    styleHistPad2->GetXaxis()->SetTitle("centrality percentile");  
    styleHistPad2->GetYaxis()->SetTitle("v_{n}/#varepsilon_{n}");   
    styleHistPad2->Draw();
    TLatex *b = new TLatex(0.95,0.93,"(b)");
    b->SetNDC();
    b->SetTextSize(0.05);
    b->Draw();
    TLegend *myLegendPad2 = new TLegend(0.6,0.64,0.86,0.94);
    myLegendSetUp(myLegendPad2,0.04); 
    myLegendPad2->AddEntry(Ratio_v22_e2CGC,"v_{2}{2, |#Delta#eta| > 1}/#varepsilon_{2}^{CGC}{2}","p");
    myLegendPad2->AddEntry(Ratio_v32_e2CGC,"v_{3}{2, |#Delta#eta| > 1}/#varepsilon_{3}^{CGC}{2}","p");
    myLegendPad2->AddEntry(Ratio_v22_e2w,"v_{2}{2, |#Delta#eta| > 1}/#varepsilon_{2}^{W}{2}","p");
    myLegendPad2->AddEntry(Ratio_v32_e2w,"v_{3}{2, |#Delta#eta| > 1}/#varepsilon_{3}^{W}{2}","p");
    //myLegendPad2->AddEntry(Ratio_v32_e3b,"v_{3}{2}/#varepsilon_{3,Ncoll}{2}","l"); 
    //myLegendPad2->AddEntry(Ratio_v34_e2b,"v_{3}{4}/#varepsilon_{3,Ncoll}{4}","l"); 
    //myLegendPad2->AddEntry(Ratio_v34_e2w,"v_{3}{4}/#varepsilon_{3}{4}","p"); 
    //myLegendPad2->AddEntry(Ratio_v32_e3w,"v_{3}{2}/#varepsilon_{3}{2}","p"); 
    myLegendPad2->Draw("same");
    //Ratio_v34_e2b->Draw("same3"); 
    //Ratio_v34_e2b->Draw("lsameX"); 
    //Ratio_v32_e3b->Draw("same3"); 
    //Ratio_v32_e3b->Draw("lsameX"); 
    //Ratio_v34_e2w->Draw("psame"); 
    //Ratio_v32_e3w->Draw("psame"); 
    
    
    
    
    Ratio_v22_e2CGC2->Draw("same3");
    Ratio_v22_e2w2->Draw("same3");
    Ratio_v32_e2CGC2->Draw("same3");
    Ratio_v32_e2w2->Draw("same3");   
    Ratio_v32_e2CGC->Draw("lsameX");
    Ratio_v22_e2CGC->Draw("lsameX");
    Ratio_v32_e2w->Draw("lsameX");
    Ratio_v22_e2w->Draw("lsameX");   
    Ratio_v32_e2CGC->Draw("psame");
    Ratio_v22_e2CGC->Draw("psame");
    Ratio_v32_e2w->Draw("psame");
    Ratio_v22_e2w->Draw("psame"); 
    
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
    currentGraph->SetLineWidth(2);
    currentGraph->SetLineColor(currentLineColor);
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

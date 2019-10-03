void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();


// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

Bool_t ifV5 = kFALSE;
Bool_t ifV6 = kFALSE;
Bool_t ifScale = kFALSE;

void Fig1() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

    
    const Double_t sys_Run2Run1_v22Gap10ratio = 0.01;
    //const Double_t sys_Run2Run1_v22Gap10ratio[] = {0.035, 0.035, 0.021, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    const Double_t sys_Run2Run1_v24Gap10ratio = 0.012;
    const Double_t sys_Run2Run1_v32Gap10ratio = 0.019;
    const Double_t sys_Run2Run1_v42Gap10ratio = 0.017;
    
    Double_t xCross[] = {2.5,7.5, 15, 25,35,45,55,65,75,85,95};
    Double_t xCrossErr[100] = {0.};
    Double_t xCrossSys[100] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    Double_t xCent[] = {0,5,10,20,30,40,50,60,70,80};
    Double_t exxx[50]={0.};
    Double_t esysxxx[]={0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02};
    
    // hydro from JYO, arXiv:1511.06289
    Double_t v22Gap10_hydro[10]={0.};
    Double_t v22Gap10Err_hydro[10]={0.};
    Double_t v32Gap10_hydro[10]={0.};
    Double_t v32Gap10Err_hydro[10]={0.};
    
    Double_t v22Gap10max_hydro[]={0.027339, 0.044692, 0.064734, 0.084598, 0.097374, 0.103430, 0.102025, 0.094148, 0.078079};
    Double_t v22Gap10min_hydro[]={0.026442, 0.043232, 0.062584, 0.081605, 0.093560, 0.098823, 0.096770, 0.088499, 0.072614};
    
    
    Double_t v32Gap10max_hydro[]={0.021506, 0.024880, 0.028030, 0.031591, 0.033745, 0.034303};
    Double_t v32Gap10min_hydro[]={0.019993, 0.023166, 0.026086, 0.029221, 0.030830, 0.030764};
    
    Double_t vn2Gap10Err_hydro[10]={0.};
    
    
    for (Int_t j=0; j<9; j++) {
        v22Gap10_hydro[j] = 0.5*(v22Gap10max_hydro[j]+v22Gap10min_hydro[j]);
        v22Gap10Err_hydro[j] = 0.5*(v22Gap10max_hydro[j]-v22Gap10min_hydro[j]);
    }
    
    TGraphErrors *Graphv22Gap10_hydro = new TGraphErrors(9,xCross,v22Gap10_hydro,xCrossErr,v22Gap10Err_hydro);
    Graphv22Gap10_hydro->SetMarkerStyle(24);
    Graphv22Gap10_hydro->SetMarkerColor(kRed);
    Graphv22Gap10_hydro->SetMarkerSize(1.);
    Graphv22Gap10_hydro->SetLineColor(kRed+2);
    Graphv22Gap10_hydro->SetLineWidth(1);
    Graphv22Gap10_hydro->SetFillStyle(0);
    Graphv22Gap10_hydro->SetFillColor(kRed+2);
    
    TGraphErrors *Graphv22Gap10Min_hydro = new TGraphErrors(9,xCross,v22Gap10min_hydro,xCrossErr,vn2Gap10Err_hydro);
    Graphv22Gap10Min_hydro->SetMarkerStyle(24);
    Graphv22Gap10Min_hydro->SetMarkerColor(kRed);
    Graphv22Gap10Min_hydro->SetMarkerSize(1.);
    Graphv22Gap10Min_hydro->SetLineColor(kRed+2);
    Graphv22Gap10Min_hydro->SetLineStyle(1);
    Graphv22Gap10Min_hydro->SetLineWidth(1);
    Graphv22Gap10Min_hydro->SetFillStyle(0);
    Graphv22Gap10Min_hydro->SetFillColor(kRed+2);
    
    TGraphErrors *Graphv22Gap10Max_hydro = new TGraphErrors(9,xCross,v22Gap10max_hydro,xCrossErr,vn2Gap10Err_hydro);
    Graphv22Gap10Max_hydro->SetMarkerStyle(24);
    Graphv22Gap10Max_hydro->SetMarkerColor(kRed);
    Graphv22Gap10Max_hydro->SetMarkerSize(1.);
    Graphv22Gap10Max_hydro->SetLineColor(kRed+2);
    Graphv22Gap10Max_hydro->SetLineStyle(1);
    Graphv22Gap10Max_hydro->SetLineWidth(1);
    Graphv22Gap10Max_hydro->SetFillStyle(0);
    Graphv22Gap10Max_hydro->SetFillColor(kRed+2);
    
    
    for (Int_t j=0; j<6; j++) {
        v32Gap10_hydro[j] = 0.5*(v32Gap10max_hydro[j]+v32Gap10min_hydro[j]);
        v32Gap10Err_hydro[j] = 0.5*(v32Gap10max_hydro[j]-v32Gap10min_hydro[j]);
    }
    TGraphErrors *Graphv32Gap10_hydro = new TGraphErrors(6,xCross,v32Gap10_hydro,xCrossErr,v32Gap10Err_hydro);
    Graphv32Gap10_hydro->SetMarkerStyle(24);
    Graphv32Gap10_hydro->SetMarkerColor(kBlue);
    Graphv32Gap10_hydro->SetMarkerSize(1.);
    Graphv32Gap10_hydro->SetLineColor(kBlue);
    Graphv32Gap10_hydro->SetLineWidth(1);
    Graphv32Gap10_hydro->SetFillStyle(0);
    Graphv32Gap10_hydro->SetFillColor(kBlue);
    
    TGraphErrors *Graphv32Gap10Min_hydro = new TGraphErrors(6,xCross,v32Gap10min_hydro,xCrossErr,vn2Gap10Err_hydro);
    Graphv32Gap10Min_hydro->SetMarkerStyle(24);
    Graphv32Gap10Min_hydro->SetMarkerColor(kBlue);
    Graphv32Gap10Min_hydro->SetMarkerSize(1.);
    Graphv32Gap10Min_hydro->SetLineColor(kBlue+2);
    Graphv32Gap10Min_hydro->SetLineStyle(1);
    Graphv32Gap10Min_hydro->SetLineWidth(1);
    Graphv32Gap10Min_hydro->SetFillStyle(0);
    Graphv32Gap10Min_hydro->SetFillColor(kBlue+2);
    
    TGraphErrors *Graphv32Gap10Max_hydro = new TGraphErrors(6,xCross,v32Gap10max_hydro,xCrossErr,vn2Gap10Err_hydro);
    Graphv32Gap10Max_hydro->SetMarkerStyle(24);
    Graphv32Gap10Max_hydro->SetMarkerColor(kBlue);
    Graphv32Gap10Max_hydro->SetMarkerSize(1.);
    Graphv32Gap10Max_hydro->SetLineColor(kBlue+2);
    Graphv32Gap10Max_hydro->SetLineStyle(1);
    Graphv32Gap10Max_hydro->SetLineWidth(1);
    Graphv32Gap10Max_hydro->SetFillStyle(0);
    Graphv32Gap10Max_hydro->SetFillColor(kBlue+2);
    
    Double_t xHydro[]={5, 15, 30, 50};
    Double_t xErrHydro[4]={0.};
    
    Double_t ratiov2_etaS020[] = {1.013, 1.024, 1.037, 1.035};
    Double_t ratiov2Err_etaS020[] = {0.015, 0.011, 0.007, 0.008};
    
    Double_t ratiov3_etaS020[] = {1.046, 1.064, 1.051, 1.051};
    Double_t ratiov3Err_etaS020[] = {0.015, 0.016, 0.011, 0.011};
    
    Double_t ratiov4_etaS020[] = {1.057, 1.106, 1.116, 1.127};
    Double_t ratiov4Err_etaS020[] = {0.016, 0.016, 0.012, 0.013};
    
    TGraphErrors *gRatioV22_etaS020 = new TGraphErrors(4,xHydro,ratiov2_etaS020,xErrHydro,ratiov2Err_etaS020);
    gRatioV22_etaS020->SetMarkerStyle(3);
    gRatioV22_etaS020->SetMarkerColor(kBlack);
    gRatioV22_etaS020->SetLineColor(0);
    gRatioV22_etaS020->SetMarkerSize(1.5);
    gRatioV22_etaS020->SetFillStyle(3008);
    gRatioV22_etaS020->SetFillColor(kRed-9);
    
    TGraphErrors *gRatioV32_etaS020 = new TGraphErrors(4,xHydro,ratiov3_etaS020,xErrHydro,ratiov3Err_etaS020);
    gRatioV32_etaS020->SetMarkerStyle(3);
    gRatioV32_etaS020->SetMarkerColor(kBlack);
    gRatioV32_etaS020->SetLineColor(kBlack);
    gRatioV32_etaS020->SetMarkerSize(1.5);
    gRatioV32_etaS020->SetFillStyle(3008);
    gRatioV32_etaS020->SetFillColor(kBlue-9);

    TGraphErrors *gRatioV42_etaS020 = new TGraphErrors(4,xHydro,ratiov4_etaS020,xErrHydro,ratiov4Err_etaS020);
    gRatioV42_etaS020->SetMarkerStyle(3);
    gRatioV42_etaS020->SetMarkerColor(kBlack);
    gRatioV42_etaS020->SetLineColor(kBlack);
    gRatioV42_etaS020->SetMarkerSize(1.5);
    gRatioV42_etaS020->SetFillStyle(3008);
    gRatioV42_etaS020->SetFillColor(kGreen-3);
    
    Double_t ratiov2_para1[] = {0.999, 1.010, 1.020, 1.019};
    Double_t ratiov2Err_para1[] = {0.015, 0.011, 0.007, 0.008};
    
    Double_t ratiov3_para1[] = {1.018, 1.034, 1.013, 1.042};
    Double_t ratiov3Err_para1[] = {0.015, 0.015, 0.011, 0.011};
    
    Double_t ratiov4_para1[] = {1.011, 1.057, 1.065, 1.092};
    Double_t ratiov4Err_para1[] = {0.015, 0.016, 0.011, 0.013};
    
    TGraphErrors *gRatioV22_para1 = new TGraphErrors(4,xHydro,ratiov2_para1,xErrHydro,ratiov2Err_para1);
    gRatioV22_para1->SetMarkerStyle(3);
    gRatioV22_para1->SetMarkerColor(kBlack);
    gRatioV22_para1->SetLineColor(0);
    gRatioV22_para1->SetMarkerSize(1.5);
    gRatioV22_para1->SetFillStyle(1001);
    gRatioV22_para1->SetFillColor(kRed-4);
    
    TGraphErrors *gRatioV32_para1 = new TGraphErrors(4,xHydro,ratiov3_para1,xErrHydro,ratiov3Err_para1);
    gRatioV32_para1->SetMarkerStyle(3);
    gRatioV32_para1->SetMarkerColor(kBlack);
    gRatioV32_para1->SetLineColor(kBlack);
    gRatioV32_para1->SetMarkerSize(1.5);
    gRatioV32_para1->SetFillStyle(1001);
    gRatioV32_para1->SetFillColor(kAzure+1);
    
    TGraphErrors *gRatioV42_para1 = new TGraphErrors(4,xHydro,ratiov4_para1,xErrHydro,ratiov4Err_para1);
    gRatioV42_para1->SetMarkerStyle(3);
    gRatioV42_para1->SetMarkerColor(kBlack);
    gRatioV42_para1->SetLineColor(kBlack);
    gRatioV42_para1->SetMarkerSize(1.5);
    gRatioV42_para1->SetFillStyle(1001);
    gRatioV42_para1->SetFillColor(kGreen-6);
    

    /*
    Double_t sysv22Gap10Run1_percent = 0.023;
    Double_t sysv24Run1_percent = 0.023;
    Double_t sysv32Gap10Run1_percent = 0.04;
    Double_t sysv42Gap10Run1_percent = 0.04;
    */

    Double_t v22Gap10Run1[50]={0.0272372, 0.0445713, 0.0640285, 0.083136, 0.0946093, 0.0994069, 0.0975893, 0.0900348, 0.0765231};
    Double_t v22Gap10Run1Err[50]={0.000123635, 0.000170732, 0.000146143, 0.000102203, 0.000155694, 0.000226115, 0.000388579, 0.00046476, 0.00202559};
    Double_t v22Gap10Run1Sys[50]={0.000680931, 0.00111428, 0.00160071, 0.00133018, 0.00151375, 0.00159051, 0.00156143, 0.00144056, 0.00122437};
    Double_t v22Gap10Run1CombErr[50]={0.000692064, 0.00112729, 0.00160737, 0.0013341, 0.00152173, 0.0016065, 0.00160905, 0.00151367, 0.00236688};
 
    Double_t v24Run1[50]={10.0143959, 0.0353041, 0.0567444, 0.074641, 0.0839708, 0.0851394, 0.0796972, 0.0689696, 0.0496503};
    Double_t v24Run1Err[50]={0.00127766, 9.79813e-05, 8.67241e-05, 0.000113496, 0.000139736, 0.000273059, 0.000533405, 0.00229005, 0.0122815};
    Double_t v24Run1Sys[50]={0.000331105, 0.000811995, 0.00130512, 0.00171674, 0.00193133, 0.00195821, 0.00183304, 0.0015863, 0.00114196};
    Double_t v24Run1CombErr[50]={0.00131987, 0.000817885, 0.001308, 0.00172049, 0.00193638, 0.00197715, 0.00190907, 0.0027858, 0.0123345};
    
    Double_t v32Gap10Run1[50]={0.0201462, 0.0232611, 0.0262213, 0.0288393, 0.030313, 0.0302248, 0.02806, 0.0233048, 0.0125242};
    Double_t v32Gap10Run1Err[50]={0.000191662, 0.000166292, 0.000202462, 0.000203758, 0.000255159, 0.000358623, 0.000641299, 0.00177173, 0.0107354};
    Double_t v32Gap10Run1Sys[50]={0.000805847, 0.000930443, 0.00104885, 0.00115357, 0.00121252, 0.00120899, 0.0011224, 0.000932192, 0.000500968};
    Double_t v32Gap10Run1CombErr[50]={0.000828326, 0.000945187, 0.00106821, 0.00117143, 0.00123908, 0.00126106, 0.00129269, 0.002002, 0.0107471};
    
    Double_t v42Gap10Run1[50]={0.00972377, 0.0113322, 0.0124995, 0.0138978, 0.0153994, 0.0151404, 0.0152379, 0.0107406, 0.00192339};
    Double_t v42Gap10Run1Err[50]={0.000289778, 0.000464912, 0.00026506, 0.00046642, 0.000566913, 0.0008171, 0.000962628, 0.00346496, 0.0636476};
    Double_t v42Gap10Run1Sys[50]={0.000388951, 0.000453288, 0.000499982, 0.000555911, 0.000615976, 0.000605617, 0.000609517, 0.000429625, 7.69358e-05};
    Double_t v42Gap10Run1CombErr[50]={0.00048503, 0.000649317, 0.000565896, 0.000725661, 0.000837148, 0.00101707, 0.00113937, 0.0034915, 0.0636477};
    
    TGraphErrors *graphv22Gap10Run1 = new TGraphErrors(9, xCross,v22Gap10Run1,exxx,v22Gap10Run1Err);
    graphv22Gap10Run1->SetMarkerStyle(33);
    graphv22Gap10Run1->SetMarkerColor(kMagenta+1);
    graphv22Gap10Run1->SetMarkerSize(1.5);
    
    TGraphErrors *graphv22Gap10Run1sys = new TGraphErrors(9, xCross,v22Gap10Run1,esysxxx,v22Gap10Run1CombErr);
    graphv22Gap10Run1sys->SetMarkerStyle(25);
    graphv22Gap10Run1sys->SetMarkerColor(kRed);
    graphv22Gap10Run1sys->SetMarkerSize(0.9);
    graphv22Gap10Run1sys->SetLineColor(kRed);
    graphv22Gap10Run1sys->SetFillColor(kRed-4);
    graphv22Gap10Run1sys->SetFillStyle(3001);
    ShiftAlongXaxis(graphv22Gap10Run1sys,-1.2);
    
    TGraphErrors *graphv24Run1sys = new TGraphErrors(8, xCross,v24Run1,esysxxx,v24Run1CombErr);
    graphv24Run1sys->SetMarkerStyle(28);
    graphv24Run1sys->SetMarkerColor(kGray+2);
    graphv24Run1sys->SetMarkerSize(1.2);
    graphv24Run1sys->SetLineColor(kGray+2);
    graphv24Run1sys->SetFillColor(kGray+2);
    graphv24Run1sys->SetFillStyle(3001);
    ShiftAlongXaxis(graphv24Run1sys,-2);
    
    TGraphErrors *graphv32Gap10Run1 = new TGraphErrors(6, xCross,v32Gap10Run1,exxx,v32Gap10Run1Err);
    graphv32Gap10Run1->SetMarkerStyle(24);
    graphv32Gap10Run1->SetMarkerColor(kBlue);
    graphv32Gap10Run1->SetMarkerSize(1.);
    
    TGraphErrors *graphv32Gap10Run1sys = new TGraphErrors(6, xCross,v32Gap10Run1,esysxxx,v32Gap10Run1CombErr);
    graphv32Gap10Run1sys->SetMarkerStyle(24);
    graphv32Gap10Run1sys->SetMarkerColor(kBlue);
    graphv32Gap10Run1sys->SetMarkerSize(0.9);
    graphv32Gap10Run1sys->SetLineColor(kBlue);
    graphv32Gap10Run1sys->SetFillColor(kBlue-8);
    graphv32Gap10Run1sys->SetFillStyle(3001);
    ShiftAlongXaxis(graphv32Gap10Run1sys,-1.2);
    
    TGraphErrors *graphv42Gap10Run1 = new TGraphErrors(6, xCross,v42Gap10Run1,exxx,v42Gap10Run1Err);
    graphv42Gap10Run1->SetMarkerStyle(27);
    graphv42Gap10Run1->SetMarkerColor(kGreen);
    graphv42Gap10Run1->SetMarkerSize(1.);
    
    TGraphErrors *graphv42Gap10Run1sys = new TGraphErrors(6, xCross,v42Gap10Run1,esysxxx,v42Gap10Run1CombErr);
    graphv42Gap10Run1sys->SetMarkerStyle(27);
    graphv42Gap10Run1sys->SetMarkerColor(kGreen+3);
    graphv42Gap10Run1sys->SetMarkerSize(1.5);
    graphv42Gap10Run1sys->SetLineColor(kGreen+3);
    graphv42Gap10Run1sys->SetFillColor(kGreen-8);
    graphv42Gap10Run1sys->SetFillStyle(3001);
    ShiftAlongXaxis(graphv42Gap10Run1sys,-1.2);
    
    
    
    // RUN2 data
    Double_t v22Gap10Run2[10]={0.0283859, 0.0456604, 0.0655068, 0.0870721, 0.099105, 0.104143, 0.10286, 0.0974591, 0.0888104};
    Double_t v22Gap10Run2Err[10]={0.000570139, 0.000643862, 0.000373923, 0.000444072, 0.000553776, 0.000732208, 0.00107352, 0.00186236, 0.00437842};
    Double_t v22Gap10Run2Sys[10]={0.000425788, 0.000684906, 0.000982603, 0.00130608, 0.00148658, 0.00156214, 0.00154291, 0.00146189, 0.00133216};
    Double_t v22Gap10Run2CombErr[10]={0.000711585, 0.000940029, 0.00105135, 0.00137951, 0.00158637, 0.00172523, 0.00187963, 0.00236759, 0.0045766};
     
    Double_t v32Gap10Run2[10]={0.0206723, 0.0231991, 0.0279915, 0.0309315, 0.0331377, 0.0323443, 0.0270494, 0.0276913};
    Double_t v32Gap10Run2Err[10]={0.000629474, 0.000797443, 0.000466369, 0.000619115, 0.000842254, 0.00136072, 0.00283039, 0.00546999};
    Double_t v32Gap10Run2Sys[10]={0.000950925, 0.00106716, 0.00128761, 0.00142285, 0.00152434, 0.00148784, 0.00124427, 0.0012738};
    Double_t v32Gap10Run2CombErr[10]={0.00114039, 0.0013322, 0.00136947, 0.00155171, 0.00174155, 0.00201624, 0.00309182, 0.00561635};
     
    Double_t v42Gap10Run2[10]={0.0114652, 0.0129068, 0.0138326, 0.0156838, 0.017018, 0.0162893, 0.0177113, 0.0198688, 0.00916289};
    Double_t v42Gap10Run2Err[10]={0.000949648, 0.00120259, 0.000788619, 0.00104461, 0.00147443, 0.00245469, 0.00406993, 0.0073516, 0.0383903};
    Double_t v42Gap10Run2Sys[10]={0.000550329, 0.000619527, 0.000663966, 0.000752823, 0.000816866, 0.000781885, 0.000850141, 0.000953704, 0.000439819};
    Double_t v42Gap10Run2CombErr[10]={0.00109759, 0.00135279, 0.00103091, 0.00128761, 0.00168559, 0.00257621, 0.00415777, 0.0074132, 0.0383928};
     
    Double_t v24Run2[10]={-10, 0.0366982, 0.0586746, 0.0779774, 0.0867288, 0.0890264, 0.0835125, 0.0806163, 0.0464366};
    Double_t v24Run2Err[10]={0.00616787, 0.00258022, 0.000331448, 0.00138072, 0.00257394, 0.00275397, 0.00207766, 0.00961209, 0.0835755};
    Double_t v24Run2Sys[10]={0.000303348, 0.000623869, 0.000997468, 0.00132562, 0.00147439, 0.00151345, 0.00141971, 0.00137048, 0.000789422};
    Double_t v24Run2CombErr[10]={0.00617532, 0.00265457, 0.00105109, 0.00191407, 0.00296631, 0.00314243, 0.0025164, 0.0097093, 0.0835792};
     
    Double_t v26Run2[10]={-1, 0.0361417, 0.0581481, 0.0774404, 0.086501, 0.0882722, 0.0827871, 0.0685936};
    Double_t v26Run2Err[10]={0.00827159, 0.00312685, 0.00031263, 0.00125004, 0.00244129, 0.00239164, 0.00455578, 0.0183562};
    Double_t v26Run2Sys[10]={0, 0.000614409, 0.000988518, 0.00131649, 0.00147052, 0.00150063, 0.00140738, 0.00116609};
    Double_t v26Run2CombErr[10]={0, 0.00318664, 0.00103678, 0.00181542, 0.00284997, 0.00282345, 0.00476821, 0.0183932};
    
    Double_t v28Run2[10]={-1, 0.0364682, 0.0580561, 0.0773885, 0.0866066, 0.0885098, 0.0826735, 0.0646657};
    Double_t v28Run2Err[10]={0, 0.00297485, 0.000345201, 0.00127119, 0.00241142, 0.00223079, 0.00610433, 0.039863};
    Double_t v28Run2Sys[10]={0, 0.00182341, 0.0029028, 0.00386943, 0.00433033, 0.00442549, 0.00413367, 0.00323328};
    Double_t v28Run2CombErr[10]={0, 0.0034892, 0.00292326, 0.00407288, 0.00495648, 0.00495595, 0.00737225, 0.0399939};
 
    TGraphErrors *graphv22Gap10 = new TGraphErrors(9,xCross,v22Gap10Run2,xCrossErr,v22Gap10Run2CombErr);
    graphv22Gap10->SetMarkerStyle(kFullSquare);
    graphv22Gap10->SetMarkerColor(kRed+1);
    graphv22Gap10->SetMarkerSize(1.);
    graphv22Gap10->SetLineColor(kRed+1);
    //ShiftAlongXaxis(graphv22Gap10,1.0);
    
    TGraphErrors *graphv32Gap10 = new TGraphErrors(8,xCross,v32Gap10Run2,xCrossErr,v32Gap10Run2CombErr);
    graphv32Gap10->SetMarkerStyle(20);
    graphv32Gap10->SetMarkerColor(kBlue+1);
    graphv32Gap10->SetMarkerSize(1.1);
    graphv32Gap10->SetMarkerColor(kBlue+1);
    //ShiftAlongXaxis(graphv32Gap10,-1.0);
    
    TGraphErrors *graphv42Gap10 = new TGraphErrors(8,xCross,v42Gap10Run2,xCrossErr,v42Gap10Run2CombErr);
    graphv42Gap10->SetMarkerStyle(33);
    graphv42Gap10->SetMarkerColor(kGreen+3);
    graphv42Gap10->SetMarkerSize(1.6);
    graphv42Gap10->SetLineColor(kGreen+3);
    //ShiftAlongXaxis(graphv42Gap10,1.5);
    
    TGraphErrors *graph_v24 = new TGraphErrors(8,xCross,v24Run2,xCrossErr,v24Run2CombErr);
    graph_v24->SetMarkerStyle(34);
    graph_v24->SetMarkerColor(kGray+2);
    graph_v24->SetLineColor(kGray+2);
    graph_v24->SetMarkerSize(1.2);
    ShiftAlongXaxis(graph_v24,-1);
    
    TGraphErrors *graph_v26 = new TGraphErrors(8,xCross,v26Run2,xCrossErr,v26Run2CombErr);
    graph_v26->SetMarkerStyle(33);
    graph_v26->SetMarkerColor(kOrange+1);
    graph_v26->SetMarkerSize(1.6);
    graph_v26->SetLineColor(kOrange+1);
    //ShiftAlongXaxis(graph_v26,-1);
    
    TGraphErrors *graph_v28 = new TGraphErrors(8,xCross,v28Run2,xCrossErr,v28Run2CombErr);
    graph_v28->SetMarkerStyle(3);
    graph_v28->SetMarkerColor(kBlack);
    graph_v28->SetLineColor(kBlack);
    graph_v28->SetMarkerSize(1.5);
    ShiftAlongXaxis(graph_v28,1);
    
    
    Double_t Rv22Gap10Run2[10]={1.02306, 1.01026, 1.01166, 1.03914, 1.04013, 1.04198, 1.05043, 1.07745, 1.15592};
    Double_t Rv22Gap10Run2Err[10]={0.0208325, 0.0145188, 0.00619087, 0.0055208, 0.00595475, 0.00753838, 0.0115733, 0.0222475, 0.05994};
    Double_t Rv22Gap10Run2Sys[10]={0.0102306, 0.0101026, 0.0101166, 0.0103914, 0.0104013, 0.0104198, 0.0105043, 0.0107745, 0.0115592};
    Double_t Rv22Gap10Run2CombErr[10]={0.023209, 0.0176878, 0.0118605, 0.0117669, 0.0119853, 0.0128608, 0.0156295, 0.0247193, 0.0610443};
    
    TGraphErrors *Ratiov22Gap10 = new TGraphErrors(9,xCross,Rv22Gap10Run2,xCrossErr,Rv22Gap10Run2Err);
    Ratiov22Gap10->SetMarkerStyle(kFullSquare);
    Ratiov22Gap10->SetMarkerColor(kRed+1);
    Ratiov22Gap10->SetMarkerSize(1.);
    Ratiov22Gap10->SetLineColor(kRed+1);
    
    TGraphErrors *Ratiov22Gap10sys = new TGraphErrors(9,xCross,Rv22Gap10Run2,xCrossSys,Rv22Gap10Run2Sys);
    Ratiov22Gap10sys->SetMarkerStyle(kFullSquare);
    Ratiov22Gap10sys->SetMarkerColor(kRed+1);
    Ratiov22Gap10sys->SetMarkerSize(1.);
    Ratiov22Gap10sys->SetLineColor(kRed+1);
    Ratiov22Gap10sys->SetFillColor(kRed-9);
    Ratiov22Gap10sys->SetFillStyle(3001);
    
    TGraphErrors *Ratiov22Gap10Comb = new TGraphErrors(9,xCross,Rv22Gap10Run2,xCrossErr,Rv22Gap10Run2CombErr);
    Ratiov22Gap10Comb->SetMarkerStyle(kFullSquare);
    Ratiov22Gap10Comb->SetMarkerColor(kRed+1);
    Ratiov22Gap10Comb->SetMarkerSize(1.);
    Ratiov22Gap10Comb->SetLineColor(kRed+1);
    Ratiov22Gap10Comb->SetFillColor(kRed-9);
    Ratiov22Gap10Comb->SetFillStyle(3001);
    
    Double_t Rv32Gap10Run2[10]={1.00616, 0.979296, 1.05403, 1.05889, 1.08273, 1.05802};
    Double_t Rv32Gap10Run2Err[10]={0.0322702, 0.0349908, 0.0186098, 0.0230543, 0.0293887, 0.049071};
    Double_t Rv32Gap10Run2Sys[10]={0.0191171, 0.0186066, 0.0200266, 0.0201189, 0.0205718, 0.0201024};
    Double_t Rv32Gap10Run2CombErr[10]={0.0375078, 0.0396303, 0.0273384, 0.0305986, 0.0358734, 0.0530289};
    
    
    TGraphErrors *Ratiov32Gap10 = new TGraphErrors(6,xCross,Rv32Gap10Run2,xCrossErr,Rv32Gap10Run2Err);
    Ratiov32Gap10->SetMarkerStyle(20);
    Ratiov32Gap10->SetMarkerColor(kBlue+1);
    Ratiov32Gap10->SetMarkerSize(1.);
    Ratiov32Gap10->SetLineColor(kBlue+1);
    
    TGraphErrors *Ratiov32Gap10sys = new TGraphErrors(6,xCross,Rv32Gap10Run2,xCrossSys,Rv32Gap10Run2Sys);
    Ratiov32Gap10sys->SetMarkerStyle(20);
    Ratiov32Gap10sys->SetMarkerColor(kBlue+1);
    Ratiov32Gap10sys->SetMarkerSize(1.);
    Ratiov32Gap10sys->SetLineColor(kBlue+1);
    Ratiov32Gap10sys->SetFillColor(kBlue-9);
    Ratiov32Gap10sys->SetFillStyle(3001);
    
    TGraphErrors *Ratiov32Gap10Comb = new TGraphErrors(6,xCross,Rv32Gap10Run2,xCrossErr,Rv32Gap10Run2CombErr);
    Ratiov32Gap10Comb->SetMarkerStyle(20);
    Ratiov32Gap10Comb->SetMarkerColor(kBlue+1);
    Ratiov32Gap10Comb->SetMarkerSize(1.);
    Ratiov32Gap10Comb->SetLineColor(kBlue+1);
    Ratiov32Gap10Comb->SetFillColor(kBlue-9);
    Ratiov32Gap10Comb->SetFillStyle(3001);
    
    
    Double_t Rv42Gap10Run2[10]={1.15567, 1.11174, 1.08522, 1.10595, 1.09148, 1.05477};
    Double_t Rv42Gap10Run2Err[10]={0.101823, 0.107807, 0.0634186, 0.0828282, 0.0960193, 0.170004};
    Double_t Rv42Gap10Run2Sys[10]={0.0196464, 0.0188996, 0.0184487, 0.0188011, 0.0185552, 0.0179311};
    Double_t Rv42Gap10Run2CombErr[10]={0.103701, 0.109451, 0.0660475, 0.0849352, 0.0977957, 0.170947};
    
    TGraphErrors *Ratiov42Gap10 = new TGraphErrors(6,xCross,Rv42Gap10Run2,xCrossErr,Rv42Gap10Run2Err);
    Ratiov42Gap10->SetMarkerStyle(34);
    Ratiov42Gap10->SetMarkerColor(kMagenta+1);
    Ratiov42Gap10->SetMarkerSize(1.5);
    Ratiov42Gap10->SetLineColor(kMagenta+1);
    ShiftAlongXaxis(Ratiov42Gap10,1.0);
    
    TGraphErrors *Ratiov42Gap10sys = new TGraphErrors(6,xCross,Rv42Gap10Run2,xCrossSys,Rv42Gap10Run2Sys);
    Ratiov42Gap10sys->SetMarkerStyle(34);
    Ratiov42Gap10sys->SetMarkerColor(kMagenta+1);
    Ratiov42Gap10sys->SetMarkerSize(1.5);
    Ratiov42Gap10sys->SetLineColor(kMagenta+1);
    Ratiov42Gap10sys->SetFillColor(kMagenta-9);
    Ratiov42Gap10sys->SetFillStyle(3001);
    ShiftAlongXaxis(Ratiov42Gap10sys,1.0);
    
    TGraphErrors *Ratiov42Gap10Comb = new TGraphErrors(6,xCross,Rv42Gap10Run2,xCrossErr,Rv42Gap10Run2CombErr);
    Ratiov42Gap10Comb->SetMarkerStyle(33);
    Ratiov42Gap10Comb->SetMarkerColor(kGreen+3);
    Ratiov42Gap10Comb->SetMarkerSize(1.5);
    Ratiov42Gap10Comb->SetLineColor(kGreen+3);
    Ratiov42Gap10Comb->SetFillColor(kGreen-9);
    Ratiov42Gap10Comb->SetFillStyle(3001);
    ShiftAlongXaxis(Ratiov42Gap10Comb,1.0);
    
    Double_t Rv24Run2[10]={-1, 1.02123, 1.02079, 1.03541, 1.0261, 1.04079, 1.04489, 1.16396};
    Double_t Rv24Run2Err[10]={0, 0.07194, 0.00637813, 0.0184073, 0.0304852, 0.0323918, 0.0280399, 0.140858};
    Double_t Rv24Run2Sys[10]={0, 0.0122548, 0.0122495, 0.0124249, 0.0123132, 0.0124894, 0.0125387, 0.0139676};
    Double_t Rv24Run2CombErr[10]={0, 0.0729764, 0.0138105, 0.0222083, 0.0328779, 0.0347162, 0.0307157, 0.141549};
    
    TGraphErrors *Ratiov24 = new TGraphErrors(9,xCross,Rv24Run2,xCrossErr,Rv24Run2Err);
    Ratiov24->SetMarkerStyle(34);
    Ratiov24->SetMarkerColor(kGray+2);
    Ratiov24->SetMarkerSize(1.5);
    Ratiov24->SetLineColor(kGray+2);
    ShiftAlongXaxis(Ratiov24,-1.0);
    
    TGraphErrors *Ratiov24sys = new TGraphErrors(9,xCross,Rv24Run2,xCrossSys,Rv24Run2Sys);
    Ratiov24sys->SetMarkerStyle(34);
    Ratiov24sys->SetMarkerColor(kGray+2);
    Ratiov24sys->SetMarkerSize(1.);
    Ratiov24sys->SetLineColor(kGray+2);
    Ratiov24sys->SetFillColor(kGray-8);
    Ratiov24sys->SetFillStyle(1001);
    ShiftAlongXaxis(Ratiov24sys,-1.0);
    
    TGraphErrors *Ratiov24Comb = new TGraphErrors(9,xCross,Rv24Run2,xCrossErr,Rv24Run2CombErr);
    Ratiov24Comb->SetMarkerStyle(34);
    Ratiov24Comb->SetMarkerColor(kGray+2);
    Ratiov24Comb->SetMarkerSize(1.2);
    Ratiov24Comb->SetLineColor(kGray+2);
    Ratiov24Comb->SetFillColor(kGray-8);
    Ratiov24Comb->SetFillStyle(1001);
    ShiftAlongXaxis(Ratiov24Comb,-1.0);
    
    
    // ***********************************************
    //          Durham data from here
    // ***********************************************
    
    //cout << "RUN1"<<endl;
    printf("centrality  v_{2}{2,|#Delta#eta|>1} (RUN1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<9;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v22Gap10Run1[p],v22Gap10Run1Err[p], v22Gap10Run1Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  v_{2}{4} (RUN1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<8;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v24Run1[p],v24Run1Err[p], v24Run1Sys[p]);
    }
    cout<<endl;
    
    
    printf("centrality  v_{3}{2,|#Delta#eta|>1} (RUN1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<6;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v32Gap10Run1[p],v32Gap10Run1Err[p], v32Gap10Run1Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  v_{4}{2,|#Delta#eta|>1} (RUN1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<6;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v42Gap10Run1[p],v42Gap10Run1Err[p], v42Gap10Run1Sys[p]);
    }
    cout<<endl;
    
    
    //cout << "Run2"<<endl;
    printf("centrality  v_{2}{2,|#Delta#eta|>1} (Run2)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<9;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v22Gap10Run2[p],v22Gap10Run2Err[p], v22Gap10Run2Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  v_{2}{4} (Run2)    +-stat.err  +-syst.err \n");
    for(Int_t p=1;p<8;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v24Run2[p],v24Run2Err[p], v24Run2Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  v_{2}{6} (Run2)    +-stat.err  +-syst.err \n");
    for(Int_t p=1;p<8;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v26Run2[p],v26Run2Err[p], v26Run2Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  v_{2}{8} (Run2)    +-stat.err  +-syst.err \n");
    for(Int_t p=1;p<8;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v28Run2[p],v28Run2Err[p], v28Run2Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  v_{3}{2,|#Delta#eta|>1} (Run2)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<9;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v32Gap10Run2[p],v32Gap10Run2Err[p], v32Gap10Run2Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  v_{4}{2,|#Delta#eta|>1} (Run2)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<9;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],v42Gap10Run2[p],v42Gap10Run2Err[p], v42Gap10Run2Sys[p]);
    }
    cout<<endl;
    
    
    printf("centrality  ratio:v_{2}{2,|#Delta#eta|>1} (Run2/Run1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<9;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],Rv22Gap10Run2[p],Rv22Gap10Run2Err[p], Rv22Gap10Run2Sys[p]);
    }
    cout<<endl;
    

    printf("centrality  ratio:v_{2}{4} (Run2/Run1)    +-stat.err  +-syst.err \n");
    for(Int_t p=1;p<8;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],Rv24Run2[p],Rv24Run2Err[p], Rv24Run2Sys[p]);
    }
    cout<<endl;
    
    
    printf("centrality  ratio:v_{3}{2,|#Delta#eta|>1} (Run2/Run1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<6;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],Rv32Gap10Run2[p],Rv32Gap10Run2Err[p], Rv32Gap10Run2Sys[p]);
    }
    cout<<endl;
    
    printf("centrality  ratio:v_{4}{2,|#Delta#eta|>1} (Run2/Run1)    +-stat.err  +-syst.err \n");
    for(Int_t p=0;p<6;p++)
    {
        printf("%f  %f  +-%f  +-%f\n",xCross[p],Rv42Gap10Run2[p],Rv42Gap10Run2Err[p], Rv42Gap10Run2Sys[p]);
    }
    cout<<endl;
    
    TLine *t2 = new TLine(0, 1, 80, 1);
    t2->SetLineWidth(2);
    t2->SetLineColor(kGray+2);
    
    TLine *t1 = new TLine(0, 0, 80, 0);
    t1->SetLineWidth(2);
    t1->SetLineColor(kGray+2);
    
    TLegend *l0 = new TLegend(0.6,0.3,0.85,0.5,"This Analysis");
    l0->SetFillColor(0);
    l0->SetBorderSize(1);
    l0->SetLineColor(0);
    //l0->SetHeader("This Analysis");
    l0->SetTextSize(0.035);
    l0->SetTextFont(42);
    l0->AddEntry(graph_v24,"v_{2}{4}","p");
    l0->AddEntry(graph_v26,"v_{2}{6}","p");
    l0->AddEntry(graph_v28,"v_{2}{8}","p");
    //l0->AddEntry(graph62211,"#LT v_{6} v_{4}^{2} v_{2}^{2} cos (6#Psi_{6} - 2#Psi_{4} - 2#Psi_{2}) #GT","p");
    
    TLatex *a = new TLatex(0.18,0.92,"ALICE Pb-Pb");
    a->SetNDC();
    a->SetTextSize(0.05);
    
    TLatex *b = new TLatex(0.7,0.92,"Hydrodynamics");
    b->SetNDC();
    b->SetTextSize(0.05);
    
    TLatex *a1 = new TLatex(0.9,0.1,"(a)");
    a1->SetNDC();
    a1->SetTextSize(0.07);
    
    TLatex *a2 = new TLatex(0.9,0.1,"(b)");
    a2->SetNDC();
    a2->SetTextSize(0.2);
    
    TLatex *a3 = new TLatex(0.9,0.5,"(c)");
    a3->SetNDC();
    a3->SetTextSize(0.135);
    
    TLegend *l1 = new TLegend(0.18,0.53,0.3,0.90);
    l1->SetFillColor(0);
    l1->SetBorderSize(1);
    l1->SetLineColor(0);
    //l1->SetHeader("ALICE Pb-Pb #sqrt{s_{NN}} = 5 TeV");
    l1->SetHeader("5.02 TeV");
    //l1->SetHeader("This Analysis");
    l1->SetTextSize(0.05);
    l1->SetTextFont(42);
    l1->AddEntry(graphv22Gap10,"#it{v}_{2 }{2, |#Delta#eta|>1}","p");
    //l1->AddEntry(GrV2ZDC,"v_{2}{#Psi_{RP}}","p");
    l1->AddEntry(graphv32Gap10,"#it{v}_{3 }{2, |#Delta#eta|>1}","p");
    l1->AddEntry(graphv42Gap10,"#it{v}_{4 }{2, |#Delta#eta|>1}","p");
    l1->AddEntry(graph_v24,"#it{v}_{2 }{4}","p");
    l1->AddEntry(graph_v26,"#it{v}_{2 }{6}","p");
    l1->AddEntry(graph_v28,"#it{v}_{2 }{8}","p");
    //l1->AddEntry(graphv32Gap10,"v_{3 }{2, |#Delta#eta|>1}","p");
    //l1->AddEntry(graphv42Gap10,"v_{4 }{2, |#Delta#eta|>1}","p");
    //l1->AddEntry(graphv5SP,"v_{5 }{2, |#Delta#eta|>1}","p");
    //l1->AddEntry(graphv6SP,"v_{6 }{2, |#Delta#eta|>1}","p");
    
    
    TLegend *l2 = new TLegend(0.42,0.65,0.65,0.90);
    l2->SetFillColor(0);
    l2->SetBorderSize(1);
    l2->SetLineColor(0);
    l2->SetHeader("2.76 TeV");
    l2->SetTextSize(0.05);
    l2->SetTextFont(42);
    l2->AddEntry(graphv22Gap10Run1sys,"#it{v}_{2 }{2, |#Delta#eta|>1}","p");
    //l2->AddEntry(SP_TPCrefMultTPConlyAll2,"#it{v}_{2 }{2, |#Delta#eta|>1}","f");
    //l2->AddEntry(SP_v3_corr2,"#it{v}_{3 }{2, |#Delta#eta|>1}","f");
    //l2->AddEntry(SP_v4_corr2,"#it{v}_{4 }{2, |#Delta#eta|>1}","f");
    l2->AddEntry(graphv32Gap10Run1sys,"#it{v}_{3 }{2, |#Delta#eta|>1}","p");
    l2->AddEntry(graphv42Gap10Run1sys,"#it{v}_{4 }{2, |#Delta#eta|>1}","p");
    l2->AddEntry(graphv24Run1sys,"#it{v}_{2 }{4}","p");
    
    TLegend *lJYO = new TLegend(0.7,0.76,0.88,0.90);
    lJYO->SetFillColor(0);
    lJYO->SetBorderSize(1);
    lJYO->SetLineColor(0);
    lJYO->SetHeader("5.02 TeV, Ref.[27]");
    lJYO->SetTextSize(0.05);
    lJYO->SetTextFont(42);
    //lJYO->AddEntry(Graphv22Gap10_hydro,"#it{v}_{2 }{2, |#Delta#eta|>1}","f");
    //lJYO->AddEntry(Graphv32Gap10_hydro,"#it{v}_{3 }{2, |#Delta#eta|>1}","f");
    lJYO->AddEntry(Graphv22Gap10_hydro,"#it{v}_{2 }{2, |#Delta#eta|>1}","");
    lJYO->AddEntry(Graphv32Gap10_hydro,"#it{v}_{3 }{2, |#Delta#eta|>1}","");
    
    
    TLegend *myLegend = new TLegend(0.3,0.68,0.55,0.88,"");
    //myLegendSetUp(myLegend,0.04);
    myLegend->SetFillColor(0);
    myLegend->SetBorderSize(1);
    myLegend->SetLineColor(0);
    myLegend->SetTextSize(0.035);
    myLegend->SetTextFont(42);
    //myLegend->AddEntry(GrV22,"v_{2}{2}","p");
    //myLegend->AddEntry(GrQaQb_V0_EtaGap00,"v_{2}{2} (#left|#Delta#eta#right| > 0)","p");
    //myLegend->AddEntry(GrSP,"v_{2}{2} (#left|#Delta#eta#right| > 1)","pf");
    //myLegend->AddEntry(GrV24,"v_{2}{4}","fp");
    //myLegend->AddEntry(GrV26,"v_{2}{6}","fp");
    //myLegend->AddEntry(GrV28,"v_{2}{8}","fp");
    //myLegend->AddEntry(graphv32Gap10,"v_{3 }{2, |#Delta#eta|>1}","p");
    //myLegend->AddEntry(graphv42Gap10,"v_{4 }{2, |#Delta#eta|>1}","p");
    if (ifV5) myLegend->AddEntry(graphv5SP,"#it{v}_{5 }{2, |#Delta#eta|>1}","p");
    if (ifV6) myLegend->AddEntry(graphv6SP,"#it{v}_{6 }{2, |#Delta#eta|>1}","p");

   
    
    TF1* fitv22Gap10 = new TF1("fitv22Gap10", "[0]", 0, 50);
    fitv22Gap10->SetLineStyle(7);
    fitv22Gap10->SetLineColor(kRed+1);
    
    TF1* fitv24 = new TF1("fitv24", "[0]", 0, 50);
    fitv24->SetLineStyle(2);
    fitv24->SetLineColor(kGreen+3);
    
    
    TF1* fitv32Gap10 = new TF1("fitv32Gap10", "[0]", 0, 50);
    fitv32Gap10->SetLineStyle(4);
    fitv32Gap10->SetLineColor(kBlue+1);
    
    TF1* fitv42Gap10 = new TF1("fitv42Gap10", "[0]", 0, 50);
    fitv42Gap10->SetLineStyle(3);
    fitv42Gap10->SetLineColor(kMagenta+1);
    
    TCanvas *cVn = new TCanvas("cVn", "cVn", 500, 600);
    gStyle->SetOptStat(0);
    cVn->Range(0,0,1,1);
    cVn->SetFillColor(0);
    cVn->SetBorderMode(0);
    cVn->SetBorderSize(2);
    cVn->SetLeftMargin(0.03333334);
    cVn->SetRightMargin(0.03125);
    cVn->SetTopMargin(0.04761905);
    cVn->SetBottomMargin(0.05);
    cVn->SetFrameBorderMode(0);
    
    TPad *PadVn_0 = new TPad("PadVn_0", "PadVn_0",0.,0.44,0.97,0.97);
    PadVn_0->Draw();
    PadVn_0->cd();
    PadVn_0->Range(0.,0.06, 6.,0.4);
    PadVn_0->SetFillColor(0);
    PadVn_0->SetBorderMode(0);
    PadVn_0->SetBorderSize(2);
    PadVn_0->SetLeftMargin(0.15);
    PadVn_0->SetRightMargin(0.02);
    PadVn_0->SetTopMargin(0.01);
    PadVn_0->SetBottomMargin(0);
    PadVn_0->SetFrameBorderMode(0);
    PadVn_0->SetFrameBorderMode(0);
    
    TH1D *framevn_0 = new TH1D("framevn_0","",80,0,80);
    framevn_0->SetMinimum(0.001);
    framevn_0->SetMaximum(0.17);
    framevn_0->SetStats(0);
    framevn_0->SetLineStyle(2);
    framevn_0->GetXaxis()->SetTitle("Centrality percentile");
    framevn_0->GetYaxis()->SetTitle(" #it{v}_{n}  ");
    framevn_0->GetXaxis()->SetNdivisions(508);
    framevn_0->GetYaxis()->SetNdivisions(505);
    framevn_0->GetYaxis()->SetLabelSize(0.07);
    framevn_0->GetYaxis()->SetTitleSize(0.07);
    framevn_0->GetYaxis()->SetTitleOffset(1.05);
    framevn_0->Draw("");
    
   
    // hydro
    /*
    Graphv22Gap10Min_hydro->Draw("lsame");
    Graphv22Gap10Max_hydro->Draw("lsame");
    Graphv32Gap10Min_hydro->Draw("lsame");
    Graphv32Gap10Max_hydro->Draw("lsame");
     */
    Graphv22Gap10_hydro->Draw("same3");
    Graphv32Gap10_hydro->Draw("same3");
    // RUN1
    graphv22Gap10Run1sys->Draw("pzsame");
    graphv24Run1sys->Draw("pzsame");
    graphv32Gap10Run1sys->Draw("pzsame");
    graphv42Gap10Run1sys->Draw("pzsame");
    // RUN2
    graph_v24->Draw("pzsame");
    graph_v26->Draw("pzsame");
    graph_v28->Draw("pzsame");
    graphv22Gap10->Draw("pzsame");
    graphv32Gap10->Draw("pzsame");
    graphv42Gap10->Draw("pzsame");
    
   
    
    
    t1->Draw("same");
    
    
    myLegend->Draw("same");
    
    
    l1->Draw("same");
    l2->Draw("same");
    lJYO->Draw("same");
    a->Draw("same");
    b->Draw("same");
    a1->Draw("same");
    
    
    TLine *thyv2max = new TLine(53, 0.145, 57, 0.145);
    thyv2max->SetLineColor(kRed+2);
    thyv2max->SetLineWidth(1);
    thyv2max->Draw("same");
    
    TLine *thyv2min = new TLine(53, 0.142, 57, 0.142);
    thyv2min->SetLineColor(kRed+2);
    thyv2min->SetLineWidth(1);
    thyv2min->Draw("same");
    
    TLine *tthyv2max = new TLine(53, 0.145, 53, 0.142);
    tthyv2max->SetLineColor(kRed+2);
    tthyv2max->SetLineWidth(1);
    tthyv2max->Draw("same");
    
    TLine *tthyv2min = new TLine(57, 0.145, 57, 0.142);
    tthyv2min->SetLineColor(kRed+2);
    tthyv2min->SetLineWidth(1);
    tthyv2min->Draw("same");
    
    TLine *thyv2min = new TLine(53, 0.142, 57, 0.142);
    thyv2min->SetLineColor(kRed+2);
    thyv2min->SetLineWidth(1);
    thyv2min->Draw("same");
    
    TLine *thyv3max = new TLine(53, 0.136, 57, 0.136);
    thyv3max->SetLineColor(kBlue+2);
    thyv3max->SetLineWidth(1);
    thyv3max->Draw("same");
    
    TLine *thyv3min = new TLine(53, 0.133, 57, 0.133);
    thyv3min->SetLineColor(kBlue+2);
    thyv3min->SetLineWidth(1);
    thyv3min->Draw("same");
    
    TLine *tthyv3max = new TLine(57, 0.133, 57, 0.136);
    tthyv3max->SetLineColor(kBlue+2);
    tthyv3max->SetLineWidth(1);
    tthyv3max->Draw("same");
    
    TLine *tthyv3min = new TLine(53, 0.133, 53, 0.136);
    tthyv3min->SetLineColor(kBlue+2);
    tthyv3min->SetLineWidth(1);
    tthyv3min->Draw("same");

 
   
    PadVn_0->Modified();
    cVn->cd();
    
    
    TLegend *myLegend22 = new TLegend(0.3,0.48,0.55,0.88,"Hydrodynamics, Ref.[25]");
    //myLegendSetUp(myLegend,0.04);
    myLegend22->SetFillColor(0);
    myLegend22->SetBorderSize(1);
    myLegend22->SetLineColor(0);
    myLegend22->SetTextSize(0.16);
    myLegend22->SetTextFont(42);
    myLegend22->AddEntry(gRatioV22_para1,"#eta/s(T), param1","f");
    myLegend22->AddEntry(gRatioV22_etaS020,"#eta/s = 0.20","f");
    
    TPad *PadVn_1 = new TPad("PadVn_1", "PadVn_1",0.,0.27,0.97,0.44);
    PadVn_1->Draw();
    PadVn_1->cd();
    PadVn_1->Range(0.,0.06, 6.,0.4);
    PadVn_1->SetFillColor(0);
    PadVn_1->SetBorderMode(0);
    PadVn_1->SetBorderSize(2);
    PadVn_1->SetLeftMargin(0.15);
    PadVn_1->SetRightMargin(0.02);
    PadVn_1->SetTopMargin(0.0);
    PadVn_1->SetBottomMargin(0);
    PadVn_1->SetFrameBorderMode(0);
    PadVn_1->SetFrameBorderMode(0);
    
    TH1D *framevn_1 = new TH1D("framevn_1","",80,0,80);
    framevn_1->SetMinimum(0.94);
    framevn_1->SetMaximum(1.24);
    framevn_1->SetStats(0);
    framevn_1->SetLineStyle(2);
    framevn_1->GetXaxis()->SetTitle("Centrality percentile");
    framevn_1->GetYaxis()->SetTitle("Ratio  ");
    framevn_1->GetXaxis()->SetNdivisions(508);
    framevn_1->GetYaxis()->SetNdivisions(503);
    framevn_1->GetYaxis()->SetLabelSize(0.22);
    
    framevn_1->GetXaxis()->SetTickLength(0.1);
    framevn_1->GetYaxis()->SetTickLength(0.03);
    framevn_1->GetYaxis()->SetTitleOffset(0.25);
    framevn_1->GetYaxis()->SetTitleSize(0.23);
    framevn_1->Draw("");
    //t2->Draw("");
    
    myLegend22->Draw("");
    
    gRatioV22_para1->Draw("3same");
    gRatioV22_etaS020->Draw("3same");
    
    Ratiov22Gap10Comb->Draw("pzsame");
    Ratiov24Comb->Draw("pzsame");
    
    //Ratiov22Gap10Comb->Fit("fitv22Gap10","R");
    //Ratiov24Comb->Fit("fitv24","R");
    
    //Ratiov22Gap10->Draw("pzsame");
    //Ratiov24->Draw("pzsame");
   
    a2->Draw("same");
    
    PadVn_1->Modified();
    cVn->cd();
    
    
    TPad *PadVn_2 = new TPad("PadVn_2", "PadVn_2",0.0,0.00,0.97,0.27);
    PadVn_2->Draw();
    PadVn_2->cd();
    PadVn_2->Range(0.,0.06, 6.,0.4);
    PadVn_2->SetFillColor(0);
    PadVn_2->SetBorderMode(0);
    PadVn_2->SetBorderSize(2);
    PadVn_2->SetLeftMargin(0.15);
    PadVn_2->SetRightMargin(0.02);
    PadVn_2->SetTopMargin(0.0);
    PadVn_2->SetBottomMargin(0.4);
    PadVn_2->SetFrameBorderMode(0);
    PadVn_2->SetFrameBorderMode(0);
    
    TH1D *framevn_2 = new TH1D("framevn_2","",80,0,80);
    framevn_2->SetMinimum(0.94);
    framevn_2->SetMaximum(1.24);
    framevn_2->SetStats(0);
    framevn_2->SetLineStyle(2);
    framevn_2->GetXaxis()->SetTitle("Centrality percentile");
    framevn_2->GetYaxis()->SetTitle("Ratio  ");
    framevn_2->GetXaxis()->SetNdivisions(508);
    framevn_2->GetYaxis()->SetNdivisions(503);
    framevn_2->GetYaxis()->SetLabelSize(0.13);
    framevn_2->GetXaxis()->SetLabelSize(0.13);
    framevn_2->GetXaxis()->SetTitleSize(0.14);
    framevn_2->GetYaxis()->SetTitleSize(0.14);
    framevn_2->GetXaxis()->SetTickLength(0.07);
    framevn_2->GetYaxis()->SetTickLength(0.03);
    framevn_2->GetYaxis()->SetTitleOffset(0.43);
    framevn_2->GetXaxis()->SetTitleOffset(0.95);
    framevn_2->Draw("");
    //t2->Draw("");
    
    gRatioV32_para1->Draw("3same");
    gRatioV42_para1->Draw("3same");
    
    gRatioV32_etaS020->Draw("3same");
    gRatioV42_etaS020->Draw("3same");

    Ratiov32Gap10Comb->Draw("pzsame");
    Ratiov42Gap10Comb->Draw("pzsame");
    
    //Ratiov32Gap10Comb->Fit("fitv32Gap10","R");
    //Ratiov42Gap10Comb->Fit("fitv42Gap10","R");
    //Ratiov32Gap10->Draw("pzsame");
    
    
    //Ratiov42Gap10->Fit("fitv42Gap10","R");
    //Ratiov42Gap10->Draw("pzsame");
    
    a3->Draw("same");
}

//________________________________
void LoadLibs() {
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
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



void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
                  int currentMarkerStyle=20, int currentMarkerColor=0,
                  int currentLineStyle=1, int currentLineColor=0, int currentFillColor=0, int currentFillStyle=0)
{
    currentGraph->SetMarkerSize(currentMarkerSize);
    currentGraph->SetMarkerStyle(currentMarkerStyle);
    currentGraph->SetMarkerColor(currentMarkerColor);
    currentGraph->SetLineStyle(currentLineStyle);
    currentGraph->SetLineColor(currentLineColor);
    currentGraph->SetFillColor(currentFillColor);
    currentGraph->SetFillStyle(currentFillStyle);
    return;
}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}


void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}

void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
			  hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  

}

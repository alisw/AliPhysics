#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TMath.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#endif

void    GetRatE(double x,double xe, double y,double ye, double &rat, double &rate);
TGraphErrors* CreateGraph(const char* nm, int np, const double *x,const double *xe, const double* y,const double *ye,Bool_t rat=kFALSE,int col=kBlack,int mrk=20,double sz=1.2,int lst=1);

double Npart_Geom[11] =  {383.851,331.921,262.756,188.639,130.721 ,86.5384,54.0822,31.1483,16.0198,7.75875,3.8324};
double Npart_rms_Geom[11] = {16.2332,17.7241,26.9529,22.2198,19.1009,15.9985,13.4516,10.7482,7.93099,4.90297,2.52505};
//
double Npart_Geom_4b[4]     = {396.624,371.028,344.477,319.018};
double Npart_rms_Geom_4b[4] = {8.84252,11.016,12.1732,12.4974};
//

//double v0SelNov21[11] = 0,79,247,577,1185,2155,3565,5527,8203,12167,15073,19889
double Npart_V0CSelNov21[11] ={380.692,328.323,259.018,185.05,128.121,84.9204,52.6411,29.79,15.1951,6.73934,2.95315};
double Npart2_V0CSelNov21[11]={379.112,323.171,253.417,183.054,129.116,86.9088,54.2474,30.2884,15.1553,6.7094,2.95315};
double Npart_rms_V0CSelNov21[11] = { 18.4857,19.5242,27.0422,21.1484,16.5675,12.8486,9.64889,6.81264,4.5665,2.72584,1.2575};

//double spdSelNov21[11] = 0,21,69,165,345,637,1069,1683,2511,3721,4585,6213}
double Npart_spdSelNov21[11] ={378.888,324.787,254.892,183.689,128.448,85.0047,52.3777,28.8656,14.3408,6.70305,3.16337};
double Npart2_spdSelNov21[11]={379.973,328.414,259.244,185.094,128.309,84.9333,52.7138,29.9245,15.3195,6.98915,3.16337};
double Npart_rms_spdSelNov21[11] = { 21.8302, 25.7908,30.7147, 24.5131,19.7898,15.9204,12.1924,8.81696,6.18204,3.78075,1.76665 };

//double spdSel[11] = {0,29,85,191,385,687,1133,1743,2567,3765,4611,6213}; 
double Npart_spdSel[11] ={379.393,328.651,259.304,185.004,128.265,84.8281,52.7315,29.0842,15.0519,6.876,3.14467};
double Npart2_spdSel[11]={378.367,322.458,248.503,176.048,121.769,80.3259,49.7426,28.1977,14.5535,6.88646,3.19064};
double Npart_rms_spdSel[11] = { 21.8302, 25.7908,30.7147, 24.5131,19.7898,15.9204,12.1924,8.81696,6.18204,3.78075,1.76665 };
//
//double v0Sel[11]  = {0,107,297,659,1305,2301,3747,5715,8361,12307,15153,19889 };
double Npart_V0CSel[11]  ={381.040,328.781,259.314,185.197,128.160,84.8531,52.6495,29.7634,15.0721,6.67109,2.93958};
double Npart2_V0CSel[11] ={375.596,313.697,240.935,170.22,118.96,79.1409,49.3726,28.2845,14.6765,6.62758,2.91953};
double Npart_rms_V0CSel[11] = {18.4857,19.5242,27.0422,21.1484,16.5675,12.8486,9.64889,6.81264,4.5665,2.72584,1.2575};
//
//double v0rSel[11]  = {0,32.5,89.5,198.5,396.5,700.5,1140.5,1744.5,2562.5,3767.5,4647.5,6367.5};
double Npart_V0CRSel[11]  ={379.114,329.655,260.318,185.996,130.454,87.0497,54.7531,31.4739,15.9883,7.28226,3.23659};
double Npart2_V0CRSel[11] ={376.795,319.268,244.685,172.492,120.471,80.0706,49.9794,28.5653,14.7274,7.13181,3.27423};
double Npart_rms_V0CRSel[11] = {18.4857,19.5242,27.0422,21.1484,16.5675,12.8486,9.64889,6.81264,4.5665,2.72584,1.2575};
//
//double tpcSel[11]  = {0,13,37,85,171,307,507,783,1157,1685,2055,2787}; // these are  tracks 
double Npart_TPCSel[11]  ={378.667,329.234,259.17,185.129,128.221,84.7075,52.6338,29.7764,15.1032,6.98597,3.22604}; // these are  tracks 
double Npart2_TPCSel[11] ={378.513,325.472,251.739,178.55,123.373,81.187,50.2066,28.4005,14.466,6.92418,3.30427};
double Npart_rms_TPCSel[11] = {18.4857,19.5242,27.0422,21.1484,16.5675,12.8486,9.64889,6.81264,4.5665,2.72584,1.2575}; // wrong
//

//double spdSelNov21[11] = 0,21,69,165,345,637,1069,1683,2511,3721,4585,6213}
double Npart_spdSelNov21_4b[4] ={396.624,371.02,344.477,319.018};
double Npart_rms_spdSelNov21_4b[11] = {8.84252,11.016,12.1732,12.4974};
//

//----------------------------------------------------------------------------------
double mltV0Ch11Nov27[] = {1590.09,1272.03,951.40,640.01,419.10,256.36,146.50,75.01,34.15,13.10,3.37};
double mltV0Ch11eNov27[] = {2.29,1.98,1.17,0.94,0.75,0.56,0.41,0.28,0.19,0.11,0.05};

// new baseline with ZDC TDC cleanup
double mltV0Ch11Nov29[]  = {1600.98,1293.63,966.21,649.30,426.09,260.65,149.44,76.23,34.75,13.28,3.38};
double mltV0Ch11eNov29[] = {2.34,1.96,1.17,0.93,0.74,0.55,0.41,0.28,0.19,0.11,0.05};

// new baseline with ZDC TDC cleanup
double mltV0Ch11Nov29_4b[4]  = {1693.65,1512.59,1365.12,1222.52};
double mltV0Ch11Nov29_4be[4] = {3.69,3.04,2.91,2.66};


//ZDC vs V0C, LHC10h_000137161_ZDCV0_24Nov10RS
// result from setting ZDC&V0 to 80% in the same V0 as V01D
////double mltZDCV0Ch8Nov24RS[11]={1599.53,1278.17,941.59,625.66,402.07,240.92,130.04,56.54,29.69,15.04,8.12};
////double mltZDCV0Ch8eNov24RS[11]={2.31,1.89,1.13,0.89,0.69,0.51,0.35,0.22,0.16,0.12,0.11};
// result from setting ZDC&V0 to have the same excess at 0 as V01D
double mltZDCV0Ch8Nov24RS[11]={1606.24,1294.40,964.34,651.15,427.81,264.23,149.39,72.52,34.75,18.79,8.92};
double mltZDCV0Ch8eNov24RS[11]={2.34,1.93,1.16,0.92,0.74,0.55,0.40,0.27,0.18,0.14,0.12};

//----------------------------------------------------------------------------------
//V0C, LHC10h8_000137161_V0C
double mltV0Ch8Nov21[11]={1606.14,1293.78,963.88,650.18,427.15,263.20,151.27,77.95,36.07,13.72,4.19};
double mltV0Ch8eNov21[11]={2.34,1.93,1.16,0.92,0.74,0.55,0.40,0.28,0.19,0.11,0.06};


//-----------------------------------
//SPD, LHC10h8_000137161_SPD
double mltSPDh8Nov21[11]={1608.81,1293.03,962.20,647.64,424.17,259.92,149.15,75.83,34.56,12.68,3.40};
double mltSPDh8eNov21[11]={2.31,1.96,1.15,0.91,0.72,0.54,0.40,0.28,0.19,0.11,0.06};


//SPD, LHC10h6_000137161_SPD (off-vertex primaries removed)
double mltSPDh6[11]={1603.99,1298.08,972.26,661.06,441.62,276.34,163.10,86.02,40.75,16.28,4.13};
double mltSPDh6e[11]={2.54,2.11,1.26,1.01,0.80,0.61,0.46,0.33,0.23,0.14,0.06};

//-----------------------------------
//SPD, LHC10h8_000137161_SPD
double mltSPDh8[11]={1613.30,1304.45,977.46,665.32,442.60,277.46,163.26,86.15,40.97,16.36,4.17};
double mltSPDh8e[11]={2.35,1.99,1.18,0.94,0.76,0.57,0.43,0.31,0.21,0.13,0.06};

//-----------------------
//V0C, LHC10h8_000137161_V0C
double mltV0Ch8[11]={1610.21,1304.62,977.99,667.33,444.72,279.27,164.76,87.79,42.34,17.53,4.89};
double mltV0Ch8e[11] = {2.38,1.97,1.18,0.95,0.76,0.57,0.43,0.31,0.21,0.13,0.06};

//-----------------------
//V0CR, LHC10h8_000137161_V0CR
double mltV0CRh8[11]={1611.13,1304.76,994.16,678.31,449.72,282.53,166.06,87.97,42.40,17.49,4.99};
double mltV0CRh8e[11]={2.38,1.98,1.19,0.97,0.77,0.58,0.43,0.31,0.21,0.13,0.06};

//-----------------------
//TPC, LHC10h8_000137161_TPC
double mltTPCh8[11]={1616.65,1304.32,957.72,656.67,435.59,273.39,160.58,86.09,40.96,16.58,4.16};
double mltTPCh8e[11]={2.44,2.01,1.15,0.92,0.74,0.56,0.42,0.31,0.21,0.13,0.06};

//--------------------------
//137366, SPD, LHC10h8_000137366_SPD
double mlt366SPDh8[11]={1574.05,1278.00,955.98,651.66,432.23,272.50,158.64,84.24,39.71,16.16,4.01};
double mlt366SPDh8e[11]={1.96,1.69,0.98,0.81,0.62,0.48,0.36,0.25,0.17,0.10,0.05};

//--------------------------
//SPD, LHC10h2_000137161_SPD  DPMJET
double  mltSPDdpmjh2[11]={1624.52,1304.84,976.87,664.33,443.03,277.41,163.68,85.95,40.76,16.58,4.25};
double  mltSPDdpmjh2e[11]={2.12,2.46,1.53,1.29,1.06,0.81,0.62,0.43,0.29,0.19,0.09};

//-------------------------------------------------------------------------
// top 4x2.5% bins in SPD 
double mltSPDh8_4b[4] = {1702.89,1516.20,1362.55,1221.55}; // w.means: 1587.4, 1259.0
double mltSPDh8_4be[4] = {3.55,3.06,2.91,2.65};
//double mltSPDh8_4b[4] = {1681.25,1480.93,1322.75,1188.11}; // w.means: 1587.4, 1259.0
//double mltSPDh8_4be[4] = {3.51,3.01,2.83,2.59};

//-------------------------------------------------------------------------
// top 4x2.5% bins in V0
double mltV0Ch8_4b[4]  = {1692.60,1514.45,1363.53,1221.21};
double mltV0Ch8_4be[4] = {3.69,3.02,2.85,2.61};
//double mltSPDh8_4b[4] = {1681.25,1480.93,1322.75,1188.11}; // w.means: 1587.4, 1259.0
//double mltSPDh8_4be[4] = {3.51,3.01,2.83,2.59};



TGraphErrors *MV0Ch11Nov27A=0,*RV0Ch11Nov27A=0;

TGraphErrors *MV0Ch11Nov29=0,*RV0Ch11Nov29=0;

TGraphErrors *MzdcV0CNov24=0,*RzdcV0CNov24=0;
//
TGraphErrors *Mspdh8Nov21=0,*Mspdh8aNov21=0,*Rspdh8Nov21=0,*Rspdh8aNov21=0;
TGraphErrors *MV0Ch8Nov21=0,*MV0Ch8aNov21=0,*RV0Ch8Nov21=0,*RV0Ch8aNov21=0;

TGraphErrors *Mspdh6=0,*Mspdh6a=0,*Rspdh6=0,*Rspdh6a=0;
TGraphErrors *Mspdh8=0,*Mspdh8a=0,*Rspdh8=0,*Rspdh8a=0;
TGraphErrors *MV0Ch8=0,*MV0Ch8a=0,*RV0Ch8=0,*RV0Ch8a=0;
TGraphErrors *MV0CRh8=0,*MV0CRh8a=0,*RV0CRh8=0,*RV0CRh8a=0;
TGraphErrors *Mtpch8=0,*Mtpch8a=0,*Rtpch8=0,*Rtpch8a=0;
TGraphErrors *Mspd366h8=0,*Mspd366h8a=0,*Rspd366h8=0,*Rspd366h8a=0;
TGraphErrors *MspdDPMJh2=0,*MspdDPMJh2a=0,*RspdDPMJh2=0,*RspdDPMJh2a=0;

TGraphErrors *Mspd4b=0,*Rspd4b=0;
TGraphErrors *MV0C4b=0,*RV0C4b=0;

TGraphErrors *MV0C29Nov4b=0,*RV0C29Nov4b=0;

TCanvas *cnvMult=0;
TCanvas *cnvMultR=0;
TH1F *frMlt,*frMltR=0;

void DrawMult()
{
  gStyle->SetOptStat(0);

  /*
  Mspd4b       =  CreateGraph("Mspd4b" ,4,Npart_spdSelNov21_4b,Npart_rms_spdSelNov21_4b, mltSPDh8_4b,mltSPDh8_4be ,kFALSE, kGreen, 20,1.2,1);
  Rspd4b       =  CreateGraph("Rspd4b" ,4,Npart_spdSelNov21_4b,Npart_rms_spdSelNov21_4b, mltSPDh8_4b,mltSPDh8_4be ,kTRUE,  kGreen, 20,1.2,1);
  //
  Mspdh8Nov21  = CreateGraph("Mspdh8Nov21" ,11,Npart_spdSelNov21 ,Npart_rms_spdSelNov21,mltSPDh8Nov21,mltSPDh8eNov21, kFALSE, kRed, 20,1.2,1);
  Mspdh8aNov21 = CreateGraph("Mspdh8aNov21",11,Npart2_spdSelNov21,Npart_rms_spdSelNov21,mltSPDh8Nov21,mltSPDh8eNov21, kFALSE, kRed, 24,1.42,2);
  Rspdh8Nov21  = CreateGraph("Rspdh8Nov21" ,11,Npart_spdSelNov21 ,Npart_rms_spdSelNov21,mltSPDh8Nov21,mltSPDh8eNov21, kTRUE,  kRed, 20,1.2,1);
  Rspdh8aNov21 = CreateGraph("Rspdh8aNov21",11,Npart2_spdSelNov21,Npart_rms_spdSelNov21,mltSPDh8Nov21,mltSPDh8eNov21, kTRUE,  kRed, 24,1.42,2);
  //
  MV0Ch8Nov21  = CreateGraph("MV0Ch8Nov21" ,11,Npart_V0CSelNov21 ,Npart_rms_V0CSelNov21,mltV0Ch8Nov21,mltV0Ch8eNov21, kFALSE, kBlue, 20,1.2, 1);
  MV0Ch8aNov21 = CreateGraph("MV0Ch8aNov21",11,Npart2_V0CSelNov21,Npart_rms_V0CSelNov21,mltV0Ch8Nov21,mltV0Ch8eNov21, kFALSE, kBlue, 24,1.42, 2);
  RV0Ch8Nov21  = CreateGraph("RV0Ch8Nov21" ,11,Npart_V0CSelNov21 ,Npart_rms_V0CSelNov21,mltV0Ch8Nov21,mltV0Ch8eNov21, kTRUE,  kBlue, 20,1.2, 1);
  RV0Ch8aNov21 = CreateGraph("RV0Ch8aNov21",11,Npart2_V0CSelNov21,Npart_rms_V0CSelNov21,mltV0Ch8Nov21,mltV0Ch8eNov21, kTRUE,  kBlue, 24,1.42, 2);
  //
  Mtpch8  = CreateGraph("Mtpch8" ,11,Npart_TPCSel ,Npart_rms_TPCSel,mltTPCh8,mltTPCh8e, kFALSE, kBlack, 23, 1.2, 1);
  Mtpch8a = CreateGraph("Mtpch8a",11,Npart2_TPCSel,Npart_rms_TPCSel,mltTPCh8,mltTPCh8e, kFALSE, kBlack, 32, 1.2, 2);
  Rtpch8  = CreateGraph("Rtpch8" ,11,Npart_TPCSel ,Npart_rms_TPCSel,mltTPCh8,mltTPCh8e, kTRUE,  kBlack, 23, 1.2, 1);
  Rtpch8a = CreateGraph("Rtpch8a",11,Npart2_TPCSel,Npart_rms_TPCSel,mltTPCh8,mltTPCh8e, kTRUE,  kBlack, 32, 1.2, 2);
  //  //
  */
  Mspd4b       =  CreateGraph("Mspd4b" ,4,Npart_Geom_4b,Npart_rms_Geom_4b, mltSPDh8_4b,mltSPDh8_4be ,kFALSE, kGreen, 20,1.2,1);
  Rspd4b       =  CreateGraph("Rspd4b" ,4,Npart_Geom_4b,Npart_rms_Geom_4b, mltSPDh8_4b,mltSPDh8_4be ,kTRUE,  kGreen, 20,1.2,1);
  //
  MV0C4b       =  CreateGraph("MV0C4b" ,4,Npart_Geom_4b,Npart_rms_Geom_4b, mltV0Ch8_4b,mltV0Ch8_4be ,kFALSE, kMagenta, 21,1.2,1);
  RV0C4b       =  CreateGraph("RV0C4b" ,4,Npart_Geom_4b,Npart_rms_Geom_4b, mltV0Ch8_4b,mltV0Ch8_4be ,kTRUE,  kMagenta, 21,1.2,1);
  //
  MV0C29Nov4b  =  CreateGraph("MV0C29Nov4b" ,4,Npart_Geom_4b,Npart_rms_Geom_4b, mltV0Ch11Nov29_4b,mltV0Ch11Nov29_4be ,kFALSE, kGreen+1, 26,1.2,1);
  RV0C29Nov4b  =  CreateGraph("RV0C29Nov4b" ,4,Npart_Geom_4b,Npart_rms_Geom_4b, mltV0Ch11Nov29_4b,mltV0Ch11Nov29_4be ,kTRUE,  kGreen+1, 26,1.2,1);
  //
  MzdcV0CNov24 = CreateGraph("MzdcV0CNov24" ,11,Npart_Geom ,Npart_rms_Geom,mltZDCV0Ch8Nov24RS,mltZDCV0Ch8eNov24RS, kFALSE, kBlack, 25,1.4,1);
  RzdcV0CNov24 = CreateGraph("RzdcV0CNov24" ,11,Npart_Geom ,Npart_rms_Geom,mltZDCV0Ch8Nov24RS,mltZDCV0Ch8eNov24RS, kTRUE,  kBlack, 25,1.4,1);
  //
  Mspdh8Nov21  = CreateGraph("Mspdh8Nov21" ,11,Npart_Geom ,Npart_rms_Geom,mltSPDh8Nov21,mltSPDh8eNov21, kFALSE, kRed, 20,1.4,1);
  Rspdh8Nov21  = CreateGraph("Rspdh8Nov21" ,11,Npart_Geom ,Npart_rms_Geom,mltSPDh8Nov21,mltSPDh8eNov21, kTRUE,  kRed, 20,1.4,1);
  //
  MV0Ch8Nov21  = CreateGraph("MV0Ch8Nov21" ,11,Npart_Geom ,Npart_rms_Geom,mltV0Ch8Nov21,mltV0Ch8eNov21, kFALSE, kBlue, 20,1.2, 1);
  RV0Ch8Nov21  = CreateGraph("RV0Ch8Nov21" ,11,Npart_Geom ,Npart_rms_Geom,mltV0Ch8Nov21,mltV0Ch8eNov21, kTRUE,  kBlue, 20,1.2, 1);
  //
  //
  MV0Ch11Nov27A  = CreateGraph("MV0Ch8Nov27A" ,11,Npart_Geom ,Npart_rms_Geom,mltV0Ch11Nov27,mltV0Ch11eNov27, kFALSE, kCyan, 20,1.2, 1);
  RV0Ch11Nov27A  = CreateGraph("RV0Ch8Nov27A" ,11,Npart_Geom ,Npart_rms_Geom,mltV0Ch11Nov27,mltV0Ch11eNov27, kTRUE,  kCyan, 20,1.2, 1);
  //
  MV0Ch11Nov29  = CreateGraph("MV0Ch8Nov29" ,11,Npart_Geom ,Npart_rms_Geom,mltV0Ch11Nov29,mltV0Ch11eNov29, kFALSE, kGreen+1, 22,1., 1);
  RV0Ch11Nov29  = CreateGraph("RV0Ch8Nov29" ,11,Npart_Geom ,Npart_rms_Geom,mltV0Ch11Nov29,mltV0Ch11eNov29, kTRUE,  kGreen+1, 22,1., 1);
  //
  Mtpch8  = CreateGraph("Mtpch8" ,11,Npart_Geom ,Npart_rms_Geom,mltTPCh8,mltTPCh8e, kFALSE, kBlack, 23, 1.2, 1);
  Rtpch8  = CreateGraph("Rtpch8" ,11,Npart_Geom ,Npart_rms_Geom,mltTPCh8,mltTPCh8e, kTRUE,  kBlack, 23, 1.2, 1);
  //
  Mspdh6  = CreateGraph("Mspdh6" ,11,Npart_spdSel ,Npart_rms_spdSel,mltSPDh6,mltSPDh6e, kFALSE, kRed, 29, 1.0,1);
  Mspdh6a = CreateGraph("Mspdh6a",11,Npart2_spdSel,Npart_rms_spdSel,mltSPDh6,mltSPDh6e, kFALSE, kRed, 30, 1.0,2);
  Rspdh6  = CreateGraph("Rspdh6" ,11,Npart_spdSel ,Npart_rms_spdSel,mltSPDh6,mltSPDh6e, kTRUE,  kRed, 29, 1.0,1);
  Rspdh6a = CreateGraph("Rspdh6a",11,Npart2_spdSel,Npart_rms_spdSel,mltSPDh6,mltSPDh6e, kTRUE,  kRed, 30, 1.0,2);
  //
  Mspdh8  = CreateGraph("Mspdh8" ,11,Npart_spdSel ,Npart_rms_spdSel,mltSPDh8,mltSPDh8e, kFALSE, kMagenta, 20,1.02,1);
  Mspdh8a = CreateGraph("Mspdh8a",11,Npart2_spdSel,Npart_rms_spdSel,mltSPDh8,mltSPDh8e, kFALSE, kMagenta, 24,1.02,2);
  Rspdh8  = CreateGraph("Rspdh8" ,11,Npart_spdSel ,Npart_rms_spdSel,mltSPDh8,mltSPDh8e, kTRUE,  kMagenta, 20,1.02,1);
  Rspdh8a = CreateGraph("Rspdh8a",11,Npart2_spdSel,Npart_rms_spdSel,mltSPDh8,mltSPDh8e, kTRUE,  kMagenta, 24,1.02,2);
  //
  MV0Ch8  = CreateGraph("MV0Ch8" ,11,Npart_V0CSel ,Npart_rms_V0CSel,mltV0Ch8,mltV0Ch8e, kFALSE, kCyan, 21,1.02, 1);
  MV0Ch8a = CreateGraph("MV0Ch8a",11,Npart2_V0CSel,Npart_rms_V0CSel,mltV0Ch8,mltV0Ch8e, kFALSE, kCyan, 25,1.02, 2);
  RV0Ch8  = CreateGraph("RV0Ch8" ,11,Npart_V0CSel ,Npart_rms_V0CSel,mltV0Ch8,mltV0Ch8e, kTRUE,  kCyan, 21,1.02, 1);
  RV0Ch8a = CreateGraph("RV0Ch8a",11,Npart2_V0CSel,Npart_rms_V0CSel,mltV0Ch8,mltV0Ch8e, kTRUE,  kCyan, 25,1.02, 2);
  //
  MV0CRh8  = CreateGraph("MV0CRh8" ,11,Npart_V0CRSel ,Npart_rms_V0CRSel,mltV0CRh8,mltV0CRh8e, kFALSE, kGreen, 22,1.2, 1);
  MV0CRh8a = CreateGraph("MV0CRh8a",11,Npart2_V0CRSel,Npart_rms_V0CRSel,mltV0CRh8,mltV0CRh8e, kFALSE, kGreen, 26,1.2, 2);
  RV0CRh8  = CreateGraph("RV0CRh8" ,11,Npart_V0CRSel ,Npart_rms_V0CRSel,mltV0CRh8,mltV0CRh8e, kTRUE,  kGreen, 22,1.2, 1);
  RV0CRh8a = CreateGraph("RV0CRh8a",11,Npart2_V0CRSel,Npart_rms_V0CRSel,mltV0CRh8,mltV0CRh8e, kTRUE,  kGreen, 26,1.2, 2);
  //
  Mspd366h8  = CreateGraph("Mspd366h8" ,11,Npart_spdSel ,Npart_rms_spdSel,mlt366SPDh8,mlt366SPDh8e, kFALSE, kRed, 27, 1.0, 1);
  Mspd366h8a = CreateGraph("Mspd366h8a",11,Npart2_spdSel,Npart_rms_spdSel,mlt366SPDh8,mlt366SPDh8e, kFALSE, kRed, 33, 1.0, 2);
  Rspd366h8  = CreateGraph("Rspd366h8" ,11,Npart_spdSel ,Npart_rms_spdSel,mlt366SPDh8,mlt366SPDh8e, kTRUE,  kRed, 27, 1.0, 1);
  Rspd366h8a = CreateGraph("Rspd366h8a",11,Npart2_spdSel,Npart_rms_spdSel,mlt366SPDh8,mlt366SPDh8e, kTRUE,  kRed, 33, 1.0, 2);
  //
  MspdDPMJh2  = CreateGraph("MspdDPMJh2" ,11,Npart_spdSel ,Npart_rms_spdSel,mltSPDdpmjh2,mltSPDdpmjh2e, kFALSE, kRed, 28, 1.0, 1);
  MspdDPMJh2a = CreateGraph("MspdDPMJh2a",11,Npart2_spdSel,Npart_rms_spdSel,mltSPDdpmjh2,mltSPDdpmjh2e, kFALSE, kRed, 34, 1.0, 2);
  RspdDPMJh2  = CreateGraph("RspdDPMJh2" ,11,Npart_spdSel ,Npart_rms_spdSel,mltSPDdpmjh2,mltSPDdpmjh2e, kTRUE,  kRed, 28, 1.0, 1);
  RspdDPMJh2a = CreateGraph("RspdDPMJh2a",11,Npart2_spdSel,Npart_rms_spdSel,mltSPDdpmjh2,mltSPDdpmjh2e, kTRUE,  kRed, 34, 1.0, 2);
  //
  TLegend *lg = new TLegend(0.4,0.2,0.65,0.5);
  lg->SetFillColor(kWhite);
  lg->SetHeader("N_{part} Glauber");
  //  lg->AddEntry(Rspdh8,"SPD corr","pl");
  //  lg->AddEntry(RV0Ch8,"V0 corr","pl");
  //  lg->AddEntry(RV0CRh8,"V0 corr resc","pl");
  //  lg->AddEntry(Rtpch8,"TPC","pl");
  //
  //  lg->AddEntry(Rspdh6,"SPD corr, old hijing","pl");
  //  lg->AddEntry(Rspd366h8,"SPD corr, 137366","pl");
  //  lg->AddEntry(RspdDPMJh2,"SPD corr, dpmjet","pl");
  //
  //  lg->AddEntry(Rspdh8Nov21,"SPD corr, 21/11/10","pl");
  lg->AddEntry(RV0Ch8Nov21,"V0  corr, 21/11/10","pl");
  //  lg->AddEntry(RzdcV0CNov24,"ZDC & V0 corr","pl");
  //  lg->AddEntry(Rspd4b,"SPD corr, 4#times2.5%","pl");
  lg->AddEntry(RV0C4b,"V0 corr, 4#times2.5%","pl");
  //  lg->AddEntry(RV0Ch11Nov27A,"V0 corr, 5+ 27/11/10","pl");
  lg->AddEntry(RV0Ch11Nov29,"V0 corr, 5+ 29/11/10","pl");
  lg->AddEntry(RV0C29Nov4b,"V0 corr, 5+ 29/11/10, 4#times2.5%","pl");
  //
  TLegend *lga = new TLegend(0.65,0.2,0.8,0.5);
  lga->SetFillColor(kWhite);
  lga->SetHeader("N_{part} data");
  lga->AddEntry(Rspdh8a," ","pl");
  lga->AddEntry(RV0Ch8a," ","pl");
  //  lga->AddEntry(RV0CRh8a, " ","pl");
  //  lga->AddEntry(Rtpch8a," ","pl");
  //
  //  lga->AddEntry(Rspdh6a," ","pl");
  //  lga->AddEntry(Rspd366h8a," ","pl");
  //  lga->AddEntry(RspdDPMJh2a," ","pl");
  //
  lga->AddEntry(Rspdh8aNov21,"SPD corr, 21/11/10","pl");
  lga->AddEntry(RV0Ch8aNov21,"V0  corr, 21/11/10","pl");
  //


  frMlt = new TH1F("dndeta","dN/d#eta",86,1,430);
  frMlt->SetMinimum(1e-3);
  frMlt->SetMaximum(2000);
  frMlt->GetXaxis()->SetTitle("N_{part}");
  frMlt->GetYaxis()->SetTitle("dN/d#eta_{|#eta|<0.5}");
  
  frMltR = new TH1F("dndetaperpart","dN/d#eta / part.pair",86,1,430);
  frMltR->SetMinimum(1e-3);
  frMltR->SetMaximum(10);
  frMltR->GetXaxis()->SetTitle("N_{part}");
  frMltR->GetYaxis()->SetTitle("dN/d#eta_{|#eta|<0.5} / part.pair");
  
  cnvMult  = new TCanvas("cnvMult","",10,10,1000,700);
  cnvMultR = new TCanvas("cnvMultR","",20,20,1000,700);
  //
  cnvMult->cd();
  gPad->Clear();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  frMlt->Draw();
  frMlt->GetYaxis()->SetTitleOffset(1.6);
  frMlt->GetYaxis()->SetLabelSize(0.03);
  frMlt->GetXaxis()->SetTitleOffset(1.6);
  frMlt->GetXaxis()->SetLabelSize(0.03);
  gPad->SetGrid(1,1);
  //
  //  Mspdh6->Draw("p");
  //  Mspdh6a->Draw("p");
  //
  //  Mspd366h8->Draw("p");
  //  Mspd366h8a->Draw("p");
  //
  //  MspdDPMJh2->Draw("p");
  //  MspdDPMJh2a->Draw("p");
  //
  //  Mspdh8->Draw("p");
  //  Mspdh8a->Draw("p");
  //
  //  Mspd4b->Draw("p");
  MV0C4b->Draw("p");
  MV0C29Nov4b->Draw("p");
  //
  //  MV0Ch8->Draw("p");
  //  MV0Ch8a->Draw("p");
  //
  //  MV0CRh8->Draw("p");
  //  MV0CRh8a->Draw("p");
  //
  //  Mtpch8->Draw("p");
  //  Mtpch8a->Draw("p");
  //
  //
  //  Mspdh8Nov21->Draw("p");
  //  Mspdh8aNov21->Draw("p");
  //
  MV0Ch8Nov21->Draw("p");
  //  MV0Ch8aNov21->Draw("p");
  //
  //  MzdcV0CNov24->Draw("p");
  //

  //  MV0Ch11Nov27A->Draw("p");
  MV0Ch11Nov29->Draw("p");

  lg->Draw();
  //  lga->Draw();
  //----------------------------------------
  //
  cnvMultR->cd();
  gPad->Clear();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  frMltR->Draw();
  frMltR->GetYaxis()->SetTitleOffset(1.6);
  frMltR->GetYaxis()->SetLabelSize(0.03);
  frMltR->GetXaxis()->SetTitleOffset(1.6);
  frMltR->GetXaxis()->SetLabelSize(0.03);
  gPad->SetGrid(1,1);
  //
  //  /*
  //  Rspdh6->Draw("p");
  //  Rspdh6a->Draw("p");
  //
  //  Rspd366h8->Draw("p");
  //  Rspd366h8a->Draw("p");
  //
  //  RspdDPMJh2->Draw("p");
  //  RspdDPMJh2a->Draw("p");
  //
  //  */
  //  Rspdh8->Draw("p");
  //  Rspdh8a->Draw("p");
  //
  //  Rspd4b->Draw("p");
  RV0C4b->Draw("p");
  //
  //  RV0Ch8->Draw("p");
  //  RV0Ch8a->Draw("p");
  //
  //  RV0CRh8->Draw("p");
  //  RV0CRh8a->Draw("p");
  //
  //  Rtpch8->Draw("p");
  //  Rtpch8a->Draw("p");
  //
  //  Rspdh8Nov21->Draw("p");
  //  Rspdh8aNov21->Draw("p");
  //
  RV0Ch8Nov21->Draw("p");
  //  RV0Ch8aNov21->Draw("p");
  //  RzdcV0CNov24->Draw("p");

  //  RV0Ch11Nov27A->Draw("p");
  RV0Ch11Nov29->Draw("p");
  RV0C29Nov4b->Draw("p");
  //
  lg->Draw();
  //  lga->Draw();
  //
}



TGraphErrors* CreateGraph(const char* nm, int np, const double *x,const double *xe, const double* y,const double *ye,Bool_t rat,int col,int mrk,double sz,int lst)
{
  TGraphErrors* gr = new TGraphErrors(np);
  gr->SetName(nm);
  gr->SetTitle(nm);
  gr->SetMarkerStyle(mrk);
  gr->SetMarkerSize(sz);
  gr->SetMarkerColor(col);
  gr->SetLineColor(col);
  gr->SetLineStyle(lst);
  int ip=0;
  for (int i=0;i<np;i++) {
    double px  = x[i];
    double pxe = xe[i];
    double py = y[i];
    double pye = ye[i];
    if (rat) {
      double rt,rte;
      GetRatE(py,pye, px/2, 0, rt,rte);
      py = rt;
      pye = rte;
    }
    gr->SetPoint(ip, px,py);
    gr->SetPointError(ip, pxe,pye);
    ip++;
  }
  //
  return gr;
}




//____________________________________________________________________________
void GetRatE(double x,double xe, double y,double ye, double &rat, double &rate)
{
  rat = 0; rate = 0;
  if (TMath::Abs(y)<1e-16 || TMath::Abs(x)<1e-16) return;
  rat = x/y;
  rate = rat*TMath::Sqrt( xe*xe/(x*x) + ye*ye/(y*y));
  //  printf("RAT %e %e / %e %e -> %e %e\n",x,xe,y,ye,rat,rate);
}

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 *                                                                        *
 * Author: Aditya Nath Mishra (amishra@cern.ch) Wigner RCP, Budapest      *
 *                                                                        *
 **************************************************************************/
/* This macro produces:  pT spectra in different multiplicity and Delta phi bins
New histograms for Dphi correlations added
   last update: 12/02/2021
*/


#include "AliAnalysisTaskUeSpectraDphi.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TFile.h>


// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <TTreeStream.h>
#include <AliHeader.h>
#include <AliAnalysisUtils.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliMultSelectionTask.h>
#include <AliAODInputHandler.h>
#include <AliAODHandler.h>
#include <AliAODVertex.h>
#include <AliAODTrack.h>
#include <AliAODPid.h>
#include <AliDataFile.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>

#include <iostream>
using namespace std;

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
ClassImp(AliAnalysisTaskUeSpectraDphi)

//_____________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::AliAnalysisTaskUeSpectraDphi():
AliAnalysisTaskSE(),
  fESD(0x0),
  fEventCuts(0x0),
  fMCEvent(0x0),
  fMCStack(0x0),
  fAnalysisMC(kFALSE),
  fTrackFilter(0x0),
  fTrackFilterDCA(0x0),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),
  fEtaCut(0.8),
  fRun(-999),
  fEventId(-999),
  fTriggeredEventMB(-999),
  fListOfObjects(0x0),
  fisPS(kFALSE),
  fisMCvtxInZcut(0x0),
  fisTracklet(kFALSE),
  fVtxBeforeCuts(0x0),
  fVtxAfterCuts(0x0),
  ptvstrackletsvsdcaData(0x0),
  ptvstrackletsvsdcaPrim(0x0),
  ptvstrackletsvsdcaDecs(0x0),
  ptvstrackletsvsdcaMatl(0x0),
  ptvstrackletsvsdcacentralData(0x0),
  ptvstrackletsvsdcacentralPrim(0x0),
  ptvstrackletsvsdcacentralDecs(0x0),
  ptvstrackletsvsdcacentralMatl(0x0),
  fdcaxy(-999),
  fdcaz(-999),
  hSelEv(0x0),
  hINEL0(0x0),
  hEvMultSel(0x0),
  hEvDphiSel(0x0),
  hPS(0x0),
  hVtxPS(0x0),
  hpT(0x0),
  hEta(0x0),
  hPhi(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  hMultvsDphiOA(0x0),
  hMultvsDphiSA(0x0),
  hMultvspTvsDphi(0x0),
  hMultvspTvsDphiWLP(0x0),
  isINEL0True(kFALSE),
  isINEL0Rec(kFALSE),
  hINEL0MCTrig(0x0),
  hINEL0MCTrue(0x0),
  hPS_MC(0x0),
  hVtxPS_MC(0x0),
  hpTMCTrue(0x0),
  hEtaMCTrue(0x0),
  hPhiMCTrue(0x0),
  hPtLMCTrue(0x0),
  hPhiLMCTrue(0x0),
  hEtaLMCTrue(0x0),
  hDphiMCTrue(0x0),
  hpTvsDphiOAMCTrue(0x0),
  hpTvsDphiSAMCTrue(0x0),
  hMultvsDphiOAMCTrue(0x0),
  hMultvsDphiSAMCTrue(0x0),
  hMultvspTvsDphiMCTrue(0x0),
  hMultvspTvsDphiWLPMCTrue(0x0),
  sigLossTrueINEL0(0x0),
  sigLossTrigINEL0(0x0),
  primaries(0x0),
  secondaries(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),
  fv0mpercentile(-999),
  hNchTSLeftvsNchTSRightvsDphi(0x0),
  hDPhiNchTSGT12(0x0),
  hDPhiNchTSLT8(0x0),
  hDPhiNchTSGT12LT12(0x0),
  hDPhiNchTSGT12LT8(0x0),
  hDPhiNchTSGT12GT0(0x0),
  hMultTSNchTSGT12(0x0),
  hMultTSNchTSLT8(0x0),
  hMultTSNchTSGT12GT0(0x0),
  hMultTSNchTSGT12LT12(0x0),
  hMultTSDNchTSGT12LT8(0x0)
{
  for(Int_t i=0; i<10; i++){
    hNchTSLeft[i]=0;
    hNchTSRight[i]=0;
    hSumptTSLeft[i]=0;
    hSumptTSRight[i]=0;
    hNchTSLeftvsNchTSRight[i]=0;
    hSumptTSLeftvsSumptTSRight[i]=0;
    hNchTSLeftvsNchTSRight[i]=0;
    hSumptTSLeftvsSumptTSRight[i]=0;
  }

  // Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::AliAnalysisTaskUeSpectraDphi(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fEventCuts(0x0),
  fMCEvent(0x0),
  fMCStack(0x0),
  fAnalysisMC(kFALSE),
  fTrackFilter(0x0),
  fTrackFilterDCA(0x0),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),
  fEtaCut(0.8),
  fRun(-999),
  fEventId(-999),
  fTriggeredEventMB(-999),
  fListOfObjects(0x0),
  fisPS(kFALSE),
  fisMCvtxInZcut(0x0),
  fisTracklet(kFALSE),
  fVtxBeforeCuts(0x0),
  fVtxAfterCuts(0x0),
  ptvstrackletsvsdcaData(0x0),
  ptvstrackletsvsdcaPrim(0x0),
  ptvstrackletsvsdcaDecs(0x0),
  ptvstrackletsvsdcaMatl(0x0),
  ptvstrackletsvsdcacentralData(0x0),
  ptvstrackletsvsdcacentralPrim(0x0),
  ptvstrackletsvsdcacentralDecs(0x0),
  ptvstrackletsvsdcacentralMatl(0x0),
  fdcaxy(-999),
  fdcaz(-999),
  hSelEv(0x0),
  hINEL0(0x0),
  hEvMultSel(0x0),
  hEvDphiSel(0x0),
  hpT(0x0),
  hEta(0x0),
  hPhi(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hMultvspTvsDphiWLP(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  hMultvsDphiOA(0x0),
  hMultvsDphiSA(0x0),
  hMultvspTvsDphi(0x0),
  isINEL0True(kFALSE),
  isINEL0Rec(kFALSE),
  hINEL0MCTrig(0x0),
  hINEL0MCTrue(0x0),
  hPS_MC(0x0),
  hVtxPS_MC(0x0),
  hpTMCTrue(0x0),
  hEtaMCTrue(0x0),
  hPhiMCTrue(0x0),
  hPtLMCTrue(0x0),
  hPhiLMCTrue(0x0),
  hEtaLMCTrue(0x0),
  hDphiMCTrue(0x0),
  hpTvsDphiOAMCTrue(0x0),
  hpTvsDphiSAMCTrue(0x0),
  hMultvsDphiOAMCTrue(0x0),
  hMultvsDphiSAMCTrue(0x0),
  hMultvspTvsDphiMCTrue(0x0),
  hMultvspTvsDphiWLPMCTrue(0x0),
  sigLossTrueINEL0(0x0),
  sigLossTrigINEL0(0x0),
  primaries(0x0),
  secondaries(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),
  fv0mpercentile(-999),
  hNchTSLeftvsNchTSRightvsDphi(0x0),
  hDPhiNchTSGT12(0x0),
  hDPhiNchTSLT8(0x0),
  hDPhiNchTSGT12LT12(0x0),
  hDPhiNchTSGT12LT8(0x0),
  hDPhiNchTSGT12GT0(0x0),
  hMultTSNchTSGT12(0x0),
  hMultTSNchTSLT8(0x0),
  hMultTSNchTSGT12GT0(0x0),
  hMultTSNchTSGT12LT12(0x0),
  hMultTSDNchTSGT12LT8(0x0)
{
  for(Int_t i=0; i<10; i++){
     hNchTSLeft[i]=0;
     hNchTSRight[i]=0;
     hSumptTSLeft[i]=0;
     hSumptTSRight[i]=0;
     hNchTSLeftvsNchTSRight[i]=0;
     hSumptTSLeftvsSumptTSRight[i]=0;
  }

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUeSpectraDphi::Exit(const char *msg) {

  Printf("%s", msg);
  return;
}


//_____________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::~AliAnalysisTaskUeSpectraDphi()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fListOfObjects;
    fListOfObjects = 0x0;
  }

}

//______________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::UserCreateOutputObjects()
{
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested
  // We also create the random generator here so it might get different seeds...

  fRandom = new TRandom(0); // 0 means random seed

  const Int_t nPtBins      = 79;
  Double_t PtBins[nPtBins+1] = {
    0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
    0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,
    1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
    4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20,22.0,
    24.0,26.0,28.0,32.0,36.0,42.0,50.0,60.0,80.0,100.0,130.0,
    160.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,700.0,800.0,1000.0};

  const Int_t nMultBins = 300;
  const Double_t MultBins[nMultBins+1] = {
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300};

  const Int_t nDphiBins = 180;
  const Double_t DphiBins[nDphiBins+1] = {
    0,0.01745,0.0349,0.05235,0.0698,0.08725,0.1047,0.12215,0.1396,0.15705,0.1745,0.19195,0.2094,0.22685,0.2443,0.26175,
    0.2792,0.29665,0.3141,0.33155,0.349,0.36645,0.3839,0.40135,0.4188,0.43625,0.4537,0.47115,0.4886,0.50605,0.5235,0.54095,
    0.5584,0.57585,0.5933,0.61075,0.6282,0.64565,0.6631,0.68055,0.698,0.71545,0.7329,0.75035,0.7678,0.78525,0.8027,0.82015,
    0.8376,0.85505,0.8725,0.88995,0.9074,0.92485,0.9423,0.95975,0.9772,0.99465,1.0121,1.02955,1.047,1.06445,1.0819,1.09935,
    1.1168,1.13425,1.1517,1.16915,1.1866,1.20405,1.2215,1.23895,1.2564,1.27385,1.2913,1.30875,1.3262,1.34365,1.3611,1.37855,
    1.396,1.41345,1.4309,1.44835,1.4658,1.48325,1.5007,1.51815,1.5356,1.55305,1.5705,1.58795,1.6054,1.62285,1.6403,1.65775,
    1.6752,1.69265,1.7101,1.72755,1.745,1.76245,1.7799,1.79735,1.8148,1.83225,1.8497,1.86715,1.8846,1.90205,1.9195,1.93695,
    1.9544,1.97185,1.9893,2.00675,2.0242,2.04165,2.0591,2.07655,2.094,2.11145,2.1289,2.14635,2.1638,2.18125,2.1987,2.21615,
    2.2336,2.25105,2.2685,2.28595,2.3034,2.32085,2.3383,2.35575,2.3732,2.39065,2.4081,2.42555,2.443,2.46045,2.4779,2.49535,
    2.5128,2.53025,2.5477,2.56515,2.5826,2.60005,2.6175,2.63495,2.6524,2.66985,2.6873,2.70475,2.7222,2.73965,2.7571,2.77455,
    2.792,2.80945,2.8269,2.84435,2.8618,2.87925,2.8967,2.91415,2.9316,2.94905,2.9665,2.98395,3.0014,3.01885,3.0363,3.05375,
    3.0712,3.08865,3.1061,3.12355,3.15
  };

  const Int_t nBinsDCAxy = 121;
  Double_t binsDCAxy[] = {-3.025,-2.975,-2.925,-2.875,-2.825,-2.775,-2.725,-2.675,-2.625,-2.575,-2.525,-2.475,-2.425,-2.375,-2.325,-2.275,-2.225,
			  -2.175,-2.125,-2.075,-2.025,-1.975,-1.925,-1.875,-1.825,-1.775,-1.725,-1.675,-1.625,-1.575,-1.525,-1.475,
			  -1.425,-1.375,-1.325,-1.275,-1.225,-1.175,-1.125,-1.075,-1.025,-0.975,-0.925,-0.875,-0.825,-0.775,-0.725,
			  -0.675,-0.625,-0.575,-0.525,-0.475,-0.425,-0.375,-0.325,-0.275,-0.225,-0.175,-0.125,-0.075,-0.025,0.025,
			  0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,
			  0.825,0.875,0.925,0.975,1.025,1.075,1.125,1.175,1.225,1.275,1.325,1.375,1.425,1.475,1.525,
			  1.575,1.625,1.675,1.725,1.775,1.825,1.875,1.925,1.975,2.025,2.075,2.125,2.175,2.225,2.275,
			  2.325,2.375,2.425,2.475,2.525,2.575,2.625,2.675,2.725,2.775,2.825,2.875,2.925,2.975,3.025};

  const Int_t nBinsDCAxyCentral = 461;
  Double_t binsDCAxyCentral[]= {-0.2305,-0.2295,-0.2285,-0.2275,-0.2265,-0.2255,-0.2245,-0.2235,-0.2225,-0.2215,-0.2205,-0.2195,-0.2185,-0.2175,-0.2165,-0.2155,-0.2145,
				-0.2135,-0.2125,-0.2115,-0.2105,-0.2095,-0.2085,-0.2075,-0.2065,-0.2055,-0.2045,-0.2035,-0.2025,-0.2015,-0.2005,-0.1995,
				-0.1985,-0.1975,-0.1965,-0.1955,-0.1945,-0.1935,-0.1925,-0.1915,-0.1905,-0.1895,-0.1885,-0.1875,-0.1865,-0.1855,-0.1845,
				-0.1835,-0.1825,-0.1815,-0.1805,-0.1795,-0.1785,-0.1775,-0.1765,-0.1755,-0.1745,-0.1735,-0.1725,-0.1715,-0.1705,-0.1695,
				-0.1685,-0.1675,-0.1665,-0.1655,-0.1645,-0.1635,-0.1625,-0.1615,-0.1605,-0.1595,-0.1585,-0.1575,-0.1565,-0.1555,-0.1545,
				-0.1535,-0.1525,-0.1515,-0.1505,-0.1495,-0.1485,-0.1475,-0.1465,-0.1455,-0.1445,-0.1435,-0.1425,-0.1415,-0.1405,-0.1395,
				-0.1385,-0.1375,-0.1365,-0.1355,-0.1345,-0.1335,-0.1325,-0.1315,-0.1305,-0.1295,-0.1285,-0.1275,-0.1265,-0.1255,-0.1245,
				-0.1235,-0.1225,-0.1215,-0.1205,-0.1195,-0.1185,-0.1175,-0.1165,-0.1155,-0.1145,-0.1135,-0.1125,-0.1115,-0.1105,-0.1095,
				-0.1085,-0.1075,-0.1065,-0.1055,-0.1045,-0.1035,-0.1025,-0.1015,-0.1005,-0.0995,-0.0985,-0.0975,-0.0965,-0.0955,-0.0945,
				-0.0935,-0.0925,-0.0915,-0.0905,-0.0895,-0.0885,-0.0875,-0.0865,-0.0855,-0.0845,-0.0835,-0.0825,-0.0815,-0.0805,-0.0795,
				-0.0785,-0.0775,-0.0765,-0.0755,-0.0745,-0.0735,-0.0725,-0.0715,-0.0705,-0.0695,-0.0685,-0.0675,-0.0665,-0.0655,-0.0645,
				-0.0635,-0.0625,-0.0615,-0.0605,-0.0595,-0.0585,-0.0575,-0.0565,-0.0555,-0.0545,-0.0535,-0.0525,-0.0515,-0.0505,-0.0495,
				-0.0485,-0.0475,-0.0465,-0.0455,-0.0445,-0.0435,-0.0425,-0.0415,-0.0405,-0.0395,-0.0385,-0.0375,-0.0365,-0.0355,-0.0345,
				-0.0335,-0.0325,-0.0315,-0.0305,-0.0295,-0.0285,-0.0275,-0.0265,-0.0255,-0.0245,-0.0235,-0.0225,-0.0215,-0.0205,-0.0195,
				-0.0185,-0.0175,-0.0165,-0.0155,-0.0145,-0.0135,-0.0125,-0.0115,-0.0105,-0.0095,-0.0085,-0.0075,-0.0065,-0.0055,-0.0045,
				-0.0035,-0.0025,-0.0015,-0.0005,0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095,0.0105,
				0.0115,0.0125,0.0135,0.0145,0.0155,0.0165,0.0175,0.0185,0.0195,0.0205,0.0215,0.0225,0.0235,0.0245,0.0255,
				0.0265,0.0275,0.0285,0.0295,0.0305,0.0315,0.0325,0.0335,0.0345,0.0355,0.0365,0.0375,0.0385,0.0395,0.0405,
				0.0415,0.0425,0.0435,0.0445,0.0455,0.0465,0.0475,0.0485,0.0495,0.0505,0.0515,0.0525,0.0535,0.0545,0.0555,
				0.0565,0.0575,0.0585,0.0595,0.0605,0.0615,0.0625,0.0635,0.0645,0.0655,0.0665,0.0675,0.0685,0.0695,0.0705,
				0.0715,0.0725,0.0735,0.0745,0.0755,0.0765,0.0775,0.0785,0.0795,0.0805,0.0815,0.0825,0.0835,0.0845,0.0855,
				0.0865,0.0875,0.0885,0.0895,0.0905,0.0915,0.0925,0.0935,0.0945,0.0955,0.0965,0.0975,0.0985,0.0995,0.1005,
				0.1015,0.1025,0.1035,0.1045,0.1055,0.1065,0.1075,0.1085,0.1095,0.1105,0.1115,0.1125,0.1135,0.1145,0.1155,
				0.1165,0.1175,0.1185,0.1195,0.1205,0.1215,0.1225,0.1235,0.1245,0.1255,0.1265,0.1275,0.1285,0.1295,0.1305,
				0.1315,0.1325,0.1335,0.1345,0.1355,0.1365,0.1375,0.1385,0.1395,0.1405,0.1415,0.1425,0.1435,0.1445,0.1455,
				0.1465,0.1475,0.1485,0.1495,0.1505,0.1515,0.1525,0.1535,0.1545,0.1555,0.1565,0.1575,0.1585,0.1595,0.1605,
				0.1615,0.1625,0.1635,0.1645,0.1655,0.1665,0.1675,0.1685,0.1695,0.1705,0.1715,0.1725,0.1735,0.1745,0.1755,
				0.1765,0.1775,0.1785,0.1795,0.1805,0.1815,0.1825,0.1835,0.1845,0.1855,0.1865,0.1875,0.1885,0.1895,0.1905,
				0.1915,0.1925,0.1935,0.1945,0.1955,0.1965,0.1975,0.1985,0.1995,0.2005,0.2015,0.2025,0.2035,0.2045,0.2055,
				0.2065,0.2075,0.2085,0.2095,0.2105,0.2115,0.2125,0.2135,0.2145,0.2155,0.2165,0.2175,0.2185,0.2195,0.2205,
				0.2215,0.2225,0.2235,0.2245,0.2255,0.2265,0.2275,0.2285,0.2295,0.2305};

  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested

  // Definition of trackcuts
  if(!fTrackFilter){
    fTrackFilter = new AliAnalysisFilter("trackFilter");
    SetTrackCuts(fTrackFilter);
  }

  if(!fTrackFilterDCA){
    fTrackFilterDCA = new AliAnalysisFilter("fTrackFilterDCA");
    SetTrackCutsDCA(fTrackFilterDCA);
  }
  // ----------------------------------------------------------------------------------------------------- //

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  //
  // Histograms
  //

  hSelEv = 0;
  hSelEv = new TH1D("hSelEv","Number of events; Events; Counts", 10, 0, 10);
  fListOfObjects->Add(hSelEv);

  hEvMultSel = new TH1D("hEvMultSel","Events in Mult Bins; Events; Counts", 10, 0, 10);
  fListOfObjects->Add(hEvMultSel);

  hEvDphiSel = new TH1D("hEvDphiSel","Events in Dphi Sel; Events; Counts", 10, 0, 10);
  fListOfObjects->Add(hEvDphiSel);

  fVtxBeforeCuts = 0;
  fVtxBeforeCuts = new TH1D("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 60, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);

  fVtxAfterCuts = 0;
  fVtxAfterCuts = new TH1D("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 60, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);

  hINEL0 = 0;
  hINEL0 = new TH1D("hINEL0","Number of events; Events; Counts", 10, 0, 10);
  fListOfObjects->Add(hINEL0);

  hRefMult08 = 0;
  hRefMult08 = new TH1D("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);N_{ch};count",nMultBins,MultBins);
  fListOfObjects->Add(hRefMult08);

  hV0Mmult = 0;
  hV0Mmult = new TH1D("hV0Mmult","V0M ;V0M percentile;count",110,0,110);
  fListOfObjects->Add(hV0Mmult);

  if (fAnalysisMC){
    hINEL0MCTrig = 0;
    hINEL0MCTrig = new TH1D("hINEL0MCTrig","Number of events; Events; Counts", 10, 0, 10);
    fListOfObjects->Add(hINEL0MCTrig);

    hINEL0MCTrue = 0;
    hINEL0MCTrue = new TH1D("hINEL0MCTrue","Number of events; Events; Counts", 10, 0, 10);
    fListOfObjects->Add(hINEL0MCTrue);

    hPS_MC = new TH1D("hPS_MC","",5,0,5);
    fListOfObjects->Add(hPS_MC);

    hVtxPS_MC = new TH1D("hVtxPS_MC","",5,0,5);
    fListOfObjects->Add(hVtxPS_MC);

    hpTMCTrue = new TH1D("hpTMCTrue","",nPtBins,PtBins);
    fListOfObjects->Add(hpTMCTrue);

    hEtaMCTrue = new TH1D("hEtaMCTrue","; #eta^{leading};counts",20,-1,1);
    fListOfObjects->Add(hEtaMCTrue);

    hPhiMCTrue = new TH1D("hPhiMCTrue", ";#phi (rad); count", 64,0,2.0*TMath::Pi());
    fListOfObjects->Add(hPhiMCTrue);

    hPtLMCTrue = 0;
    hPtLMCTrue = new TH1D("hPtLMCTrue",";#it{p}_{T}^{leading} (GeV/#it{c});counts",nPtBins,PtBins);
    fListOfObjects->Add(hPtLMCTrue);

    hEtaLMCTrue = 0;
    hEtaLMCTrue = new TH1D("hEtaLMCTrue","; #eta^{leading};counts",20,-1,1);
    fListOfObjects->Add(hEtaLMCTrue);

    hPhiLMCTrue = 0;
    hPhiLMCTrue = new TH1D("hPhiLMCTrue","; #phi^{leading} (rad);counts",64,0,2.0*TMath::Pi());
    fListOfObjects->Add(hPhiLMCTrue);

    hDphiMCTrue = 0;
    hDphiMCTrue = new TH1D("hDphiMCTrue","",64,-2*TMath::Pi(),2*TMath::Pi());
    fListOfObjects->Add(hDphiMCTrue);

    hpTvsDphiOAMCTrue = 0;
    hpTvsDphiOAMCTrue= new TH2D("hpTvsDphiOAMCTrue","p_{T} vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",360,-3.15,3.15,nPtBins,PtBins);
    fListOfObjects->Add(hpTvsDphiOAMCTrue);

    hpTvsDphiSAMCTrue = 0;
    hpTvsDphiSAMCTrue = new TH2D("hpTvsDphiSAMCTrue","p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",nDphiBins,DphiBins,nPtBins,PtBins);
    fListOfObjects->Add(hpTvsDphiSAMCTrue);

    hMultvsDphiOAMCTrue = 0;
    hMultvsDphiOAMCTrue= new TH2D("hMultvsDphiOAMCTrue","Multiplicity vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);Multiplicity",360,-3.15,3.15,nMultBins,MultBins);
    fListOfObjects->Add(hMultvsDphiOAMCTrue);

    hMultvsDphiSAMCTrue = 0;
    hMultvsDphiSAMCTrue = new TH2D("hMultvsDphiSAMCTrue","Multiplicity vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity",nDphiBins,DphiBins,nMultBins,MultBins);
    fListOfObjects->Add(hMultvsDphiSAMCTrue);

    hMultvspTvsDphiMCTrue = 0;
    hMultvspTvsDphiMCTrue = new TH3D("hMultvspTvsDphiMCTrue","Multiplicity vs p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
    fListOfObjects->Add(hMultvspTvsDphiMCTrue);

    hMultvspTvsDphiWLPMCTrue = 0;
    hMultvspTvsDphiWLPMCTrue = new TH3D("hMultvspTvsDphiWLPMCTrue","Multiplicity vs p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
    fListOfObjects->Add(hMultvspTvsDphiWLPMCTrue);

    sigLossTrueINEL0 = new TH3D("sigLossTrueINEL0","Multiplicity vs p_{T} vs #Delta#phi ;#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
    fListOfObjects->Add(sigLossTrueINEL0);

    sigLossTrigINEL0 = new TH3D("sigLossTrigINEL0","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
    fListOfObjects->Add(sigLossTrigINEL0);

    primaries = new TH3D("primaries","Multiplicity vs p_{T} vs #Delta#phi ;#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
    fListOfObjects->Add(primaries);

    secondaries = new TH3D("secondaries","Multiplicity vs p_{T} vs #Delta#phi ;#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
    fListOfObjects->Add(secondaries);

    // HISTOS FOR FEED_DOWN CORRECTION ------------------------------------------------------------------------
    ptvstrackletsvsdcaPrim = new TH3D("ptvstrackletsvsdcaPrim","pt vs tracklets vs dca Primaries",nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvstrackletsvsdcaPrim);

    ptvstrackletsvsdcaDecs = new TH3D("ptvstrackletsvsdcaDecs","pt vs tracklets vs dca Decays",nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvstrackletsvsdcaDecs);

    ptvstrackletsvsdcaMatl = new TH3D("ptvstrackletsvsdcaMatl","pt vs tracklets vs dca Material",nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvstrackletsvsdcaMatl);

    ptvstrackletsvsdcacentralPrim = new TH3D("ptvstrackletsvsdcacentralPrim","pt vs tracklets vs dca central Primaries",nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvstrackletsvsdcacentralPrim);

    ptvstrackletsvsdcacentralDecs = new TH3D("ptvstrackletsvsdcacentralDecs","pt vs tracklets vs dca central Decays",nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvstrackletsvsdcacentralDecs);

    ptvstrackletsvsdcacentralMatl = new TH3D("ptvstrackletsvsdcacentralMatl","pt vs tracklets vs dca central Material",nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvstrackletsvsdcacentralMatl);
  }

  else // for the data
    {

      hPS = new TH1D("hPS","",5,0,5);
      fListOfObjects->Add(hPS);

      hVtxPS = new TH1D("hVtxPS","",5,0,5);
      fListOfObjects->Add(hVtxPS);

      // histos for FEED_DOWN CORRECTION ----------------------------------------------------------------------------
      ptvstrackletsvsdcaData = new TH3D("ptvstrackletsvsdcaData","pt vs tracklets vs dca Data", nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxy,binsDCAxy);
      fListOfObjects->Add(ptvstrackletsvsdcaData);

      ptvstrackletsvsdcacentralData = new TH3D("ptvstrackletsvsdcacentralData","pt vs tracklets vs dca central Data",nPtBins,PtBins,nMultBins,MultBins,nBinsDCAxyCentral,binsDCAxyCentral);
      fListOfObjects->Add(ptvstrackletsvsdcacentralData);
    }


  hpT = new TH1D("hpT","",nPtBins,PtBins);
  fListOfObjects->Add(hpT);

  hEta = new TH1D("hEta","; #eta^{leading};counts",20,-1,1);
  fListOfObjects->Add(hEta);

  hPhi = new TH1D("hPhi", ";#phi (rad); count", 64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hPhi);

  hPtL = 0;
  hPtL = new TH1D("hPtL",";#it{p}_{T}^{leading} (GeV/#it{c});counts",nPtBins,PtBins);
  fListOfObjects->Add(hPtL);

  hEtaL = 0;
  hEtaL = new TH1D("hEtaL","; #eta^{leading};counts",20,-1,1);
  fListOfObjects->Add(hEtaL);

  hPhiL = 0;
  hPhiL = new TH1D("hPhiL","; #phi^{leading} (rad);counts",64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hPhiL);

  hDphi = 0;
  hDphi = new TH1D("hDphi","",64,-2*TMath::Pi(),2*TMath::Pi());
  fListOfObjects->Add(hDphi);

  hDPhiNchTSGT12        = new TH1D("hDPhiNchTSGT12",";#Delta#phi (rad.);N_{ch} (|#eta| < 0.8)",64,-2.0*TMath::Pi(),2.0*TMath::Pi());
  hDPhiNchTSLT8         = new TH1D("hDPhiNchTSLT8",";#Delta#phi (rad.);N_{ch} (|#eta| < 0.8)",64,-2.0*TMath::Pi(),2.0*TMath::Pi());
  hDPhiNchTSGT12GT0     = new TH1D("hDPhiNchTSGT12GT0",";#Delta#phi (rad.);N_{ch} (|#eta| < 0.8)",64,-2.0*TMath::Pi(),2.0*TMath::Pi());
  hDPhiNchTSGT12LT12    = new TH1D("hDPhiNchTSGT12LT12",";#Delta#phi (rad.);N_{ch} (|#eta| < 0.8)",64,-2.0*TMath::Pi(),2.0*TMath::Pi());
  hDPhiNchTSGT12LT8     = new TH1D("hDPhiNchTSGT12LT8",";#Delta#phi (rad.);N_{ch} (|#eta| < 0.8)",64,-2.0*TMath::Pi(),2.0*TMath::Pi());

  fListOfObjects->Add(hDPhiNchTSGT12);
  fListOfObjects->Add(hDPhiNchTSLT8);
  fListOfObjects->Add(hDPhiNchTSGT12GT0);
  fListOfObjects->Add(hDPhiNchTSGT12LT12);
  fListOfObjects->Add(hDPhiNchTSGT12LT8);
  
  hMultTSNchTSGT12      = new TH2D("hMultTSNchTSGT12",";N_{ch}^{TS-left} (|#eta| < 0.8);N_{ch}^{TS-right} (|#eta| < 0.8)",nMultBins,MultBins,nMultBins,MultBins);
  hMultTSNchTSLT8       = new TH2D("hMultTSNchTSLT8",";N_{ch}^{TS-left} (|#eta| < 0.8);N_{ch}^{TS-right} (|#eta| < 0.8)",nMultBins,MultBins,nMultBins,MultBins);
  hMultTSNchTSGT12GT0   = new TH2D("hMultTSNchTSGT12GT0",";N_{ch}^{TS-left} (|#eta| < 0.8);N_{ch}^{TS-right} (|#eta| < 0.8)",nMultBins,MultBins,nMultBins,MultBins); 
  hMultTSNchTSGT12LT12  = new TH2D("hMultTSNchTSGT12LT12",";N_{ch}^{TS-left} (|#eta| < 0.8);N_{ch}^{TS-right} (|#eta| < 0.8)",nMultBins,MultBins,nMultBins,MultBins);
  hMultTSDNchTSGT12LT8  = new TH2D("hMultTSDNchTSGT12LT8",";N_{ch}^{TS-left} (|#eta| < 0.8);N_{ch}^{TS-right} (|#eta| < 0.8)",nMultBins,MultBins,nMultBins,MultBins);
  
  fListOfObjects->Add(hMultTSNchTSGT12);
  fListOfObjects->Add(hMultTSNchTSLT8);
  fListOfObjects->Add(hMultTSNchTSGT12GT0);
  fListOfObjects->Add(hMultTSNchTSGT12LT12);
  fListOfObjects->Add(hMultTSDNchTSGT12LT8);
    
  hpTvsDphiOA = 0;
  hpTvsDphiOA = new TH2D("hpTvsDphiOA","p_{T} vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",360,-3.15,3.15,nPtBins,PtBins);
  fListOfObjects->Add(hpTvsDphiOA);
  
  hpTvsDphiSA = 0;
  hpTvsDphiSA = new TH2D("hpTvsDphiSA","p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",nDphiBins,DphiBins,nPtBins,PtBins);
  fListOfObjects->Add(hpTvsDphiSA);
  
  hMultvsDphiOA = 0;
  hMultvsDphiOA= new TH2D("hMultvsDphiOA","Multiplicity vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);Multiplicity",360,-3.15,3.15,nMultBins,MultBins);
  fListOfObjects->Add(hMultvsDphiOA);

  hMultvsDphiSA = 0;
  hMultvsDphiSA = new TH2D("hMultvsDphiSA","Multiplicity vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity",nDphiBins,DphiBins,nMultBins,MultBins);
  fListOfObjects->Add(hMultvsDphiSA);

  hMultvspTvsDphi = 0;
  hMultvspTvsDphi = new TH3D("hMultvspTvsDphi","Multiplicity vs p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
  fListOfObjects->Add(hMultvspTvsDphi);

  hMultvspTvsDphiWLP = 0;
  hMultvspTvsDphiWLP = new TH3D("hMultvspTvsDphiWLP","Multiplicity vs p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
  fListOfObjects->Add(hMultvspTvsDphiWLP);

  hNchTSLeftvsNchTSRightvsDphi= new TH3D("hNchTSLeftvsNchTSRightvsDphi",";N_{ch}^{TS-left} (|#eta| < 0.8);N_{ch}^{TS-right} (|#eta| < 0.8);#Delta#phi",200,0,200,200,0,200,64,-2.0*TMath::Pi(),-2.0*TMath::Pi());
  fListOfObjects->Add(hNchTSLeftvsNchTSRightvsDphi);
  
  for(Int_t i=0; i<10; i++){
    hNchTSLeft[i]               = new TH1D(Form("hNchTSLeft%d",i),";N_{ch}^{TS-left} (|#eta| < 0.8);",nMultBins,MultBins);
    hNchTSRight[i]              = new TH1D(Form("hNchTSRight%d",i),";N_{ch}^{TS-right} (|#eta| < 0.8);",nMultBins,MultBins);
    hSumptTSLeft[i]             = new TH1D(Form("hSumptTSLeft%d",i),";#sump_{T}^{TS-left} (GeV/c);",nPtBins,PtBins);
    hSumptTSRight[i]            = new TH1D(Form("hSumptTSRight%d",i),";#sump_{T}^{TS-right} (GeV/c);",nPtBins,PtBins);  
    hNchTSLeftvsNchTSRight[i]   = new TH2D(Form("hNchTSLeftvsNchTSRight%d",i),";N_{ch}^{TS-left};N_{ch}^{TS-right}",nMultBins,MultBins,nMultBins,MultBins);
    hSumptTSLeftvsSumptTSRight[i]= new TH2D(Form("hSumptTSLeftvsSumptTSRight%d",i),";#sump_{T}{TS-left};#sump_{T}^{TS-right}",nPtBins,PtBins,nPtBins,PtBins);

    fListOfObjects->Add(hNchTSLeft[i]);
    fListOfObjects->Add(hNchTSRight[i]);
    fListOfObjects->Add(hSumptTSLeft[i]);
    fListOfObjects->Add(hSumptTSRight[i]);
    fListOfObjects->Add(hNchTSLeftvsNchTSRight[i]);
    fListOfObjects->Add(hSumptTSLeftvsSumptTSRight[i]);
  }

  // fEventCuts.AddQAplotsToList(fListOfObjects);
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::UserExec(Option_t *)
{

  // -----------------------------------------------------
  //			 InputEvent
  // -----------------------------------------------------

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  // -----------------------------------------------------
  //			 E S D
  // -----------------------------------------------------
  fESD = dynamic_cast<AliESDEvent*>(event);

  if(!fESD){
    Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }

  // -----------------------------------------------------
  //			 MC
  // -----------------------------------------------------

  if (fAnalysisMC){
    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent)  {
      Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }

    fMCStack = fMCEvent->Stack();
    if(!fMCStack){
      Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }


    AliHeader *header = fMCEvent->Header();
    if(!header) {AliDebug( AliLog::kError , "Header not avaible" ); return; }
  }


  /************ BEGINNING OF EVENT SELECTION *******************/
  // Get trigger decision
  fTriggeredEventMB = 0; //init
  if (!fAnalysisMC){  // for data check event selection as well
    if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit )  fTriggeredEventMB = 1;  //event triggered as minimum bias
  }
  else {
    if (ftrigBit)  fTriggeredEventMB = 1;
  }

  Bool_t SPDvsClustersBG = kFALSE;

  AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
  if (!AnalysisUtils)
    {
      cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
      return;
    }
  else SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG

  Bool_t isNotPileUp = !fESD->IsPileupFromSPD(5,0.8);
  Bool_t IncompleteDAQ = fESD->IsIncompleteDAQ(); // we want is not incomplete DAQ

  // vertex
  const AliESDVertex * vertex    =    fESD->GetPrimaryVertex(); // tracks vertex, if not -> spd vertex, if not TPC vertex
  Bool_t isVtxGood = vertex->GetStatus() && selectVertex2015pp( fESD ,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm
  Double_t vertex_z = vertex->GetZ();
  Bool_t isVtxInZCut = (TMath::Abs(vertex_z) <= fVtxCut); // Zvtx in +- 10

  // Implement INEL>0
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  Bool_t isINEL0 = kFALSE;
  for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i){ if (TMath::Abs(mult->GetEta(i)) < 1.) isINEL0 = kTRUE;}

  /********** IS PHYSICS SELECTION FLAG ****************************/
  fisPS = fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp;

  // recontructed INEL > 0 is PS + vtx + Zvtx inside +-10 ------
  isINEL0Rec = kFALSE;
  if ( isINEL0 && fisPS && isVtxGood && isVtxInZCut) isINEL0Rec = kTRUE;

  if (fisPS) fVtxBeforeCuts->Fill(vertex_z);              // VZ hack
  if (isINEL0Rec) fVtxAfterCuts->Fill(vertex_z);
  hSelEv->Fill(0); // all events
  if(fTriggeredEventMB) hSelEv->Fill(1); // triggered events
  if (fTriggeredEventMB && !IncompleteDAQ ) hSelEv->Fill(2); // trigger + IsIncompleteDAQ
  if (fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG) hSelEv->Fill(3); // trigger + IsIncompleteDAQ + BG rejection
  if(fPileUpRej)
    {
      if(fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp) //
	hSelEv->Fill(4); // trigger + IsIncompleteDAQ + BG rejection + PileUp
    }
  if (fisPS)hSelEv->Fill(5); //PS: trigger + IsIncompleteDAQ + BG rejection + PileUp + 1 Tracklet in eta +-1
  if (fisPS && isVtxGood) hSelEv->Fill(6); //PS + GetPrimaryVertex
  if (isINEL0Rec) hSelEv->Fill(7); //PS + GetPrimaryVertex + isVtxInZCut

  // -------------------------------------- multiplcity estimators section ------------------------------------------ //
  ftrackmult08 = -999;
  fv0mpercentile = -999;

  //ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //reference
  ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);     //tracklets
  //ftrackmult08 = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD); //Combined estimator

  hRefMult08->Fill(ftrackmult08);

  fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
  if (!fMultSelection)
    cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
  else
    fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  hV0Mmult->Fill(fv0mpercentile);

  // cout<<"------- V0M mult ==  "<<fv0mpercentile<<"--------"<<endl;

  // ------------------------------------------ end of mult estimators -------------------------------------------------//
  if (fAnalysisMC){ // analysis for generated MC
    const AliVVertex *vertexMC = (AliVVertex*) fMCEvent->GetPrimaryVertex();
    fisMCvtxInZcut     = (TMath::Abs(vertexMC->GetZ()) <= fVtxCut);   // ZMCvtx in +- 10
    isINEL0True = isMCEventTrueINEL0(fMCEvent);

    // for trigger efficiency PS / INEL > 0 true
    if ( fisPS ){
      hINEL0MCTrig->Fill(0); // for trigger efficiency Trig/True
      hPS_MC->Fill(0);       // for missing vtx correction hPS_MC/hVtxPS_MC
    }
    if (fisPS && isVtxGood) hVtxPS_MC->Fill(0);
    if (isINEL0True) hINEL0MCTrue->Fill(0);
    AnalyzeMC(fMCEvent); // analysis for MC

  }
  else
    {
      if (isINEL0Rec) hINEL0->Fill(0);
      // Two histos for missing vtx correction
      if ( fisPS ) hPS->Fill(0);
      if ( fisPS && isVtxGood ) hVtxPS->Fill(0);
    }

  hEvDphiSel ->Fill(0);
  if (isINEL0Rec) {
    AnalyzeESD(fESD);
    AnalyzeESDforDCA(fESD);
  }
  //  cout<<"hello!!!"<<endl;

  // Post output data.
  PostData(1, fListOfObjects);

}
//________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::AnalyzeMC(AliMCEvent* fMCEvent){

  // selection on leading particle
  Double_t pt_leading    = 0;
  Double_t eta_leading    = 0;
  Double_t phi_leading    = 0;
  Int_t    i_leading = 0;
  Int_t mult = 0;

  for ( int iT = 0 ; iT < fMCStack->GetNtrack(); iT++ ){ // loop over TRUE MC

    TParticle *mcParticle = fMCStack->Particle(iT);

    if (!mcParticle){
      cout<<"no mcParticle"<<endl;
      continue;
    }

    if(!(fMCStack->IsPhysicalPrimary(iT))) continue;

    Double_t eta = mcParticle->Eta();
    Double_t pt = mcParticle->Pt();
    Double_t phi = mcParticle->Phi();

    int partPDG = TMath::Abs(mcParticle->GetPdgCode());
    if ( TMath::Abs(eta) > fEtaCut ) continue;
    if ( pt < 0.15 ) continue;
    //if (isINEL0Rec) {
    if ( isINEL0True && fisMCvtxInZcut){
      if ( !(TMath::Abs(mcParticle->GetPDG()->Charge() ) == 3) ) continue;
      mult++;

      if(pt>pt_leading){
	pt_leading      = pt;
	eta_leading     = eta;
	phi_leading     = phi;
	i_leading = iT;
      }
      hpTMCTrue->Fill(pt);
      hEtaMCTrue->Fill(eta);
      hPhiMCTrue->Fill(phi);
    }
  }// end loop over tracks

  hPtLMCTrue->Fill(pt_leading);
  hEtaLMCTrue->Fill(eta_leading);
  hPhiLMCTrue->Fill(phi_leading);

  for ( int iT = 0 ; iT < fMCStack->GetNtrack(); iT++ ){ // loop over TRUE MC

    TParticle *mcParticle = fMCStack->Particle(iT);

    if (!mcParticle){
      cout<<"no mcParticle"<<endl;
      continue;
    }

    if(!fMCStack->IsPhysicalPrimary(iT)) continue;

    Double_t eta = mcParticle->Eta();
    Double_t pt = mcParticle->Pt();
    Double_t phi = mcParticle->Phi();

    int partPDG = TMath::Abs(mcParticle->GetPdgCode());
    if ( TMath::Abs(eta) > fEtaCut ) continue;
    if ( pt < 0.15 ) continue;

    Double_t DPhiOA = DeltaPhiOA( phi, phi_leading );
    Double_t DPhiSA = TMath::Abs(DPhiOA);
    //if (isINEL0Rec) {
      if ( isINEL0True && fisMCvtxInZcut ){
      if (!(TMath::Abs(mcParticle->GetPDG()->Charge() ) == 3) ) continue;
      hDphiMCTrue->Fill(DPhiOA);
      hpTvsDphiOAMCTrue->Fill(DPhiOA,pt);
      hpTvsDphiSAMCTrue->Fill(DPhiSA,pt);
      hMultvsDphiOAMCTrue->Fill(DPhiOA,mult);
      hMultvsDphiSAMCTrue->Fill(DPhiSA,mult);
      hMultvspTvsDphiMCTrue->Fill(DPhiSA,mult,pt);
      if(DPhiSA > 0) hMultvspTvsDphiWLPMCTrue->Fill(DPhiSA,mult,pt);
    }
    if ( !( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) ) continue;
    if ( isINEL0Rec ) sigLossTrigINEL0->Fill(DPhiSA,mult,pt);
    if ( isINEL0True && fisMCvtxInZcut) sigLossTrueINEL0->Fill(DPhiSA,mult,pt);
  }// end loop over tracks
}
//________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::AnalyzeESD(AliESDEvent* fESD){

  Int_t MultBinsMin[9] = {0, 0,11,21,31,41,51,61,91};
  Int_t MultBinsMax[9] = {500,10,20,30,40,50,60,500,500};

  hINEL0->Fill(1);

  fRun  = fESD->GetRunNumber();
  fEventId = 0;
  if(fESD->GetHeader()) fEventId = GetEventIdAsLong(fESD->GetHeader());

  // selection on leading particle
  Double_t pt_leading    = 0;
  Double_t eta_leading    = 0;
  Double_t phi_leading    = 0;
  Int_t    i_leading = 0;
  Int_t mult=0, nchTSRight=0, nchTSLeft=0, nchTS = 0;
  Double_t sumpt=0., sumptTSRight=0., sumptTSLeft=0., sumptTS=0.;

  if (fAnalysisMC){ // for reconstructed MC
    for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

      AliESDtrack* esdTrack = fESD->GetTrack(i);
      if(!esdTrack) continue;
      Double_t eta      = esdTrack->Eta();
      Double_t phi      = esdTrack->Phi();
      Double_t pt       = esdTrack->Pt();

      if(TMath::Abs(eta) > fEtaCut) continue;
      //quality cuts, standard 2015 track cuts
      if(!fTrackFilter->IsSelected(esdTrack)) continue;
      if(pt<0.15) continue;
      mult++;

      Int_t mcLabel = TMath::Abs(esdTrack->GetLabel());
      TParticle *mcParticle = fMCEvent->GetTrack(mcLabel)->Particle();

      if(!mcParticle) {
	printf("ERROR: mcParticle not available-------\n");	\
	continue;
      }

      eta = mcParticle->Eta(); // generated eta and pT used intead of recontructed
      pt = mcParticle->Pt();
      phi = mcParticle->Phi();
      if(TMath::Abs(eta) > fEtaCut) continue;
      if(pt<0.15) continue;

      Int_t partPDG = TMath::Abs(mcParticle->GetPdgCode());
      if ((TMath::Abs(mcParticle->GetPDG()->Charge()) == 3)){ // only for charged particles
	if(pt>pt_leading){
	  pt_leading      = pt;
	  eta_leading     = eta;
	  phi_leading     = phi;
	  i_leading = i;
	}

	hpT->Fill(pt);
	hEta->Fill(eta);
	hPhi->Fill(phi);
      }
    }// end loop over tracks

    hPtL->Fill(pt_leading);
    hEtaL->Fill(eta_leading);
    hPhiL->Fill(phi_leading);

    for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

      AliESDtrack* esdTrack = fESD->GetTrack(i);
      Double_t eta      = esdTrack->Eta();
      Double_t phi      = esdTrack->Phi();
      Double_t pt       = esdTrack->Pt();

      if(TMath::Abs(eta) > fEtaCut) continue;
      //quality cuts, standard 2015 track cuts
      if(!fTrackFilter->IsSelected(esdTrack)) continue;
      if(pt<0.15) continue;

      Int_t mcLabel = TMath::Abs(esdTrack->GetLabel());
      TParticle *mcParticle = fMCEvent->GetTrack(mcLabel)->Particle();

      if(!mcParticle) {
	printf("ERROR: mcParticle not available-------\n");	\
	continue;
      }

      eta = mcParticle->Eta(); // generated eta and pT used intead of recontructed
      pt = mcParticle->Pt();
      phi = mcParticle->Phi();
      if(TMath::Abs(eta) > fEtaCut) continue;
      if(pt<0.15) continue;

      Int_t partPDG = TMath::Abs(mcParticle->GetPdgCode());
      if ((TMath::Abs(mcParticle->GetPDG()->Charge()) == 3)){ // only for charged particles

	Double_t DPhiOA = DeltaPhiOA( phi, phi_leading );
	Double_t DPhiSA = TMath::Abs(DPhiOA);

	primaries->Fill(DPhiSA,mult,pt);
	if (!(fMCEvent->IsPhysicalPrimary(mcLabel))) secondaries->Fill(DPhiSA,mult,pt); // secondary particles
	else{
	  hDphi->Fill(DPhiOA);
	  hpTvsDphiOA->Fill(DPhiOA,pt);
	  hpTvsDphiSA->Fill(DPhiSA,pt);
	  hMultvsDphiOA->Fill(DPhiOA,mult);
	  hMultvsDphiSA->Fill(DPhiSA,mult);
	  hMultvspTvsDphi->Fill(DPhiSA,mult,pt);
	  if(DPhiSA > 0) hMultvspTvsDphiWLP->Fill(DPhiSA,mult,pt);
	}
      }
    }// end loop over tracks
  }

  else{//data
    for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

      AliESDtrack* esdTrack = fESD->GetTrack(i);

      Double_t eta      = esdTrack->Eta();
      Double_t phi      = esdTrack->Phi();
      Double_t pt       = esdTrack->Pt();

      if(TMath::Abs(eta) > fEtaCut) continue;
      //quality cuts, standard 2015 track cuts
      if(!fTrackFilter->IsSelected(esdTrack)) continue;
      if(pt<0.15) continue;
      mult++;
      sumpt+=pt;

      Double_t DPhi   = DeltaPhi( phi, phi_leading );
      if( DPhi <= -TMath::Pi()/3.0 || DPhi >= 4*TMath::Pi()/3.0){sumptTSLeft+=pt; nchTSLeft++;}
      if( DPhi >=  TMath::Pi()/3.0 && DPhi <= 2*TMath::Pi()/3.0){sumptTSRight+=pt; nchTSRight++;}

      if(pt>pt_leading){
	pt_leading      = pt;
	eta_leading     = eta;
	phi_leading     = phi;
	i_leading = i;
      }

      hDphi->Fill(DPhi);
      hpT->Fill(pt);
      hEta->Fill(eta);
      hPhi->Fill(phi);

    }// end loop over tracks

    hPtL->Fill(pt_leading);
    hEtaL->Fill(eta_leading);
    hPhiL->Fill(phi_leading);
  

    for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

      AliESDtrack* esdTrack = fESD->GetTrack(i);
      Double_t eta      = esdTrack->Eta();
      Double_t phi      = esdTrack->Phi();
      Double_t pt       = esdTrack->Pt();

   
      if(TMath::Abs(eta) > fEtaCut) continue;
      //quality cuts, standard 2015 track cuts
      if(!fTrackFilter->IsSelected(esdTrack)) continue;
      if(pt<0.15) continue;

      Double_t DPhi   = DeltaPhi( phi, phi_leading );
      Double_t DPhiOA = DeltaPhiOA( phi, phi_leading );
      Double_t DPhiSA = TMath::Abs(DPhiOA);

      hNchTSLeftvsNchTSRightvsDphi->Fill(nchTSLeft,nchTSRight,DPhi);
      hpTvsDphiOA->Fill(DPhiOA,pt);
      hpTvsDphiSA->Fill(DPhiSA,pt);
      hMultvsDphiOA->Fill(DPhiOA,mult);
      hMultvsDphiSA->Fill(DPhiSA,mult);
      hMultvspTvsDphi->Fill(DPhiSA,mult,pt);
      if(DPhiSA > 0) hMultvspTvsDphiWLP->Fill(DPhiSA,mult,pt);

      hEvDphiSel ->Fill(1);
      
      if(mult<60) continue;
      hEvDphiSel ->Fill(2);
      if(nchTSLeft >= 12  && nchTSRight >= 12) {
	hEvDphiSel ->Fill(3);
	Double_t DeltaphiNchTSGT12 = DeltaPhi( phi, phi_leading );
	hDPhiNchTSGT12->Fill(DeltaphiNchTSGT12);
	hMultTSNchTSGT12->Fill(nchTSLeft,nchTSRight);
      }
      if(nchTSLeft <=  8  && nchTSRight <=  8) {
	hEvDphiSel ->Fill(4);
	Double_t DeltaphiNchTSLT8 = DeltaPhi( phi, phi_leading );
	hDPhiNchTSLT8->Fill(DeltaphiNchTSLT8);
	hMultTSNchTSLT8->Fill(nchTSLeft,nchTSRight);
      }
      
      if((nchTSLeft >= 12  && nchTSRight <= 12) || (nchTSRight >= 12  && nchTSLeft <= 12))  {
	hEvDphiSel ->Fill(5);
	Double_t DPhiNchTSGT12LT12 = DeltaPhi( phi, phi_leading );
	hDPhiNchTSGT12LT12->Fill(DPhiNchTSGT12LT12);
	hMultTSNchTSGT12LT12->Fill(nchTSLeft,nchTSRight);
      }
      
      if((nchTSLeft >= 12  && nchTSRight >= 0) || (nchTSRight >= 12  && nchTSLeft >=0))  {
	hEvDphiSel ->Fill(6);
	Double_t DPhiNchTSGT12GT0 = DeltaPhi( phi, phi_leading );
	hDPhiNchTSGT12GT0->Fill(DPhiNchTSGT12GT0);
	hMultTSNchTSGT12GT0->Fill(nchTSLeft,nchTSRight);
      }
      
      if((nchTSLeft >= 12  && nchTSRight <= 8) || (nchTSRight >= 12  && nchTSLeft <= 8))  {
	hEvDphiSel ->Fill(7);
	Double_t DPhiNchTSGT12LT8 = DeltaPhi( phi, phi_leading );
	hDPhiNchTSGT12LT8->Fill(DPhiNchTSGT12LT8);
	hMultTSDNchTSGT12LT8->Fill(nchTSLeft,nchTSRight);
      }
    }// end loop over tracks
    
    for(Int_t imult=0;imult<9;imult++) {
      if(mult >=MultBinsMin[imult] && mult <=MultBinsMax[imult]) {
	hEvMultSel->Fill(imult);
	hNchTSLeft[imult]->Fill(nchTSLeft);
	hNchTSRight[imult]->Fill(nchTSRight);
	hSumptTSLeft[imult]->Fill(sumptTSLeft);    
	hSumptTSRight[imult]->Fill(sumptTSRight);
	hNchTSLeftvsNchTSRight[imult]->Fill(nchTSLeft,nchTSRight);
	hSumptTSLeftvsSumptTSRight[imult]->Fill(sumptTSLeft,sumptTSRight);
      }
    }
  }
}
//______________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::AnalyzeESDforDCA(AliESDEvent* fESD)
{
  hINEL0->Fill(2);

  fRun  = fESD->GetRunNumber();
  fEventId = 0;
  if(fESD->GetHeader()) fEventId = GetEventIdAsLong(fESD->GetHeader());
  const Int_t nESDTracks = fESD->GetNumberOfTracks();

  if (fAnalysisMC){
    for ( int iT = 0 ; iT < nESDTracks ; iT++ ){

      AliESDtrack* esdTrack = fESD->GetTrack(iT);
      if(!esdTrack) continue;

      //track cuts
      UInt_t selectDebug = 0;
      if (fTrackFilterDCA){
	selectDebug = fTrackFilterDCA->IsSelected(esdTrack);
	if (!selectDebug) continue;
      }

      Double_t eta = esdTrack->Eta();
      Double_t pt  = esdTrack->Pt();

      if ( TMath::Abs(eta) > fEtaCut) continue;
      if ( pt < 0.15) continue;

      Int_t mcLabel = TMath::Abs(esdTrack->GetLabel());

      fMCStack = fMCEvent->Stack();
      if(!fMCStack){
	cout<<"------------No Ali Stack------------------"<<endl;
	return;
      }

      TParticle *mcParticle = fMCStack->Particle(mcLabel);
      if(!mcParticle) {printf("----------------ERROR: mcParticle not available------------------\n"); continue;}

      eta = mcParticle->Eta();
      pt = mcParticle->Pt();

      if ( TMath::Abs(eta) > fEtaCut) continue;
      if ( pt < 0.15) continue;

      if (!(TMath::Abs(fMCStack->Particle(mcLabel)->GetPDG()->Charge()) == 3) ) continue;
      esdTrack->GetImpactParameters(fdcaxy,fdcaz);

      if (!fMCStack->IsPhysicalPrimary(mcLabel)){
	Int_t index = fMCStack->Particles()->IndexOf(mcParticle);
	if ( fMCStack->IsSecondaryFromWeakDecay(index)){
	  ptvstrackletsvsdcaDecs->Fill(pt,ftrackmult08,fdcaxy);
	  ptvstrackletsvsdcacentralDecs->Fill(pt,ftrackmult08,fdcaxy);
	}

	if ( fMCStack->IsSecondaryFromMaterial(index) )
	  {
	    ptvstrackletsvsdcaMatl->Fill(pt,ftrackmult08,fdcaxy);
	    ptvstrackletsvsdcacentralMatl->Fill(pt,ftrackmult08,fdcaxy);
	  }
	continue;
      }

      ptvstrackletsvsdcaPrim->Fill(pt,ftrackmult08,fdcaxy);
      ptvstrackletsvsdcacentralPrim->Fill(pt,ftrackmult08,fdcaxy);
    }

  }
  else{
    for(Int_t iT = 0; iT < nESDTracks; iT++){

      AliESDtrack* esdTrack = fESD->GetTrack(iT);

      //only golden track cuts
      UInt_t selectDebug = 0;
      if (fTrackFilterDCA) {
	selectDebug = fTrackFilterDCA->IsSelected(esdTrack);
	if (!selectDebug) continue;
      }

      Double_t eta  = esdTrack->Eta();
      Double_t pt   = esdTrack->Pt();

      if( TMath::Abs(eta) > fEtaCut )continue;
      if( pt < 0.15 )continue;

      esdTrack->GetImpactParameters(fdcaxy,fdcaz);
      ptvstrackletsvsdcaData->Fill(pt,ftrackmult08,fdcaxy);
      ptvstrackletsvsdcacentralData->Fill(pt,ftrackmult08,fdcaxy);
    }//end of track loop
  }

  hINEL0->Fill(3);
}
//______________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraDphi::selectVertex2015pp(AliESDEvent *esd,
							    Bool_t checkSPDres, //enable check on vtx resolution
							    Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
							    Bool_t checkProximity) //apply cut on relative position of spd and trk verteces
{

  if (!esd) return kFALSE;

  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();

  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;

  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD) return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
      if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
    }
  }
  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraDphi::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskUeSpectraDphi::GetEventIdAsLong(AliVHeader* header)
{
  // To have a unique id for each event in a run!
  // Modified from AliRawReader.h
  return ((ULong64_t)header->GetBunchCrossNumber()+
	  (ULong64_t)header->GetOrbitNumber()*3564+
	  (ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

//____________________________________________________________
void AliAnalysisTaskUeSpectraDphi::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;

  // standar parameters ------------------- //
  double maxdcaz = 2.;
  double minratiocrossrowstpcover = 0.8;
  double maxfraclusterstpcshared = 0.4;
  double maxchi2perclustertpc = 4.0;
  double maxchi2perclusterits = 36.;
  double geowidth = 3.;
  double geolenght = 130.;
  //double mincrossedrows = 120.0;
  //double chi2tpcconstrainedglobal = 36.;

  // TPC
  esdTrackCuts->SetCutGeoNcrNcl(geowidth,geolenght,1.5,0.85,0.7);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minratiocrossrowstpcover);
  esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perclustertpc);
  esdTrackCuts->SetMaxFractionSharedTPCClusters(maxfraclusterstpcshared);

  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMaxChi2PerClusterITS(maxchi2perclusterits);

  // primary selection
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMaxDCAToVertexZ(maxdcaz);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // (7*(------))

  fTrackFilter->AddCuts(esdTrackCuts);
}

//____________________________________________________________
void AliAnalysisTaskUeSpectraDphi::SetTrackCutsDCA(AliAnalysisFilter* fTrackFilter){
  // track cuts for Feed-Down correction ---------------------------------------------------------------- //
  AliESDtrackCuts* esdTrackCutsDCA = new AliESDtrackCuts;

  Double_t maxdcaz = 2.;
  Double_t minratiocrossrowstpcover = 0.8;
  Double_t maxfraclusterstpcshared = 0.4;
  Double_t maxchi2perclustertpc = 4.0;
  Double_t maxchi2perclusterits = 36.;
  Double_t geowidth = 3.;
  Double_t geolenght = 130.;

  // TPC

  esdTrackCutsDCA->SetRequireTPCRefit(kTRUE);
  esdTrackCutsDCA->SetMinRatioCrossedRowsOverFindableClustersTPC(minratiocrossrowstpcover);
  esdTrackCutsDCA->SetMaxChi2PerClusterTPC(maxchi2perclustertpc);
  esdTrackCutsDCA->SetMaxFractionSharedTPCClusters(maxfraclusterstpcshared);

  // ITS
  esdTrackCutsDCA->SetRequireITSRefit(kTRUE);
  esdTrackCutsDCA->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCutsDCA->SetMaxChi2PerClusterITS(maxchi2perclusterits);

  // primary selection
  esdTrackCutsDCA->SetDCAToVertex2D(kFALSE);
  esdTrackCutsDCA->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsDCA->SetMaxDCAToVertexZ(maxdcaz);
  esdTrackCutsDCA->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsDCA->SetCutGeoNcrNcl(geowidth,geolenght,1.5,0.85,0.7); // Affects more or less 10% !!
  fTrackFilterDCA->AddCuts(esdTrackCutsDCA);

  fTrackFilterDCA->AddCuts(esdTrackCutsDCA);

}

//____________________________________________________________
Bool_t AliAnalysisTaskUeSpectraDphi::isMCEventTrueINEL0(AliMCEvent* fMCEvent)
{
  Bool_t isINEL0 = kFALSE;
  for ( int iT = 0 ; iT < fMCStack->GetNtrack(); iT++ ) // loop over TRUE MC
    {
      TParticle *mcParticle = fMCStack->Particle(iT);

      if (!mcParticle)
	{
	  cout<<"no mcParticle"<<endl;
	  continue;
	}

      if(!fMCStack->IsPhysicalPrimary(iT))
	continue;

      if(!(mcParticle->Pt()>0.0))
	continue;

      Double_t eta = mcParticle->Eta();
      if ( TMath::Abs(eta) > 1.0 )
	continue;

      if ( !( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) )
	continue;

      isINEL0 = kTRUE;
      break;

    }

  return isINEL0;

}

Double_t AliAnalysisTaskUeSpectraDphi::DeltaPhiOA(Double_t phia, Double_t phib, Double_t range)
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();

  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi > range)        dphi -= 2*pi;
  else if (dphi < -range)  dphi += 2*pi;

  return dphi;
}

Double_t AliAnalysisTaskUeSpectraDphi::DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin, Double_t rangeMax)
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();

  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;

  return dphi;
}



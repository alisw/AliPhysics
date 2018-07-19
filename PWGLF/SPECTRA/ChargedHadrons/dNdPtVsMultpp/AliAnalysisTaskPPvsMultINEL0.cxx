/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
 * provided "as is" without express or implied warranty.                  * **************************************************************************/

//Authors: Sergio Iga ,sigabuit@cern.ch


#include "AliAnalysisTaskPPvsMultINEL0.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include "AliCentrality.h" 
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 
#include <AliAODHeader.h>


#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"

//#include "AliTransverseEventShape.h"

// STL includes
#include <iostream>
using namespace std;

ClassImp(AliAnalysisTaskPPvsMultINEL0)
	//_____________________________________________________________________________
	//AliAnalysisTaskPPvsMultINEL0::AliAnalysisTaskQAHighPtDeDx(const char *name):
	AliAnalysisTaskPPvsMultINEL0::AliAnalysisTaskPPvsMultINEL0():
	AliAnalysisTaskSE(),
	fMCEvent(0x0),
	fMCStack(0x0),
	fESD(0x0),
	fAOD(0x0),
	fAnalysisMC(kFALSE),
	fTrackFilterDCA(0x0),
	fAnalysisType("ESD"),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
	fPPVsMultUtils(0x0),
	fMultSelection(0x0),
	fVtxCut(10.0),  
	fEtaCut(0.9),  
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtx_SPD(-9999),
	fRun(-999),
	fEventId(-999),
	fdcaxy(-999),
	fdcaz(-999),
	fisPS(kFALSE),
	fisTracklet(kFALSE),
	fisMCvtxInZcut(kFALSE),
	fListOfObjects(0),
	fEvents(0x0),
	fVtxBeforeCuts(0x0), 
	fVtxAfterCuts(0x0),
	fn1(0x0),
	fNchargedTrue(0x0),      
	ftrackmult08(-999),   
	fv0mpercentile(-999),
	isINEL0Rec(kFALSE),
	isINEL0True(kFALSE),
	refMultvsV0mINEL0(0x0),
	refMultvsV0mTrigINEL0(0x0),
	refMultvsV0mTrueINEL0(0x0),
	ptvstrackletsvsdcaData(0x0),
	ptvstrackletsvsdcaPrim(0x0),
	ptvstrackletsvsdcaDecs(0x0),
	ptvstrackletsvsdcaMatl(0x0),
	ptvsv0mvsdcaData(0x0),
	ptvsv0mvsdcaPrim(0x0),
	ptvsv0mvsdcaDecs(0x0),
	ptvsv0mvsdcaMatl(0x0),
	ptvstrackletsvsdcacentralData(0x0),
	ptvstrackletsvsdcacentralPrim(0x0),
	ptvstrackletsvsdcacentralDecs(0x0),
	ptvstrackletsvsdcacentralMatl(0x0),
	ptvsv0mvsdcacentralData(0x0),
	ptvsv0mvsdcacentralPrim(0x0),
	ptvsv0mvsdcacentralDecs(0x0),
	ptvsv0mvsdcacentralMatl(0x0),
	effcomputationGen(0x0),
	sigLossTrueINEL0(0x0),
	sigLossTrigINEL0(0x0),
	nchtruevsrefmult08(0x0),
	effcomputationGenPi(0x0),
	effcomputationGenK(0x0),
	effcomputationGenP(0x0),
	effcomputationGenSm(0x0),
	effcomputationGenSp(0x0),
	effcomputationGenO(0x0),
	effcomputationGenXi(0x0),
	effcomputationGenL(0x0),
	effcomputationGenRest(0x0),  
	effcomputationRecPi(0x0),
	effcomputationRecK(0x0),
	effcomputationRecP(0x0),
	effcomputationRecSm(0x0),
	effcomputationRecSp(0x0),
	effcomputationRecO(0x0),
	effcomputationRecXi(0x0),
	effcomputationRecL(0x0),
	effcomputationRecRest(0x0),
	fPS_MC(0x0),
	fVtxPS_MC(0x0),
	fPS(0x0),
	fVtxPS(0x0),
	refMultvsZvtx(0x0),
	v0mPercentileQA(0x0)

{
	for ( int i = 0 ; i < 18 ; i++ )
	{
	  secondaries[i] = 0;
	  primariesTrackFilter[i] = 0;
	  ptvsv0m[i] = 0; 
	  ptvstracklets[i] = 0;
	  effcomputationRec[i] = 0;
  	  fTrackFilter[i] = 0;

	}

}


AliAnalysisTaskPPvsMultINEL0::AliAnalysisTaskPPvsMultINEL0(const char *name):
	AliAnalysisTaskSE(name),
	fMCEvent(0x0),
	fMCStack(0x0),
	fESD(0x0),
	fAOD(0x0),
	fAnalysisMC(kFALSE),
	fTrackFilterDCA(0x0),
	fAnalysisType("ESD"),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
	fPPVsMultUtils(0x0),
	fMultSelection(0x0),
	fVtxCut(10.0),  
	fEtaCut(0.9),  
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtx_SPD(-9999),
	fRun(-999),
	fEventId(-999),
	fdcaxy(-999),
	fdcaz(-999),
	fisPS(kFALSE),
	fisTracklet(kFALSE),
	fisMCvtxInZcut(kFALSE),
	fListOfObjects(0),
	fEvents(0x0), 
	fVtxBeforeCuts(0x0), 
	fVtxAfterCuts(0x0),
	fn1(0x0),
	fNchargedTrue(0x0),      
	ftrackmult08(-999),   
	fv0mpercentile(-999),
	isINEL0Rec(kFALSE),
	isINEL0True(kFALSE),
	refMultvsV0mINEL0(0x0),
	refMultvsV0mTrigINEL0(0x0),
	refMultvsV0mTrueINEL0(0x0),
	ptvstrackletsvsdcaData(0x0),
	ptvstrackletsvsdcaPrim(0x0),
	ptvstrackletsvsdcaDecs(0x0),
	ptvstrackletsvsdcaMatl(0x0),
	ptvsv0mvsdcaData(0x0),
	ptvsv0mvsdcaPrim(0x0),
	ptvsv0mvsdcaDecs(0x0),
	ptvsv0mvsdcaMatl(0x0),
	ptvstrackletsvsdcacentralData(0x0),
	ptvstrackletsvsdcacentralPrim(0x0),
	ptvstrackletsvsdcacentralDecs(0x0),
	ptvstrackletsvsdcacentralMatl(0x0),
	ptvsv0mvsdcacentralData(0x0),
	ptvsv0mvsdcacentralPrim(0x0),
	ptvsv0mvsdcacentralDecs(0x0),
	ptvsv0mvsdcacentralMatl(0x0),
	effcomputationGen(0x0),
	sigLossTrueINEL0(0x0),
	sigLossTrigINEL0(0x0),
	nchtruevsrefmult08(0x0),
	effcomputationGenPi(0x0),
	effcomputationGenK(0x0),
	effcomputationGenP(0x0),
	effcomputationGenSm(0x0),
	effcomputationGenSp(0x0),
	effcomputationGenO(0x0),
	effcomputationGenXi(0x0),
	effcomputationGenL(0x0),
	effcomputationGenRest(0x0),  
	effcomputationRecPi(0x0),
	effcomputationRecK(0x0),
	effcomputationRecP(0x0),
	effcomputationRecSm(0x0),
	effcomputationRecSp(0x0),
	effcomputationRecO(0x0),
	effcomputationRecXi(0x0),
	effcomputationRecL(0x0),
	effcomputationRecRest(0x0),
	fPS_MC(0x0),
	fVtxPS_MC(0x0),
	fPS(0x0),
	fVtxPS(0x0),
	refMultvsZvtx(0x0),
	v0mPercentileQA(0x0)
	
{
	for ( int i = 0 ; i < 18 ; i++ )
	{
	  secondaries[i] = 0;
	  primariesTrackFilter[i] = 0;
	  ptvsv0m[i] = 0; 
	  ptvstracklets[i] = 0;
	  fTrackFilter[i] = 0;
	  effcomputationRec[i] = 0;
	}	


	DefineOutput(1, TList::Class());//esto es nuevo

}




AliAnalysisTaskPPvsMultINEL0::~AliAnalysisTaskPPvsMultINEL0()
{
  //
  // Destructor
  //

  if (fListOfObjects) {
    delete fListOfObjects;
    fListOfObjects = 0x0;
  }
  
  if (fPPVsMultUtils)
  {
    delete fPPVsMultUtils;
    fPPVsMultUtils = 0x0;
  }
  
  if (fMultSelection)
  {
    delete fMultSelection;
    fMultSelection = 0x0;
  }

}

//______________________________________________________________________________
void AliAnalysisTaskPPvsMultINEL0::UserCreateOutputObjects()
{ 
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 
  // We also create the random generator here so it might get different seeds...

  fRandom = new TRandom(0); // 0 means random seed


  //OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  //
  // Histograms
  //  
  fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 6, 0, 6);
  fListOfObjects->Add(fEvents);

  fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 12, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 12, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);
  

  fn1=new TH1F("fn1","fn1",10,0,10);
  fListOfObjects->Add(fn1);
  

  const Int_t nPtBins = 69;
  Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
	  0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
	  1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
	  2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
	  4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
	  11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
	  26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 60.0,200 };
	  
	const int nBinsDCAxy = 121;
	double binsDCAxy[] = {-3.025,-2.975,-2.925,-2.875,-2.825,-2.775,-2.725,-2.675,-2.625,-2.575,-2.525,-2.475,-2.425,-2.375,-2.325,-2.275,-2.225,
	-2.175,-2.125,-2.075,-2.025,-1.975,-1.925,-1.875,-1.825,-1.775,-1.725,-1.675,-1.625,-1.575,-1.525,-1.475,
	-1.425,-1.375,-1.325,-1.275,-1.225,-1.175,-1.125,-1.075,-1.025,-0.975,-0.925,-0.875,-0.825,-0.775,-0.725,
	-0.675,-0.625,-0.575,-0.525,-0.475,-0.425,-0.375,-0.325,-0.275,-0.225,-0.175,-0.125,-0.075,-0.025,0.025,
	0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,
	0.825,0.875,0.925,0.975,1.025,1.075,1.125,1.175,1.225,1.275,1.325,1.375,1.425,1.475,1.525,
	1.575,1.625,1.675,1.725,1.775,1.825,1.875,1.925,1.975,2.025,2.075,2.125,2.175,2.225,2.275,
	2.325,2.375,2.425,2.475,2.525,2.575,2.625,2.675,2.725,2.775,2.825,2.875,2.925,2.975,3.025};
	
	const int nBinsDCAxyCentral = 461;
	
double binsDCAxyCentral[]= {-0.2305,-0.2295,-0.2285,-0.2275,-0.2265,-0.2255,-0.2245,-0.2235,-0.2225,-0.2215,-0.2205,-0.2195,-0.2185,-0.2175,-0.2165,-0.2155,-0.2145,
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
	  
	  
  const Int_t multBins = 200;
  const Double_t xmultBins[201] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
				  20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,
				  40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,
				  60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
				  80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,
				  100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,
				  120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,
				  140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
				  160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,
				  180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,
				  200};

				    
  Double_t V0Mbins[11] ={0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

  
  AliESDtrackCuts* esdTrackCutsRun2[18] = {0};
	
  for ( int iTc = 0 ; iTc < 18 ; iTc++ )
  {
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
	// ------------------------------------- //
	
	// variations of the track cuts -------- //
    if ( iTc == 1) maxdcaz = 1.0;
    if ( iTc == 2) maxdcaz = 5.0;
    if ( iTc == 5) minratiocrossrowstpcover = 0.7;
    if ( iTc == 6) minratiocrossrowstpcover = 0.9;
    if ( iTc == 7) maxfraclusterstpcshared = 0.2;
    if ( iTc == 8) maxfraclusterstpcshared = 1.0;
    if ( iTc == 9) maxchi2perclustertpc = 3.0;
    if ( iTc == 10) maxchi2perclustertpc = 5.0;
    if ( iTc == 11) maxchi2perclusterits = 25.0;
    if ( iTc == 12) maxchi2perclusterits = 49.0;
    /*
    if ( iTc == 14) mincrossedrows = 100.0;
    if ( iTc == 15) mincrossedrows = 130.0;
    if ( iTc == 16) chi2tpcconstrainedglobal = 25.0;
    if ( iTc == 17) chi2tpcconstrainedglobal = 49.0;  
    */
    
    if ( iTc == 14) geowidth = 2.0;
    if ( iTc == 15) geowidth = 4.0;
    if ( iTc == 16) geolenght = 120.0;
    if ( iTc == 17) geolenght = 140.0;
	
	// variations of the track cuts -------- //
    
    
    fTrackFilter[iTc] = new AliAnalysisFilter(Form("fTrackFilter%d",iTc));
    esdTrackCutsRun2[iTc] = new AliESDtrackCuts(Form("esdTrackCutsRun2%d",iTc));
    
    // TPC
    
    //esdTrackCutsRun2[iTc]->SetMinNCrossedRowsTPC(mincrossedrows); // not anymore, GeCuts do the job
	esdTrackCutsRun2[iTc]->SetCutGeoNcrNcl(geowidth,geolenght,1.5,0.85,0.7);
    esdTrackCutsRun2[iTc]->SetRequireTPCRefit(kTRUE);
    esdTrackCutsRun2[iTc]->SetMinRatioCrossedRowsOverFindableClustersTPC(minratiocrossrowstpcover); 
    esdTrackCutsRun2[iTc]->SetMaxChi2PerClusterTPC(maxchi2perclustertpc);
    esdTrackCutsRun2[iTc]->SetMaxFractionSharedTPCClusters(maxfraclusterstpcshared); 
    
    // ITS
    esdTrackCutsRun2[iTc]->SetRequireITSRefit(kTRUE);
    if ( iTc != 13 )
      esdTrackCutsRun2[iTc]->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    
    esdTrackCutsRun2[iTc]->SetMaxChi2PerClusterITS(maxchi2perclusterits);

    // primary selection
    esdTrackCutsRun2[iTc]->SetDCAToVertex2D(kFALSE);
    esdTrackCutsRun2[iTc]->SetRequireSigmaToVertex(kFALSE);
    esdTrackCutsRun2[iTc]->SetMaxDCAToVertexZ(maxdcaz); 
    esdTrackCutsRun2[iTc]->SetAcceptKinkDaughters(kFALSE);
    
    if ( iTc == 3 )
      esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("4*(0.0026+0.0050/pt^1.01)");
    else if ( iTc == 4 )
      esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("10*(0.0026+0.0050/pt^1.01)");
    else
      esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // (7*(------))
    
    //esdTrackCutsRun2[iTc]->SetMaxChi2TPCConstrainedGlobal(chi2tpcconstrainedglobal); // not used 
    
    fTrackFilter[iTc]->AddCuts(esdTrackCutsRun2[iTc]);
  }
  
  // ---------------------------------------------------------------------------------------------------- //

  // track cuts for Feed-Down correction ---------------------------------------------------------------- //
  
  AliESDtrackCuts* esdTrackCutsDCA = 0;
  
  double maxdcaz = 2.; 
  double minratiocrossrowstpcover = 0.8;
  double maxfraclusterstpcshared = 0.4; 
  double maxchi2perclustertpc = 4.0; 
  double maxchi2perclusterits = 36.; 
  double geowidth = 3.;
  double geolenght = 130.;
  
  fTrackFilterDCA = new AliAnalysisFilter("fTrackFilterDCA");
  esdTrackCutsDCA = new AliESDtrackCuts("esdTrackCutsDCA");
  
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
  
  // ----------------------------------------------------------------------------------------------------- //
    
  if ( fAnalysisMC )
  {
	// generated pT
    effcomputationGen = new TH3F("effcomputationGen","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGen);
	
	// histos for signal loss
    
    sigLossTrueINEL0 = new TH3F("sigLossTrueINEL0","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(sigLossTrueINEL0);
    
    sigLossTrigINEL0 = new TH3F("sigLossTrigINEL0","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(sigLossTrigINEL0);
	
	// histos for event loss (trigger eff)
    
    refMultvsV0mTrigINEL0 = new TH2F("refMultvsV0mTrigINEL0","RefMult N_{ch} vs v0m percentile;#N_{ch}^{RefMult};V0M^{Rec} percentile",
			      multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(refMultvsV0mTrigINEL0);
    
    refMultvsV0mTrueINEL0 = new TH2F("refMultvsV0mTrueINEL0","RefMult N_{ch} vs v0m percentile;#N_{ch}^{RefMult};V0M^{Rec} percentile",
			      multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(refMultvsV0mTrueINEL0);
    
	// nch vs tracklets (RefMult08)
    
    nchtruevsrefmult08 = new TH2F("nchtruevsrefmult08","N_{ch}^{True} vs RefMult N_{ch};N_{ch}^{True};N_{ch}^{Rec}",
			      multBins,xmultBins,multBins,xmultBins);
    fListOfObjects->Add(nchtruevsrefmult08);
	
	// histos for Vtx correction MC ( not used, only to see the difference with data)
    
    fPS_MC = new TH2F("fPS_MC","RefMult N_{ch} vs v0m percentile;#N_{ch}^{RefMult};V0M^{Rec} percentile",
			      multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(fPS_MC);
    
    fVtxPS_MC = new TH2F("fVtxPS_MC","RefMult N_{ch} vs v0m percentile;#N_{ch}^{RefMult};V0M^{Rec} percentile",
			      multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(fVtxPS_MC);
    
    // HISTOS FOR FEED_DOWN CORRECTION ------------------------------------------------------------------------
    
    //RefMult
    ptvstrackletsvsdcaPrim = new TH3F("ptvstrackletsvsdcaPrim","pt vs tracklets vs dca Primaries",
				      nPtBins,xBins,200,xmultBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvstrackletsvsdcaPrim);
    
    ptvstrackletsvsdcaDecs = new TH3F("ptvstrackletsvsdcaDecs","pt vs tracklets vs dca Decays",
				      nPtBins,xBins,200,xmultBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvstrackletsvsdcaDecs);
    
    ptvstrackletsvsdcaMatl = new TH3F("ptvstrackletsvsdcaMatl","pt vs tracklets vs dca Material",
				      nPtBins,xBins,200,xmultBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvstrackletsvsdcaMatl);
    // V0M
    ptvsv0mvsdcaPrim = new TH3F("ptvsv0mvsdcaPrim","pt vs v0m vs dca Primaries",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsv0mvsdcaPrim);
    
    ptvsv0mvsdcaDecs = new TH3F("ptvsv0mvsdcaDecs","pt vs v0m vs dca Decays",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsv0mvsdcaDecs);
    
    ptvsv0mvsdcaMatl = new TH3F("ptvsv0mvsdcaMatl","pt vs v0m vs dca Material",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsv0mvsdcaMatl);
	
	
	
	
	//RefMult
    ptvstrackletsvsdcacentralPrim = new TH3F("ptvstrackletsvsdcacentralPrim","pt vs tracklets vs dca central Primaries",
				      nPtBins,xBins,200,xmultBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvstrackletsvsdcacentralPrim);
    
    ptvstrackletsvsdcacentralDecs = new TH3F("ptvstrackletsvsdcacentralDecs","pt vs tracklets vs dca central Decays",
				      nPtBins,xBins,200,xmultBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvstrackletsvsdcacentralDecs);
    
    ptvstrackletsvsdcacentralMatl = new TH3F("ptvstrackletsvsdcacentralMatl","pt vs tracklets vs dca central Material",
				      nPtBins,xBins,200,xmultBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvstrackletsvsdcacentralMatl);
    // V0M
    ptvsv0mvsdcacentralPrim = new TH3F("ptvsv0mvsdcacentralPrim","pt vs v0m vs dca central Primaries",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsv0mvsdcacentralPrim);
    
    ptvsv0mvsdcacentralDecs = new TH3F("ptvsv0mvsdcacentralDecs","pt vs v0m vs dca central Decays",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsv0mvsdcacentralDecs);
    
    ptvsv0mvsdcacentralMatl = new TH3F("ptvsv0mvsdcacentralMatl","pt vs v0m vs dca central Material",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsv0mvsdcacentralMatl);
  
    // ------------------------------------------------------------------------------------------------------------------
    
    
	// histos for particle composition
    
    //pions
    effcomputationGenPi = new TH3F("effcomputationGenPi","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenPi);
    //kaons
    effcomputationGenK = new TH3F("effcomputationGenK","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenK);
    //protons
    effcomputationGenP = new TH3F("effcomputationGenP","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenP);
    //Sigma-
    effcomputationGenSm = new TH3F("effcomputationGenSm","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenSm);
    //sigma+
    effcomputationGenSp = new TH3F("effcomputationGenSp","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenSp);
    //omega-
    effcomputationGenO = new TH3F("effcomputationGenO","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenO);
    //xi-
    effcomputationGenXi = new TH3F("effcomputationGenXi","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenXi);
    //Lambda
    effcomputationGenL = new TH3F("effcomputationGenL","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenL);
    //Others
    effcomputationGenRest = new TH3F("effcomputationGenRest","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationGenRest);
    
    
    //pions
    effcomputationRecPi = new TH3F("effcomputationRecPi","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecPi);
    //kaons
    effcomputationRecK = new TH3F("effcomputationRecK","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecK);
    //protons
    effcomputationRecP = new TH3F("effcomputationRecP","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecP);
    //Sigma-
    effcomputationRecSm = new TH3F("effcomputationRecSm","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecSm);
    //sigma+
    effcomputationRecSp = new TH3F("effcomputationRecSp","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecSp);
    //omega-
    effcomputationRecO = new TH3F("effcomputationRecO","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecO);
    //xi-
    effcomputationRecXi = new TH3F("effcomputationRecXi","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecXi);
    //Lambda
    effcomputationRecL = new TH3F("effcomputationRecL","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecL);
    //Others
    effcomputationRecRest = new TH3F("effcomputationRecRest","p_{T}^{gen} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Gen};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(effcomputationRecRest);
  
    for ( int i = 0 ; i < 18 ; i++ ) // rec pT, secondaries, and "primaries" (all) that passes track fileters, i == 0  standard
    {	  

      secondaries[i] = new TH3F(Form("secondaries%d",i),"p_{T}^{rec} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
      fListOfObjects->Add(secondaries[i]);
      
      effcomputationRec[i] = new TH3F(Form("effcomputationRec%d",i),"p_{T}^{rec} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
      fListOfObjects->Add(effcomputationRec[i]);

      primariesTrackFilter[i] = new TH3F(Form("primariesTrackFilter%d",i),"p_{T}^{rec} vs True N_{ch} vs v0m percentile;#it{p}_{T} (GeV/c);N_{ch}^{Rec};V0M^{Rec} percentile",
			      nPtBins,xBins,multBins,xmultBins,10,V0Mbins);
      fListOfObjects->Add(primariesTrackFilter[i]);
    }
    
    
  }
  else // for the data
  {
    // histos for FEED_DOWN CORRECTION ----------------------------------------------------------------------------
    
    //RefMult
    ptvstrackletsvsdcaData = new TH3F("ptvstrackletsvsdcaData","pt vs tracklets vs dca Data",
				      nPtBins,xBins,multBins,xmultBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvstrackletsvsdcaData);
    //V0M
    ptvsv0mvsdcaData = new TH3F("ptvsv0mvsdcaData","pt vs v0m vs dca Data",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsv0mvsdcaData);
	
	
	
	
	//RefMult
    ptvstrackletsvsdcacentralData = new TH3F("ptvstrackletsvsdcacentralData","pt vs tracklets vs dca central Data",
				      nPtBins,xBins,multBins,xmultBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvstrackletsvsdcacentralData);
    //V0M
    ptvsv0mvsdcacentralData = new TH3F("ptvsv0mvsdcacentralData","pt vs v0m vs dca central Data",
				      nPtBins,xBins,10,V0Mbins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsv0mvsdcacentralData);
    
    
    
    // --------------------------------------------------------------------------------------------------------------
    
    // event counter, tracklets vs V0M, projections (XY) are the event counter for each mult class
    
    refMultvsV0mINEL0 = new TH2F("refMultvsV0mINEL0","RefMult N_{ch} vs v0m percentile;#N_{ch}^{RefMult};V0M^{Rec} percentile",
			      multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(refMultvsV0mINEL0);
	
	// histos for vtx corrections, data driven
    
    fPS = new TH2F("fPS","RefMult N_{ch} vs v0m percentile;#N_{ch}^{RefMult};V0M^{Rec} percentile",
			      multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(fPS);
    
    fVtxPS = new TH2F("fVtxPS","RefMult N_{ch} vs v0m percentile;#N_{ch}^{RefMult};V0M^{Rec} percentile",
			      multBins,xmultBins,10,V0Mbins);
    fListOfObjects->Add(fVtxPS);
  
    for ( int i = 0 ; i < 18 ; i++ )
    {
      // Histos for data -----------------------------------------------------------------------------------------------------------------
    
      ptvsv0m[i] = new TH2F(Form("ptvsv0m%d",i),"p_{T} vs V0M percentile; #it{p}_{T} (GeV/c);V0M percentile",nPtBins,xBins,10,V0Mbins);
      fListOfObjects->Add(ptvsv0m[i]);
      
      ptvstracklets[i] = new TH2F(Form("ptvstracklets%d",i),"p_{T} vs Tracklets at |#eta| < 0.8; #it{p}_{T} (GeV/c);Tracklets_{|#eta|<0.8}",
			    nPtBins,xBins,multBins,xmultBins);
      fListOfObjects->Add(ptvstracklets[i]);

    // ----------------------------------------------------------------------------------------------------------------------------------
    }
  }
  
	// QA mult plots

	refMultvsZvtx = new TH2F("refMultvsZvtx"," ; Vtx_z ; Tracklets ; Entries",40,-10,10,multBins,xmultBins);
	fListOfObjects->Add(refMultvsZvtx);

	v0mPercentileQA = new TH1F("v0mPercentileQA","V0M percentile ; v0mPercentile ; Entries",100,0,100);
	fListOfObjects->Add(v0mPercentileQA);


  // Post output data.
  PostData(1,fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskPPvsMultINEL0::UserExec(Option_t *) 
{
	// Main loop
	//
	// First we make sure that we have valid input(s)!
	//

	AliVEvent *event = InputEvent();
	if (!event) 
	{
		Error("UserExec", "Could not retrieve event");
		return;
	}

	if (fAnalysisMC)
	{
		fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
		if (!fMCEvent)
		{ 
			cout<<"Could not retrieve MC event"<<endl;
			return;
		}

		AliHeader *header = fMCEvent->Header();
		if(!header) {AliDebug( AliLog::kError , "Header not avaible" ); return; }

		fMCStack = fMCEvent->Stack();
		if(!fMCStack){ cout<<"No Ali Stack"<<endl;return; }
			
	}


	if (fAnalysisType == "ESD")
	{
		fESD = dynamic_cast<AliESDEvent*>(event);
		if(!fESD)
		{
			Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}    
	}
	else
	{
		fAOD = dynamic_cast<AliAODEvent*>(event);
		if(!fAOD)
		{
			Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}    
	}

	/************ BEGGINING OF EVENT SELECTION ************************************************************************************/

	// Get trigger decision
	fTriggeredEventMB = 0; //init    
	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigBit )
	{
		fTriggeredEventMB = 1;  //event triggered as minimum bias
	}

	// check for spd vs clusters background 
	Bool_t SPDvsClustersBG = kFALSE;

	AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
	if (!AnalysisUtils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
		return;
	}
	else
		SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG
		
	// is pile up ? We want isNotPileUp
	Bool_t isNotPileUp = AliPPVsMultUtils::IsNotPileupSPDInMultBins( fESD ); // my analysis 
	//Bool_t isNotPileUp = !fESD->IsPileupFromSPD(5,0.8);                    // GSI analysis
	//Bool_t isNotPileUp = !fESD->IsPileupFromSPD();                         // Gyula's note 

	Bool_t IncompleteDAQ = fESD->IsIncompleteDAQ(); // we want is not incomplete DAQ

	// -------------------------------------- multiplcity estimators section ------------------------------------------ //

	ftrackmult08 = -999;
	fv0mpercentile = -999;

	//ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //reference
	ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);     //tracklets
	//ftrackmult08 = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD); //Combined estimator

	fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
	if (!fMultSelection)
		cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
	else
	fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");



	// ------------------------------------------ end of mult estimators -------------------------------------------------//

	// vertex bussines 
	//const AliESDVertex * spdVertex =    fESD->GetPrimaryVertexSPD();
	const AliESDVertex * vertex    =    fESD->GetPrimaryVertex(); // tracks vertex, if not -> spd vertex, if not TPC vertex

	Bool_t isVtxGood = vertex->GetStatus() &&
				selectVertex2015pp( fESD ,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm

	double vertex_z = vertex->GetZ();
	Bool_t isVtxInZCut = (TMath::Abs(vertex_z)   <= fVtxCut); // Zvtx in +- 10

	// physics selection		     
	fisTracklet = (AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1);

	/****************** IS PHYSICS SELECTION FLAG *************************************************/
	fisPS = fTriggeredEventMB && 
			!IncompleteDAQ && 
			!SPDvsClustersBG && 
			isNotPileUp &&
			fisTracklet; // *tracklet now included in the physics selection !!!
		

	// recontructed INEL > 0 is PS + vtx + Zvtx inside +-10 ------
	if ( fisPS && isVtxGood && isVtxInZCut)
		isINEL0Rec = kTRUE;
	else
		isINEL0Rec = kFALSE;

	if (fisPS)
		fVtxBeforeCuts->Fill(vertex_z);
	if (isINEL0Rec)
		fVtxAfterCuts->Fill(vertex_z);

	fn1->Fill(0); // all events

	if(fTriggeredEventMB)  
		fn1->Fill(1); // triggered events

	if ( fTriggeredEventMB && !IncompleteDAQ )  
		fn1->Fill(2); // trigger + IsIncompleteDAQ

	if (fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG)   
		fn1->Fill(3); // trigger + IsIncompleteDAQ + BG rejection

	if(fPileUpRej)
	{
		if(fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp) //     
			fn1->Fill(4); // trigger + IsIncompleteDAQ + BG rejection + PileUp
	}
	if (fisPS)  
		fn1->Fill(5); //PS: trigger + IsIncompleteDAQ + BG rejection + PileUp + 1 Tracklet in eta +-1

	if (fisPS && isVtxGood)
		fn1->Fill(6); //PS + GetPrimaryVertex

	if (isINEL0Rec)
		fn1->Fill(7); //PS + GetPrimaryVertex + isVtxInZCut


	if (fAnalysisMC) // analysis for generated MC
	{ 
		const AliVVertex *vertexMC = (AliVVertex*) fMCEvent->GetPrimaryVertex();
		fisMCvtxInZcut     = (TMath::Abs(vertexMC->GetZ()) <= fVtxCut);   // ZMCvtx in +- 10
		isINEL0True = isMCEventTrueINEL0(fMCEvent);

		// for trigger efficiency PS / INEL > 0 true
		if ( fisPS )
		{
			refMultvsV0mTrigINEL0->Fill(ftrackmult08,fv0mpercentile);
		}
		if (isINEL0True)
		{
			refMultvsV0mTrueINEL0->Fill(ftrackmult08,fv0mpercentile);
		}

		// MC (not used) next two histos for missing vtx correction
		if ( fisPS ) 
			fPS_MC->Fill(ftrackmult08,fv0mpercentile);

		if (fisPS && isVtxGood)
			fVtxPS_MC->Fill(ftrackmult08,fv0mpercentile);

		AnalyzeMC(fMCEvent);
	}
	else
	{
		if (isINEL0Rec)
			refMultvsV0mINEL0->Fill(ftrackmult08,fv0mpercentile);

		// Two histos for missing vtx correction
		if ( fisPS )
		{
			fPS->Fill(ftrackmult08,fv0mpercentile);
		}   
		if ( fisPS && isVtxGood )
		{
			fVtxPS->Fill(ftrackmult08,fv0mpercentile);
		}
	}


	if( (fAnalysisType == "ESD") && isINEL0Rec ) 
	{
		v0mPercentileQA->Fill(fv0mpercentile);
		refMultvsZvtx->Fill(vertex_z,ftrackmult08);
		AnalyzeESD(fESD);
		AnalyzeESDforDCA(fESD);

	}

	// Post output data.
	PostData(1, fListOfObjects);
}

//________________________________________________________________________

void AliAnalysisTaskPPvsMultINEL0::AnalyzeMC(AliMCEvent* fMCEvent)
{
	fEvents->Fill(0);

	fNchargedTrue = 0;

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

		Double_t eta = mcParticle->Eta();
		Double_t pt = mcParticle->Pt();
		int partPDG = TMath::Abs(mcParticle->GetPdgCode());
		if ( TMath::Abs(eta) > fEtaCut ) continue;

		if ( ( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) )
			fNchargedTrue++;

		if ( pt < 0.15 ) continue;    



		if (isINEL0Rec) 
		{
			if ( ( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) )
			{
				effcomputationGen->Fill(pt,ftrackmult08,fv0mpercentile); //inclusive
				if (partPDG == 211)  effcomputationGenPi->Fill(pt,ftrackmult08,fv0mpercentile); //pions
				else if (partPDG == 321)  effcomputationGenK->Fill(pt,ftrackmult08,fv0mpercentile);//kaons
				else if (partPDG == 2212) effcomputationGenP->Fill(pt,ftrackmult08,fv0mpercentile); //protons
				else if (partPDG == 3112) effcomputationGenSm->Fill(pt,ftrackmult08,fv0mpercentile);//sigma-
				else if (partPDG == 3222) effcomputationGenSp->Fill(pt,ftrackmult08,fv0mpercentile);//sigma+
				else if (partPDG == 3334) effcomputationGenO->Fill(pt,ftrackmult08,fv0mpercentile);//omega-
				else if (partPDG == 3312) effcomputationGenXi->Fill(pt,ftrackmult08,fv0mpercentile);//Xi-
				else effcomputationGenRest->Fill(pt,ftrackmult08,fv0mpercentile);//rest
			}
			else
			{
				if (partPDG == 3122) effcomputationGenL->Fill(pt,ftrackmult08,fv0mpercentile); //Lambda 
			}
		}

		if ( !( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) )
			continue;

		if ( isINEL0Rec ) // fisPS is wrong
		{
			sigLossTrigINEL0->Fill(pt,ftrackmult08,fv0mpercentile);
		}
		if ( isINEL0True && fisMCvtxInZcut)
		{
			sigLossTrueINEL0->Fill(pt,ftrackmult08,fv0mpercentile);
		}

	}

	if (isINEL0Rec)
		nchtruevsrefmult08->Fill(fNchargedTrue,ftrackmult08);

	fEvents->Fill(1);

}


//________________________________________________________________________

void AliAnalysisTaskPPvsMultINEL0::AnalyzeESD(AliESDEvent* esdEvent)
{

  
	fRun  = esdEvent->GetRunNumber();
	fEventId = 0;
	if(esdEvent->GetHeader())
		fEventId = GetEventIdAsLong(esdEvent->GetHeader());

	fEvents->Fill(2);

	const Int_t nESDTracks = esdEvent->GetNumberOfTracks();


	if (fAnalysisMC) // for recontructed montecarlo
	{
		for ( int iT = 0 ; iT < nESDTracks ; iT++ )
		{
			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);
			if(!esdTrack) continue;
			
			Double_t eta = esdTrack->Eta();
			Double_t pt  = esdTrack->Pt();

			if ( TMath::Abs(eta) > fEtaCut) continue;
			if ( pt < 0.15) continue;
			
			Int_t mcLabel = TMath::Abs(esdTrack->GetLabel());

			TParticle *mcParticle = fMCStack->Particle(mcLabel);
			if(!mcParticle) {printf("----------------ERROR: mcParticle not available------------------\n"); continue;}

			eta = mcParticle->Eta(); // generated eta and pT used intead of recontructed
			pt = mcParticle->Pt();
			
			if ( TMath::Abs(eta) > fEtaCut) continue;
			if ( pt < 0.15) continue;
			
			int partPDG = TMath::Abs(mcParticle->GetPdgCode());
			
			if ((TMath::Abs(mcParticle->GetPDG()->Charge()) == 3) ) // only for charged particles
			{      
				for ( int i = 0 ; i < 18 ; i++ ) // run over all the track cuts, i == 0 standard 
				{ 
					
					UInt_t selectDebug = 0;
					if (fTrackFilter[i])
					{
						selectDebug = fTrackFilter[i]->IsSelected(esdTrack);
						if (!selectDebug)
							continue;
					}
					
					primariesTrackFilter[i]->Fill(pt,ftrackmult08,fv0mpercentile); // "primaries" filtered with track filter
					
					if (!fMCStack->IsPhysicalPrimary(mcLabel)) // secondary particles
					{
						secondaries[i]->Fill(pt,ftrackmult08,fv0mpercentile);
					}
					else
					{
						effcomputationRec[i]->Fill(pt,ftrackmult08,fv0mpercentile); // rec pT for efficiency true primary particles
						if (i == 0) // for particle composition, only standar track cuts requires (i == 0)
						{
							if (partPDG == 211)  effcomputationRecPi->Fill(pt,ftrackmult08,fv0mpercentile); //pions
							else if (partPDG == 321)  effcomputationRecK->Fill(pt,ftrackmult08,fv0mpercentile);//kaons
							else if (partPDG == 2212) effcomputationRecP->Fill(pt,ftrackmult08,fv0mpercentile); //protons
							else if (partPDG == 3112) effcomputationRecSm->Fill(pt,ftrackmult08,fv0mpercentile);//sigma-
							else if (partPDG == 3222) effcomputationRecSp->Fill(pt,ftrackmult08,fv0mpercentile);//sigma+
							else if (partPDG == 3334) effcomputationRecO->Fill(pt,ftrackmult08,fv0mpercentile);//omega-
							else if (partPDG == 3312) effcomputationRecXi->Fill(pt,ftrackmult08,fv0mpercentile);//Xi-
							else effcomputationRecRest->Fill(pt,ftrackmult08,fv0mpercentile);//Rest
						}
					}	
				}
			}
			else // only for neutral particles (searching Lambdas)
			{
				UInt_t selectDebug = 0;
				if (fTrackFilter[0])
				{
					selectDebug = fTrackFilter[0]->IsSelected(esdTrack);
					if (!selectDebug)
						continue;
				}
				
				if (!fMCStack->IsPhysicalPrimary(mcLabel)) continue;
				
				if (partPDG == 3122) effcomputationRecL->Fill(pt,ftrackmult08,fv0mpercentile); //Lambda	

			}
		}
	}
	else // for real data
	{  
		for(Int_t iT = 0; iT < nESDTracks; iT++)
		{

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);
			
			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
			continue;
			if( pt < 0.15 )
			continue;
			
			for ( int i = 0 ; i < 18 ; i++ )
			{
				UInt_t selectDebug = 0;
				if (fTrackFilter[i])
				{
					selectDebug = fTrackFilter[i]->IsSelected(esdTrack);
					if (!selectDebug)
						continue;
				}	
				ptvstracklets[i]->Fill(pt,ftrackmult08);
				ptvsv0m[i]->Fill(pt,fv0mpercentile);
			}
		}//end of track loop
	}  
	fEvents->Fill(3);  
}


void AliAnalysisTaskPPvsMultINEL0::AnalyzeESDforDCA(AliESDEvent* esdEvent)
{

	fEvents->Fill(4);

	fRun  = esdEvent->GetRunNumber();
	fEventId = 0;
	if(esdEvent->GetHeader())
	fEventId = GetEventIdAsLong(esdEvent->GetHeader());


	const Int_t nESDTracks = esdEvent->GetNumberOfTracks();

	if (fAnalysisMC)
	{

		for ( int iT = 0 ; iT < nESDTracks ; iT++ )
		{

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);
			if(!esdTrack) continue;

			//track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterDCA)
			{
				selectDebug = fTrackFilterDCA->IsSelected(esdTrack);
				if (!selectDebug) 
					continue;
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
			
			
			if (!fMCStack->IsPhysicalPrimary(mcLabel)) 
			{
				Int_t index = fMCStack->Particles()->IndexOf(mcParticle);			  

				if ( fMCStack->IsSecondaryFromWeakDecay(index) )
				{
					ptvstrackletsvsdcaDecs->Fill(pt,ftrackmult08,fdcaxy);
					ptvsv0mvsdcaDecs->Fill(pt,fv0mpercentile,fdcaxy);
					
					ptvstrackletsvsdcacentralDecs->Fill(pt,ftrackmult08,fdcaxy);
					ptvsv0mvsdcacentralDecs->Fill(pt,fv0mpercentile,fdcaxy);
				}

				if ( fMCStack->IsSecondaryFromMaterial(index) )
				{
					ptvstrackletsvsdcaMatl->Fill(pt,ftrackmult08,fdcaxy);
					ptvsv0mvsdcaMatl->Fill(pt,fv0mpercentile,fdcaxy);
					
					ptvstrackletsvsdcacentralMatl->Fill(pt,ftrackmult08,fdcaxy);
					ptvsv0mvsdcacentralMatl->Fill(pt,fv0mpercentile,fdcaxy);
				}

				continue;			  
			}
			
			ptvstrackletsvsdcaPrim->Fill(pt,ftrackmult08,fdcaxy);
			ptvsv0mvsdcaPrim->Fill(pt,fv0mpercentile,fdcaxy);
			
			ptvstrackletsvsdcacentralPrim->Fill(pt,ftrackmult08,fdcaxy);
			ptvsv0mvsdcacentralPrim->Fill(pt,fv0mpercentile,fdcaxy);
			
		}

	}
	else
	{
		for(Int_t iT = 0; iT < nESDTracks; iT++)
		{

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			//only golden track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterDCA) {
				selectDebug = fTrackFilterDCA->IsSelected(esdTrack);
				if (!selectDebug)
					continue;
			}

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			if( pt < 0.15 )
				continue;

			esdTrack->GetImpactParameters(fdcaxy,fdcaz);

			ptvstrackletsvsdcaData->Fill(pt,ftrackmult08,fdcaxy);
			ptvsv0mvsdcaData->Fill(pt,fv0mpercentile,fdcaxy);
			
			ptvstrackletsvsdcacentralData->Fill(pt,ftrackmult08,fdcaxy);
			ptvsv0mvsdcacentralData->Fill(pt,fv0mpercentile,fdcaxy);

		}//end of track loop
	}

	fEvents->Fill(5);

}



//________________________________________________________________________
void AliAnalysisTaskPPvsMultINEL0::AnalyzeAOD(AliAODEvent* aodEvent)
{
  fRun  = aodEvent->GetRunNumber();
  fEventId = 0;
  if(aodEvent->GetHeader())
	  fEventId = GetEventIdAsLong(aodEvent->GetHeader());

  //Int_t     ftrackmult08 = 0; // no pt cuts
  //ftrackmult08 =  (dynamic_cast<AliAODHeader*>(aodEvent->GetHeader()))->GetRefMultiplicityComb08();

  //UInt_t    time      = 0; // Missing AOD info? aodEvent->GetTimeStamp();

  //Int_t     trackmult = 0; // no pt cuts
  //Int_t     nadded    = 0;

  Bool_t isNotPileUp = aodEvent->IsPileupFromSPD();
  if(fPileUpRej)
	  if(isNotPileUp)
		  return;


  if(fTriggeredEventMB) {// Only MC case can we have not triggered events

	  // accepted event

	  Int_t nAODTracks = aodEvent->GetNumberOfTracks();
	  //First loop to get ptleading
	  Double_t ptL = 0;

	  for(Int_t iT = 0; iT < nAODTracks; iT++) {

		  //AliAODTrack* aodTrack = AODevent->GetTrack(iT);
		  AliVTrack   *trackTmp = (AliVTrack *)aodEvent->GetTrack(iT);
		  AliAODTrack * aodTrack  = dynamic_cast<AliAODTrack *>(trackTmp);

		  if (fTrackFilter) {     
			  // "Global track RAA analysis QM2011 + Chi2ITS<36"; bit 1024
			  if(!aodTrack->TestFilterBit(1024))
				  continue;
		  }

		  Double_t eta  = aodTrack->Eta();
		  Double_t pt   = aodTrack->Pt();

		  if( TMath::Abs(eta) > fEtaCut )
			  continue;
		  //ptmin cut, only > 150 Mev/c
		  if( pt < 0.15 )
			  continue;
		  //extracting the ptleading
		  if( pt > ptL )
			  ptL = pt;

	  }//end of track loop

	  


  } // end if triggered

}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskPPvsMultINEL0::GetEventIdAsLong(AliVHeader* header)
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

Bool_t AliAnalysisTaskPPvsMultINEL0::selectVertex2015pp(AliESDEvent *esd,
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

Bool_t AliAnalysisTaskPPvsMultINEL0::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskPPvsMultINEL0::isMCEventTrueINEL0(AliMCEvent* fMCEvent)
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

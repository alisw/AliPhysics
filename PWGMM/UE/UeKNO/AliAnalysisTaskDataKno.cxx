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
 * provided "as is" without express or implied warranty.     
 *                                                                        *
 * Authors: Sushanta Tripathy (Sushanta.Tripathy@cern.ch)                 *
 *          Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)               *
 **************************************************************************/

/* AliAnaysisTaskDataKno source code
 * The analysis task produce all the histos needed for MC closure test studies
 * Results include only the KNO properties
 */

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;


#include <Riostream.h>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TList.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
//#include "AliV0vertexer.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include <AliHeader.h>
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultVariable.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
//#include "AliGenPythiaEventHeader.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include <AliESDVertex.h>
#include <AliMultiplicity.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TBits.h>
#include <AliAnalysisFilter.h>

using std::cout;
using std::endl;

#include "AliAnalysisTaskDataKno.h"

const Char_t * NameReg_1_data[3]={"NS","AS","TS"};

const Int_t nchNbins_data = 100;
Double_t nchbins_1_data[nchNbins_data+1]={-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5,48.5,49.5,50.5,51.5,52.5,53.5,54.5,55.5,56.5,57.5,58.5,59.5,60.5,61.5,62.5,63.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5,71.5,72.5,73.5,74.5,75.5,76.5,77.5,78.5,79.5,80.5,81.5,82.5,83.5,84.5,85.5,86.5,87.5,88.5,89.5,90.5,91.5,92.5,93.5,94.5,95.5,96.5,97.5,98.5,99.5};

Double_t nchbins_2_data[nchNbins_data+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
				21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
				39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,
				57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
				75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,
				93,94,95,96,97,98,99,100};

const Int_t nchNbins_data_3=10;
Double_t nchbins_3_data[nchNbins_data_3+1]={0,10,20,30,40,50,60,70,80,90,100};


const Int_t nDeltabins_data = 64;
Double_t Deltabins_data[nDeltabins_data+1]={-1.0472, -0.957204, -0.867211, -0.777217, -0.687223, -0.59723, -0.507236, -0.417243, -0.327249, -0.237256, -0.147262, -0.0572686, 0.0327249, 0.122718, 0.212712, 0.302706, 0.392699, 0.482693, 0.572686, 0.66268, 0.752673, 0.842667, 0.93266, 1.02265, 1.11265, 1.20264, 1.29263, 1.38263, 1.47262, 1.56262, 1.65261, 1.7426, 1.8326, 1.92259, 2.01258, 2.10258, 2.19257, 2.28256, 2.37256, 2.46255, 2.55254, 2.64254, 2.73253, 2.82252, 2.91252, 3.00251, 3.09251, 3.1825, 3.27249, 3.36249, 3.45248, 3.54247, 3.63247, 3.72246, 3.81245, 3.90245, 3.99244, 4.08243, 4.17243, 4.26242, 4.35241, 4.44241, 4.5324, 4.6224, 4.71239};

const Int_t ptNbins_data = 36;
Double_t ptbins1_1_data[ptNbins_data+1] = {
  0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0, 1.25, 1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0, 7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  30.0,  40.0,  50.0
};

const int nBinsDCAxy_data = 121;
  double binsDCAxy_data[nBinsDCAxy_data+1] = {-3.025,-2.975,-2.925,-2.875,-2.825,-2.775,-2.725,-2.675,-2.625,-2.575,-2.525,-2.475,-2.425,-2.375,-2.325,-2.275,-2.225, -2.175,-2.125,-2.075,-2.025,-1.975,-1.925,-1.875,-1.825,-1.775,-1.725,-1.675,-1.625,-1.575,-1.525,-1.475,-1.425,-1.375,-1.325,-1.275,-1.225,-1.175,-1.125,-1.075,-1.025,-0.975,-0.925,-0.875,-0.825,-0.775,-0.725,-0.675,-0.625,-0.575,-0.525,-0.475,-0.425,-0.375,-0.325,-0.275,-0.225,-0.175,-0.125,-0.075,-0.025,0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975,1.025,1.075,1.125,1.175,1.225,1.275,1.325,1.375,1.425,1.475,1.525,1.575,1.625,1.675,1.725,1.775,1.825,1.875,1.925,1.975,2.025,2.075,2.125,2.175,2.225,2.275,2.325,2.375,2.425,2.475,2.525,2.575,2.625,2.675,2.725,2.775,2.825,2.875,2.925,2.975,3.025};
	
  const int nBinsDCAxy_dataCentral = 461;
	
  double binsDCAxyCentral_data[nBinsDCAxy_dataCentral+1]= {-0.2305,-0.2295,-0.2285,-0.2275,-0.2265,-0.2255,-0.2245,-0.2235,-0.2225,-0.2215,-0.2205,-0.2195,-0.2185,-0.2175,-0.2165,-0.2155,-0.2145,-0.2135,-0.2125,-0.2115,-0.2105,-0.2095,-0.2085,-0.2075,-0.2065,-0.2055,-0.2045,-0.2035,-0.2025,-0.2015,-0.2005,-0.1995,-0.1985,-0.1975,-0.1965,-0.1955,-0.1945,-0.1935,-0.1925,-0.1915,-0.1905,-0.1895,-0.1885,-0.1875,-0.1865,-0.1855,-0.1845,-0.1835,-0.1825,-0.1815,-0.1805,-0.1795,-0.1785,-0.1775,-0.1765,-0.1755,-0.1745,-0.1735,-0.1725,-0.1715,-0.1705,-0.1695,-0.1685,-0.1675,-0.1665,-0.1655,-0.1645,-0.1635,-0.1625,-0.1615,-0.1605,-0.1595,-0.1585,-0.1575,-0.1565,-0.1555,-0.1545,-0.1535,-0.1525,-0.1515,-0.1505,-0.1495,-0.1485,-0.1475,-0.1465,-0.1455,-0.1445,-0.1435,-0.1425,-0.1415,-0.1405,-0.1395,-0.1385,-0.1375,-0.1365,-0.1355,-0.1345,-0.1335,-0.1325,-0.1315,-0.1305,-0.1295,-0.1285,-0.1275,-0.1265,-0.1255,-0.124,-0.1235,-0.1225,-0.1215,-0.1205,-0.1195,-0.1185,-0.1175,-0.1165,-0.1155,-0.1145,-0.1135,-0.1125,-0.1115,-0.1105,-0.1095,-0.1085,-0.1075,-0.1065,-0.1055,-0.1045,-0.1035,-0.1025,-0.1015,-0.1005,-0.0995,-0.0985,-0.0975,-0.0965,-0.0955,-0.0945,-0.0935,-0.0925,-0.0915,-0.0905,-0.0895,-0.0885,-0.0875,-0.0865,-0.0855,-0.0845,-0.0835,-0.0825,-0.0815,-0.0805,-0.0795,-0.0785,-0.0775,-0.0765,-0.0755,-0.0745,-0.0735,-0.0725,-0.0715,-0.0705,-0.0695,-0.0685,-0.0675,-0.0665,-0.0655,-0.0645,-0.0635,-0.0625,-0.0615,-0.0605,-0.0595,-0.0585,-0.0575,-0.0565,-0.0555,-0.0545,-0.0535,-0.0525,-0.0515,-0.0505,-0.0495,-0.0485,-0.0475,-0.0465,-0.0455,-0.0445,-0.0435,-0.0425,-0.0415,-0.0405,-0.0395,-0.0385,-0.0375,-0.0365,-0.0355,-0.0345,-0.0335,-0.0325,-0.0315,-0.0305,-0.0295,-0.0285,-0.0275,-0.0265,-0.0255,-0.0245,-0.0235,-0.0225,-0.0215,-0.0205,-0.0195,-0.0185,-0.0175,-0.0165,-0.0155,-0.0145,-0.0135,-0.0125,-0.0115,-0.0105,-0.0095,-0.0085,-0.0075,-0.0065,-0.0055,-0.0045,-0.0035,-0.0025,-0.0015,-0.0005,0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095,0.0105,0.0115,0.0125,0.0135,0.0145,0.0155,0.0165,0.0175,0.0185,0.0195,0.0205,0.0215,0.0225,0.0235,0.0245,0.0255,0.0265,0.0275,0.0285,0.0295,0.0305,0.0315,0.0325,0.0335,0.0345,0.0355,0.0365,0.0375,0.0385,0.0395,0.0405,0.0415,0.0425,0.0435,0.0445,0.0455,0.0465,0.0475,0.0485,0.0495,0.0505,0.0515,0.0525,0.0535,0.0545,0.0555,0.0565,0.0575,0.0585,0.0595,0.0605,0.0615,0.0625,0.0635,0.0645,0.0655,0.0665,0.0675,0.0685,0.0695,0.0705,0.0715,0.0725,0.0735,0.0745,0.0755,0.0765,0.0775,0.0785,0.0795,0.0805,0.0815,0.0825,0.0835,0.0845,0.0855,0.0865,0.0875,0.0885,0.0895,0.0905,0.0915,0.0925,0.0935,0.0945,0.0955,0.0965,0.0975,0.0985,0.0995,0.1005,0.1015,0.1025,0.1035,0.1045,0.1055,0.1065,0.1075,0.1085,0.1095,0.1105,0.1115,0.1125,0.1135,0.1145,0.1155,0.1165,0.1175,0.1185,0.1195,0.1205,0.1215,0.1225,0.1235,0.1245,0.1255,0.1265,0.1275,0.1285,0.1295,0.1305,0.1315,0.1325,0.1335,0.1345,0.1355,0.1365,0.1375,0.1385,0.1395,0.1405,0.1415,0.1425,0.1435,0.1445,0.1455,0.1465,0.1475,0.1485,0.1495,0.1505,0.1515,0.1525,0.1535,0.1545,0.1555,0.1565,0.1575,0.1585,0.1595,0.1605,0.1615,0.1625,0.1635,0.1645,0.1655,0.1665,0.1675,0.1685,0.1695,0.1705,0.1715,0.1725,0.1735,0.1745,0.1755,0.1765,0.1775,0.1785,0.1795,0.1805,0.1815,0.1825,0.1835,0.1845,0.1855,0.1865,0.1875,0.1885,0.1895,0.1905,0.1915,0.1925,0.1935,0.1945,0.1955,0.1965,0.1975,0.1985,0.1995,0.2005,0.2015,0.2025,0.2035,0.2045,0.2055,0.2065,0.2075,0.2085,0.2095,0.2105,0.2115,0.2125,0.2135,0.2145,0.2155,0.2165,0.2175,0.2185,0.2195,0.2205,0.2215,0.2225,0.2235,0.2245,0.2255,0.2265,0.2275,0.2285,0.2295,0.2305};


const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
Float_t MultV0M_data, MultRef_data;
class AliAnalysisTaskDataKno;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskDataKno) // classimp: necessary for root

AliAnalysisTaskDataKno::AliAnalysisTaskDataKno() : AliAnalysisTaskSE(),
  fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0), fUseMC(kFALSE), fIsMCclosure(kFALSE), fIspPb(kFALSE), fIsTPConly(kTRUE), fTPCclustersVar1(kFALSE), fTPCclustersVar2(kFALSE), fNcrVar1(kFALSE), fNcrVar2(kFALSE), fGeoTPCVar1(kFALSE), fGeoTPCVar2(kFALSE), fGeoTPCVar3(kFALSE), fGeoTPCVar4(kFALSE), fChisqTPCVar1(kFALSE), fChisqTPCVar2(kFALSE), fChisqITSVar1(kFALSE), fChisqITSVar2(kFALSE), fChisqITSmTPCVar1(kFALSE), fChisqITSmTPCVar2(kFALSE), fDcazVar1(kFALSE), fDcazVar2(kFALSE), fSPDreqVar1(kFALSE), fLeadingTrackFilter(0x0), fTrackFilter(0x0),fTrackFilterwoDCA(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5),fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fV0Mmin(0.0),fV0Mmax(100.0), fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0), fRecLeadPt(0), fRecLeadIn(0),ftrackmult08(0), fv0mpercentile(0), fv0mpercentilebefvtx(0), fdcaxy(-999), fdcaz(-999), fMultSelection(0x0), fMultSelectionbefvtx(0x0), hNchTSGen(0), hNchGen(0), hNchTSRec(0), hNchData(0), hNchTSData(0),hNchRec(0), hCounter(0),hRefMult08(0), hV0Mmult(0), hV0Mmultbefvtx(0), hRefMultvsV0Mmult(0),hV0MmultvsUE(0),hRefmultvsUE(0), hITSclustersvsUE(0),hITSclustersvsNch(0), hPtVsV0MData(0), hDphiVsUEData(0), hDphiVsNchData(0), hPTVsDCAData(0), hPTVsDCAcentData(0)
{
  for(Int_t i=0;i<3;++i){
    hPtVsUEData[i]=0;
    hPtVsNchData[i]=0;
    hPhiGen[i]=0;
    hPhiRec[i]=0;
  }


  // default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskDataKno::AliAnalysisTaskDataKno(const char* name) : AliAnalysisTaskSE(name),
								   fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0), fUseMC(kFALSE), fIsMCclosure(kFALSE), fIspPb(kFALSE),fIsTPConly(kTRUE), fTPCclustersVar1(kFALSE), fTPCclustersVar2(kFALSE), fNcrVar1(kFALSE), fNcrVar2(kFALSE), fGeoTPCVar1(kFALSE), fGeoTPCVar2(kFALSE), fGeoTPCVar3(kFALSE), fGeoTPCVar4(kFALSE), fChisqTPCVar1(kFALSE), fChisqTPCVar2(kFALSE), fChisqITSVar1(kFALSE), fChisqITSVar2(kFALSE), fChisqITSmTPCVar1(kFALSE), fChisqITSmTPCVar2(kFALSE), fDcazVar1(kFALSE), fDcazVar2(kFALSE),fSPDreqVar1(kFALSE), fLeadingTrackFilter(0x0), fTrackFilter(0x0),fTrackFilterwoDCA(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5),fLeadPtCutMin(5.0), fLeadPtCutMax(40.0),fV0Mmin(0.0),fV0Mmax(100.0),  fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0), fRecLeadPt(0), fRecLeadIn(0),ftrackmult08(0), fv0mpercentile(0), fv0mpercentilebefvtx(0), fdcaxy(-999), fdcaz(-999),fMultSelection(0x0), fMultSelectionbefvtx(0x0), hNchTSGen(0), hNchGen(0), hNchTSRec(0), hNchData(0), hNchTSData(0), hNchRec(0), hCounter(0),hRefMult08(0),hV0Mmult(0), hV0Mmultbefvtx(0), hRefMultvsV0Mmult(0), hV0MmultvsUE(0),hRefmultvsUE(0),  hITSclustersvsUE(0), hITSclustersvsNch(0), hPtVsV0MData(0), hDphiVsUEData(0), hDphiVsNchData(0), hPTVsDCAData(0), hPTVsDCAcentData(0)
{
  for(Int_t i=0;i<3;++i){
    hPtVsUEData[i]=0;
    hPtVsNchData[i]=0;
    hPhiGen[i]=0;
    hPhiRec[i]=0;
  }

  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it, does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskDataKno::~AliAnalysisTaskDataKno()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    fOutputList = 0x0;
  }

}
//_____________________________________________________________________________
void AliAnalysisTaskDataKno::UserCreateOutputObjects()
{

  // parametrization of efficiency
  // fCuts *** leading particle ***
  if(!fLeadingTrackFilter){
    fLeadingTrackFilter = new AliAnalysisFilter("trackFilter2015");
    AliESDtrackCuts * fCuts1 = new AliESDtrackCuts();
    fCuts1->SetAcceptKinkDaughters(kFALSE);// Default
    fCuts1->SetRequireTPCRefit(kTRUE);// Default
    fCuts1->SetRequireITSRefit(kTRUE);// Default
    fCuts1->SetDCAToVertex2D(kFALSE);// Default
    fCuts1->SetRequireSigmaToVertex(kFALSE);// Default
    fCuts1->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");// Default

    if (fTPCclustersVar1) {fCuts1->SetMaxFractionSharedTPCClusters(0.2);}
    else if (fTPCclustersVar2) {fCuts1->SetMaxFractionSharedTPCClusters(1.);}
    else {fCuts1->SetMaxFractionSharedTPCClusters(0.4);}// Default

    if (fNcrVar1) {fCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}//
    else if (fNcrVar2) {fCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}//
    else {fCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);}// Default

    if (fGeoTPCVar1) {fCuts1->SetCutGeoNcrNcl(2., 130., 1.5, 0.85, 0.7);}//
    else if (fGeoTPCVar2) {fCuts1->SetCutGeoNcrNcl(4., 130., 1.5, 0.85, 0.7);}//
    else if (fGeoTPCVar3) {fCuts1->SetCutGeoNcrNcl(3., 120., 1.5, 0.85, 0.7);}//
    else if (fGeoTPCVar4) {fCuts1->SetCutGeoNcrNcl(3., 140., 1.5, 0.85, 0.7);}//
    else {fCuts1->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);}// Default

    if (fChisqTPCVar1) {fCuts1->SetMaxChi2PerClusterTPC(3);}
    else if (fChisqTPCVar2) {fCuts1->SetMaxChi2PerClusterTPC(5);}
    else {fCuts1->SetMaxChi2PerClusterTPC(4);}// Default

    if (!fSPDreqVar1) {fCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);}// Default    

    if (fChisqITSmTPCVar1) {fCuts1->SetMaxChi2TPCConstrainedGlobal(25);}//
    else if (fChisqITSmTPCVar2) {fCuts1->SetMaxChi2TPCConstrainedGlobal(49);}//
    else {fCuts1->SetMaxChi2TPCConstrainedGlobal(36);}// Default

    if (fDcazVar1) {fCuts1->SetMaxDCAToVertexZ(1);} // DCAz = 1 cm
    else if (fDcazVar2) {fCuts1->SetMaxDCAToVertexZ(5);} // DCAz = 5 cm
    else {fCuts1->SetMaxDCAToVertexZ(2);}// Default

    if (fChisqITSVar1) {fCuts1->SetMaxChi2PerClusterITS(25);}//
    else if (fChisqITSVar2) {fCuts1->SetMaxChi2PerClusterITS(49);}//
    else {fCuts1->SetMaxChi2PerClusterITS(36);}// Default

    fLeadingTrackFilter->AddCuts(fCuts1);
  }

  ///track quality =====
  // TPC ***  multiplicity in transverse side ***
  if(!fTrackFilter){
    if(fIsTPConly) //Default option
      {
    fTrackFilter = new AliAnalysisFilter("trackFilterTPConly");
    AliESDtrackCuts * fCuts2 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fCuts2->SetRequireTPCRefit(kTRUE);
    fCuts2->SetRequireITSRefit(kTRUE);
    fCuts2->SetEtaRange(-0.8,0.8);
    fTrackFilter->AddCuts(fCuts2);
      }
    else //For systematic uncertainties (same cut for NchTS and tracks)
      {
	fTrackFilter = new AliAnalysisFilter("trackFilter2015");
	AliESDtrackCuts * fCuts2_1 = new AliESDtrackCuts();
	fCuts2_1->SetMaxFractionSharedTPCClusters(0.4);//
	fCuts2_1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
	fCuts2_1->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
	fCuts2_1->SetMaxChi2PerClusterTPC(4);//
	fCuts2_1->SetAcceptKinkDaughters(kFALSE);//
	fCuts2_1->SetRequireTPCRefit(kTRUE);//
	fCuts2_1->SetRequireITSRefit(kTRUE);//
	fCuts2_1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);//
	fCuts2_1->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
	fCuts2_1->SetMaxChi2TPCConstrainedGlobal(36);//
	fCuts2_1->SetMaxDCAToVertexZ(2);//
	fCuts2_1->SetDCAToVertex2D(kFALSE);//
	fCuts2_1->SetRequireSigmaToVertex(kFALSE);//
	fCuts2_1->SetMaxChi2PerClusterITS(36);//
	fCuts2_1->SetEtaRange(-0.8,0.8);
	fTrackFilter->AddCuts(fCuts2_1);
      }
  }

  //track cuts to find contamination via DCA distribution
  if (!fTrackFilterwoDCA){

    fTrackFilterwoDCA = new AliAnalysisFilter("trackFilter2015");
    AliESDtrackCuts * fCuts3 = new AliESDtrackCuts();
    fCuts3->SetMaxFractionSharedTPCClusters(0.4);//
    fCuts3->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
    fCuts3->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
    fCuts3->SetMaxChi2PerClusterTPC(4);//
    fCuts3->SetAcceptKinkDaughters(kFALSE);//
    fCuts3->SetRequireTPCRefit(kTRUE);//
    fCuts3->SetRequireITSRefit(kTRUE);//
    fCuts3->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
	 			     AliESDtrackCuts::kAny);//
    //fCuts3->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
    //fCuts3->SetMaxChi2TPCConstrainedGlobal(36);//
    fCuts3->SetMaxDCAToVertexZ(2);//
    fCuts3->SetDCAToVertex2D(kFALSE);//
    fCuts3->SetRequireSigmaToVertex(kFALSE);//
    fCuts3->SetMaxChi2PerClusterITS(36);//
    fTrackFilterwoDCA->AddCuts(fCuts3);
  }
    
  // create output objects

  OpenFile(1);
  fOutputList = new TList();          // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written  to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects and will delete them if requested

  hNchTSRec = new TH1D("hNchTSRec","",3000,-0.5,2999.5);
  fOutputList->Add(hNchTSRec);

  hNchTSData = new TH1D("hNchTSData","",3000,-0.5,2999.5); 
  fOutputList->Add(hNchTSData);

  hNchRec = new TH1D("hNchRec","",3000,-0.5,2999.5);
  fOutputList->Add(hNchTSRec);

  hNchData = new TH1D("hNchData","",3000,-0.5,2999.5); 
  fOutputList->Add(hNchData);

  hRefMult08 = 0;
  hRefMult08 = new TH1D("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);N_{ch};count",3000,-0.5,2999.5);   
  fOutputList->Add(hRefMult08);

  hV0Mmult = 0;
  hV0Mmult = new TH1D("hV0Mmult","V0M ;V0M percentile;count",100,0,100);   
  fOutputList->Add(hV0Mmult);

  hV0Mmultbefvtx = 0;
  hV0Mmultbefvtx = new TH1D("hV0Mmultbefvtx","V0M ;V0M percentile bef. vtx;count",100,0,100);   
  fOutputList->Add(hV0Mmultbefvtx);

  hRefMultvsV0Mmult = 0;
  hRefMultvsV0Mmult = new TH2D("hRefMultvsV0Mmult","N_{ch} vs V0M percentile;N_{ch}; v0M percentile",3000,-0.5,2999.5,10,0,100);
  fOutputList->Add(hRefMultvsV0Mmult);

  hV0MmultvsUE = 0;
  hV0MmultvsUE = new TH2D("hV0MmultvsUE","V0M percentile vs NchTS; v0M percentile;N_{ch}^{TS}",100,0,100,3000,-0.5,2999.5);
  fOutputList->Add(hV0MmultvsUE);

  hRefmultvsUE = 0;
  hRefmultvsUE = new TH2D("hRefmultvsUE","Ref Mult vs NchTS; Ref. mult.;N_{ch}^{TS}",3000,-0.5,2999.5,3000,-0.5,2999.5);
  fOutputList->Add(hRefmultvsUE);

  hITSclustersvsUE = 0;
  hITSclustersvsUE = new TH2D("hITSclustersvsUE","ITSclusters vs NchTS; ITSclusters;N_{ch}^{TS}",3000,-0.5,2999.5,3000,-0.5,2999.5);
  fOutputList->Add(hITSclustersvsUE);

  hITSclustersvsNch = 0;
  hITSclustersvsNch = new TH2D("hITSclustersvsNch","ITSclusters vs Nch; ITSclusters;N_{ch}",3000,-0.5,2999.5,3000,-0.5,2999.5);
  fOutputList->Add(hITSclustersvsNch);

  for(Int_t i=0;i<3;++i){
    hPhiRec[i]= new TH1D(Form("hPhiRec_%s",NameReg_1_data[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
    fOutputList->Add(hPhiRec[i]);
  }

   hPtVsV0MData = new TH2D("hPtVsV0MData","data pT vs V0M",100,0.,100.,ptNbins_data,ptbins1_1_data);
   fOutputList->Add(hPtVsV0MData);
     
  for(Int_t i=0;i<3;++i){

    hPtVsUEData[i] = new TH2D(Form("hPtVsUEData_%s",NameReg_1_data[i]),"data pT vs nch_transverse",3000,-0.5,2999.5,ptNbins_data,ptbins1_1_data);
    fOutputList->Add(hPtVsUEData[i]);

    hPtVsNchData[i] = new TH2D(Form("hPtVsNchData_%s",NameReg_1_data[i]),"data pT vs nch_transverse",3000,-0.5,2999.5,ptNbins_data,ptbins1_1_data);
    fOutputList->Add(hPtVsNchData[i]);

  }

  hDphiVsUEData = new TH2D("hDphiVsUEData","Delta phi vs nch_transverse",3000,-0.5,2999.5,nDeltabins_data,Deltabins_data);
  fOutputList->Add(hDphiVsUEData);

  hDphiVsNchData = new TH2D("hDphiVsNchData","Delta phi vs nch_transverse",3000,-0.5,2999.5,nDeltabins_data,Deltabins_data);
  fOutputList->Add(hDphiVsNchData);


  //hV0MVsUEvsRef = 0;
  //hV0MVsUEvsRef = new TH3D("hV0MVsUEvsRef","nch_transverse vs ref multiplicity vs V0M",nTSBins_1_data,TSBins_1,nchNbins_data_3,nchbins_3_data,nchNbins_data_3,nchbins_3_data);
  //fOutputList->Add(hV0MVsUEvsRef);

  hPTVsDCAData = new TH2D("hPTVsDCAData","pT vs dcaxy",ptNbins_data,ptbins1_1_data,nBinsDCAxy_data,binsDCAxy_data);
  fOutputList->Add(hPTVsDCAData);

  hPTVsDCAcentData = new TH2D("hPTVsDCAcentData","pT vs dcaxy",ptNbins_data,ptbins1_1_data,nBinsDCAxy_dataCentral,binsDCAxyCentral_data);
  fOutputList->Add(hPTVsDCAcentData);

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
	
}
//_____________________________________________________________________________
void AliAnalysisTaskDataKno::UserExec(Option_t *)
{

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  fESD = dynamic_cast<AliESDEvent*>(event);

  if(!fESD){
    Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }

  if (fUseMC) {

    //      E S D
    fMC = dynamic_cast<AliMCEvent*>(MCEvent());
    if(!fMC){
      Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
    fMCStack = fMC->Stack();
  }


  //hCounter->Fill(0);
  // Before trigger selection
  GetLeadingObject(kFALSE);// leading particle at rec level
  // Now we get the leading particle pT
	


  // Trigger selection
  UInt_t fSelectMask= fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
  if(!isINT7selected)
    return;
  //hCounter->Fill(1);

  // Good events
  if (!fEventCuts.AcceptEvent(event)) {
    PostData(1, fOutputList);
    return;
  }
  
  fMultSelectionbefvtx = (AliMultSelection*) fESD->FindListObject("MultSelection");
  if (fIspPb) {fv0mpercentilebefvtx = fMultSelectionbefvtx->GetMultiplicityPercentile("V0A");}
  else {fv0mpercentilebefvtx = fMultSelectionbefvtx->GetMultiplicityPercentile("V0M");}

  hV0Mmultbefvtx->Fill(fv0mpercentilebefvtx);

   // Good vertex
  Bool_t hasRecVertex = kFALSE;
  hasRecVertex=HasRecVertex();
  if(!hasRecVertex)return;

  // Multiplicity Estimation
  ftrackmult08 = -999;
  fv0mpercentile = -999;
	
  ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //tracklets
  hRefMult08->Fill(ftrackmult08);

  fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
  //if (!fMultSelection)
  //cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
  if (fIspPb) {fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0A");}
  else {fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");}
	  
  hV0Mmult->Fill(fv0mpercentile);

  //cout<<"------- V0M mult ==  "<<fv0mpercentile<<"--------"<<endl;
	
  hRefMultvsV0Mmult->Fill(ftrackmult08,fv0mpercentile);
	
  //analysis
  if (fv0mpercentile>fV0Mmin && fv0mpercentile<=fV0Mmax)
    {
	  GetMB();
	  // KNO scaling
	  if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
	    GetMultiplicityDistributionsData();

    }

  PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file

}


//______________________________________________________________________________
void AliAnalysisTaskDataKno::Terminate(Option_t *)
{

}


//_____________________________________________________________________________
void AliAnalysisTaskDataKno::GetLeadingObject(Bool_t isMC) {

  Double_t flPt = 0;// leading pT
  Double_t flPhi = 0;
  Int_t flIndex = 0;

  if(isMC){
    for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

      AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
      if (!particle) continue;

      if (!fMC->IsPhysicalPrimary(i)) continue;  //  Particles produced including products of strong and electromagnetic decays but excluding feed-down from weak decays of strange particles like Ks,Lambda etc)
      if (particle->Charge() == 0) continue;
      if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
      if( particle->Pt() < fPtMin)continue;

      if (flPt<particle->Pt()){
	flPt = particle->Pt();
	flPhi = particle->Phi();
	flIndex = i;
      }
    }

    fGenLeadPhi = flPhi;
    fGenLeadPt  = flPt;
    fGenLeadIn  = flIndex;
  }
  else{

    Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
    for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

      AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

      if(!track) continue;

      if(!fLeadingTrackFilter->IsSelected(track))
	continue;

      if(TMath::Abs(track->Eta()) > fEtaCut)
	continue;

      if( track->Pt() < fPtMin)continue;

      if (flPt<track->Pt()){
	flPt  = track->Pt();
	flPhi = track->Phi();
	flIndex = i;
      }

    } 
    fRecLeadPhi = flPhi;
    fRecLeadPt  = flPt;
    fRecLeadIn  = flIndex;
  }


}

void AliAnalysisTaskDataKno::GetMB(){
  Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
  for(Int_t i=0; i < iTracks; i++) {

    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

    if(!track) continue;

    if(!fTrackFilter->IsSelected(track))
      continue;

    if( track->Pt() < 0.15)continue;

    if(TMath::Abs(track->Eta()) > fEtaCut)
      continue;

    hPtVsV0MData->Fill(fv0mpercentile,track->Pt());

  }
  Int_t multTSrec=0;
  Int_t multrec=0;
  
   for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks


    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

    if(!track) continue;

    if(!fTrackFilter->IsSelected(track))
      continue;

    if(TMath::Abs(track->Eta()) > fEtaCut)
      continue;

    if( track->Pt() < 0.15)continue;

    Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
    multrec++;
    // definition of the topological regions
    if(TMath::Abs(DPhi)<pi/3.0){// near side
      continue;
    }
    else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
      continue;
    }
    else{// transverse side
      multTSrec++;
    }

  }

   for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

    if(!track) continue;

    if(!fLeadingTrackFilter->IsSelected(track))
      continue;

    if(TMath::Abs(track->Eta()) > fEtaCut)
      continue;

    if( track->Pt() < 0.15)continue;

    hITSclustersvsUE->Fill(track->GetITSNcls(),multTSrec);
    hITSclustersvsNch->Fill(track->GetITSNcls(),multrec);
    
   }
}
//______________________________________________________________
void AliAnalysisTaskDataKno::GetMultiplicityDistributionsData(){
  
  Int_t multTSrec=0;
  Int_t multrec=0;

  //if(fv0mpercentile>fV0Mmin && fv0mpercentile<=fV0Mmax)
  //{	
  Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
  for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

    if(i==fRecLeadIn)
      continue;

    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

    if(!track) continue;

    if(!fTrackFilter->IsSelected(track))
      continue;

    if(TMath::Abs(track->Eta()) > fEtaCut)
      continue;

    if( track->Pt() < fPtMin)continue;

    Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

    // definition of the topological regions
    if(TMath::Abs(DPhi)<pi/3.0){// near side
      continue;
    }
    else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
      continue;
    }
    else{// transverse side
      multTSrec++;
    }

    multrec++;

  }
  hNchTSData->Fill(multTSrec);
  hNchData->Fill(multrec);
	
  hV0MmultvsUE->Fill(fv0mpercentile,multTSrec);
  hRefmultvsUE->Fill(ftrackmult08,multTSrec);

  // Filling rec pT vs UE (for pT I use 2015 track cuts, UE uses TPC-only)
  for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

    if(i==fRecLeadIn)
      continue;

    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

    if(!track) continue;

    if(!fLeadingTrackFilter->IsSelected(track))
      continue;

    if(TMath::Abs(track->Eta()) > fEtaCut)
      continue;

    if( track->Pt() < fPtMin)continue;

    Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
    

    // definition of the topological regions
    if(TMath::Abs(DPhi)<pi/3.0){// near side
      hPtVsUEData[0]->Fill(multTSrec,track->Pt());
      hPtVsNchData[0]->Fill(multrec,track->Pt());
    }
    else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
      hPtVsUEData[1]->Fill(multTSrec,track->Pt());
      hPtVsNchData[1]->Fill(multrec,track->Pt());			
    }
    else{// transverse side
      hPtVsUEData[2]->Fill(multTSrec,track->Pt());
      hPtVsNchData[2]->Fill(multrec,track->Pt());
    }

    hDphiVsUEData->Fill(multTSrec,DPhi);
    hDphiVsNchData->Fill(multrec,DPhi);
      
  }
  //}
  for(Int_t i=0; i < iTracks; i++) {

     if(i==fRecLeadIn)
      continue;
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)
    if(!track) continue;
    if(!fTrackFilterwoDCA->IsSelected(track))
      continue;
    if(TMath::Abs(track->Eta()) > fEtaCut)
      continue;
    if( track->Pt() < fPtMin)continue;

    track->GetImpactParameters(fdcaxy,fdcaz);

    hPTVsDCAData->Fill(track->Pt(),fdcaxy);
    hPTVsDCAcentData->Fill(track->Pt(),fdcaxy);

  }
  

}
//____________________________________________________________

Double_t AliAnalysisTaskDataKno::DeltaPhi(Double_t phia, Double_t phib,
					Double_t rangeMin, Double_t rangeMax)
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
Bool_t AliAnalysisTaskDataKno::HasRecVertex(){


  float fMaxDeltaSpdTrackAbsolute = 0.5f;
  float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  float fMaxResolutionSPDvertex = 0.25f;
  float fMaxDispersionSPDvertex = 1.e14f;

  Bool_t fRequireTrackVertex = true;
  unsigned long fFlag;
  fFlag =BIT(AliEventCuts::kNoCuts);

  const AliVVertex* vtTrc = fESD->GetPrimaryVertex();
  bool isTrackV = true;
  if(vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ()) isTrackV=false;
  const AliVVertex* vtSPD = fESD->GetPrimaryVertexSPD();


  if (vtSPD->GetNContributors() > 0) fFlag |= BIT(AliEventCuts::kVertexSPD);

  if (vtTrc->GetNContributors() > 1 && isTrackV ) fFlag |= BIT(AliEventCuts::kVertexTracks);

  if (((fFlag & BIT(AliEventCuts::kVertexTracks)) ||  !fRequireTrackVertex) && (fFlag & BIT(AliEventCuts::kVertexSPD))) fFlag |= BIT(AliEventCuts::kVertex);

  const AliVVertex* &vtx = bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
  AliVVertex   *fPrimaryVertex = const_cast<AliVVertex*>(vtx);
  if(!fPrimaryVertex)return kFALSE;

  /// Vertex quality cuts
  double covTrc[6],covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = bool(fFlag & AliEventCuts::kVertexSPD) && bool(fFlag & AliEventCuts::kVertexTracks) ? vtTrc->GetZ() - vtSPD->GetZ() : 0.; /// If one of the two vertices is not available this cut is always passed.
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex* vtSPDESD = dynamic_cast<const AliESDVertex*>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if (
      (TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute && nsigTot <= fMaxDeltaSpdTrackNsigmaSPD && nsigTrc <= fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
      (!vtSPD->IsFromVertexerZ() || TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
      (!vtSPD->IsFromVertexerZ() || vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for run1, only for ESD
      ) // quality cut on vertexer SPD z
    fFlag |= BIT(AliEventCuts::kVertexQuality);  

  Bool_t hasVtx = (TESTBIT(fFlag,AliEventCuts::kVertex))&&(TESTBIT(fFlag,AliEventCuts::kVertexQuality));

  return hasVtx;

}


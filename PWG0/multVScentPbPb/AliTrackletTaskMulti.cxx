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
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
// Class AliTrackletTaskMulti                                            //
// Analysis task to produce data and MC histos needed for tracklets      //
// dNdEta extraction in multiple bins in one go                          //
// Author:  ruben.shahoyan@cern.ch                                       //
///////////////////////////////////////////////////////////////////////////
/*
  Important parameters to set:
  1) make sure to initialize correct geometry in UserCreateOutputObjects
  2) The cut on signal selection variable (delta, dphi ...) should be decided beforehand
...
*/

#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TH3F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TNtuple.h"
#include "TObjArray.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliESDEvent.h"  
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
#include "../ANALYSIS/EventMixing/AliMixEventInputHandler.h"
#include "AliCDBPath.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliMagF.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliRunLoader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "../ITS/AliITSRecPoint.h"
#include "../ITS/AliITSgeomTGeo.h"
#include "../ITS/AliITSMultReconstructor.h" 

#include "AliLog.h"

#include "AliPhysicsSelection.h"
#include "AliESDCentrality.h" 
#include "AliTrackletTaskMulti.h"
#include "AliITSMultRecBg.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliESDtrackCuts.h"

ClassImp(AliTrackletTaskMulti)

// centrality definition with SPD2corr
const Float_t  AliTrackletTaskMulti::fgkCentBinDefSPD2[] = 
{0,21,69,165,345,637,1069,1683,2511,3721,4585,6500}; // Alberica 2010-11-21 at 22:17 
//corresponding npart: {379.973,328.414,259.244,185.094,128.309,84.9333,52.7138,29.9245,15.3195,6.98915,3.16337};

// special test with splitting upper 10% to 4 equidistant bins
//  {3721,4125,4584,5102,6500};


//  {0,29,85,191,385,687,1133,1743,2567,3765,4611,6500}; // before 21/11/10 
//corresponding npart: {3.14467,6.876,15.0519,29.0842,52.7315,84.8281,128.265,185.004,259.304,328.651,379.393};


// centrality definition with V0 Rescaled
const Float_t AliTrackletTaskMulti::fgkCentBinDefV0CR[] = 
{0,32.5,89.5,198.5,396.5,700.5,1140.5,1744.5,2562.5,3767.5,4647.5,6500.};
// on V0 corrected and rescaled
// corresponding npart:  {3.23659,7.28226,15.9883,31.4739,54.7531,87.0497,130.454,185.996,260.318,329.655,379.114}


const Float_t AliTrackletTaskMulti::fgkCentBinDefV0[] = 
  {12191,13529,15079,16815,21000}; // special test with splitting upper 10% to 4 equidistant bins for new baseline with ZDC timing cleanup, 11/29/2010 01:10:37 PM

//  {0,79,239,559,1165,2135,3555,5525,8213,12191,15079,21000}; // new baseline with ZDC timing cleanup, 11/29/2010 01:10:37 PM

  //  {12165,13527,15079,16801,21000}; // special test with splitting upper 10% to 4 equidistant bins for Alberica 2010-11-21 at 22:17 
//  {0,79,247,577,1185,2155,3565,5527,8203,12167,15073,21000}; // Alberica 2010-11-21 at 22:17 
// corresponding npart: {379.112,323.171,253.417,183.054,129.116,86.9088,54.2474,30.2884,15.1553,6.7094,2.95315}

//  {0,107,297,659,1305,2301,3747,5715,8361,12307,15153,21000}; // before 21/11/10 on V0 corrected for non-linearity
// corresponding npart: {2.93958,6.67109,15.0721,29.7634,52.6495,84.8531,128.160,185.197,259.314,328.781,381.040};

//{0,105,291,644,1280,2269,3728,5704,8366,12111,14706,19540}; // after  18/11/10, V0 not corrected
//{0,124.5,274.5,574.5,1224.5,2174.5,3624.5,5574.5,8274.5,12024.5,14674.5,20000}; // before 18/11/10


// centrality selection with TPC tracls only
const Float_t AliTrackletTaskMulti::fgkCentBinDefTrTPC[] = 
  {0,13,37,85,171,307,507,783,1157,1685,2055,2900};
//corresponding npart:  {3.22604,6.98597,15.1032,29.7764,52.6338,84.7075,128.221,185.129,259.17,329.234,378.667};

//------------------------- 2D boundaries ---------------------------------------
const Float_t AliTrackletTaskMulti::fgkCentBinDefZDCV0X[] =
  {0.,150.0,310.0,610.0,1190.0,2170.0,3570.0,5530.0,8210.0,12150.0,15070.0,21000.};
const Float_t AliTrackletTaskMulti::fgkCentBinDefZDCV0Y[] =
  {0.0,1224.6,1999.4,2665.2,3106.0,3312.4,3343.0,3238.1,2758.1,1697.8,978.7,0.00};
const Float_t AliTrackletTaskMulti::fgkCentBinDefZDCV0S[] =
  {-1.,-0.1225,-0.3001,-0.9252,-2.0376,-10.1909,62.0849,9.1496,4.1895,3.6557,4.7364,1.};

/*
// 80% of ZDC&V0 set to 247 as 80% of 1D V0
const Float_t AliTrackletTaskMulti::fgkCentBinDefZDCV0X[] =
  {0.,130.,247.,490.,1010.0,1950.,3310.,5270.,7930.,11930.,14930.,21000.};
const Float_t AliTrackletTaskMulti::fgkCentBinDefZDCV0Y[] =
  {0.00,1061.35,1732.80,2459.06,3017.64,3288.97,3346.46,3264.99,2822.64,1758.14,1008.61,0.00};
const Float_t AliTrackletTaskMulti::fgkCentBinDefZDCV0S[] =
  {-1.,0.1225,-0.3001,-0.5631,-2.0376,-7.3697,125.8992,10.2539,4.4412,3.6288,4.6181,1.};
*/

const Float_t *fkCentBinDef    = 0;
const Float_t *fkCentBinDef2DY = 0;
const Float_t *fkCentBinDef2DS = 0;

const char*  AliTrackletTaskMulti::fgkPDGNames[] = {
"#pi^{+}",
"p",
"K^{+}",
"K^{*+}",
"e^{-}",
"#mu^{-}",
"#rho^{+}",
"D^{+}",
"D^{*+}",
"D_{s}^{+}",
"D_{s}^{*+}",
"#Delta^{-}",
"#Delta^{+}",
"#Delta^{++}",
"#Sigma^{-}",
"#Sigma^{+}",
"#Sigma^{*-}",
"#Sigma^{*+}",
"#Sigma^{*+}_{c}",
"#Sigma^{*++}_{c}",
"#Xi^{-}",
"#Xi^{*-}",
"#Lambda^{+}_{c}",
"n",
"#Delta^{0}",
"#gamma",
"K^{0}_{S}",
"K^{0}_{L}",
"K^{0}",
"K^{*}",
"#eta",
"#pi^{0}",
"#rho^{0}",
"#varphi",
"#eta'",
"#omega",
"#Lambda",
"#Sigma^{0}",
"#Sigma^{*0}_{c}",
"#Sigma^{*0}",
"D^{0}",
"D^{*0}",
"#Xi_{0}",
"#Xi^{*0}",
"#Xi^{0}_{c}",
"#Xi^{*0}_{c}",
"Nuclei",
"Others"
};

const int AliTrackletTaskMulti::fgkPDGCodes[] = {
  211,
 2212, 
  321, 
  323, 
   11, 
   13, 
  213, 
  411, 
  413, 
  431, 
  433, 
 1114, 
 2214, 
 2224, 
 3112, 
 3222, 
 3114, 
 3224, 
 4214, 
 4224, 
 3312, 
 3314, 
 4122, 
 2112, 
 2114, 
   22, 
  310, 
  130, 
  311, 
  313, 
  221, 
  111, 
  113, 
  333, 
  331, 
  223, 
 3122, 
 3212, 
 4114, 
 3214, 
  421, 
  423, 
 3322, 
 3324, 
 4132, 
 4314
// nuclei
// unknown
};

//________________________________________________________________________
/*//Default constructor
AliTrackletTaskMulti::AliTrackletTaskMulti(const char *name)
  : AliAnalysisTaskSE(name),
*/  
//________________________________________________________________________
AliTrackletTaskMulti::AliTrackletTaskMulti(const char *name) 
  : AliAnalysisTaskSE(name), 
//
  fOutput(0), 
//
  fDoNormalReco(kFALSE),
  fDoInjection(kFALSE),
  fDoRotation(kFALSE),
  fDoMixing(kFALSE),
  //
  fUseMC(kFALSE),
  fCheckReconstructables(kFALSE),
//
  fHistosTrData(0),
  fHistosTrInj(0),
  fHistosTrRot(0),
  fHistosTrMix(0),
//
  fHistosTrPrim(0),
  fHistosTrSec(0),
  fHistosTrComb(0),
  fHistosTrCombU(0),
//
  fHistosTrRcblPrim(0),
  fHistosTrRcblSec(0),
  fHistosCustom(0),
//
  fEtaCut(3.0),
  fZVertexMin(-20),
  fZVertexMax( 20),
//
  fScaleDTBySin2T(kFALSE),
  fCutOnDThetaX(kFALSE),
  fNStdDev(1.),
  fDPhiWindow(0.08),
  fDThetaWindow(0.025),
  fDPhiShift(0.0045),
  fPhiOverlapCut(0.005),
  fZetaOverlap(0.05),
  fPhiRot(0.),
  fInjScale(1.),
  fRemoveOverlaps(kFALSE),
//
  fDPhiSCut(0.06),
  fNStdCut(1.),
  fMCV0Scale(0.7520),
//
  fMultReco(0),
  fRPTree(0),
  fRPTreeMix(0),
  fStack(0),
  fMCEvent(0),
  fTrackCuts(0),  
  //
  fNPart(0),
  fNBColl(0),
  fCurrCentBin(-1),
  fNCentBins(0),
  fUseCentralityVar(kCentV0)
  /*
  ,
  fTrigger(AliTriggerAnalysis::kAcceptAll),
  fMCCentralityBin(AliAnalysisTaskSPDdNdEta::kall),
  fCentrLowLim(0),
  fCentrUpLim(0),
  fCentrEst("")
  */
{
  // Constructor

  DefineOutput(1, TList::Class());
  //
  SetScaleDThetaBySin2T();
  SetNStdDev();
  SetPhiWindow();
  SetThetaWindow();
  SetPhiShift();
  SetPhiOverlapCut();
  SetZetaOverlapCut();
  SetPhiRot();
  SetRemoveOverlaps();
  //
}

//________________________________________________________________________
AliTrackletTaskMulti::~AliTrackletTaskMulti()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {  //RRR
    printf("Deleteing output\n");
    delete fOutput;
    fOutput = 0;
  }
  //
  delete fMultReco;
  delete fTrackCuts;
  //
  delete fHistosTrData;
  delete fHistosTrPrim;
  delete fHistosTrSec;
  delete fHistosTrComb;
  delete fHistosTrCombU;
  delete fHistosTrInj;
  delete fHistosTrRot;
  delete fHistosTrMix;
  delete fHistosTrRcblPrim;
  delete fHistosTrRcblSec;
  delete fHistosCustom;
  //
}

//________________________________________________________________________
void AliTrackletTaskMulti::UserCreateOutputObjects() 
{
  //
  fOutput = new TList();
  fOutput->SetOwner(); 
  //
  if      (fUseCentralityVar == kCentV0)   {
    fNCentBins = sizeof(fgkCentBinDefV0)/sizeof(Float_t) - 1;
    fkCentBinDef = fgkCentBinDefV0;
  }
  else if (fUseCentralityVar == kCentV0CR)   {
    fNCentBins = sizeof(fgkCentBinDefV0CR)/sizeof(Float_t) - 1;
    fkCentBinDef = fgkCentBinDefV0CR;
  }
  else if (fUseCentralityVar == kCentSPD2) {
    fNCentBins = sizeof(fgkCentBinDefSPD2)/sizeof(Float_t) - 1;
    fkCentBinDef = fgkCentBinDefSPD2;
  }
  else if (fUseCentralityVar == kCentTrTPC) {
    fNCentBins = sizeof(fgkCentBinDefTrTPC)/sizeof(Float_t) - 1;
    fkCentBinDef = fgkCentBinDefTrTPC;
  }
  else if (fUseCentralityVar == kCentZDCV0) {
    fNCentBins = sizeof(fgkCentBinDefZDCV0X)/sizeof(Float_t) - 1;
    fkCentBinDef    = fgkCentBinDefZDCV0X;
    fkCentBinDef2DY = fgkCentBinDefZDCV0Y;
    fkCentBinDef2DS = fgkCentBinDefZDCV0S;    
  }
  else {
    AliFatal(Form("Unknown cenrality parameter %d",fUseCentralityVar));
  }  
  //
  AliCDBManager *man = AliCDBManager::Instance();
  if (fUseMC) {
    Bool_t newGeom = kTRUE;
    man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
    if (newGeom) {
      // new geom
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data");
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130844,6,-1)) AliFatal("Failed to misalign geometry");
    }
    else {
      // old geom
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data");
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130845,5,-1)) AliFatal("Failed to misalign geometry");
    }
  }
  else {
    man->SetDefaultStorage("raw://"); man->SetRun(137045);
    AliCDBEntry*  obj = man->Get("GRP/Geometry/Data");
    AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
    if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",137045,8,-1)) AliFatal("Failed to misalign geometry");
  }
  //
  fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  //
  // Create histograms
  //---------------------------------------------Standard histos per tracklet type--->>
  UInt_t hPattern = 0xffffffff;
  fHistosTrData                      = BookHistosSet("TrData",hPattern);
  if (GetDoInjection()) fHistosTrInj = BookHistosSet("TrInj",hPattern);
  if (GetDoRotation())  fHistosTrRot = BookHistosSet("TrRot",hPattern);
  if (GetDoMixing())    fHistosTrMix = BookHistosSet("TrMix",hPattern);
  if (fUseMC) {
    fHistosTrPrim  = BookHistosSet("TrPrim",hPattern);
    fHistosTrSec   = BookHistosSet("TrSec",hPattern);
    fHistosTrComb  = BookHistosSet("TrComb",hPattern);
    fHistosTrCombU = BookHistosSet("TrCombU",hPattern);
    if (fCheckReconstructables) {
      fHistosTrRcblPrim = BookHistosSet("TrRcblPrim",hPattern);
      fHistosTrRcblSec  = BookHistosSet("TrRcblSec",hPattern);
    }
  }
  //---------------------------------------------Standard histos per tracklet type---<<
  //
  //---------------------------------------------Custom Histos----------------------->>
  // put here any non standard histos
  fHistosCustom = BookCustomHistos();
  //
  //---------------------------------------------Custom Histos-----------------------<<
  int nhist = fOutput->GetEntries();
  for (int i=0;i<nhist;i++) {
    TObject* hst = fOutput->At(i);
    if (!hst || !(hst->InheritsFrom(TH1::Class()))) continue;
    ((TH1*)hst)->Sumw2();
  }
  //
  PostData(1, fOutput);
  //
}

//________________________________________________________________________
void AliTrackletTaskMulti::UserExec(Option_t *) 
{
  // Main loop
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  fRPTree = fRPTreeMix = 0;
  AliESDInputHandlerRP *handRP = (AliESDInputHandlerRP*)anMan->GetInputEventHandler();
  if (!handRP) { printf("No RP handler\n"); return; }
  AliESDEvent *esd  = handRP->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }
  //
  // do we need to initialize the field?
  AliMagF* field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!field && !esd->InitMagneticField()) {printf("Failed to initialize the B field\n");return;}
  //
  /* // RS to be clarified
  // Trigger selection
  static AliTriggerAnalysis* triggerAnalysis = 0; 
  Bool_t eventTriggered = triggerAnalysis->IsTriggerFired(esd, fTrigger);
  if (!eventTriggered) {printf("No trigger\n"); return;}
  //
  // Centrality selection 
  Bool_t eventInCentralityBin = kFALSE;
  // Centrality selection
  AliESDCentrality *centrality = esd->GetCentrality();
  if (fCentrEst=="") eventInCentralityBin = kTRUE;
  else {
    if(!centrality) {
      AliError("Centrality object not available"); 
    }  else {
      if (centrality->IsEventInCentralityClass(fCentrLowLim,fCentrUpLim,fCentrEst.Data())) eventInCentralityBin = kTRUE;
    }
  }
  */
  //  if (!fUseMC && !ZDCTimeTrigger(esd)) return;
  //
  const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
  if (vtxESD->GetNContributors()<1) return;
  if (vtxESD->GetDispersion()>0.04) return;
  if (vtxESD->GetZRes()>0.25) return;
  const AliMultiplicity* multESD = esd->GetMultiplicity();
  const AliESDVertex* vtxESDTPC = esd->GetPrimaryVertexTPC();
  float nSPD1 = multESD->GetNumberOfITSClusters(0);
  float nSPD2 = multESD->GetNumberOfITSClusters(1);
  //
    /*
  if (vtxESDTPC->GetNContributors()<1 ||
      vtxESDTPC->GetNContributors()<(-10.+0.25*nSPD1)) return;
    */
  //
  TH1* hstat = (TH1*)fHistosCustom->UncheckedAt(kHStat);
  //
  hstat->Fill(kEvTot); // RS
  //
  Double_t esdvtx[3];
  vtxESD->GetXYZ(esdvtx);
  for (int i=3;i--;) fESDVtx[i] = esdvtx[i];
  //
  float vtxf[3] = {vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ()};
  //
  //------------------------------------------------------
  // ZDC cut
  AliESDZDC *esdZDC = esd->GetESDZDC();
  // --- ZDC offline trigger ---
  // Returns if ZDC triggered, based on TDC information
  Bool_t tdc[32] = {kFALSE};
  for(Int_t itdc=0; itdc<32; itdc++){
    for(Int_t i=0; i<4; i++){
      if (0.025*esdZDC->GetZDCTDCData(itdc, i) != 0){
	tdc[itdc] = kTRUE;
      }
    }
  }
  Bool_t zdcNA = tdc[12];
  Bool_t zdcNC = tdc[10];
  Bool_t zdcPA = tdc[13];
  Bool_t zdcPC = tdc[11];
  //
  Bool_t zdcA= ((zdcPA) || (zdcNA));
  Bool_t zdcC= ((zdcPC) || (zdcNC));
  //  if (!fUseMC && !(zdcA&&zdcC)) return;
  if (!fUseMC && !(zdcNA&&zdcNC)) return;
  float zdcEnergy = esdZDC->GetZDCN1Energy() + esdZDC->GetZDCP1Energy() + esdZDC->GetZDCN2Energy()+ esdZDC->GetZDCP2Energy();
  zdcEnergy /= 8.;
  //
  //-----------------------------------------------------
  Float_t multV0A=0,multV0C=0,multV0=0,multV0Corr=0,multV0CorrResc=0;
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  if (esdV0) {
    multV0A = esdV0->GetMTotV0A();
    multV0C = esdV0->GetMTotV0C();
  }
  if (fUseMC) {
    multV0A *= fMCV0Scale;
    multV0C *= fMCV0Scale;    
  }
  multV0   = multV0A + multV0C;
  if (!fUseMC) multV0Corr = GetCorrV0(esd, multV0CorrResc);
  else {
    multV0Corr = multV0;
    multV0CorrResc = multV0;
  }
  //
  // correct nSPD2 for Zv dependence
  float nSPD2Corr = GetCorrSPD2(nSPD2,esdvtx[2]);
  //
  float multTPC = fTrackCuts->GetReferenceMultiplicity(esd,kTRUE);
  //
  // registed Ntracklets and ZVertex of the event
  ((TH1*)fHistosCustom->UncheckedAt(kHZVtxNoSel))->Fill(esdvtx[2]);
  ((TH1*)fHistosCustom->UncheckedAt(kHNTrackletsNoSel))->Fill(multESD->GetNumberOfTracklets());      
  ((TH1*)fHistosCustom->UncheckedAt(kHNClSPD1NoSel))->Fill(nSPD1);
  ((TH1*)fHistosCustom->UncheckedAt(kHNClSPD2NoSel))->Fill(nSPD2);
  ((TH1*)fHistosCustom->UncheckedAt(kHV0NoSel))->Fill(multV0Corr);
  ((TH1*)fHistosCustom->UncheckedAt(kHV0CCNoSel))->Fill(multV0CorrResc);
  ((TH1*)fHistosCustom->UncheckedAt(kHMultTPCNoSel))->Fill(multTPC);
  //
  //  ((TH2F*)fHistosCustom->UncheckedAt(kHV0NClSPD2NoSel))->Fill(multV0A+multV0C,multESD->GetNumberOfITSClusters(1));
  //
  //  printf("ESD vertex! %f %f %f, %d contributors\n",esdvtx[0],esdvtx[1],esdvtx[2],vtxESD->GetNContributors());

  if(vtxf[2] < fZVertexMin || vtxf[2] > fZVertexMax) return;
  ((TH2F*)fHistosCustom->UncheckedAt(kHV0NClSPD2NoSel))->Fill(multV0Corr,nSPD2);
  ((TH2F*)fHistosCustom->UncheckedAt(kHV0CCNClSPD2NoSel))->Fill(multV0CorrResc,nSPD2);
  //
  ///  double mltTst = fUseMC ?  multESD->GetNumberOfITSClusters(1) : multV0A+multV0C;  
  //  double mltTst = multESD->GetNumberOfITSClusters(1); //RRR
  float mltTst = -1;
  float mltTst2 = -1;
  if      (fUseCentralityVar == kCentSPD2) mltTst = nSPD2Corr;
  else if (fUseCentralityVar == kCentV0)   mltTst = multV0Corr; // Cvetan's corrected V0
  else if (fUseCentralityVar == kCentV0CR) mltTst = multV0CorrResc; // Cvetan's rescaled V0
  else if (fUseCentralityVar == kCentTrTPC)   mltTst = multTPC; // TPC Mult
  else if (fUseCentralityVar == kCentZDCV0) {
    mltTst = multV0Corr;
    mltTst2 = zdcEnergy;
  }
  else AliFatal(Form("Unknown cenrality parameter %d",fUseCentralityVar));
  //
  fCurrCentBin = GetCentralityBin(mltTst,mltTst2);
  if (fCurrCentBin<0) {
    printf("Reject: %.1f : V0:%.1f V0Cor:%.1f V0CR:%.1f SPD2c:%.1f\n",mltTst, multV0,multV0Corr,multV0CorrResc,nSPD2Corr);
    return;
  }
  //
  ((TH1*)fHistosCustom->UncheckedAt(kHStatCent))->Fill(mltTst);
  ((TH1*)fHistosCustom->UncheckedAt(kHStatCentBin))->Fill(fCurrCentBin);
  printf("Bin %d (mlt=%f) Multiplicity from ESD:\n",fCurrCentBin,mltTst);
  //  multESD->Print();
  //
  AliMCEventHandler* eventHandler = 0;
  fMCEvent = 0;
  fStack = 0;
  //
  if (fUseMC) {
    eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
    if (!eventHandler) { printf("ERROR: Could not retrieve MC event handler\n"); return; }
    fMCEvent = eventHandler->MCEvent();
    if (!fMCEvent) { printf("ERROR: Could not retrieve MC event\n"); return; }
    fStack = fMCEvent->Stack();
    if (!fStack) { printf("Stack not available\n"); return; }
  }
  //
  fRPTree = handRP->GetTreeR("ITS");
  if (!fRPTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  //
  // =============================================================================>>>
  // MC Generator info
  AliGenEventHeader* mcGenH = 0;
  fNPart  = 0;
  fNBColl = 0;
  if (fUseMC) {
    mcGenH = fMCEvent->GenEventHeader();
    if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())) {
      AliGenHijingEventHeader* hHijing = (AliGenHijingEventHeader*)mcGenH;
      fNPart  = (hHijing->ProjectileParticipants()+hHijing->TargetParticipants())/2.;
      fNBColl = hHijing->NN()+hHijing->NNw()+hHijing->NwN()+hHijing->NwNw();
    }
    else if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
      AliGenDPMjetEventHeader* hDpmJet = (AliGenDPMjetEventHeader*)mcGenH;
      fNPart  = (hDpmJet->ProjectileParticipants()+hDpmJet->TargetParticipants())/2.;
      fNBColl = hDpmJet->NN()+hDpmJet->NNw()+hDpmJet->NwN()+hDpmJet->NwNw();
    }
    else {} // unknown generator
  }
  //
  // register Ntracklets and ZVertex of the event
  ((TH2*)fHistosCustom->UncheckedAt(kHZVtx))->Fill(esdvtx[2],fCurrCentBin);
  ((TH2*)fHistosCustom->UncheckedAt(kHNTracklets))->Fill(multESD->GetNumberOfTracklets(),fCurrCentBin);
  //
  if (fUseMC) FillMCPrimaries();
  // fill N clusters
  ((TH2*)fHistosCustom->UncheckedAt(kHNClSPD1))->Fill(nSPD1,fCurrCentBin);
  ((TH2*)fHistosCustom->UncheckedAt(kHNClSPD2))->Fill(nSPD2,fCurrCentBin);
  ((TH2*)fHistosCustom->UncheckedAt(kHV0))->Fill(multV0Corr,fCurrCentBin);
  ((TH2*)fHistosCustom->UncheckedAt(kHV0CC))->Fill(multV0CorrResc,fCurrCentBin);
  ((TH2*)fHistosCustom->UncheckedAt(kHMultTPC))->Fill(multTPC,fCurrCentBin);
  //
  // normal reconstruction
  hstat->Fill(kBinEntries+kEvProcData + kEntriesPerBin*fCurrCentBin);
  //
  if (GetDoNormalReco() || GetDoInjection()) { // for the injection the normal reco should be done
    InitMultReco();
    fMultReco->Run(fRPTree, vtxf);
    printf("Multiplicity Reconstructed:\n");
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    if (mlt) mlt->Print();
    if (GetDoNormalReco()) FillHistos(kData,mlt);
    FillClusterInfo();
    //
  }
  if (!GetDoNormalReco()) FillHistos(kData,multESD); // fill data histos from ESD
  //
  // Injection: it must come right after the normal reco since needs its results
  if (GetDoInjection()) {
    if (!fMultReco) InitMultReco(); // in principle, not needed, the reco is created above
    fMultReco->SetRecType(AliITSMultRecBg::kBgInj);
    fMultReco->Run(fRPTree, vtxf);
    printf("Multiplicity from Injection:\n");
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    if (mlt) mlt->Print();
    hstat->Fill(kBinEntries + kEvProcInj + kEntriesPerBin*fCurrCentBin);
    FillHistos(kBgInj,mlt);
  }
  //
  // Rotation
  if (GetDoRotation()) {
    InitMultReco();
    fMultReco->SetRecType(AliITSMultRecBg::kBgRot);
    fMultReco->SetPhiRotationAngle(fPhiRot);
    fMultReco->Run(fRPTree, vtxf);
    printf("Multiplicity from Rotation:\n");
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    if (mlt) mlt->Print();
    hstat->Fill(kBinEntries + kEvProcRot + kEntriesPerBin*fCurrCentBin);
    FillHistos(kBgRot,mlt);
  }
  //
  if (GetDoMixing()) {
    AliMixEventInputHandler* handToMix = (AliMixEventInputHandler*)handRP->MixingHandler();
    if (!handToMix) { printf("No Mixing handler\n"); return; }
    handToMix->GetEntry();
    if(handToMix->MixedEventNumber()<1) {printf("Mixing: No enough events in pool\n"); return;}
    AliESDInputHandlerRP* handRPMix = (AliESDInputHandlerRP*) handToMix->InputEventHandler(0);

    if (!handRPMix) { printf("No Mixing RP handler\n"); return; }
    fRPTreeMix = handRPMix->GetTreeR("ITS");
    if (!fRPTreeMix) { AliError(" Invalid ITS cluster tree of the 2nd event!\n"); return; }
    //
    AliESDEvent *esdMix = handRPMix->GetEvent();
    const AliESDVertex* vtxESDmix = esdMix->GetVertex();
    ((TH2*)fHistosCustom->UncheckedAt(kHZVtxMixDiff))->Fill(vtxESDmix->GetZ()-esdvtx[2],fCurrCentBin);
    ((TH2*)fHistosCustom->UncheckedAt(kHNTrMixDiff) )->
      Fill(esdMix->GetMultiplicity()->GetNumberOfTracklets() - multESD->GetNumberOfTracklets(),fCurrCentBin);
    //
    InitMultReco();
    fMultReco->SetRecType(AliITSMultRecBg::kBgMix);
    fMultReco->Run(fRPTree, vtxf,fRPTreeMix);
    printf("Multiplicity from Mixing:\n");
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    if (mlt) mlt->Print();
    hstat->Fill(kBinEntries + kEvProcMix + kEntriesPerBin*fCurrCentBin);
    FillHistos(kBgMix,mlt);
    //
  }
  // =============================================================================<<<
  //
  delete fMultReco; 
  fMultReco = 0;
  //
}      

//________________________________________________________________________
Float_t AliTrackletTaskMulti::GetCorrSPD2(Float_t spd2raw,Float_t zv) const
{
  //renormalize N spd2 clusters at given Zv to acceptance at Zv=0
  const double pars[] = {8.10030e-01,-2.80364e-03,-7.19504e-04};
  zv -= pars[0];
  float corr = 1 + zv*(pars[1] + zv*pars[2]);
  return corr>0 ? spd2raw/corr : -1;
}

//________________________________________________________________________
Float_t AliTrackletTaskMulti::GetCorrV0(const AliESDEvent* esd, float &v0CorrResc) const
{
  // correct V0 non-linearity, prepare a version rescaled to SPD2 corr
  const Double_t par0[64] = { 6.71e-02 , 6.86e-02 , 7.06e-02 , 6.32e-02 , 
			      5.91e-02 , 6.07e-02 , 5.78e-02 , 5.73e-02 , 5.91e-02 , 6.22e-02 , 
			      5.90e-02 , 6.11e-02 , 5.55e-02 , 5.29e-02 , 5.19e-02 , 5.56e-02 , 
			      6.25e-02 , 7.03e-02 , 5.64e-02 , 5.81e-02 , 4.57e-02 , 5.30e-02 , 
			      5.13e-02 , 6.43e-02 , 6.27e-02 , 6.48e-02 , 6.07e-02 , 1.01e-01 , 
			      6.68e-02 , 7.16e-02 , 6.36e-02 , 5.95e-02 , 2.52e-02 , 2.82e-02 , 
			      2.56e-02 , 2.86e-02 , 2.82e-02 , 2.10e-02 , 2.13e-02 , 2.32e-02 , 
			      2.75e-02 , 4.34e-02 , 3.78e-02 , 4.52e-02 , 4.11e-02 , 3.89e-02 , 
			      4.10e-02 , 3.73e-02 , 4.51e-02 , 5.07e-02 , 5.42e-02 , 4.74e-02 , 
			      4.33e-02 , 4.44e-02 , 4.64e-02 , 3.01e-02 , 6.38e-02 , 5.26e-02 , 
			      4.99e-02 , 5.26e-02 , 5.47e-02 , 3.84e-02 , 5.00e-02 , 5.20e-02 };
  const Double_t par1[64] = { -6.68e-05 , -7.78e-05 , -6.88e-05 , -5.92e-05 , 
			      -2.43e-05 , -3.54e-05 , -2.91e-05 , -1.99e-05 , -1.40e-05 , -4.01e-05 , 
			      -2.29e-05 , -3.68e-05 , -2.53e-05 , -2.44e-06 , -9.22e-06 , -1.51e-05 , 
			      -2.80e-05 , -2.34e-05 , -1.72e-05 , -1.81e-05 , -1.29e-05 , -2.65e-05 , 
			      -1.61e-05 , -2.86e-05 , -1.74e-05 , -4.23e-05 , -3.41e-05 , -1.05e-04 , 
			      -2.76e-05 , -4.71e-05 , -3.06e-05 , -2.32e-05 , -1.55e-06 , 2.15e-05 , 
			      1.40e-05 , 2.16e-05 , 1.21e-05 , 3.05e-06 , 1.67e-05 , -3.84e-06 , 
			      3.09e-06 , 1.50e-05 , 3.47e-06 , 4.87e-06 , -3.71e-07 , -1.75e-06 , 
			      -1.80e-06 , 9.99e-06 , -6.46e-06 , -4.91e-06 , 1.33e-05 , -2.52e-07 , 
			      -3.85e-06 , 4.94e-06 , -2.48e-07 , -1.20e-05 , 2.07e-06 , 6.12e-06 , 
			      -1.18e-06 , 4.54e-06 , -1.54e-05 , -1.25e-05 , 1.46e-06 , -6.67e-06 };
  const Double_t par2[64] = { 1.29e-08 , 1.51e-08 , 1.43e-08 , 1.11e-08 , 
			      5.04e-09 , 6.99e-09 , 5.58e-09 , 4.15e-09 , 4.00e-09 , 8.22e-09 , 
			      4.97e-09 , 7.66e-09 , 4.91e-09 , 1.10e-09 , 2.64e-09 , 3.64e-09 , 
			      5.76e-09 , 5.46e-09 , 3.38e-09 , 3.47e-09 , 2.43e-09 , 4.13e-09 , 
			      2.80e-09 , 5.80e-09 , 3.86e-09 , 7.46e-09 , 5.98e-09 , 2.58e-08 , 
			      5.50e-09 , 8.72e-09 , 5.23e-09 , 4.37e-09 , 2.33e-09 , -6.01e-10 , 
			      3.99e-11 , -2.02e-10 , 7.67e-10 , 2.03e-09 , 1.17e-10 , 2.56e-09 , 
			      1.16e-09 , -4.75e-10 , 1.28e-09 , 1.23e-09 , 1.62e-09 , 1.61e-09 , 
			      1.93e-09 , 2.97e-10 , 2.21e-09 , 2.16e-09 , 5.22e-10 , 1.03e-09 , 
			      1.56e-09 , 5.00e-10 , 1.01e-09 , 2.93e-09 , 1.05e-09 , 9.96e-11 , 
			      1.21e-09 , 7.45e-10 , 3.07e-09 , 2.31e-09 , 6.70e-10 , 1.89e-09 };
  //
  Float_t multCorr = 0;
  Float_t multCorr2 = 0;
  Float_t multChCorr[64];
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  for(Int_t i = 0; i < 64; ++i) {
    Double_t b = (esdV0->GetMultiplicity(i)*par1[i]-par0[i]);
    Double_t s = (b*b-4.*par2[i]*esdV0->GetMultiplicity(i)*esdV0->GetMultiplicity(i));
    Double_t n;
    if (s<0) {
      printf("FPE %d %.2f %.2f %.2e\n",i,esdV0->GetMultiplicity(i),b,(b*b-4.*par2[i]*esdV0->GetMultiplicity(i)*esdV0->GetMultiplicity(i)));
      n = -b;
    }
    else {
      n = (-b + TMath::Sqrt(s));
    }
    multChCorr[i] = 2.*esdV0->GetMultiplicity(i)/n*par0[i];
    multCorr += multChCorr[i];
    multCorr2 += (multChCorr[i]/par0[i]/64.);
  }
  v0CorrResc =  multCorr2;
  return multCorr;
}


//________________________________________________________________________
void AliTrackletTaskMulti::Terminate(Option_t *) 
{
  Printf("Terminating...");
  TH1* hstat;
  TList *lst = dynamic_cast<TList*>(GetOutputData(1));
  printf("Term: %p %p %p\n",fOutput,lst,fHistosCustom);
  if (lst && (hstat=(TH1*)lst->FindObject("hStat"))) {
    Info("Terminate","Registering used settings");
    // fill used settings
    hstat->Fill(kOneUnit,1.);    
    hstat->Fill(kDPhi,fDPhiWindow);
    hstat->Fill(kDTht,fDThetaWindow);
    hstat->Fill(kNStd,fNStdDev);
    hstat->Fill(kPhiShift,fDPhiShift);
    hstat->Fill(kThtS2,fScaleDTBySin2T);  
    hstat->Fill(kThtCW,fCutOnDThetaX);  
    hstat->Fill(kPhiOvl,fPhiOverlapCut);
    hstat->Fill(kZEtaOvl,fZetaOverlap);
    hstat->Fill(kNoOvl,fRemoveOverlaps);
    //
    hstat->Fill(kPhiRot,fPhiRot);
    hstat->Fill(kInjScl,fInjScale);
    hstat->Fill(kEtaCut,fEtaCut);
    hstat->Fill(kZVMin,fZVertexMin);
    hstat->Fill(kZVMax,fZVertexMax);
    //
    hstat->Fill(kDPiSCut,fDPhiSCut);
    hstat->Fill(kNStdCut,fNStdCut);    
    hstat->Fill(kMCV0Scale, fMCV0Scale);
    //
  }
  //
  //  AliAnalysisTaskSE::Terminate();
}


//_________________________________________________________________________
void AliTrackletTaskMulti::InitMultReco()
{
  // create mult reconstructor
  if (fMultReco) delete fMultReco;
  fMultReco = new AliITSMultRecBg();
  fMultReco->SetCreateClustersCopy(kTRUE);
  fMultReco->SetScaleDThetaBySin2T(fScaleDTBySin2T);
  fMultReco->SetNStdDev(fNStdDev);
  fMultReco->SetPhiWindow( fDPhiWindow );
  fMultReco->SetThetaWindow( fDThetaWindow );
  fMultReco->SetPhiShift( fDPhiShift );
  fMultReco->SetRemoveClustersFromOverlaps(fRemoveOverlaps);
  fMultReco->SetPhiOverlapCut(fPhiOverlapCut);
  fMultReco->SetZetaOverlapCut(fZetaOverlap);
  fMultReco->SetHistOn(kFALSE); 
  fMultReco->SetRecType( AliITSMultRecBg::kData );
}

//_________________________________________________________________________
TObjArray* AliTrackletTaskMulti::BookCustomHistos()
{
  // book custom histos, not related to specific tracklet type
  TObjArray* histos = new TObjArray();
  TH1F* hstat;
  //
  // ------------ job parameters, statistics ------------------------------>>>
  int nbs = kBinEntries + fNCentBins*kEntriesPerBin;
  hstat = new TH1F("hStat","Run statistics",nbs,0.5,nbs+0.5);
  //
  hstat->GetXaxis()->SetBinLabel(kEvTot, "Ev.Tot");
  hstat->GetXaxis()->SetBinLabel(kOneUnit,"ScaleMerge");
  hstat->GetXaxis()->SetBinLabel(kNWorkers,"Workers");
  //
  hstat->GetXaxis()->SetBinLabel(kDPhi,  "#Delta#varphi");
  hstat->GetXaxis()->SetBinLabel(kDTht,  "#Delta#theta");
  hstat->GetXaxis()->SetBinLabel(kNStd,  "N.std");
  hstat->GetXaxis()->SetBinLabel(kPhiShift,"#delta#varphi");
  hstat->GetXaxis()->SetBinLabel(kThtS2,"scale #Delta#theta");
  hstat->GetXaxis()->SetBinLabel(kPhiOvl,"#varpho_{Ovl}");
  hstat->GetXaxis()->SetBinLabel(kZEtaOvl,"#z_{Ovl}");
  hstat->GetXaxis()->SetBinLabel(kNoOvl, "rem.ovl");
  //
  hstat->GetXaxis()->SetBinLabel(kPhiRot,"#varphi_{rot}");
  hstat->GetXaxis()->SetBinLabel(kInjScl,"inj");
  hstat->GetXaxis()->SetBinLabel(kEtaCut,"#eta cut");
  hstat->GetXaxis()->SetBinLabel(kZVMin,"ZV_{min} cut");
  hstat->GetXaxis()->SetBinLabel(kZVMax,"ZV_{max} cut");
  //
  hstat->GetXaxis()->SetBinLabel(kDPiSCut,"#Delta#varphi-#delta_{#phi} cut");
  hstat->GetXaxis()->SetBinLabel(kNStdCut,"#Delta cut");
  //
  hstat->GetXaxis()->SetBinLabel(kMCV0Scale,"MC V0 scale");
  //
  for (int i=0;i<fNCentBins;i++) {
    TString bnt = "b"; bnt+= i;
    int offs = kBinEntries + i*kEntriesPerBin;
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcData, bnt+" Ev.ProcData");
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcInj,  bnt+" Ev.ProcInj");
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcRot,  bnt+" Ev.ProcRot");
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcMix,  bnt+" Ev.ProcMix");
    //
  }
  //
  hstat->Fill(kNWorkers);
  //  
  AddHisto(histos,hstat,kHStat);
  //
  // ------------------------ events per centrality bin ----------------------
  TH1D* hCentAx = new TH1D("EvCentr","Events per centrality",fNCentBins,fkCentBinDef);
  hCentAx->GetXaxis()->SetTitle("Centrality parameter");
  AddHisto(histos,hCentAx,kHStatCent);
  //
  TH1D* hCentBin = new TH1D("EvCentrBin","Events per centrality bin",fNCentBins,-0.5,fNCentBins-0.5);
  hCentBin->GetXaxis()->SetTitle("Centrality Bin");
  AddHisto(histos,hCentBin,kHStatCentBin);
  //  
  // ------------ job parameters, statistics ------------------------------<<<
  //
  double etaMn=-3,etaMx=3;
  double zMn=-30, zMx=30;  
  int nEtaBins = int((etaMx-etaMn)/0.1);
  if (nEtaBins<1) nEtaBins = 1;
  //
  int nZVBins = int(zMx-zMn);
  if (nZVBins<1) nZVBins = 1;
  //
  // Z vertex distribution for events before selection
  TH1F* hzvns = new  TH1F("zvNoSel","Z vertex before selection",nZVBins,zMn,zMx);
  hzvns->GetXaxis()->SetTitle("Zvertex");
  AddHisto(histos,hzvns,kHZVtxNoSel);
  //
  int nbmltSPD2 = 700;
  double maxmltSPD2 = 7000;
  int nbmltV0 = 1000;
  double maxmltV0 = 20000;
  // N tracklets for processed events
  TH1F* hntns = new  TH1F("NtrackletsNoSel","N Tracklets Before Selection",nbmltSPD2,0,maxmltSPD2);
  hntns->GetXaxis()->SetTitle("N tracklets");
  AddHisto(histos,hntns,kHNTrackletsNoSel);
  //
  // N SPD1 clusters
  TH1F* hncl1ns = new  TH1F("NClustersSPD1NoSel","N Clusters on SPD1 Before Selection",nbmltSPD2,0,maxmltSPD2);
  hncl1ns->GetXaxis()->SetTitle("N Clus SPD1");
  AddHisto(histos,hncl1ns,kHNClSPD1NoSel);
  //
  // N SPD2 clusters
  TH1F* hncl2ns = new  TH1F("NClustersSPD2NoSel","N Clusters on SPD2 Before Selection",nbmltSPD2,0,maxmltSPD2);
  hncl2ns->GetXaxis()->SetTitle("N Clus SPD2");
  AddHisto(histos,hncl2ns,kHNClSPD2NoSel);
  //
  // V0 
  TH1F* hnV0ns = new  TH1F("V0NoSel","V0 signal Before Selection",nbmltV0,0,maxmltV0);
  hnV0ns->GetXaxis()->SetTitle("V0 signal");
  AddHisto(histos,hnV0ns,kHV0NoSel);
  //
  // V0 Corr
  TH1F* hnV0CCns = new  TH1F("V0CorrNoSel","V0 Corr signal Before Selection",nbmltSPD2,0,maxmltSPD2); //!!! Same scale as SPD2
  hnV0ns->GetXaxis()->SetTitle("V0 Corr signal");
  AddHisto(histos,hnV0CCns,kHV0CCNoSel);
  //
  // V0 
  //TH2F* hnV0SPD2ns = new  TH2F("V0NDP2NoSel","NSPD2 vs V0 signal Before Selection",2500,0,20000,1400,0,maxmltSPD2);
  TH2F* hnV0SPD2ns = new  TH2F("V0NDP2NoMltSel","NSPD2 vs V0 signal Before Mlt Selection",100,0,maxmltV0,100,0,maxmltSPD2);
  hnV0SPD2ns->GetXaxis()->SetTitle("V0 signal");
  hnV0SPD2ns->GetYaxis()->SetTitle("N Clus SPD2 ");
  AddHisto(histos,hnV0SPD2ns,kHV0NClSPD2NoSel);
  //
  // V0 corr
  //TH2F* hnV0SPD2ns = new  TH2F("V0SPD2NoSel","NSPD2 vs V0 signal Before Selection",2500,0,20000,1400,0,maxmltSPD2);
  TH2F* hnV0CCSPD2ns = new  TH2F("V0CCSPD2NoMltSel","NSPD2 vs V0 Corr signal Before Mlt Selection",100,0,maxmltV0,100,0,maxmltSPD2);
  hnV0CCSPD2ns->GetXaxis()->SetTitle("V0 Corr signal");
  hnV0CCSPD2ns->GetYaxis()->SetTitle("N Clus SPD2 ");
  AddHisto(histos,hnV0CCSPD2ns,kHV0CCNClSPD2NoSel);
  //
  // TPC ref mult before selection
  TH1F* hntpc = new  TH1F("TPCMultNoSel","TPC Multipliplicity Before Selection",300,0,3000);
  hntpc->GetXaxis()->SetTitle("N TPC tracks");
  AddHisto(histos,hntpc,kHMultTPCNoSel);
  //

  TH2F* hzv = new  TH2F("zv","Z vertex after Selection per Cent.Bin",nZVBins,zMn,zMx, fNCentBins, -0.5,fNCentBins-0.5);
  hzv->GetXaxis()->SetTitle("Zvertex");
  hzv->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hzv,kHZVtx);
  //
  // N tracklets for processed events
  TH2F* hnt = new  TH2F("Ntracklets","N Tracklets per Cent.Bin",nbmltSPD2,0,maxmltSPD2,  fNCentBins, -0.5,fNCentBins-0.5);
  hnt->GetXaxis()->SetTitle("N tracklets");
  hnt->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hnt,kHNTracklets);
  //
  // N SPD1 clusters
  TH2F* hncl1 = new  TH2F("NClustersSPD1","N Clusters on SPD1 per Cent.Bin",nbmltSPD2,0,maxmltSPD2,  fNCentBins, -0.5,fNCentBins-0.5);
  hncl1->GetXaxis()->SetTitle("N Clus SPD1");
  hncl1->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hncl1,kHNClSPD1);
  //
  // N SPD2 clusters
  TH2F* hncl2 = new  TH2F("NClustersSPD2","N Clusters on SPD2 per Cent.Bin",nbmltSPD2,0,maxmltSPD2, fNCentBins, -0.5,fNCentBins-0.5);
  hncl2->GetXaxis()->SetTitle("N Clus SPD2");
  hncl2->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hncl2,kHNClSPD2);
  //
  // V0 
  TH2F* hnV0 = new  TH2F("V0","V0 signalper Cent.Bin ",nbmltV0,0,maxmltV0, fNCentBins, -0.5,fNCentBins-0.5);
  hnV0->GetXaxis()->SetTitle("V0 signal");
  hnV0->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hnV0,kHV0);
  //
  // V0 corr 
  TH2F* hnV0CC = new  TH2F("V0Corr","V0 Corr signal per Cent.Bin ",nbmltSPD2,0,maxmltSPD2, fNCentBins, -0.5,fNCentBins-0.5);
  hnV0CC->GetXaxis()->SetTitle("V0 Corr signal");
  hnV0CC->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hnV0CC,kHV0CC);
  //
  // TPC 
  TH2F* hnTPC = new  TH2F("TPCMult","TPC Mult Cent.Bin ",300,0,3000, fNCentBins, -0.5,fNCentBins-0.5);
  hnTPC->GetXaxis()->SetTitle("TPC Mult");
  hnTPC->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hnTPC,kHMultTPC);
  //
  //----------------------------------------------------------------------
  int nEtaBinsS = int(2*fEtaCut/0.1);
  if (nEtaBinsS<1) nEtaBins = 1;
  //
  int nZVBinsS = int(fZVertexMax-fZVertexMin);
  if (nZVBinsS<1) nZVBinsS = 1;

  if (fUseMC) {
    // Z vertex vs Eta distribution for primaries
    char buffn[100],bufft[500];
    for (int ib=0;ib<fNCentBins;ib++) {
      sprintf(buffn,"b%d_zvEtaPrimMC",ib);
      sprintf(bufft,"bin%d Zvertex vs #eta PrimMC",ib);
      TH2F* hzvetap = new  TH2F(buffn,bufft, nEtaBinsS,-fEtaCut,fEtaCut,nZVBinsS,fZVertexMin,fZVertexMax);
      hzvetap->GetXaxis()->SetTitle("#eta");
      hzvetap->GetYaxis()->SetTitle("Zvertex");
      AddHisto(histos,hzvetap,kHZVEtaPrimMC+ib);
    }
    //
    // <n> primaries according to MC generator
    TH1F* hnprimM = new  TH1F("nPrimMean","<N> primaries",fNCentBins, -0.5,fNCentBins-0.5);
    hnprimM->GetXaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hnprimM,kHNPrimMeanMC);
    //
    // <n> primaries per part.pair according to MC generator
    TH1F* hnprim2partM = new  TH1F("nPrim2Part","<N> primaries per part.pair",fNCentBins, -0.5,fNCentBins-0.5);
    hnprim2partM->GetXaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hnprim2partM,kHNPrim2PartMC);
    //
    // <n> primaries per part.pair vs npart.pair according to MC generator
    TH2F* hnprim2partNp = new  TH2F("nPrim2Part_vs_NPart","<N> primaries per part.pair vs N part.pairs",105,0,210,200,0,40);
    hnprim2partNp->GetXaxis()->SetTitle("N.part.pairs");
    hnprim2partNp->GetYaxis()->SetTitle("N.prim/N.part.pairs");
    AddHisto(histos,hnprim2partNp,kHNPrim2PartNpMC);
    //
    // <n> primaries per b.coll vs npart.pair according to MC generator
    TH2F* hnprim2BCollNp = new  TH2F("nPrim2BColl_vs_NPart","<N> primaries per bin.coll vs N part.pairs",105,0,210,200,0,40);
    hnprim2BCollNp->GetXaxis()->SetTitle("N.part.pairs");
    hnprim2BCollNp->GetYaxis()->SetTitle("N.prim/N.bin.coll.");
    AddHisto(histos,hnprim2BCollNp,kHNPrim2BCollNpMC);
    //
    // <n> primaries per bin.coll. according to MC generator
    TH1F* hnprim2BCollM = new  TH1F("nPrim2BColl","<N> primaries per bin.coll",fNCentBins, -0.5,fNCentBins-0.5);
    hnprim2BCollM->GetXaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hnprim2BCollM,kHNPrim2BCollMC);
    //
    // n participants according to MC generator
    TH2F* hnpart = new  TH2F("nPart","N participant pairs",210,0,210,fNCentBins, -0.5,fNCentBins-0.5);
    hnpart->GetXaxis()->SetTitle("N part. pairs");
    hnpart->GetYaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hnpart,kHNPartMC);
    //
    // <n> participants according to MC generator
    TH1F* hnpartM = new  TH1F("nPartMean","<N> participant pairs",fNCentBins, -0.5,fNCentBins-0.5);
    hnpartM->GetXaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hnpartM,kHNPartMeanMC);
    //
    // n bin coll. according to MC generator
    TH2F* hnbcoll = new  TH2F("nBColl","N bin. coll",2000,0,2000,fNCentBins, -0.5,fNCentBins-0.5);
    hnbcoll->GetXaxis()->SetTitle("N bin. coll");
    hnbcoll->GetYaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hnbcoll,kHNBCollMC);
    //
    // <n> bin col according to MC generator
    TH1F* hnbcollM = new  TH1F("nBCollMean","<N> bin.colls",fNCentBins, -0.5,fNCentBins-0.5);
    hnbcollM->GetXaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hnbcollM,kHNBCollMeanMC);
    //    
  }
  //
  if (GetDoMixing()) {
    //
    // Difference in Z vertex for mixed events
    TH2F* hzdiff = new TH2F("MixSPDVertexDiff","SPD #Delta Z Vertex distribution per mult bin ",100,-5,5, fNCentBins, -0.5,fNCentBins-0.5);
    hzdiff->GetXaxis()->SetTitle("#Delta Z Vertex [cm]");
    hzdiff->GetYaxis()->SetTitle(Form("Entries / %1.2f [cm] per mult bin",10./100.));
    AddHisto(histos,hzdiff,kHZVtxMixDiff);
    //
    // Difference in N tracklets for mixed events
    TH2F* hntdiff = new TH2F("MixNTrackletsDiff"," SPD tracklets Diff ",200,-1000,1000, fNCentBins, -0.5,fNCentBins-0.5);
    hntdiff->GetXaxis()->SetTitle("# tracklet diff");
    AddHisto(histos,hntdiff,kHNTrMixDiff);
  }
  // 
  // --------------------------------------------------
  if (fUseMC) {
    int npdg = sizeof(fgkPDGNames)/sizeof(char*);
    TH2F* hpdgP = new TH2F("pdgPrim","primary PDG",npdg,0,npdg,fNCentBins, -0.5,fNCentBins-0.5);
    AddHisto(histos,hpdgP,kHPrimPDG);
    TH2F* hpdgS = new TH2F("pdgSec","secondary PDG",npdg,0,npdg,fNCentBins, -0.5,fNCentBins-0.5);
    AddHisto(histos,hpdgS,kHSecPDG);
    TH2F* hpdgPP = new TH2F("pdgPrimPar","primary parent PDG ",npdg,0,npdg,fNCentBins, -0.5,fNCentBins-0.5);
    AddHisto(histos,hpdgPP,kHPrimParPDG);
    TH2F* hpdgSP = new TH2F("pdgSecPar","secondary parent PDG",npdg,0,npdg,fNCentBins, -0.5,fNCentBins-0.5);
    AddHisto(histos,hpdgSP,kHSecParPDG);
    for (int i=0;i<npdg;i++) {
      hpdgP->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
      hpdgS->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
      hpdgPP->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
      hpdgSP->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
    }
  }
  //
  // -------------------------------------------------
  TH2F* hclinf=0;
  hclinf = new TH2F("cl0InfoUsed","#phi vs Z of used clusters, Lr0",60,-15,15, 80,0,2*TMath::Pi());
  AddHisto(histos,hclinf,kHClUsedInfoL0);
  hclinf = new TH2F("cl1InfoUsed","#phi vs Z of used clusters, Lr1",60,-15,15, 2*80,0,2*TMath::Pi());
  AddHisto(histos,hclinf,kHClUsedInfoL1);
  hclinf = new TH2F("cl0InfoAll","#phi vs Z of all clusters, Lr0",60,-15,15, 80,0,2*TMath::Pi());
  AddHisto(histos,hclinf,kHClAllInfoL0);
  hclinf = new TH2F("cl1InfoAll","#phi vs Z of all clusters, Lr1",60,-15,15, 2*80,0,2*TMath::Pi());
  AddHisto(histos,hclinf,kHClAllInfoL1);
  //
  // -------------------------------------------------
  histos->SetOwner(kFALSE);
  //
  return histos;
}

//_________________________________________________________________________
TObjArray* AliTrackletTaskMulti::BookHistosSet(const char* pref, UInt_t selHistos) 
{
  // book standard set of histos attaching the pref in front of the name/title
  //
  const int kNDPhiBins = 100;
  const int kNDThtBins = 100;
  int nDistBins = int(fNStdDev)*5;
  //
  int nEtaBins = int(2*fEtaCut/0.1);
  if (nEtaBins<1) nEtaBins = 1;
  //
  int nZVBins = int(fZVertexMax-fZVertexMin);
  if (nZVBins<1) nZVBins = 1;
  float dphir = fDPhiWindow*TMath::Sqrt(fNStdDev);
  float dthtr = fDThetaWindow*TMath::Sqrt(fNStdDev);
  //
  TObjArray* histos = new TObjArray();
  TH2F* h2;
  TH1F* h1;
  char buffn[100],bufft[500];
  //
  for (int ib=0;ib<fNCentBins;ib++) {
    //
    int offs = ib*kNStandardH;
    if (selHistos & (0x1<<kHEtaZvCut) ) {
      sprintf(buffn,"b%d_%s_ZvEtaCutT",ib,pref);
      sprintf(bufft,"bin%d (%s) Zv vs Eta with tracklet cut",ib,pref);
      h2 = new TH2F(buffn,bufft,nEtaBins,-fEtaCut,fEtaCut, nZVBins, fZVertexMin,fZVertexMax);
      h2->GetXaxis()->SetTitle("#eta");
      h2->GetYaxis()->SetTitle("Zv");
      AddHisto(histos,h2,offs+kHEtaZvCut);
    }
    //
    if (selHistos & (0x1<<kHDPhiDTheta) ) {
      sprintf(buffn,"b%d_%s_dPhidTheta",ib,pref);
      sprintf(bufft,"bin%d (%s) #Delta#theta vs #Delta#varphi",ib,pref);
      h2 = new TH2F(buffn,bufft,kNDPhiBins,-dphir,dphir,kNDThtBins,-dthtr,dthtr);
      h2->GetXaxis()->SetTitle("#Delta#varphi [rad]");
      h2->GetYaxis()->SetTitle("#Delta#theta [rad]");
      AddHisto(histos,h2,offs+kHDPhiDTheta);
    }
    //
    if (selHistos & (0x1<<kHDPhiSDThetaX) ) {
      sprintf(buffn,"b%d_%s_dPhiSdThetaX",ib,pref);
      sprintf(bufft,"bin%d (%s) #Delta#theta%s vs #Delta#varphi-#delta_{#varphi}",ib,pref,fScaleDTBySin2T ? "/sin^{2}(#theta)":"");
      h2 = new TH2F(buffn,bufft,kNDPhiBins,-dphir,dphir,kNDThtBins,-dthtr,dthtr);
      h2->GetXaxis()->SetTitle("#Delta#varphi-#delta_{#varphi} [rad]");
      sprintf(bufft,"#Delta#theta%s",fScaleDTBySin2T ? "/sin^{2}(#theta)":"");
      h2->GetYaxis()->SetTitle(bufft);
      AddHisto(histos,h2,offs+kHDPhiSDThetaX);
    }
    //
    if (selHistos & (0x1<<kHWDist) ) {
      sprintf(buffn,"b%d_%s_WDist",ib,pref);
      sprintf(bufft,"bin%d #Delta=[(#Delta#varphi-#delta_{#varphi})/#sigma#varphi]^{2}+"
	      "[#Delta#theta%s/#sigma#theta]^{2}",ib,fScaleDTBySin2T ? "*sin^{-2}(#theta)":"");
      h1 = new TH1F(buffn,bufft,nDistBins,0,fNStdDev);
      sprintf(bufft,"#Delta=[(#Delta#varphi-#delta_{#varphi})/#sigma#varphi]^{2}+"
	      "[#Delta#theta%s/#sigma#theta]^{2}",fScaleDTBySin2T ? "*sin^{-2}(#theta)":"");
      h1->GetXaxis()->SetTitle(bufft);
      AddHisto(histos,h1,offs+kHWDist);
    }
    //
  }
  //
  histos->SetOwner(kFALSE);
  return histos;
}

//_________________________________________________________________________
void AliTrackletTaskMulti::AddHisto(TObjArray* histos, TObject* h, Int_t at)
{
  // add single histo to the set
  if (at>=0) histos->AddAtAndExpand(h,at);
  else       histos->Add(h);
  fOutput->Add(h);
}

//_________________________________________________________________________
void AliTrackletTaskMulti::FillHistos(Int_t type, const AliMultiplicity* mlt)
{
  // fill histos of given type
  if (!mlt) return;
  //
  TObjArray* histos = 0;
  if      (type == kData)  histos = fHistosTrData;
  else if (type == kBgInj) histos = fHistosTrInj;
  else if (type == kBgRot) histos = fHistosTrRot;
  else if (type == kBgMix) histos = fHistosTrMix;
  //
  Bool_t fillMC = (type==kData) && fUseMC && fStack;
  //
  //
  //---------------------------------------- CHECK ------------------------------>>>
  TArrayF vtxMC;
  AliGenHijingEventHeader* pyHeader = 0;
  //
  if (fUseMC) {
    pyHeader = (AliGenHijingEventHeader*) fMCEvent->GenEventHeader();//header->GenEventHeader();
    pyHeader->PrimaryVertex(vtxMC);
  }
  //---------------------------------------- CHECK ------------------------------<<<
  //
  if (!histos) return;
  int ntr = mlt->GetNumberOfTracklets();
  for (int itr=ntr;itr--;) {
    //
    //---------------------------------------- CHECK ------------------------------>>>
    /*
    if (fUseMC) {
      Bool_t reject = kFALSE;
      while(1) {
	int lab0 = mlt->GetLabel(itr,0);
	int lab1 = mlt->GetLabel(itr,1);
	if (lab0!=lab1) break;
	if (!fStack->IsPhysicalPrimary(lab0)) break;
	//
	TParticle* part = fStack->Particle(lab0);
	Float_t dz = part->Vz() - vtxMC[2];
	if (TMath::Abs(dz)<1e-6) break;
	reject = kTRUE; 
	break;
      }
      if (reject) continue;
    }
    */
    //---------------------------------------- CHECK ------------------------------<<<
    //
    double theta  = mlt->GetTheta(itr);
    double eta    = -TMath::Log(TMath::Tan(theta/2));
    if (TMath::Abs(eta)>fEtaCut) continue;
    //
    double dtheta = mlt->GetDeltaTheta(itr);
    double dThetaX = dtheta;
    if (fScaleDTBySin2T) {
      double sint   =  TMath::Sin(theta);
      dThetaX /= (sint*sint);
    }
    if (fCutOnDThetaX && TMath::Abs(dThetaX)>fDThetaWindow) continue;
    //    double phi    = mlt->GetPhi(itr);
    double dphi   = mlt->GetDeltaPhi(itr);
    double dist   = mlt->CalcDist(itr);
    //
    FillHistosSet(histos,eta,/*phi,theta,*/dphi,dtheta,dThetaX,dist);
    // special handling for mc info
    if (fillMC && fStack) {
      int lab0 = mlt->GetLabel(itr,0);
      int lab1 = mlt->GetLabel(itr,1);
      int typeMC = 2; // comb.bg.
      if (lab0 == lab1)	typeMC = fStack->IsPhysicalPrimary(lab0) ? 0:1; // prim or sec
      if      (typeMC==0) FillHistosSet(fHistosTrPrim,eta,/*phi,theta,*/dphi,dtheta,dThetaX,dist); // primary
      else if (typeMC==1) FillHistosSet(fHistosTrSec, eta,/*phi,theta,*/dphi,dtheta,dThetaX,dist); // secondary
      else {
	FillHistosSet(fHistosTrComb,eta,/*phi,theta,*/dphi,dtheta,dThetaX,dist); // comb
	// for combinatorals fill also the uncorrelated part
	if (fMultReco) {
	  float *trl = fMultReco->GetTracklet(itr);
	  int clId0 = (int)trl[AliITSMultReconstructor::kClID1];
	  int clId1 = (int)trl[AliITSMultReconstructor::kClID2];
	  float *clLabs0 = fMultReco->GetClusterOfLayer(0,clId0) + AliITSMultReconstructor::kClMC0;
	  float *clLabs1 = fMultReco->GetClusterOfLayer(1,clId1) + AliITSMultReconstructor::kClMC0;
	  if (!HaveCommonParent(clLabs0,clLabs1)) 
	    FillHistosSet(fHistosTrCombU,eta,/*phi,theta,*/dphi,dtheta,dThetaX,dist);
	}
      } // combinatorials
      
      if (dist<fNStdCut) {
	double dphiS  = TMath::Abs(dphi) - fDPhiShift; if (dphi<0) dphiS = -dphiS;
	if (dphiS<fDPhiSCut) FillSpecies(typeMC, lab0);
      }
      if (fCheckReconstructables) CheckReconstructables();
    }
  }
  //
}

//_________________________________________________________________________
void AliTrackletTaskMulti::FillMCPrimaries()
{
  // fill all MC primaries Zv vs Eta
  if (!fStack || !fMCEvent) return;

  //---------------------------------------- CHECK ------------------------------>>>
  TArrayF vtxMC;
  AliGenHijingEventHeader* pyHeader = 0;
  //
  if (fUseMC) {
    pyHeader = (AliGenHijingEventHeader*) fMCEvent->GenEventHeader();//header->GenEventHeader();
    pyHeader->PrimaryVertex(vtxMC);
  }
  //---------------------------------------- CHECK ------------------------------<<<
  //
  int ntr = fStack->GetNtrack();
  TH2* hprimEtaZ = (TH2F*)fHistosCustom->UncheckedAt(kHZVEtaPrimMC+fCurrCentBin);
  int nprim = 0;
  for (int itr=ntr;itr--;) {
    if (!fStack->IsPhysicalPrimary(itr)) continue;
    AliMCParticle *part  = (AliMCParticle*)fMCEvent->GetTrack(itr);
    if (!part->Charge()) continue;
    //
    //---------------------------------------- CHECK ------------------------------>>>
    /*
    Float_t dz = part->Zv() - vtxMC[2];
    if (TMath::Abs(dz)>1e-6) continue; // reject
    */
    //---------------------------------------- CHECK ------------------------------<<<
    //
    Float_t theta = part->Theta();
    if (theta<1e-6 || theta>TMath::Pi()-1e-6) continue;
    Float_t eta = part->Eta();
    if (TMath::Abs(eta)>fEtaCut) continue;
    hprimEtaZ->Fill(eta, fESDVtx[2]);
    nprim++;
  }
  //
  ((TH1*)fHistosCustom->UncheckedAt(kHNPrimMeanMC))->Fill(fCurrCentBin,nprim);
  if (fNPart>0) {
    ((TH1*)fHistosCustom->UncheckedAt(kHNPrim2PartMC))->Fill(fCurrCentBin,nprim/fNPart);
    ((TH2*)fHistosCustom->UncheckedAt(kHNPrim2PartNpMC))->Fill(fNPart,nprim/fNPart);
    ((TH2*)fHistosCustom->UncheckedAt(kHNPartMC))->Fill(fNPart,fCurrCentBin);
    ((TH1*)fHistosCustom->UncheckedAt(kHNPartMeanMC))->Fill(fCurrCentBin,fNPart);
  }
  if (fNBColl>0) {
    ((TH1*)fHistosCustom->UncheckedAt(kHNPrim2BCollMC))->Fill(fCurrCentBin,nprim/fNBColl);
    ((TH2*)fHistosCustom->UncheckedAt(kHNPrim2BCollNpMC))->Fill(fNPart,nprim/fNBColl);
    ((TH2*)fHistosCustom->UncheckedAt(kHNBCollMC))->Fill(fNBColl,fCurrCentBin);
    ((TH1*)fHistosCustom->UncheckedAt(kHNBCollMeanMC))->Fill(fCurrCentBin,fNBColl);
  }
  //
}

//_________________________________________________________________________
 void AliTrackletTaskMulti::FillHistosSet(TObjArray* histos, double eta,
					  //double /*phi*/,double /*theta*/,
					  double dphi,double dtheta,double dThetaX,
					  double dist) 
{
  // fill standard set of histos
  if (dist>fNStdDev) return;
  //
  double dphiS  = TMath::Abs(dphi) - fDPhiShift;
  if (dphi<0) dphiS = -dphiS;
  //
  int offs = fCurrCentBin*kNStandardH;
  //
  if (histos->UncheckedAt(offs+kHDPhiSDThetaX)) 
    ((TH2*)histos->UncheckedAt(offs+kHDPhiSDThetaX))->Fill(dphiS,dThetaX);
  //
  if (histos->UncheckedAt(offs+kHDPhiDTheta)) 
    ((TH2*)histos->UncheckedAt(offs+kHDPhiDTheta))->Fill(dphi,dtheta);
  //
  if (histos->UncheckedAt(kHWDist))
    ((TH2*)histos->UncheckedAt(offs+kHWDist))->Fill(dist);
  //
  if (dist<fNStdCut && dphiS<fDPhiSCut && histos->UncheckedAt(offs+kHEtaZvCut))
    ((TH2*)histos->UncheckedAt(offs+kHEtaZvCut))->Fill(eta,fESDVtx[2]);
  //
}
 
//__________________________________________________________________
void AliTrackletTaskMulti::FillSpecies(Int_t primsec, Int_t id)
{
  // fill PDGcode 
  TH1 *hPart=0,*hParent=0;
  if (primsec==0) {
    hPart   = (TH1*)fHistosCustom->UncheckedAt(kHPrimPDG);
    hParent = (TH1*)fHistosCustom->UncheckedAt(kHPrimParPDG);
  } 
  else if (primsec==1) {
    hPart   = (TH1*)fHistosCustom->UncheckedAt(kHSecPDG);
    hParent = (TH1*)fHistosCustom->UncheckedAt(kHSecParPDG);    
  }
  else return;
  int ntr = fStack->GetNtrack();
  TParticle* part = fStack->Particle(id);
  int pdgCode = TMath::Abs(part->GetPdgCode());
  int pdgBin = GetPdgBin(pdgCode);
  int parID = part->GetFirstMother();
  int pdgCodePar = -1;
  int pdgBinPar = -1;
  while (parID>=0 && parID<ntr) {
    part = fStack->Particle(parID);
    pdgCodePar = TMath::Abs(part->GetPdgCode());
    parID = part->GetFirstMother();
  }
  if (pdgCodePar>0) pdgBinPar = GetPdgBin(pdgCodePar);
  //
  hPart->Fill(pdgBin,fCurrCentBin);
  hParent->Fill(pdgBinPar,fCurrCentBin);
  //
}

//_________________________________________________________________________
Int_t AliTrackletTaskMulti::GetCentralityBin(Float_t multX, Float_t multY)
{
  // calculate centrality bin
  if (fUseCentralityVar<kCent2D) { // 1D bins
    for (int i=0;i<fNCentBins;i++) if (multX>=fkCentBinDef[i] && multX<fkCentBinDef[i+1]) return i;
  }
  else {
    for (int i=fNCentBins;i--;) {
      // X on boundary corresponding to multY
      float xOnLineUp = (multY-fkCentBinDef2DY[i+1])/fkCentBinDef2DS[i+1] + fkCentBinDef[i+1];
      float xOnLineLw = (multY-fkCentBinDef2DY[i]) / fkCentBinDef2DS[i]   + fkCentBinDef[i];
      if (multX>=xOnLineLw && multX<xOnLineUp) return i;
    }
  }
  return -1;
}

//_________________________________________________________________________
Int_t AliTrackletTaskMulti::GetPdgBin(Int_t pdgCode)
{
  // return my pdg bin
  int ncodes = sizeof(fgkPDGCodes)/sizeof(int);
  int pdgBin=0;
  for (pdgBin=0;pdgBin<ncodes;pdgBin++) if (pdgCode==fgkPDGCodes[pdgBin]) break;
  if (pdgBin>=ncodes) {
    if (float(pdgCode)>1e9) pdgBin = ncodes; // nuclei
    else pdgBin = ncodes+1; // unknown
  }
  return pdgBin;
}

//_________________________________________________________________________
Bool_t AliTrackletTaskMulti::HaveCommonParent(const float* clLabs0,const float* clLabs1)
{
  // do 2 clusters have common parrent
  const int kMaxPar = 50;
  static int pars[2][50];
  int npars[2]={0,0};
  const float *labs[2] = {clLabs0,clLabs1};
  int ntr = fStack->GetNtrack();
  for (int il=0;il<2;il++) {
    for (int ilb=0;ilb<3;ilb++) {
      int lbl = (int)labs[il][ilb];
      if (lbl<0 || lbl>=ntr) continue;
      //
      while (npars[il]<kMaxPar-1) {
	pars[il][ npars[il]++ ] = lbl;
	TParticle* part = fStack->Particle(lbl);
	if (!part) break;
	lbl = part->GetFirstMother();
	if (lbl<1 || lbl>=ntr) break;
      }
    }
  }
  // compare array of parents
  for (int i0=npars[0];i0--;) for (int i1=npars[1];i1--;) if (pars[0][i0]==pars[1][i1]) return kTRUE;
  return kFALSE;
}


//_________________________________________________________________________
void AliTrackletTaskMulti::CheckReconstructables()
{
  // fill reconstructable tracklets hitsos
  static TArrayI trInd;
  static TBits   isPrimArr;
  //
  if (!fMultReco || !fMultReco->IsRecoDone()) {AliInfo("To check reconstructables the reco had to be requested"); return;}
  if (!fStack) {AliInfo("MC Stack is not availalble"); return;}
  const double kPtMin = 0.05;
  //
  TClonesArray *clArr[2];
  for (int ilr=0;ilr<2;ilr++) {
    clArr[ilr] = fMultReco->GetClustersOfLayer(ilr);
    if (!clArr[ilr]) {AliInfo("Clusters are not available"); return;}
  }
  //
  int ntr = fStack->GetNtrack();
  if (!ntr) return;
  trInd.Reset();
  if (trInd.GetSize()<ntr) trInd.Set(ntr);
  isPrimArr.ResetAllBits();
  // count track wich may be reconstructable
  //
  int ntrStore=0,ntrStorePrim=0; 
  Int_t *trIndArr = trInd.GetArray();
  for (Int_t it=0; it<ntr; it++) {
    TParticle* part = fStack->Particle(it);
    if (TMath::Abs(part->Eta())>2.2)       continue;
    if (TMath::Abs(part->Pt())<kPtMin)      continue;
    if (fStack->IsPhysicalPrimary(it)) {
      isPrimArr.SetBitNumber(it);
      ntrStorePrim++;
    }
    else { // check if secondary is worth cheking
      TParticlePDG* pdgPart = part->GetPDG();
      if (TMath::Abs(pdgPart->Charge())!=3)  continue;
      if (part->R()>5.)                      continue;
    }
    trIndArr[it] = ++ntrStore;
  }
  //
  AliInfo(Form("Selected %d MC particles (%d prim) out of %d in the stack\n",ntrStore,ntrStorePrim,ntr));
  //
  const int kMaxCl=3;
  AliITSRecPoint **clIndL[2];
  clIndL[0] = new AliITSRecPoint*[kMaxCl*ntrStore]; // max 2 clusters per layer
  clIndL[1] = new AliITSRecPoint*[kMaxCl*ntrStore]; // max 2 clusters per layer
  memset(clIndL[0],0,kMaxCl*ntrStore*sizeof(AliITSRecPoint*));
  memset(clIndL[1],0,kMaxCl*ntrStore*sizeof(AliITSRecPoint*));
  //
  for (int ilr=0;ilr<2;ilr++) {
    TClonesArray *clusters = clArr[ilr];
    int ncl = clusters->GetEntriesFast();
    for (int icl=ncl;icl--;) {
      AliITSRecPoint *cl = (AliITSRecPoint*)clusters->UncheckedAt(icl);
      for (int ilb=3;ilb--;) {
	int lbl = cl->GetLabel(ilb); if (lbl<0 || lbl>=ntr) continue;
	int lblI = trIndArr[lbl];
	if (--lblI<0) continue;    // not kept
	for (int icc=0;icc<kMaxCl;icc++) if (!clIndL[ilr][lblI+icc*ntrStore]) {clIndL[ilr][lblI+ntrStore*icc] = cl; break;} // first empty one
      }
    }
  }
  //
  Float_t clusterLay[2][AliITSMultReconstructor::kClNPar];
  double trComp[6][kMaxCl*kMaxCl];
  int indQual[kMaxCl*kMaxCl];
  //
  for (int itr=ntr;itr--;) {
    int lblI = trIndArr[itr];
    if (--lblI<0) continue; // discarded
    int ntrCand = 0;        // number of tracklet candidates (including overlaps)
    for (int icl0=0;icl0<kMaxCl;icl0++) {
      AliITSRecPoint *cl0 = clIndL[0][lblI+icl0*ntrStore];
      if (!cl0 || !clIndL[1][lblI]) break;
      cl0->GetGlobalXYZ( clusterLay[0] );
      fMultReco->ClusterPos2Angles(clusterLay[0], fESDVtx);
      for (int icl1=0;icl1<kMaxCl;icl1++) {
	AliITSRecPoint *cl1 = clIndL[1][lblI+icl1*ntrStore];
	if (!cl1) break;
	cl1->GetGlobalXYZ( clusterLay[1] );
	fMultReco->ClusterPos2Angles(clusterLay[1], fESDVtx);
	trComp[AliITSMultReconstructor::kTrPhi][ntrCand]    = clusterLay[0][AliITSMultReconstructor::kClPh];
	trComp[AliITSMultReconstructor::kTrTheta][ntrCand]  = clusterLay[0][AliITSMultReconstructor::kClTh];      
	trComp[AliITSMultReconstructor::kTrDTheta][ntrCand] = clusterLay[0][AliITSMultReconstructor::kClTh] - clusterLay[1][AliITSMultReconstructor::kClTh]; 
	trComp[AliITSMultReconstructor::kTrDPhi][ntrCand]   = clusterLay[0][AliITSMultReconstructor::kClPh] - clusterLay[1][AliITSMultReconstructor::kClPh];
	trComp[AliITSMultReconstructor::kTrLab1][ntrCand]   = icl1*10 + icl0;
	double &dphi = trComp[ntrCand][3];
	if (dphi>TMath::Pi()) dphi=2.*TMath::Pi()-dphi;     // take into account boundary condition
	trComp[5][ntrCand] = fMultReco->CalcDist(trComp[AliITSMultReconstructor::kTrDPhi][ntrCand], 
						 trComp[AliITSMultReconstructor::kTrDTheta][ntrCand], 
						 trComp[AliITSMultReconstructor::kTrTheta][ntrCand]);
	ntrCand++;
      }
    }
    if (!ntrCand) continue; // no tracklets
    if (ntrCand>1) TMath::Sort(ntrCand,trComp[5],indQual,kFALSE); else indQual[0] = 0; // sort in weighted distance
    if (fRemoveOverlaps) ntrCand = 1; // select the best
    //
    // disable worst tracklet with shared cluster
    for (int itc=0;itc<ntrCand;itc++) {
      int ind = indQual[itc];
      if (trComp[AliITSMultReconstructor::kTrLab1][ind]<0) continue; // already disabled
      for (int jtc=itc+1;jtc<ntrCand;jtc++) {
	int jnd = indQual[jtc];
	if (trComp[AliITSMultReconstructor::kTrLab1][jnd]<0) continue; // already disabled
	if ( int(trComp[AliITSMultReconstructor::kTrLab1][ind])/10 == int(trComp[AliITSMultReconstructor::kTrLab1][jnd])/10 ||
	     int(trComp[AliITSMultReconstructor::kTrLab1][ind])%10 == int(trComp[AliITSMultReconstructor::kTrLab1][jnd])%10) trComp[AliITSMultReconstructor::kTrLab1][jnd] = -1;
      }
    }
    //
    // store, but forbid cluster reusing
    TObjArray* histos = isPrimArr.TestBitNumber(itr) ? fHistosTrRcblPrim : fHistosTrRcblSec;
    for (int itc=0;itc<ntrCand;itc++) {
      int ind = indQual[itc];
      if (trComp[4][ind]<0) continue; // discarded
      double eta    = -TMath::Log(TMath::Tan(trComp[AliITSMultReconstructor::kTrTheta][ind]/2));
      if (TMath::Abs(eta)>fEtaCut) continue;
      double dThetaX = trComp[AliITSMultReconstructor::kTrTheta][ind];
      if (fScaleDTBySin2T) {
	double sint   =  TMath::Sin(trComp[AliITSMultReconstructor::kTrTheta][ind]);
	dThetaX /= (sint*sint);
      }
      FillHistosSet(histos,eta,
		    //trComp[AliITSMultReconstructor::kTrPhi][ind],trComp[AliITSMultReconstructor::kTrTheta][ind],
		    trComp[AliITSMultReconstructor::kTrDPhi][ind],trComp[AliITSMultReconstructor::kTrDTheta][ind],
		    dThetaX,trComp[5][ind]);
    }
  }
  //
  delete[] clIndL[0];
  delete[] clIndL[1];
}

//_________________________________________________________________________
void AliTrackletTaskMulti::FillClusterInfo()
{
  // fill info on clusters associated to good tracklets
  if (!fMultReco) return;
  int ntr = fMultReco->GetNTracklets();
  int clID[2];
  TH2F *hclU[2] = {(TH2F*)fHistosCustom->UncheckedAt(kHClUsedInfoL0),(TH2F*)fHistosCustom->UncheckedAt(kHClUsedInfoL1)};
  TH2F *hclA[2] = {(TH2F*)fHistosCustom->UncheckedAt(kHClAllInfoL0),(TH2F*)fHistosCustom->UncheckedAt(kHClAllInfoL1)};
  for (int itr=ntr;itr--;) {
    Float_t *trc = fMultReco->GetTracklet(itr);
    if (TMath::Abs(TMath::Abs(trc[AliITSMultReconstructor::kTrDPhi])-fDPhiShift)>fDPhiSCut) continue;
    if (fMultReco->CalcDist(trc[AliITSMultReconstructor::kTrDPhi],
			    trc[AliITSMultReconstructor::kTrDTheta],
			    trc[AliITSMultReconstructor::kTrTheta]) > fNStdCut) continue;
    clID[0] = (int)trc[AliITSMultReconstructor::kClID1];
    clID[1] = (int)trc[AliITSMultReconstructor::kClID2];
    for (int il=0;il<2;il++) {
      Float_t *clinf = fMultReco->GetClusterOfLayer(il,clID[il]);
      hclU[il]->Fill( clinf[AliITSMultReconstructor::kClZ], clinf[AliITSMultReconstructor::kClPh]);
    }
  }
  //
  for (int il=0;il<2;il++) for (int ic=fMultReco->GetNClustersLayer(il);ic--;) {
      Float_t *clinf = fMultReco->GetClusterOfLayer(il,ic);
      hclA[il]->Fill( clinf[AliITSMultReconstructor::kClZ], clinf[AliITSMultReconstructor::kClPh]);
    }
  //
}


//_________________________________________________________________
Bool_t AliTrackletTaskMulti::ZDCTimeTrigger(const AliESDEvent *aEsd) const
{
  // This method implements a selection
  // based on the timing in both sides of zdcN
  // It can be used in order to eliminate
  // parasitic collisions
  Bool_t zdcAccept = kFALSE;

  AliESDZDC *esdZDC = aEsd->GetESDZDC();

  const Float_t refSum = -568.5;
  const Float_t refDelta = -2.1;
  const Float_t sigmaSum = 3.25;
  const Float_t sigmaDelta = 2.25;
  for(Int_t i = 0; i < 4; ++i) {
    if (esdZDC->GetZDCTDCData(10,i) != 0) {
      Float_t tdcC = 0.025*(esdZDC->GetZDCTDCData(10,i)-esdZDC->GetZDCTDCData(14,i)); 
      for(Int_t j = 0; j < 4; ++j) {
	if (esdZDC->GetZDCTDCData(12,j) != 0) {
	  Float_t tdcA = 0.025*(esdZDC->GetZDCTDCData(12,j)-esdZDC->GetZDCTDCData(14,j));
	  if (((tdcC-tdcA-refDelta)*(tdcC-tdcA-refDelta)/(sigmaDelta*sigmaDelta) +
	       (tdcC+tdcA-refSum)*(tdcC+tdcA-refSum)/(sigmaSum*sigmaSum))< 1.0)
	    zdcAccept = kTRUE;
	}
      }
    }
  }
  return zdcAccept;
}

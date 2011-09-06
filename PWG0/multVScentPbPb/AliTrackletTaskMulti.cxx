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
#include "AliCentrality.h"
#include "../ITS/AliITSRecPoint.h"
#include "../ITS/AliITSgeomTGeo.h"
#include "../ITS/AliITSMultReconstructor.h" 

#include "AliLog.h"

#include "AliPhysicsSelection.h"
#include "AliTrackletTaskMulti.h"
#include "AliITSMultRecBg.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliESDtrackCuts.h"

ClassImp(AliTrackletTaskMulti)

// centrality percentile (inverted: 100% - most central)
const Float_t  AliTrackletTaskMulti::fgkCentPerc[] = {0,100};//{0,5,10,20,30};
//const Float_t  AliTrackletTaskMulti::fgkCentPerc[] = {0,5,10,20,30,40};
  //{0,10,20,30,40,50,60,70,80,90,95,101};

const char* AliTrackletTaskMulti::fgCentSelName[] = {"V0M","FMD","TRK","TKL","CL0","CL1","V0MvsFMD","TKLvsV0M","ZENvsZDC"};

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
  fEtaMin(-3.0),
  fEtaMax(3.0),
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
  //
  fNPart(0),
  fNBColl(0),
  fCurrCentBin(-1),
  fNCentBins(0),
  fUseCentralityVar(kCentV0M)
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
  //
  Bool_t needGeom = GetDoNormalReco() || GetDoInjection() || GetDoRotation() || GetDoMixing();
  if (needGeom) {
    AliCDBManager *man = AliCDBManager::Instance();
    if (fUseMC) {
      Bool_t newGeom = kTRUE;
      man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
      if (newGeom) {
	// new geom
	AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",130844,8);
	AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
	if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130844,6,-1)) AliFatal("Failed to misalign geometry");
      }
      else {
	// old geom
	AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",130845,7);
	AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
	if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130845,5,-1)) AliFatal("Failed to misalign geometry");
      }
    }
    else {
      man->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB"); //man->SetRun(137366);
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",137366);
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",137366,-1,-1)) AliFatal("Failed to misalign geometry");
    }
  }
    //
  // Create histograms
  fNCentBins = sizeof(fgkCentPerc)/sizeof(Float_t)-1;
  //---------------------------------------------Standard histos per tracklet type--->>
  UInt_t hPattern = 0xffffffff;
  fHistosTrData                      = BookHistosSet("TrData",hPattern);
  hPattern &= ~(BIT(kHEtaZvSPD1));  // fill single clusters for "data" only
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
  if (fUseCentralityVar<0||fUseCentralityVar>kNCentTypes) {
    printf("Wrong centrality type %d\n",fUseCentralityVar);
    exit(1);
  }
  AliInfo(Form("Centrality type selected: %s\n",fgCentSelName[fUseCentralityVar]));
  PostData(1, fOutput);
  //
}

//________________________________________________________________________
void AliTrackletTaskMulti::UserExec(Option_t *) 
{
  // Main loop
  //
  Bool_t needRecPoints = GetDoNormalReco() || GetDoInjection() || GetDoRotation() || GetDoMixing();
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  fRPTree = fRPTreeMix = 0;
  AliESDInputHandler *handler = (AliESDInputHandler*)anMan->GetInputEventHandler();
  AliESDInputHandlerRP *handRP = 0;
  if (needRecPoints) {
    handRP = (AliESDInputHandlerRP*)handler;
    if (!handRP) { printf("No RP handler\n"); return; }
  }
  AliESDEvent *esd  = handler->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }
  //
  // do we need to initialize the field?
  AliMagF* field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!field && !esd->InitMagneticField()) {printf("Failed to initialize the B field\n");return;}
  //
  //
  TH1* hstat = (TH1*)fHistosCustom->UncheckedAt(kHStat);
  hstat->Fill(kEvTot0); // RS
  const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
  const AliMultiplicity* multESD = esd->GetMultiplicity();
  //
  if (vtxESD->GetNContributors()<1) return;
  TString vtxTyp = vtxESD->GetTitle();
  if (vtxTyp.Contains("vertexer: Z")) {
    if (vtxESD->GetDispersion()>0.04) return;
    if (vtxESD->GetZRes()>0.25) return;
  }
  //
  AliCentrality *centrality = esd->GetCentrality();  
  //
  hstat->Fill(kEvTot); // RS
  //
  Double_t esdvtx[3];
  vtxESD->GetXYZ(esdvtx);
  for (int i=3;i--;) fESDVtx[i] = esdvtx[i];
  //
  float vtxf[3] = {vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ()};
  //
  // registed Ntracklets and ZVertex of the event
  ((TH1*)fHistosCustom->UncheckedAt(kHZVtxNoSel))->Fill(esdvtx[2]);
  //
  if(vtxf[2] < fZVertexMin || vtxf[2] > fZVertexMax) return;
  //
  //  centrality->Dump();
  Float_t multV0=0;
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  if (esdV0) {
    multV0 = esdV0->GetMTotV0A()+esdV0->GetMTotV0C();
    if (fUseMC) multV0 *= fMCV0Scale;
  }
  ((TH1*)fHistosCustom->UncheckedAt(kHV0NoSel))->Fill(multV0);
  //
  float nSPD2 = multESD->GetNumberOfITSClusters(1);
  ((TH1*)fHistosCustom->UncheckedAt(kHNClSPD2NoSel))->Fill(nSPD2);
  //
  //------------------------------------------------------
  AliESDZDC *esdZDC = esd->GetESDZDC();
  float zdcEnergy=0,zemEnergy=0;
  if (esdZDC) {
    zdcEnergy = (esdZDC->GetZDCN1Energy() + esdZDC->GetZDCP1Energy() + esdZDC->GetZDCN2Energy()+ esdZDC->GetZDCP2Energy());
    zemEnergy = (esdZDC->GetZDCEMEnergy(0)+esdZDC->GetZDCEMEnergy(1))/8.; 
  }
  ((TH2*)fHistosCustom->UncheckedAt(kHZDCZEMNoSel))->Fill(zemEnergy,zdcEnergy);
  //
  Float_t centPercentile = centrality->GetCentralityPercentileUnchecked(fgCentSelName[fUseCentralityVar]);

  // temporary >>>>>>>>>>>>>>>>>>>>>>>>
  if (fUseCentralityVar==kCentZEMvsZDC) {
    float zdcEn = zdcEnergy;
    float zemEn = zemEnergy;
    Float_t slope;
    Float_t zdcPercentile;
    if (zemEn > 295.) {
      slope = (zdcEn + 15000.)/(zemEn - 295.);
      slope += 2.23660e+02;
      zdcPercentile = (TMath::ATan(slope) - 1.56664)/8.99571e-05;
      if (zdcPercentile<0) zdcPercentile = 0;
    }
    else zdcPercentile = 100;
    centPercentile = zdcPercentile;
  }
  // temporary >>>>>>>>>>>>>>>>>>>>>>>>
  
  fCurrCentBin = GetCentralityBin(centPercentile);
  //
  //  printf("CentPerc: %f : Bin %d\n",centPercentile, fCurrCentBin);
  if (fCurrCentBin<0) {
    //printf("Reject: %.1f : V0:%.1f V0Cor:%.1f V0CR:%.1f SPD2c:%.1f\n",mltTst, multV0,multV0Corr,multV0CorrResc,nSPD2Corr);
    return;
  }
  //
  ((TH1*)fHistosCustom->UncheckedAt(kHStatCentBin))->Fill(fCurrCentBin);
  ((TH1*)fHistosCustom->UncheckedAt(kHStatCent))->Fill(centPercentile);
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
  if (needRecPoints) {
    fRPTree = handRP->GetTreeR("ITS");
    if (!fRPTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  }
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
  ((TH2*)fHistosCustom->UncheckedAt(kHV0))->Fill(multV0,fCurrCentBin);
  ((TH2*)fHistosCustom->UncheckedAt(kHNClSPD2))->Fill(nSPD2,fCurrCentBin);
  ((TH3*)fHistosCustom->UncheckedAt(kHZDCZEM))->Fill(zemEnergy,zdcEnergy,fCurrCentBin);
  //
  if (fUseMC) FillMCPrimaries();
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
  if (!GetDoNormalReco()) {
    FillHistos(kData,multESD); // fill data histos from ESD
    FillClusterInfoFromMult(multESD, vtxf[2] );
  }
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
    /*
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
    */
    //
  }
  // =============================================================================<<<
  //
  if (fMultReco) delete fMultReco; 
  fMultReco = 0;
  //
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
    hstat->Fill(kCentVar,fUseCentralityVar);
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
    hstat->Fill(kEtaMin,fEtaMin);
    hstat->Fill(kEtaMax,fEtaMax);
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
  hstat->GetXaxis()->SetBinLabel(kEvTot, "Ev.Tot0");
  hstat->GetXaxis()->SetBinLabel(kEvTot, "Ev.Tot");
  hstat->GetXaxis()->SetBinLabel(kOneUnit,"ScaleMerge");
  hstat->GetXaxis()->SetBinLabel(kNWorkers,"Workers");
  //
  hstat->GetXaxis()->SetBinLabel(kCentVar,  fgCentSelName[fUseCentralityVar]);
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
  hstat->GetXaxis()->SetBinLabel(kEtaMin,"#eta_{min}");
  hstat->GetXaxis()->SetBinLabel(kEtaMax,"#eta_{max}");
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
  TH1D* hCentAx = new TH1D("EvCentr","Events per centrality",fNCentBins,fgkCentPerc);
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
  // V0 for events before selection
  int nbmltV0 = 250;
  double maxmltV0 = 25000;
  //
  TH1F* hnV0ns = new  TH1F("V0NoSel","V0 signal Before Cent. Selection",nbmltV0,0,maxmltV0);
  hnV0ns->GetXaxis()->SetTitle("V0 signal");
  AddHisto(histos,hnV0ns,kHV0NoSel);
  //
  // N SPD2 clusters
  int nbmltSPD2 = 175;
  double maxmltSPD2 = 7000;
  TH1F* hncl2ns = new  TH1F("NClustersSPD2NoSel","N Clusters on SPD2 Before Cent Selection",nbmltSPD2,0,maxmltSPD2);
  hncl2ns->GetXaxis()->SetTitle("N Clus SPD2");
  AddHisto(histos,hncl2ns,kHNClSPD2NoSel);
  //
  int nbzdc=50;
  double maxZDC=6000., maxZEM=2500.;
  TH2F* hzdczemns = new  TH2F("ZDCZEMNoSel","ZDC vs ZEM Before Cent Selection",
			      nbzdc,0,maxZEM,nbzdc,0,maxZDC);
  hzdczemns->GetXaxis()->SetTitle("ZEM");
  hzdczemns->GetXaxis()->SetTitle("ZDC");
  AddHisto(histos,hzdczemns,kHZDCZEMNoSel);
  //
  TH2F* hzv = new  TH2F("zv","Z vertex after Selection per Cent.Bin",nZVBins,zMn,zMx, fNCentBins, -0.5,fNCentBins-0.5);
  hzv->GetXaxis()->SetTitle("Zvertex");
  hzv->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hzv,kHZVtx);
  //
  // V0 
  TH2F* hnV0 = new  TH2F("V0","V0 signal per Cent.Bin ",nbmltV0,0,maxmltV0, fNCentBins, -0.5,fNCentBins-0.5);
  hnV0->GetXaxis()->SetTitle("V0 signal");
  hnV0->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hnV0,kHV0);
  //
  // N SPD2 clusters
  TH2F* hncl2 = new  TH2F("NClustersSPD2","N Clusters on SPD2 per Cent.Bin",nbmltSPD2,0,maxmltSPD2, fNCentBins, -0.5,fNCentBins-0.5);
  hncl2->GetXaxis()->SetTitle("N Clus SPD2");
  hncl2->GetYaxis()->SetTitle("Cent.Bin ID");
  AddHisto(histos,hncl2,kHNClSPD2);
  //
  // ZDCZEM
  TH3F* hzdczem = new TH3F("ZDCZEM","ZDC vs ZEM per Cent.Bin",nbzdc,0,maxZEM,nbzdc,0,maxZDC, fNCentBins, -0.5,fNCentBins-0.5);
  hzdczem->GetXaxis()->SetTitle("ZEM");
  hzdczem->GetYaxis()->SetTitle("ZDC");
  AddHisto(histos,hzdczem,kHZDCZEM);
  //
  //----------------------------------------------------------------------
  int nEtaBinsS = int((fEtaMax-fEtaMin)/0.1);
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
      TH2F* hzvetap = new  TH2F(buffn,bufft, nEtaBinsS,fEtaMin,fEtaMax,nZVBinsS,fZVertexMin,fZVertexMax);
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
  hclinf = new TH2F("cl0InfoUsed","#phi vs Z of used clusters, Lr0",64,-16,16, 80,0,2*TMath::Pi());
  AddHisto(histos,hclinf,kHClUsedInfoL0);
  hclinf = new TH2F("cl1InfoUsed","#phi vs Z of used clusters, Lr1",64,-16,16, 2*80,0,2*TMath::Pi());
  AddHisto(histos,hclinf,kHClUsedInfoL1);
  hclinf = new TH2F("cl0InfoAll","#phi vs Z of all clusters, Lr0",64,-16,16, 80,0,2*TMath::Pi());
  AddHisto(histos,hclinf,kHClAllInfoL0);
  hclinf = new TH2F("cl1InfoAll","#phi vs Z of all clusters, Lr1",64,-16,16, 2*80,0,2*TMath::Pi());
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
  int nEtaBins = int((fEtaMax-fEtaMin)/0.1);
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
      h2 = new TH2F(buffn,bufft,nEtaBins,fEtaMin,fEtaMax, nZVBins, fZVertexMin,fZVertexMax);
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
    if (selHistos & (0x1<<kHEtaZvSPD1) ) {
      sprintf(buffn,"b%d_%s_ZvEtaSPD1",ib,pref);
      sprintf(bufft,"bin%d (%s) Zv vs Eta SPD1 clusters",ib,pref);
      h2 = new TH2F(buffn,bufft,nEtaBins,fEtaMin,fEtaMax, nZVBins, fZVertexMin,fZVertexMax);
      h2->GetXaxis()->SetTitle("#eta");
      h2->GetYaxis()->SetTitle("Zv");
      AddHisto(histos,h2,offs+kHEtaZvSPD1);
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
    if (eta<fEtaMin || eta>fEtaMax) continue;
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
  //-------------------------------------------------------------TMP RS - singles ------->>>
  int offsH = fCurrCentBin*kNStandardH;
  TH2* hSingles = (TH2*)histos->UncheckedAt(offsH+kHEtaZvSPD1);
  if (hSingles) {
    int nclS = mlt->GetNumberOfSingleClusters();
    double *thtS = mlt->GetThetaSingle();
    for (int ics=nclS;ics--;) {
      double etaS = -TMath::Log(TMath::Tan(thtS[ics]/2));
      if (etaS<fEtaMin || etaS>fEtaMax) continue;
      hSingles->Fill(etaS,fESDVtx[2]);
    }
  }
  //-------------------------------------------------------------TMP RS - singles -------<<<
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
    if (eta<fEtaMin || eta>fEtaMax) continue;
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
Int_t AliTrackletTaskMulti::GetCentralityBin(Float_t perc) const
{
  // calculate centrality bin
  for (int i=0;i<fNCentBins;i++) if (perc>=fgkCentPerc[i] && perc<=fgkCentPerc[i+1]) return i;
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
  //  printf("\nHaveCommonParent \n");
  for (int il=0;il<2;il++) {
    
    for (int ilb=0;ilb<3;ilb++) {
      int lbl = (int)labs[il][ilb];
      if (lbl<2 || lbl>=ntr) continue;
      //
      while (npars[il]<kMaxPar-1) {
	TParticle* part = fStack->Particle(lbl);
	if (!part) break;
	int code = TMath::Abs(part->GetPdgCode());
	int q = (int)TMath::Abs(part->GetPDG()->Charge());
	if (code==21 || code<10 || q==1 || q==2 || q==4 ) break;
	//printf("%d/%d/%d: %4d (%d)%s|",il,ilb,npars[il],lbl,part->GetStatusCode(),part->GetName());
	pars[il][ npars[il]++ ] = lbl;
	lbl = part->GetFirstMother();
	if (lbl<1 || lbl>=ntr) break;
      }
      //      printf("\n");
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
      if (eta<fEtaMin || eta>fEtaMax) continue;
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

//_________________________________________________________________________
void AliTrackletTaskMulti::FillClusterInfoFromMult(const AliMultiplicity* mlt, double zVertex)
{
  // fill info on clusters taking them from Multiplicity object
  const double kRSPD2 = 3.9;
  TH2F *hclU = (TH2F*)fHistosCustom->UncheckedAt(kHClUsedInfoL0);
  TH2F *hclA = (TH2F*)fHistosCustom->UncheckedAt(kHClAllInfoL0);
  int ntr = mlt->GetNumberOfTracklets();
  for (int itr=ntr;itr--;) {
    Bool_t goodTracklet = kTRUE;
    if (TMath::Abs( mlt->GetDeltaPhi(itr)-fDPhiShift)>fDPhiSCut) goodTracklet = kFALSE;
    if (mlt->CalcDist(itr) > fNStdCut) goodTracklet = kFALSE;
    double phi   = mlt->GetPhi(itr);
    double z     = kRSPD2/TMath::Tan(mlt->GetTheta(itr)) + zVertex;
    if (goodTracklet) hclU->Fill(z,phi);
    hclA->Fill(z,phi);
  }
  //
  int ncl = mlt->GetNumberOfSingleClusters();
  for (int icl=ncl;icl--;) {
    double phi   = mlt->GetPhiSingle(icl);
    double z     = kRSPD2/TMath::Tan(mlt->GetThetaSingle(icl)) + zVertex;
    hclA->Fill(z,phi);
  }
  //
}


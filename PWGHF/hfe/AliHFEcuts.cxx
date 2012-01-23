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
**************************************************************************/
//
// Cut menagement class implemented by the
// ALICE Heavy Flavour Electron Group
// Interface to the correction framework
// Provides a set of standard cuts 
// 
// Authors:
//   Raphaelle Bailhache <R.Bailhache@gsi.de>
//   Markus Fasel <M.Fasel@gsi.de>
//   Markus Heide <mheide@uni-muenster.de>
//   Matus Kalisky <m.kalisky@uni-muenster.de>
//
// Overview over the 18 steps in the correction Framework
// 0. Generated Electrons
// 1. Signal Electrons
// 2. Electron in Acceptance
// ------------------------------------------------------------
// 3. Rec without cuts (MC information)
// 4. Rec Kine ITS/TPC (MC Information)
// 5. Rec Primary (MC Information)
// 6. HFE ITS (MC Information)
// 7. HFE TRD (MC Information)
// 8. PID (MC Information) 
// ............................................................
// 9. Rec without cuts(MC Information for tracks which are already registered)
// 10. Rec Kine ITS/TPC (MC Information for tracks which are already registered)
// 11. RecPrimary (MC Information for tracks which are already registered)
// 12. HFE ITS (MC Information for tracks which are already registered)
// 13. HFE TPC (MC Information for tracks which are already registered)
// 14. PID (MC Information for tracks which are already registered)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 15. Rec without cuts
// 16. Rec Kine ITS/TPC 
// 17. Rec Primary
// 18. HFE ITS
// 19. HFE TRD
// 20. PID
//
#include <TClass.h>
#include <TList.h>
#include <TObjArray.h>
#include <TString.h>

#include "AliCFAcceptanceCuts.h"
#include "AliCFCutBase.h"
#include "AliCFEventGenCuts.h"
#include "AliCFEventRecCuts.h"
#include "AliCFManager.h"
#include "AliCFParticleGenCuts.h"
#include "AliCFTrackIsPrimaryCuts.h"
#include "AliCFTrackKineCuts.h"
#include "AliCFTrackQualityCuts.h"
#include "AliESDtrack.h"
#include "AliHFEextraEventCuts.h"

#include "AliHFEcuts.h"

ClassImp(AliHFEcuts)

const Char_t *AliHFEcuts::fgkMCCutName[AliHFEcuts::kNcutStepsMCTrack] = {
  "MCGenerated",
  "MCGeneratedZOutNoPileUpCentralityFine",
  "MCGeneratedEventCut",
  "MCInAcceptance"
};

const Char_t * AliHFEcuts::fgkRecoCutName[AliHFEcuts::kNcutStepsRecTrack] = {
  "NoCuts",
  "RecKineITSTPC",
  "Primary",
  "HFEITS",
  "HFETOF",
  "HFETPC",
  "HFETRD"
};

const Char_t * AliHFEcuts::fgkDECutName[AliHFEcuts::kNcutStepsDETrack] = {
  "HFEDCA"
};

const Char_t * AliHFEcuts::fgkSecvtxCutName[AliHFEcuts::kNcutStepsSecvtxTrack] = {
  "HFESecvtx"
};

const Char_t * AliHFEcuts::fgkEventCutName[AliHFEcuts::kNcutStepsEvent] = {
  "EventStepGenerated",
  "EventStepRecNoCut",
  "EventStepRecNoPileUp",
  "EventStepRecCentralityOK",
  "EventStepZRange",
  "EventStepReconstructed"
};

const Char_t * AliHFEcuts::fgkUndefined = "Undefined";

//__________________________________________________________________
AliHFEcuts::AliHFEcuts():
  TNamed(),
  fRequirements(0),
  fTPCclusterDef(AliHFEextraCuts::kFound),
  fTPCratioDef(AliHFEextraCuts::kFoundOverFindable),
  fMinClustersTPC(0),
  fMinClustersITS(0),
  fMinTrackletsTRD(0),
  fCutITSPixel(0),
  fCheckITSLayerStatus(kTRUE),
  fMaxChi2clusterITS(-1.),
  fMaxChi2clusterTPC(0.),
  fMinClusterRatioTPC(0.),
  fVertexRangeZ(20.),
  fTOFPIDStep(kFALSE),
  fTOFMISMATCHStep(kFALSE),
  fTPCPIDCLEANUPStep(kFALSE),
  fUseMixedVertex(kTRUE),
  fIsIPSigmacut(kFALSE),
  fFractionOfSharedTPCClusters(-1.0),
  fMaxImpactParameterRpar(kFALSE),
  fHistQA(0x0),
  fCutList(0x0),
  fDebugLevel(0)
{
  //
  // Dummy Constructor
  //
  memset(fProdVtx, 0, sizeof(Double_t) * 4);
  memset(fDCAtoVtx, 0, sizeof(Double_t) * 2);
  memset(fPtRange, 0, sizeof(Double_t) * 2);
  memset(fIPCutParams, 0, sizeof(Float_t) * 4);
  memset(fSigmaToVtx, 0, sizeof(Double_t) * 3);

}

//__________________________________________________________________
AliHFEcuts::AliHFEcuts(const Char_t *name, const Char_t *title):
  TNamed(name, title),
  fRequirements(0),
  fTPCclusterDef(AliHFEextraCuts::kFound),
  fTPCratioDef(AliHFEextraCuts::kFoundOverFindable),
  fMinClustersTPC(0),
  fMinClustersITS(0),
  fMinTrackletsTRD(0),
  fCutITSPixel(0),
  fCheckITSLayerStatus(kTRUE),
  fMaxChi2clusterITS(-1.),
  fMaxChi2clusterTPC(0.),
  fMinClusterRatioTPC(0.),
  fVertexRangeZ(20.),
  fTOFPIDStep(kFALSE),
  fTOFMISMATCHStep(kFALSE),
  fTPCPIDCLEANUPStep(kFALSE),
  fUseMixedVertex(kTRUE),
  fIsIPSigmacut(kFALSE),
  fFractionOfSharedTPCClusters(-1.0),
  fMaxImpactParameterRpar(kFALSE),
  fHistQA(0x0),
  fCutList(0x0),
  fDebugLevel(0)
{
  //
  // Default Constructor
  //
  memset(fProdVtx, 0, sizeof(Double_t) * 4);
  memset(fDCAtoVtx, 0, sizeof(Double_t) * 2);
  memset(fPtRange, 0, sizeof(Double_t) * 2);
  memset(fIPCutParams, 0, sizeof(Float_t) * 4);
  memset(fSigmaToVtx, 0, sizeof(Double_t) * 3);
}

//__________________________________________________________________
AliHFEcuts::AliHFEcuts(const AliHFEcuts &c):
  TNamed(c),
  fRequirements(c.fRequirements),
  fTPCclusterDef(c.fTPCclusterDef),
  fTPCratioDef(c.fTPCratioDef),
  fMinClustersTPC(0),
  fMinClustersITS(0),
  fMinTrackletsTRD(0),
  fCutITSPixel(0),
  fCheckITSLayerStatus(0),
  fMaxChi2clusterITS(-1.),
  fMaxChi2clusterTPC(0),
  fMinClusterRatioTPC(0),
  fVertexRangeZ(20.),
  fTOFPIDStep(kFALSE),
  fTOFMISMATCHStep(kFALSE),
  fTPCPIDCLEANUPStep(kFALSE),
  fUseMixedVertex(kTRUE),
  fIsIPSigmacut(kFALSE),
  fFractionOfSharedTPCClusters(-1.0),
  fMaxImpactParameterRpar(kFALSE),
  fHistQA(0x0),
  fCutList(0x0),
  fDebugLevel(0)
{
  //
  // Copy Constructor
  //
  c.Copy(*this);
}

//__________________________________________________________________
AliHFEcuts &AliHFEcuts::operator=(const AliHFEcuts &c){
  //
  // Make assignment
  //
  if(&c != this) c.Copy(*this);
  return *this;
}

//__________________________________________________________________
void AliHFEcuts::Copy(TObject &c) const {
  //
  // Performing copy
  //
  AliHFEcuts &target = dynamic_cast<AliHFEcuts &>(c);

  target.fRequirements = fRequirements;
  target.fTPCclusterDef = fTPCclusterDef;
  target.fTPCratioDef = fTPCratioDef;
  target.fMinClustersTPC = fMinClustersTPC;
  target.fMinClustersITS = fMinClustersITS;
  target.fMinTrackletsTRD = fMinTrackletsTRD;
  target.fCutITSPixel = fCutITSPixel;
  target.fCheckITSLayerStatus = fCheckITSLayerStatus;
  target.fMaxChi2clusterITS = fMaxChi2clusterITS;
  target.fMaxChi2clusterTPC = fMaxChi2clusterTPC;
  target.fMinClusterRatioTPC = fMinClusterRatioTPC;
  target.fVertexRangeZ = fVertexRangeZ;
  target.fTOFPIDStep = fTOFPIDStep;
  target.fTOFMISMATCHStep = fTOFMISMATCHStep;
  target.fTPCPIDCLEANUPStep = fTPCPIDCLEANUPStep;
  target.fUseMixedVertex = fUseMixedVertex;
  target.fIsIPSigmacut = fIsIPSigmacut;
  target.fFractionOfSharedTPCClusters = fFractionOfSharedTPCClusters;
  target.fMaxImpactParameterRpar = fMaxImpactParameterRpar;
  target.fDebugLevel = 0;

  memcpy(target.fProdVtx, fProdVtx, sizeof(Double_t) * 4);
  memcpy(target.fDCAtoVtx, fDCAtoVtx, sizeof(Double_t) * 2);
  memcpy(target.fPtRange, fPtRange, sizeof(Double_t) *2);
  memcpy(target.fIPCutParams, fIPCutParams, sizeof(Float_t) * 4);
  memcpy(target.fSigmaToVtx, fSigmaToVtx, sizeof(Double_t) * 3);

  // Copy cut List
  if(target.fCutList){
    target.fCutList->Clear();
    delete target.fCutList;
  }
  if(fCutList){
    target.fCutList = dynamic_cast<TObjArray *>(fCutList->Clone());
    if(target.fCutList) target.fCutList->SetOwner(); // Coverity
  }
  if(target.fHistQA){
    delete target.fHistQA;
  }
  if(fHistQA){
    // If the QA list was already produced, then we create it new, loop over the cuts and connect all the histos with this list
    target.fHistQA = new TList;
    target.fHistQA->SetName(Form("%s_CutQAhistograms", GetName()));
    fHistQA->SetOwner(kTRUE);
    TIter cit(target.fCutList);
    TObjArray *clist = NULL;
    AliCFCutBase *co = NULL;
    while((clist = dynamic_cast<TObjArray *>(cit()))){
      TIter cit1(clist);
      while((co = dynamic_cast<AliCFCutBase *>(cit1()))) co->SetQAOn(target.fHistQA);
    }
  }
}

//__________________________________________________________________
AliHFEcuts::~AliHFEcuts(){
  //
  // Destruktor
  //
  if(fCutList){
    fCutList->Delete();
    delete fCutList;
  }
  fCutList = 0x0;
  if(fHistQA) delete fHistQA;
}

//__________________________________________________________________
Long64_t AliHFEcuts::Merge(const TCollection *list){
  //
  // Merge function doing nothing, just writing the object to file
  //
  if(!list) return 0;
  if(list->IsEmpty()) return 1;
  Long64_t counts = 0;
  // just count the number of objects to merge
  TIter iter(list);
  while(iter()){ counts++; }
  return counts+1;
}

//__________________________________________________________________
void AliHFEcuts::Initialize(AliCFManager *cfm){
  //
  // Initializes the cut objects from the correction framework
  // Publishes the cuts to the correction framework manager
  //
  AliDebug(2, "Called");
  const Int_t kMCOffset = kNcutStepsMCTrack;
  const Int_t kRecOffset = kNcutStepsRecTrack;
  if(fCutList)
    fCutList->Delete();
  else{
    fCutList = new TObjArray;
    fCutList->SetOwner();
  }
  if(IsQAOn()){
    fHistQA = new TList;
    fHistQA->SetName(Form("%s_CutQAhistograms", GetName()));
    fHistQA->SetOwner(kTRUE);
  }
 
  // Call all the setters for the cuts
  SetParticleGenCutList();
  SetAcceptanceCutList();
  SetRecKineITSTPCCutList();
  SetRecPrimaryCutList();
  SetHFElectronITSCuts();
  SetHFElectronTOFCuts();
  SetHFElectronTPCCuts();
  SetHFElectronTRDCuts();
  SetHFElectronDcaCuts();

  // Publish to the cuts which analysis type they are (ESD Analysis by default)
  if(IsAOD()){
    AliInfo("Setting AOD Analysis");
    TObjArray *genCuts = dynamic_cast<TObjArray *>(fCutList->FindObject("fPartGenCuts"));
    if(genCuts){
      AliCFParticleGenCuts *myGenCut = dynamic_cast<AliCFParticleGenCuts *>(genCuts->FindObject("fCutsGenMC"));
      if(myGenCut) myGenCut->SetAODMC(kTRUE);
    }
  }

  // Connect the event cuts
  SetEventCutList(kEventStepGenerated);
  SetEventCutList(kEventStepReconstructed);
  cfm->SetEventCutsList(kEventStepGenerated, dynamic_cast<TObjArray *>(fCutList->FindObject("fEvGenCuts")));
  cfm->SetEventCutsList(kEventStepReconstructed, dynamic_cast<TObjArray *>(fCutList->FindObject("fEvRecCuts")));
  
  // Connect the particle cuts
  // 1st MC
  cfm->SetParticleCutsList(kStepMCGenerated, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartGenCuts")));
  cfm->SetParticleCutsList(kStepMCInAcceptance, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartAccCuts")));
  // 2nd Reco
  cfm->SetParticleCutsList(kStepRecKineITSTPC + kMCOffset, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartRecKineITSTPCCuts")));
  cfm->SetParticleCutsList(kStepRecPrim + kMCOffset, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartPrimCuts")));
  cfm->SetParticleCutsList(kStepHFEcutsITS + kMCOffset, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartHFECutsITS")));
  cfm->SetParticleCutsList(kStepHFEcutsTOF+ kMCOffset, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartHFECutsTOF")));
  cfm->SetParticleCutsList(kStepHFEcutsTPC+ kMCOffset, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartHFECutsTPC")));
  cfm->SetParticleCutsList(kStepHFEcutsTRD + kMCOffset, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartHFECutsTRD")));
  cfm->SetParticleCutsList(kStepHFEcutsDca + kRecOffset + kMCOffset, dynamic_cast<TObjArray *>(fCutList->FindObject("fPartHFECutsDca")));

}

//__________________________________________________________________
void AliHFEcuts::Initialize(){
  // Call all the setters for the cuts
  AliDebug(2, "Called\n");
   if(fCutList)
    fCutList->Delete();
  else{
    fCutList = new TObjArray;
    fCutList->SetOwner();
  }
  if(IsQAOn()){
    fHistQA = new TList;
    fHistQA->SetName(Form("%s_CutQAhistograms", GetName()));
    fHistQA->SetOwner(kFALSE);
  }
  SetParticleGenCutList();
  SetAcceptanceCutList();
  SetRecKineITSTPCCutList();
  SetRecPrimaryCutList();
  SetHFElectronITSCuts();
  SetHFElectronTOFCuts();
  SetHFElectronTPCCuts();
  SetHFElectronTRDCuts();
  SetHFElectronDcaCuts();

  // Connect the event cuts
  SetEventCutList(kEventStepGenerated);
  SetEventCutList(kEventStepReconstructed);


}

//__________________________________________________________________
void AliHFEcuts::SetEventCutList(Int_t istep){
  // 
  // Cuts for Event Selection
  //
  AliDebug(2, "Called\n");
  TObjArray *arr = new TObjArray;
  if(istep == kEventStepGenerated){
    AliCFEventGenCuts *evGenCuts = new AliCFEventGenCuts((Char_t *)"fCutsEvGen", (Char_t *)"Event Generated cuts");
    //evGenCuts->SetNTracksCut(1);
    evGenCuts->SetRequireVtxCuts(kTRUE);
    //evGenCuts->SetVertexXCut(-1, 1);
    //evGenCuts->SetVertexYCut(-1, 1);
    evGenCuts->SetVertexZCut(-fVertexRangeZ, fVertexRangeZ);
    if(IsQAOn()) evGenCuts->SetQAOn(fHistQA);

    arr->SetName("fEvGenCuts");
    arr->AddLast(evGenCuts);
  } else {

    if(!fUseMixedVertex) {
      AliCFEventRecCuts *evRecCuts = new AliCFEventRecCuts((Char_t *)"fCutsEvRec", (Char_t *)"Event Reconstructed cuts");
      //evRecCuts->SetNTracksCut(1);
      evRecCuts->SetRequireVtxCuts(kTRUE);
      //evRecCuts->SetVertexXCut(-1, 1);
      //evRecCuts->SetVertexYCut(-1, 1);
      //evRecCuts->SetVertexZCut(-30, 30);
      evRecCuts->SetVertexZCut(-fVertexRangeZ, fVertexRangeZ);
      evRecCuts->SetVertexNContributors(1,(Int_t)1.e9);
      if(IsQAOn()) evRecCuts->SetQAOn(fHistQA);
      arr->SetName("fEvRecCuts");
      arr->AddLast(evRecCuts);
    } else {
      AliHFEextraEventCuts *evRecCuts = new AliHFEextraEventCuts((Char_t *)"fCutsEvRec", (Char_t *)"Event Reconstructed cuts");
      evRecCuts->SetRequireVtxCuts(kTRUE);
      evRecCuts->SetUseMixedVertex();
      evRecCuts->SetVertexZCut(-fVertexRangeZ, fVertexRangeZ);
      //evRecCuts->SetVertexNContributors(1,(Int_t)1.e9);
      if(IsQAOn()) evRecCuts->SetQAOn(fHistQA);
      arr->SetName("fEvRecCuts");
      arr->AddLast(evRecCuts);
    }
  }
  fCutList->AddLast(arr);
}

//__________________________________________________________________
void AliHFEcuts::SetParticleGenCutList(){
  //
  // Initialize Particle Cuts for Monte Carlo Tracks
  // Production Vertex Radius: < 3cm
  // Particle Species: Electrons
  // Eta: < 0.8
  //
  
  TObjArray *mcCuts = new TObjArray;
  mcCuts->SetName("fPartGenCuts");

  // 
  AliDebug(2, "Called\n");
  AliCFParticleGenCuts *genCuts = new AliCFParticleGenCuts("fCutsGenMC", "Particle Generation Cuts");
  genCuts->SetRequireIsCharged();
  if(IsRequirePrimary()) { 
    genCuts->SetRequireIsPrimary();
  }
  if(IsRequireProdVertex()){
    AliDebug(3, Form("Vertex Range: fProdVtx[0] %f, fProdVtx[1] %f, fProdVtx[2] %f, fProdVtx[3] %f", fProdVtx[0], fProdVtx[1], fProdVtx[2], fProdVtx[3]));
    genCuts->SetProdVtxRangeX(fProdVtx[0], fProdVtx[1]);
    genCuts->SetProdVtxRangeY(fProdVtx[2], fProdVtx[3]);
    genCuts->SetProdVtxRange2D(kTRUE);  // Use ellipse
  }
  genCuts->SetRequirePdgCode(11, kTRUE);
  if(IsQAOn()) genCuts->SetQAOn(fHistQA);

  // Add
  mcCuts->AddLast(genCuts);
  
  //
  if(IsRequireKineMCCuts()) {  
    AliCFTrackKineCuts *kineMCcuts = new AliCFTrackKineCuts((Char_t *)"fCutsKineMC", (Char_t *)"MC Kine Cuts");
    kineMCcuts->SetPtRange(fPtRange[0], fPtRange[1]);
    kineMCcuts->SetEtaRange(-0.8, 0.8);
    if(IsQAOn()) kineMCcuts->SetQAOn(fHistQA);
    mcCuts->AddLast(kineMCcuts);
  }
   
  fCutList->AddLast(mcCuts);
}

//__________________________________________________________________
void AliHFEcuts::SetAcceptanceCutList(){
  //
  // Initialize Particle (Monte Carlo) acceptance cuts
  // Min. Required Hist for Detectors:
  //          ITS [3]
  //          TPC [2]
  //          TRD [2*nTracklets]
  //          TOF [0]
  //
  AliDebug(2, "Called\n");
  AliCFAcceptanceCuts *accCuts = new AliCFAcceptanceCuts("fCutsAccMC", "MC Acceptance Cuts");
  accCuts->SetMinNHitITS(3);
  accCuts->SetMinNHitTPC(2);
  accCuts->SetMinNHitTRD(2*fMinTrackletsTRD);
  if(IsQAOn()) accCuts->SetQAOn(fHistQA);
  
  TObjArray *partAccCuts = new TObjArray();
  partAccCuts->SetName("fPartAccCuts");
  partAccCuts->AddLast(accCuts);
  fCutList->AddLast(partAccCuts);
}

//__________________________________________________________________
void AliHFEcuts::SetRecKineITSTPCCutList(){
  //
  // Track Kinematics and Quality cuts (Based on the Standard cuts from PWG0)
  //
  // ITS refit
  // Variances:
  //  y: 2
  //  z: 2
  //  sin(phi): 0.5
  //  tan(theta): 0.5
  //  1/pt: 2
  // Min. Number of Clusters:
  //  TPC: 50
  // RefitRequired:
  //  TPC
  // Chi2 per TPC cluster: 3.5
  //
  // Kinematics:
  //  Momentum Range: 100MeV - 20GeV
  //  Eta: < 0.9 (TRD-TOF acceptance)
  //
  AliDebug(2, "Called\n");
  AliCFTrackQualityCuts *trackQuality = new AliCFTrackQualityCuts((Char_t *)"fCutsQualityRec", (Char_t *)"REC Track Quality Cuts");
  trackQuality->SetMinNClusterITS(fMinClustersITS);
  trackQuality->SetMaxChi2PerClusterTPC(fMaxChi2clusterTPC);
  if(fMaxChi2clusterITS >= 0.) trackQuality->SetMaxChi2PerClusterITS(fMaxChi2clusterITS);
  trackQuality->SetStatus(AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit);
  //trackQuality->SetMaxCovDiagonalElements(2., 2., 0.5, 0.5, 2); 

  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("fCutsHFElectronGroupTPC","Extra cuts from the HFE group");
  hfecuts->SetDebugLevel(fDebugLevel);

  // Set the cut in the TPC number of clusters
  hfecuts->SetMinNClustersTPC(fMinClustersTPC, fTPCclusterDef);
  hfecuts->SetClusterRatioTPC(fMinClusterRatioTPC, fTPCratioDef);
  if(fFractionOfSharedTPCClusters > 0.0) hfecuts->SetFractionOfTPCSharedClusters(fFractionOfSharedTPCClusters); 
  
  AliCFTrackKineCuts *kineCuts = new AliCFTrackKineCuts((Char_t *)"fCutsKineRec", (Char_t *)"REC Kine Cuts");
  kineCuts->SetPtRange(fPtRange[0], fPtRange[1]);
  kineCuts->SetEtaRange(-0.8, 0.8);
  
  if(IsQAOn()){
    trackQuality->SetQAOn(fHistQA);
    hfecuts->SetQAOn(fHistQA);
    kineCuts->SetQAOn(fHistQA);
  }
  
  TObjArray *recCuts = new TObjArray;
  recCuts->SetName("fPartRecKineITSTPCCuts");
  recCuts->AddLast(trackQuality);
  recCuts->AddLast(hfecuts);
  recCuts->AddLast(kineCuts);
  fCutList->AddLast(recCuts);
}

//__________________________________________________________________
void AliHFEcuts::SetRecPrimaryCutList(){
  //
  // Primary cuts (based on standard cuts from PWG0):
  //  DCA to Vertex: 
  //    XY: 3. cm
  //    Z:  10. cm
  //  No Kink daughters
  //
  AliDebug(2, "Called\n");
  AliCFTrackIsPrimaryCuts *primaryCut = new AliCFTrackIsPrimaryCuts((Char_t *)"fCutsPrimaryCuts", (Char_t *)"REC Primary Cuts");
  if(IsRequireDCAToVertex()){
    //primaryCut->SetDCAToVertex2D(kTRUE);
    primaryCut->SetMaxDCAToVertexXY(fDCAtoVtx[0]);
    primaryCut->SetMaxDCAToVertexZ(fDCAtoVtx[1]);
  }
  if(IsRequireSigmaToVertex()){
    primaryCut->SetRequireSigmaToVertex(kTRUE);
    if(fSigmaToVtx[0]) primaryCut->SetMaxNSigmaToVertex(fSigmaToVtx[0]);
    if(fSigmaToVtx[1]) primaryCut->SetMaxSigmaDCAxy(fSigmaToVtx[1]);
    if(fSigmaToVtx[2]) primaryCut->SetMaxSigmaDCAz(fSigmaToVtx[2]);
  }
  primaryCut->SetAcceptKinkDaughters(kFALSE);
  if(IsQAOn()) primaryCut->SetQAOn(fHistQA);
  
  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("fCutsPrimaryCutsextra","Extra cuts from the HFE group");
  hfecuts->SetMaxImpactParameterRpar(fMaxImpactParameterRpar);

  TObjArray *primCuts = new TObjArray;
  primCuts->SetName("fPartPrimCuts");
  primCuts->AddLast(primaryCut);
  if(fMaxImpactParameterRpar){
    primCuts->AddLast(hfecuts);
  }
  fCutList->AddLast(primCuts);
}

//__________________________________________________________________
void AliHFEcuts::SetHFElectronITSCuts(){
  //
  // Special Cuts introduced by the HFElectron Group: ITS
  //
  AliDebug(2, "Called\n");
  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("fCutsHFElectronGroupPixels","Extra cuts from the HFE group");
  if(IsRequireITSpixel()){
    hfecuts->SetRequireITSpixel(AliHFEextraCuts::ITSPixel_t(fCutITSPixel));
    hfecuts->SetCheckITSstatus(fCheckITSLayerStatus);
  }
  
  if(IsQAOn()) hfecuts->SetQAOn(fHistQA);
  hfecuts->SetDebugLevel(fDebugLevel);

  TObjArray *hfeCuts = new TObjArray;
  hfeCuts->SetName("fPartHFECutsITS");
  hfeCuts->AddLast(hfecuts);
  fCutList->AddLast(hfeCuts);
}

//__________________________________________________________________
void AliHFEcuts::SetHFElectronTOFCuts(){
  //
  // Special Cuts introduced by the HFElectron Group: TRD
  //
  AliDebug(2, "Called\n");
  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("fCutsHFElectronGroupTOF","Extra cuts from the HFE group on TOF PID");
  if(fTOFPIDStep) hfecuts->SetTOFPID(kTRUE);
  if(fTOFMISMATCHStep) hfecuts->SetTOFMISMATCH(kTRUE);
  if(IsQAOn()) hfecuts->SetQAOn(fHistQA);
  hfecuts->SetDebugLevel(fDebugLevel);
  
  TObjArray *hfeCuts = new TObjArray;
  hfeCuts->SetName("fPartHFECutsTOF");
  hfeCuts->AddLast(hfecuts);
  fCutList->AddLast(hfeCuts);
}

//__________________________________________________________________
void AliHFEcuts::SetHFElectronTPCCuts(){
  //
  // Special Cuts introduced by the HFElectron Group: TPC
  //
  AliDebug(2, "Called\n");
  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("fCutsHFElectronGroupTPC","Extra cuts from the HFE group on TPC PID");
  if(fTPCPIDCLEANUPStep) hfecuts->SetTPCPIDCleanUp(kTRUE);
  if(IsQAOn()) hfecuts->SetQAOn(fHistQA);
  hfecuts->SetDebugLevel(fDebugLevel);
  
  TObjArray *hfeCuts = new TObjArray;
  hfeCuts->SetName("fPartHFECutsTPC");
  hfeCuts->AddLast(hfecuts);
  fCutList->AddLast(hfeCuts);
}


//__________________________________________________________________
void AliHFEcuts::SetHFElectronTRDCuts(){
  //
  // Special Cuts introduced by the HFElectron Group: TRD
  //
  AliDebug(2, "Called\n");
  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("fCutsHFElectronGroupTRD","Extra cuts from the HFE group on TRD PID");
  if(fMinTrackletsTRD > 0.) hfecuts->SetMinTrackletsTRD(fMinTrackletsTRD);
  if(IsQAOn()) hfecuts->SetQAOn(fHistQA);
  hfecuts->SetDebugLevel(fDebugLevel);
  
  TObjArray *hfeCuts = new TObjArray;
  hfeCuts->SetName("fPartHFECutsTRD");
  hfeCuts->AddLast(hfecuts);
  fCutList->AddLast(hfeCuts);
}

//__________________________________________________________________
void AliHFEcuts::SetHFElectronDcaCuts(){
  //
  // Special Cuts introduced by the HFElectron Group: minimum of impact parameter
  //
  AliDebug(2, "Called\n");
  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("fCutsHFElectronGroupDCA","Extra cuts from the HFE group");
  hfecuts->SetMinHFEImpactParamR(fIPCutParams,fIsIPSigmacut);
  if(IsQAOn()) hfecuts->SetQAOn(fHistQA);
  hfecuts->SetDebugLevel(fDebugLevel);

  TObjArray *hfeCuts = new TObjArray;
  hfeCuts->SetName("fPartHFECutsDca");
  hfeCuts->AddLast(hfecuts);
  fCutList->AddLast(hfeCuts);
}

//__________________________________________________________________
Bool_t AliHFEcuts::CheckParticleCuts(UInt_t step, TObject *o){
  //
  // Checks the cuts without using the correction framework manager
  // 
  AliDebug(2, "Called\n");
  TString stepnames[kNcutStepsMCTrack + kNcutStepsRecTrack + kNcutStepsDETrack + kNcutStepsSecvtxTrack + 1] = {"fPartGenCuts","fPartEvCutPileupZ","fPartEvCut","fPartAccCuts","fPartRecNoCuts","fPartRecKineITSTPCCuts", "fPartPrimCuts", "fPartHFECutsITS","fPartHFECutsTOF","fPartHFECutsTPC","fPartHFECutsTRD","fPartHFECutsDca", "fPartHFECutsSecvtx"};
  AliDebug(2, Form("Doing cut %s", stepnames[step].Data()));
 TObjArray *cuts = dynamic_cast<TObjArray *>(fCutList->FindObject(stepnames[step].Data()));
  if(!cuts) return kTRUE;
  TIter it(cuts);
  AliCFCutBase *mycut;
  Bool_t status = kTRUE;
  while((mycut = dynamic_cast<AliCFCutBase *>(it()))){
    status &= mycut->IsSelected(o);
  }
  return status;
}


//__________________________________________________________________
Bool_t AliHFEcuts::CheckEventCuts(const char*namestep, TObject *o){
  //
  // Checks the cuts without using the correction framework manager
  // 
  AliDebug(2, "Called\n");
  TObjArray *cuts = dynamic_cast<TObjArray *>(fCutList->FindObject(namestep));
  if(!cuts) return kTRUE;
  TIter it(cuts);
  AliCFCutBase *mycut;
  Bool_t status = kTRUE;
  while((mycut = dynamic_cast<AliCFCutBase *>(it()))){
    status &= mycut->IsSelected(o);
  }
  return status;
}

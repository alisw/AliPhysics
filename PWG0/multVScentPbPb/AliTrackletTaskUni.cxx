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

//////////////////////////////////////////////////////////////////////////
// Class AliTrackletTaskUni                                             //
// Analysis task to study performance and combinatorial background      //
// in tracklet reconstruction                                           //
// Author:  M. Nicassio (INFN Bari)                                     //
// Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it         //
//////////////////////////////////////////////////////////////////////////

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
//#include "../ANALYSIS/EventMixing/AliMixEventInputHandler.h"
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
#include "AliTrackletTaskUni.h"
#include "AliITSMultRecBg.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

ClassImp(AliTrackletTaskUni)

const char*  AliTrackletTaskUni::fgkPDGNames[] = {
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

const int AliTrackletTaskUni::fgkPDGCodes[] = {
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
AliTrackletTaskUni::AliTrackletTaskUni(const char *name)
  : AliAnalysisTaskSE(name),
*/  
//________________________________________________________________________
AliTrackletTaskUni::AliTrackletTaskUni(const char *name) 
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
  fMultCutMin(0),
  fMultCutMax(99999),
  fMCV0Scale(0.7520),
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
  fMultReco(0),
  fRPTree(0),
  fRPTreeMix(0),
  fStack(0),
  fMCEvent(0),
  fDontMerge(kFALSE)    
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
AliTrackletTaskUni::~AliTrackletTaskUni()
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
void AliTrackletTaskUni::UserCreateOutputObjects() 
{
  //

  if (fDontMerge) {
    OpenFile(1);
    printf("No merging will be done\n");
  }

  fOutput = new TList();
  fOutput->SetOwner(); 
  //
  AliCDBManager *man = AliCDBManager::Instance();
  if (fUseMC) {
    printf("Loading MC geometry\n");
    Bool_t newGeom = kTRUE;
    man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
    if (newGeom) {
      // new geom
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",130844,-1,-1);
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130844,-1,-1)) AliFatal("Failed to misalign geometry");
    }
    else {
      // old geom
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",130844,7,-1);
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130845,5,-1)) AliFatal("Failed to misalign geometry");
    }
  }
  else {
    printf("Loading Raw geometry\n");
    man->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB"); //man->SetRun(137366);
    AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",137366);
    AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
    if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",137366,-1,-1)) AliFatal("Failed to misalign geometry");
  }
  //
  printf("Geometry loaded\n");  
  // Create histograms
  //---------------------------------------------Standard histos per tracklet type--->>
  UInt_t hPattern = 0xffffffff;
  hPattern &= ~(BIT(kHEtaZvDPhiS));  // pattern of histos to fill
  fHistosTrData                      = BookHistosSet("TrData",hPattern);
  if (GetDoInjection()) fHistosTrInj = BookHistosSet("TrInj",hPattern);
  if (GetDoRotation())  fHistosTrRot = BookHistosSet("TrRot",hPattern);
  if (GetDoMixing())    fHistosTrMix = BookHistosSet("TrMix",hPattern);
  if (fUseMC) {
    hPattern = 0xffffffff;  hPattern &= ~(BIT(kHEtaZvDist)|(BIT(kHEtaZvDPhiS)));
    fHistosTrPrim  = BookHistosSet("TrPrim",hPattern);
    //
    hPattern = 0xffffffff;  hPattern &= ~(BIT(kHEtaZvDist)|(BIT(kHEtaZvDPhiS)));
    fHistosTrSec   = BookHistosSet("TrSec",hPattern);
    //
    hPattern = 0xffffffff;
    fHistosTrComb  = BookHistosSet("TrComb",hPattern);
    //
    hPattern = 0xffffffff;  hPattern &= ~(BIT(kHEtaZvDist)|(BIT(kHEtaZvDPhiS)));
    fHistosTrCombU = BookHistosSet("TrCombU",hPattern);
    //
    if (fCheckReconstructables) {
      hPattern = 0xffffffff;  hPattern &= ~(BIT(kHEtaZvDist)|(BIT(kHEtaZvDPhiS)));
      fHistosTrRcblPrim = BookHistosSet("TrRcblPrim",hPattern);
      //
      hPattern = 0xffffffff;  hPattern &= ~(BIT(kHEtaZvDist)|(BIT(kHEtaZvDPhiS)));
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
void AliTrackletTaskUni::RegisterStat() 
{
  static Bool_t done = kFALSE;
  if (done) return;
  TH1* hstat;
  Printf("Registering used settings");
  TList *lst = dynamic_cast<TList*>(GetOutputData(1));
  if (lst && (hstat=(TH1*)lst->FindObject("hStat"))) {
    // fill used settings
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
    hstat->Fill(kTrcMin,fMultCutMin);    
    hstat->Fill(kTrcMax,fMultCutMax);    
    hstat->Fill(kMCV0Scale, fMCV0Scale);
    //
    hstat->Fill(kOneUnit,1.);    
  }
  else Printf("Did not find stat histo");
  done = kTRUE;
  //
}


//________________________________________________________________________
void AliTrackletTaskUni::UserExec(Option_t *) 
{
  // Main loop
  //
  if (fDontMerge) RegisterStat();
  Bool_t needRecPoints = GetDoNormalReco() || GetDoInjection() || GetDoRotation() || GetDoMixing();

  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  fRPTree = fRPTreeMix = 0;
  AliESDInputHandler *handler = (AliESDInputHandler*)anMan->GetInputEventHandler();
  AliESDInputHandlerRP *handRP = 0;
  if (needRecPoints) {
    handRP = (AliESDInputHandlerRP*)handler;
    if (!handRP) { printf("No RP handler\n"); return; }
  }
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
  TH1F* hstat = (TH1F*)fHistosCustom->UncheckedAt(kHStat);
  //
  hstat->Fill(kEvTot0); // RS
  const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
  if (vtxESD->GetNContributors()<1) return;
  TString vtxTyp = vtxESD->GetTitle();
  if (vtxTyp.Contains("vertexer: Z")) {
    if (vtxESD->GetDispersion()>0.04) return;
    if (vtxESD->GetZRes()>0.25) return;
  }
  // pile-up rejection
  if (esd->IsPileupFromSPD(3,0.8)) {
    hstat->Fill(kEvTotPlp); // RS
    return;
  }
  //
  const AliMultiplicity* multESD = esd->GetMultiplicity();
  /*
  const AliESDVertex* vtxESDTPC = esd->GetPrimaryVertexTPC();
  if (vtxESDTPC->GetNContributors()<1 ||
      vtxESDTPC->GetNContributors()<(-10.+0.25*multESD->GetNumberOfITSClusters(0))) return;
  */
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
   //------------------------------------------------------
  AliESDZDC *esdZDC = esd->GetESDZDC();
  float zdcEnergy=0,zemEnergy=0;
  if (esdZDC) {
    zdcEnergy = (esdZDC->GetZDCN1Energy() + esdZDC->GetZDCP1Energy() + esdZDC->GetZDCN2Energy()+ esdZDC->GetZDCP2Energy());
    zemEnergy = (esdZDC->GetZDCEMEnergy(0)+esdZDC->GetZDCEMEnergy(1))/8.; 
  }
  //-----------------------------------------------------
  Float_t multV0A=0,multV0C=0;
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  if (esdV0) {
    multV0A = esdV0->GetMTotV0A();
    multV0C = esdV0->GetMTotV0C();
  }
  if (fUseMC) {
    const double v0Scale = fMCV0Scale;
    multV0A *= v0Scale;
    multV0C *= v0Scale;    
  }
  // registed Ntracklets and ZVertex of the event
  ((TH1F*)fHistosCustom->UncheckedAt(kHZVtxNoSel))->Fill(esdvtx[2]);
  ((TH1F*)fHistosCustom->UncheckedAt(kHNTrackletsNoSel))->Fill(multESD->GetNumberOfTracklets());      
  ((TH1F*)fHistosCustom->UncheckedAt(kHNClSPD1NoSel))->Fill(multESD->GetNumberOfITSClusters(0));
  ((TH1F*)fHistosCustom->UncheckedAt(kHNClSPD2NoSel))->Fill(multESD->GetNumberOfITSClusters(1));
  ((TH1F*)fHistosCustom->UncheckedAt(kHV0NoSel))->Fill(multV0A+multV0C);
  //  ((TH2F*)fHistosCustom->UncheckedAt(kHV0NClSPD2NoSel))->Fill(multV0A+multV0C,multESD->GetNumberOfITSClusters(1));
  //
  //  printf("ESD vertex! %f %f %f, %d contributors\n",esdvtx[0],esdvtx[1],esdvtx[2],vtxESD->GetNContributors());

  if(vtxf[2] < fZVertexMin || vtxf[2] > fZVertexMax) return;
  ((TH2F*)fHistosCustom->UncheckedAt(kHV0NClSPD2NoSel))->Fill(multV0A+multV0C,multESD->GetNumberOfITSClusters(1));
  //
  ///  double mltTst = fUseMC ?  multESD->GetNumberOfITSClusters(1) : multV0A+multV0C;  
  //  double mltTst = multESD->GetNumberOfITSClusters(1); //RRR
  double mltTst = multV0A+multV0C; //RRR

  if ((mltTst<fMultCutMin) || (mltTst>fMultCutMax)) return;
  //
  //  if((multESD->GetNumberOfTracklets() < fMultCutMin) || (multESD->GetNumberOfTracklets() > fMultCutMax)) return;
  //
  //  printf("Multiplicity from ESD:\n");
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
  if (needRecPoints) {
    fRPTree = handRP->GetTreeR("ITS");
    if (!fRPTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  }
  //
  //
  // registed Ntracklets and ZVertex of the event
  ((TH1F*)fHistosCustom->UncheckedAt(kHZVtx))->Fill(esdvtx[2]);
  ((TH1F*)fHistosCustom->UncheckedAt(kHNTracklets))->Fill(multESD->GetNumberOfTracklets());      
  //
  if (fUseMC) FillMCPrimaries( ((TH2F*)fHistosCustom->UncheckedAt(kHZVEtaPrimMC)) );
  // fill N clusters
  ((TH1F*)fHistosCustom->UncheckedAt(kHNClSPD1))->Fill(multESD->GetNumberOfITSClusters(0));
  ((TH1F*)fHistosCustom->UncheckedAt(kHNClSPD2))->Fill(multESD->GetNumberOfITSClusters(1));
  ((TH1F*)fHistosCustom->UncheckedAt(kHV0))->Fill(multV0A+multV0C);
  //
  // normal reconstruction
  static int prnRec = 10;
  static int prnInj = 10;
  //
  hstat->Fill(kEvProcData);
  if (GetDoNormalReco() || GetDoInjection()) { // for the injection the normal reco should be done
    InitMultReco();
    fMultReco->Run(fRPTree, vtxf);
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    if (mlt && (--prnRec)>0) {
      printf("Multiplicity Reconstructed: %d\n",prnRec);
      mlt->Print();
    }
    if (GetDoNormalReco()) FillHistos(kData,mlt);
    //
  }
  if (!GetDoNormalReco()) FillHistos(kData,multESD); // fill data histos from ESD
  //
  // Injection: it must come right after the normal reco since needs its results
  if (GetDoInjection()) {
    if (!fMultReco) InitMultReco(); // in principle, not needed, the reco is created above
    fMultReco->SetRecType(AliITSMultRecBg::kBgInj);
    fMultReco->Run(fRPTree, vtxf);    
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    if (mlt && (--prnInj)>0) {
      printf("Multiplicity from Injection: %d\n",prnInj);
      mlt->Print();
    }
    // if (mlt) mlt->Print();
    hstat->Fill(kEvProcInj);
    FillHistos(kBgInj,mlt);
  }
  //
  // Rotation
  if (GetDoRotation()) {
    InitMultReco();
    fMultReco->SetRecType(AliITSMultRecBg::kBgRot);
    fMultReco->SetPhiRotationAngle(fPhiRot);
    fMultReco->Run(fRPTree, vtxf);
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    hstat->Fill(kEvProcRot);
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
    ((TH1F*)fHistosCustom->UncheckedAt(kHZVtxMixDiff))->Fill(vtxESDmix->GetZ()-esdvtx[2]);
    ((TH1F*)fHistosCustom->UncheckedAt(kHNTrMixDiff) )->
      Fill(esdMix->GetMultiplicity()->GetNumberOfTracklets() - multESD->GetNumberOfTracklets());
    //
    InitMultReco();
    fMultReco->SetRecType(AliITSMultRecBg::kBgMix);
    fMultReco->Run(fRPTree, vtxf,fRPTreeMix);
    printf("Multiplicity from Mixing:\n");
    AliMultiplicity* mlt = fMultReco->GetMultiplicity();
    if (mlt) mlt->Print();
    hstat->Fill(kEvProcMix);
    FillHistos(kBgMix,mlt);
    */
    AliFatal("Mixing is outphased");
    //
  }
  // =============================================================================<<<
  //
  delete fMultReco; 
  fMultReco = 0;
  //
}      
//________________________________________________________________________
void AliTrackletTaskUni::Terminate(Option_t *) 
{
  Printf("Terminating...");
  RegisterStat();
  //  AliAnalysisTaskSE::Terminate();
}


//_________________________________________________________________________
void AliTrackletTaskUni::InitMultReco()
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
TObjArray* AliTrackletTaskUni::BookCustomHistos()
{
  // book custom histos, not related to specific tracklet type
  TObjArray* histos = new TObjArray();
  TH1F* hstat;
  //
  // ------------ job parameters, statistics ------------------------------>>>
  hstat = new TH1F("hStat","Run statistics",kNStatBins,0.5,kNStatBins+0.5);
  //
  for (int ib=1;ib<=kNStatBins;ib++) hstat->GetXaxis()->SetBinLabel(ib,"-"); // dummy label
  hstat->GetXaxis()->SetBinLabel(kEvTot0, "Ev.Tot0");
  hstat->GetXaxis()->SetBinLabel(kEvTot, "Ev.Tot");
  hstat->GetXaxis()->SetBinLabel(kEvTotPlp, "Ev.Tot Plp");
  //
  hstat->GetXaxis()->SetBinLabel(kEvProcData,"Ev.ProcData");
  hstat->GetXaxis()->SetBinLabel(kEvProcInj,"Ev.ProcInj");
  hstat->GetXaxis()->SetBinLabel(kEvProcRot,"Ev.ProcRot");
  hstat->GetXaxis()->SetBinLabel(kEvProcMix,"Ev.ProcMix");
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
  hstat->GetXaxis()->SetBinLabel(kEtaMin,"#eta_{min}");
  hstat->GetXaxis()->SetBinLabel(kEtaMax,"#eta_{max}");
  hstat->GetXaxis()->SetBinLabel(kZVMin,"ZV_{min} cut");
  hstat->GetXaxis()->SetBinLabel(kZVMax,"ZV_{max} cut");
  hstat->GetXaxis()->SetBinLabel(kTrcMin,"Mult_{min} cut");
  hstat->GetXaxis()->SetBinLabel(kTrcMax,"Mult_{max} cut");  
  //
  hstat->GetXaxis()->SetBinLabel(kOneUnit,"ScaleMerge");
  hstat->GetXaxis()->SetBinLabel(kNWorkers,"Workers");
  hstat->Fill(kNWorkers);
  //  
  AddHisto(histos,hstat,kHStat);
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
  TH1F* hzvns = new  TH1F("zvNoSel","Z vertex Before Selection",nZVBins,zMn,zMx);
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
  // V0 
  //TH2F* hnV0SPD2ns = new  TH2F("V0NDP2NoSel","NSPD2 vs V0 signal Before Selection",2500,0,20000,1400,0,maxmltSPD2);
  TH2F* hnV0SPD2ns = new  TH2F("V0NDP2NoMltSel","NSPD2 vs V0 signal Before Mlt Selection",100,0,maxmltV0,100,0,maxmltSPD2);
  hnV0SPD2ns->GetXaxis()->SetTitle("V0 signal");
  hnV0SPD2ns->GetYaxis()->SetTitle("N Clus SPD2 ");
  AddHisto(histos,hnV0SPD2ns,kHV0NClSPD2NoSel);
  //
  //
  // Z vertex distribution for selected events
  TH1F* hzv = new  TH1F("zv","Z vertex",nZVBins,zMn,zMx);
  hzv->GetXaxis()->SetTitle("Zvertex");
  AddHisto(histos,hzv,kHZVtx);
  //
  // N tracklets for processed events
  TH1F* hnt = new  TH1F("Ntracklets","N Tracklets",nbmltSPD2,0,maxmltSPD2);
  hnt->GetXaxis()->SetTitle("N tracklets");
  AddHisto(histos,hnt,kHNTracklets);
  //
  // N SPD1 clusters
  TH1F* hncl1 = new  TH1F("NClustersSPD1","N Clusters on SPD1",nbmltSPD2,0,maxmltSPD2);
  hncl1->GetXaxis()->SetTitle("N Clus SPD1");
  AddHisto(histos,hncl1,kHNClSPD1);
  //
  // N SPD2 clusters
  TH1F* hncl2 = new  TH1F("NClustersSPD2","N Clusters on SPD2",nbmltSPD2,0,maxmltSPD2);
  hncl2->GetXaxis()->SetTitle("N Clus SPD2");
  AddHisto(histos,hncl2,kHNClSPD2);
  //
  // V0 
  TH1F* hnV0 = new  TH1F("V0","V0 signal",nbmltV0,0,maxmltV0);
  hnV0->GetXaxis()->SetTitle("V0 signal");
  AddHisto(histos,hnV0,kHV0);
  //
  //----------------------------------------------------------------------
  int nEtaBinsS = int((fEtaMax-fEtaMin)/0.1);
  if (nEtaBinsS<1) nEtaBins = 1;
  //
  int nZVBinsS = int(fZVertexMax-fZVertexMin);
  if (nZVBinsS<1) nZVBinsS = 1;

  if (fUseMC) {
    // Z vertex vs Eta distribution for primaries
    TH2F* hzvetap = new  TH2F("zvEtaPrimMC","Z vertex vs eta PrimMC",nEtaBinsS,fEtaMin,fEtaMax,nZVBinsS,fZVertexMin,fZVertexMax);
    hzvetap->GetXaxis()->SetTitle("#eta");
    hzvetap->GetYaxis()->SetTitle("Zvertex");
    AddHisto(histos,hzvetap,kHZVEtaPrimMC);
  }
  //
  if (GetDoMixing()) {
    //
    /*
    // Difference in Z vertex for mixed events
    TH1F* hzdiff = new TH1F("MixSPDVertexDiff","SPD #Delta Z Vertex distribution ",100,-5,5);
    hzdiff->GetXaxis()->SetTitle("#Delta Z Vertex [cm]");
    hzdiff->GetYaxis()->SetTitle(Form("Entries / %1.2f [cm]",10./100.));
    AddHisto(histos,hzdiff,kHZVtxMixDiff);
    //
    // Difference in N tracklets for mixed events
    TH1F* hntdiff = new TH1F("MixNTrackletsDiff"," SPD tracklets Diff ",2000,-3000,3000);
    hntdiff->GetXaxis()->SetTitle("# tracklet diff");
    AddHisto(histos,hntdiff,kHNTrMixDiff);
    */
  }
  // 
  // --------------------------------------------------
  int nDistBins = int(fNStdDev)*2;
  int npdg = sizeof(fgkPDGNames)/sizeof(char*);
  TH2F* hpdgP = new TH2F("pdgPrim","primary PDG",npdg,0,npdg,nDistBins,0,fNStdDev);
  AddHisto(histos,hpdgP,kHPrimPDG);
  TH2F* hpdgS = new TH2F("pdgSec","secondary PDG",npdg,0,npdg,nDistBins,0,fNStdDev);
  AddHisto(histos,hpdgS,kHSecPDG);
  TH2F* hpdgPP = new TH2F("pdgPrimPar","primary parent PDG ",npdg,0,npdg,nDistBins,0,fNStdDev);
  AddHisto(histos,hpdgPP,kHPrimParPDG);
  TH2F* hpdgSP = new TH2F("pdgSecPar","secondary parent PDG",npdg,0,npdg,nDistBins,0,fNStdDev);
  AddHisto(histos,hpdgSP,kHSecParPDG);
  for (int i=0;i<npdg;i++) {
    hpdgP->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
    hpdgS->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
    hpdgPP->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
    hpdgSP->GetXaxis()->SetBinLabel(i+1,fgkPDGNames[i]);
  }
  //
  // -------------------------------------------------
  histos->SetOwner(kFALSE);

  return histos;
}

//_________________________________________________________________________
TObjArray* AliTrackletTaskUni::BookHistosSet(const char* pref, UInt_t selHistos) 
{
  // book standard set of histos attaching the pref in front of the name/title
  //
  const int kNDPhiBins = 100;
  const int kNDThtBins = 100;
  int nDistBins = int(fNStdDev)*4;

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
  THnSparseF* hsp;
  char buffn[100],bufft[500];
  //
  if (selHistos & (0x1<<kHEtaZvDist) ) {
    // sparse 3d histo for w.dist vs zv vs eta 
    sprintf(buffn,"%s_DistZvEta",pref);
    sprintf(bufft,"(%s)Weighted Dist.(#Delta) vs Zv vs Eta",pref);
    int nbnEZD[3]    = {nEtaBins,nZVBins,nDistBins};
    double xmnEZD[3] = { fEtaMin,fZVertexMin,0};
    double xmxEZD[3] = { fEtaMax,fZVertexMax,fNStdDev};
    hsp = new THnSparseF(buffn,bufft,3,nbnEZD,xmnEZD,xmxEZD);
    hsp->GetAxis(0)->SetTitle("#eta");
    hsp->GetAxis(1)->SetTitle("Zv");
    sprintf(bufft,"#Delta=[(#Delta#varphi-#delta_{#varphi})/#sigma#varphi]^{2}+"
	    "[#Delta#theta%s/#sigma#theta]^{2}",fScaleDTBySin2T ? "*sin^{-2}(#theta)":"");
    hsp->GetAxis(2)->SetTitle(bufft);
    AddHisto(histos,hsp,kHEtaZvDist);
  }
  //
  if (selHistos & (0x1<<kHEtaZvDPhiS) ) {
    // sparse 3d histo for dphi vs zv vs eta 
    sprintf(buffn,"%s_DistZvDPhiS",pref);
    sprintf(bufft,"(%s) #Delta#varphi-#delta_{#varphi} vs Zv vs Eta",pref);
    int nbnEZP[3]    = {nEtaBins,nZVBins, int(dphir*2/0.005)};
    double xmnEZP[3] = { fEtaMin,fZVertexMin,-dphir};
    double xmxEZP[3] = { fEtaMax,fZVertexMax, dphir};
    hsp = new THnSparseF(buffn,bufft,3,nbnEZP,xmnEZP,xmxEZP);
    hsp->GetAxis(0)->SetTitle("#eta");
    hsp->GetAxis(1)->SetTitle("Zv");
    hsp->GetAxis(2)->SetTitle("#Delta#varphi-#delta_{#varphi} [rad]");
    AddHisto(histos,hsp,kHEtaZvDPhiS);
  }
  //
  if (selHistos & (0x1<<kHEtaZvCut) ) {
    sprintf(buffn,"%s_ZvEtaCutT",pref);
    sprintf(bufft,"(%s) Zv vs Eta with tracklet cut",pref);
    h2 = new TH2F(buffn,bufft,nEtaBins,fEtaMin,fEtaMax, nZVBins, fZVertexMin,fZVertexMax);
    h2->GetXaxis()->SetTitle("#eta");
    h2->GetYaxis()->SetTitle("Zv");
    AddHisto(histos,h2,kHEtaZvCut);
  }
  //
  if (selHistos & (0x1<<kHDPhiDTheta) ) {
    sprintf(buffn,"%s_dPhidTheta",pref);
    sprintf(bufft,"(%s) #Delta#theta vs #Delta#varphi",pref);
    h2 = new TH2F(buffn,bufft,kNDPhiBins,-dphir,dphir,kNDThtBins,-dthtr,dthtr);
    h2->GetXaxis()->SetTitle("#Delta#varphi [rad]");
    h2->GetYaxis()->SetTitle("#Delta#theta [rad]");
    AddHisto(histos,h2,kHDPhiDTheta);
  }
  //
  if (selHistos & (0x1<<kHDPhiSDThetaX) ) {
    sprintf(buffn,"%s_dPhiSdThetaX",pref);
    sprintf(bufft,"(%s) #Delta#theta%s vs #Delta#varphi-#delta_{#varphi}",pref,fScaleDTBySin2T ? "/sin^{2}(#theta)":"");
    h2 = new TH2F(buffn,bufft,kNDPhiBins,-dphir,dphir,kNDThtBins,-dthtr,dthtr);
    h2->GetXaxis()->SetTitle("#Delta#varphi-#delta_{#varphi} [rad]");
    sprintf(bufft,"#Delta#theta%s",fScaleDTBySin2T ? "/sin^{2}(#theta)":"");
    h2->GetYaxis()->SetTitle(bufft);
    AddHisto(histos,h2,kHDPhiSDThetaX);
  }
  //
  if (selHistos & (0x1<<kHEtaDPhiS) ) {
    sprintf(buffn,"%s_EtaDPhiS",pref);
    sprintf(bufft,"(%s) #Delta#varphi-#delta_{#varphi} vs #eta",pref);
    h2 = new TH2F(buffn,bufft,nEtaBins, fEtaMin,fEtaMax,kNDPhiBins,-dphir,dphir);
    h2->GetXaxis()->SetTitle("#eta");
    h2->GetYaxis()->SetTitle("#Delta#varphi-#delta_{#varphi} [rad]");
    AddHisto(histos,h2,kHEtaDPhiS);
  }
  //
  if (selHistos & (0x1<<kHEtaDThetaX) ) {
    sprintf(buffn,"%s_EtaDThetaX",pref);
    sprintf(bufft,"(%s) #Delta#theta%s vs #eta",pref,fScaleDTBySin2T ? "/sin^{2}(#theta)":"");
    h2 = new TH2F(buffn,bufft,nEtaBins, fEtaMin,fEtaMax,kNDThtBins,-dthtr,dthtr);
    h2->GetXaxis()->SetTitle("#eta");
    sprintf(bufft,"#Delta#theta%s",fScaleDTBySin2T ? "/sin^{2}(#theta)":"");
    h2->GetYaxis()->SetTitle(bufft);
    AddHisto(histos,h2,kHEtaDThetaX);
  }
  //
  if (selHistos & (0x1<<kHEtaDist) ) {
    sprintf(buffn,"%s_EtaDist",pref);
    sprintf(bufft,"(%s) Weighted Dist.(#Delta) vs #eta",pref);
    h2 = new TH2F(buffn,bufft,nEtaBins, fEtaMin,fEtaMax,nDistBins,0,fNStdDev);
    h2->GetXaxis()->SetTitle("#eta");
    sprintf(bufft,"#Delta=[(#Delta#varphi-#delta_{#varphi})/#sigma#varphi]^{2}+"
	    "[#Delta#theta%s/#sigma#theta]^{2}",fScaleDTBySin2T ? "*sin^{-2}(#theta)":"");
    h2->GetYaxis()->SetTitle(bufft);
    AddHisto(histos,h2,kHEtaDist);
  }
  //
  if (selHistos & (0x1<<kHZvDPhiS) ) {
    sprintf(buffn,"%s_ZvDPhiS",pref);
    sprintf(bufft,"(%s) #Delta#varphi-#delta_{#varphi} vs Zv",pref);
    h2 = new TH2F(buffn,bufft, nZVBins, fZVertexMin,fZVertexMax, kNDPhiBins,-dphir,dphir);
    h2->GetXaxis()->SetTitle("Zv");
    h2->GetYaxis()->SetTitle("#Delta#varphi-#delta_{#varphi} [rad]");
    AddHisto(histos,h2,kHZvDPhiS);
  }
  //
  if (selHistos & (0x1<<kHZvDThetaX) ) {
    sprintf(buffn,"%s_ZvDThetaX",pref);
    sprintf(bufft,"(%s) #Delta#theta%s vs Zv",pref,fScaleDTBySin2T ? "/sin^{2}(#theta)":"");
    h2 = new TH2F(buffn,bufft, nZVBins, fZVertexMin,fZVertexMax, kNDThtBins,-dthtr,dthtr);
    h2->GetXaxis()->SetTitle("Zv");
    h2->GetYaxis()->SetTitle("#Delta#varphi-#delta_{#varphi} [rad]");
    AddHisto(histos,h2,kHZvDThetaX);
  }
  //
  if (selHistos & (0x1<<kHZvDist) ) {
    sprintf(buffn,"%s_ZvDist",pref);
    sprintf(bufft,"(%s) Weighted Dist.(#Delta) vs Zv",pref);
    h2 = new TH2F(buffn,bufft, nZVBins, fZVertexMin,fZVertexMax, nDistBins,0,fNStdDev);
    h2->GetXaxis()->SetTitle("Zv");
    sprintf(bufft,"#Delta=[(#Delta#varphi-#delta_{#varphi})/#sigma#varphi]^{2}+"
	    "[#Delta#theta%s/#sigma#theta]^{2}",fScaleDTBySin2T ? "*sin^{-2}(#theta)":"");
    h2->GetYaxis()->SetTitle(bufft);
    AddHisto(histos,h2,kHZvDist);
  }
  //
  histos->SetOwner(kFALSE);
  return histos;
}

//_________________________________________________________________________
void AliTrackletTaskUni::AddHisto(TObjArray* histos, TObject* h, Int_t at)
{
  // add single histo to the set
  if (at>=0) histos->AddAtAndExpand(h,at);
  else       histos->Add(h);
  fOutput->Add(h);
}

//_________________________________________________________________________
void AliTrackletTaskUni::FillHistos(Int_t type, const AliMultiplicity* mlt)
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
  //
  for (int itr=ntr;itr--;) {
    //
    //---------------------------------------- CHECK ------------------------------>>>
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
    //---------------------------------------- CHECK ------------------------------<<<
    //
    double theta  = mlt->GetTheta(itr);
    double phi    = mlt->GetPhi(itr);
    double dtheta = mlt->GetDeltaTheta(itr);
    double dphi   = mlt->GetDeltaPhi(itr);
    double dist   = mlt->CalcDist(itr);
    //
    FillHistosSet(histos,phi,theta,dphi,dtheta,dist);
    // special handling for mc info
    if (fillMC && fStack) {
      int lab0 = mlt->GetLabel(itr,0);
      int lab1 = mlt->GetLabel(itr,1);
      int typeMC = 2; // comb.bg.
      if (lab0 == lab1)	typeMC = fStack->IsPhysicalPrimary(lab0) ? 0:1; // prim or sec
      if      (typeMC==0) FillHistosSet(fHistosTrPrim,phi,theta,dphi,dtheta,dist); // primary
      else if (typeMC==1) FillHistosSet(fHistosTrSec, phi,theta,dphi,dtheta,dist); // secondary
      else {
	FillHistosSet(fHistosTrComb,phi,theta,dphi,dtheta,dist); // comb
	// for combinatorals fill also the uncorrelated part
	if (fMultReco) {
	  float *trl = fMultReco->GetTracklet(itr);
	  int clId0 = (int)trl[AliITSMultReconstructor::kClID1];
	  int clId1 = (int)trl[AliITSMultReconstructor::kClID2];
	  float *clLabs0 = fMultReco->GetClusterOfLayer(0,clId0) + AliITSMultReconstructor::kClMC0;
	  float *clLabs1 = fMultReco->GetClusterOfLayer(1,clId1) + AliITSMultReconstructor::kClMC0;
	  if (!HaveCommonParent(clLabs0,clLabs1)) 
	    FillHistosSet(fHistosTrCombU,phi,theta,dphi,dtheta,dist);
	}
      } // combinatorials
      FillSpecies(typeMC, lab0, dist);
      if (fCheckReconstructables) CheckReconstructables();
    }
  }
  //
}

//_________________________________________________________________________
void AliTrackletTaskUni::FillMCPrimaries(TH2F* hetaz)
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
  for (int itr=ntr;itr--;) {
    if (!fStack->IsPhysicalPrimary(itr)) continue;
    AliMCParticle *part  = (AliMCParticle*)fMCEvent->GetTrack(itr);
    if (!part->Charge()) continue;
    //
    //---------------------------------------- CHECK ------------------------------>>>
    Float_t dz = part->Zv() - vtxMC[2];
    if (TMath::Abs(dz)>1e-6) continue; // reject
    //---------------------------------------- CHECK ------------------------------<<<
    //

    Float_t theta = part->Theta();
    if (theta<1e-6 || theta>TMath::Pi()-1e-6) continue;
    Float_t eta = part->Eta();
    hetaz->Fill(eta, fESDVtx[2]);
  }
  //
}

//_________________________________________________________________________
void AliTrackletTaskUni::FillHistosSet(TObjArray* histos, double /*phi*/,double theta,
				       double dphi,double dtheta,double dist) 
{
  // fill standard set of histos
  if (dist>fNStdDev) return;
  //
  double dphiS  = TMath::Abs(dphi) - fDPhiShift;
  if (dphi<0) dphiS = -dphiS;
  double eta    = -TMath::Log(TMath::Tan(theta/2));
  if (eta<fEtaMin || eta>fEtaMax) return;
  //
  double dThetaX = dtheta;
  if (fScaleDTBySin2T) {
    double sint   =  TMath::Sin(theta);
    dThetaX /= (sint*sint);
  }
  //
  if (fCutOnDThetaX && TMath::Abs(dThetaX)>fDThetaWindow) return;
  //
  if (histos->UncheckedAt(kHEtaZvDist)) {
    double ezd[3] = {eta,fESDVtx[2],dist};
    ((THnSparseF*)histos->UncheckedAt(kHEtaZvDist))->Fill(ezd);
  }
  //
  if (histos->UncheckedAt(kHEtaZvDPhiS)) {
    double ezp[3] = {eta,fESDVtx[2],dphiS};
    ((THnSparseF*)histos->UncheckedAt(kHEtaZvDPhiS))->Fill(ezp);
  }
  //
  if (histos->UncheckedAt(kHDPhiDTheta)) 
    ((TH2F*)histos->UncheckedAt(kHDPhiDTheta))->Fill(dphi,dtheta);
  //
  if (histos->UncheckedAt(kHDPhiSDThetaX)) 
    ((TH2F*)histos->UncheckedAt(kHDPhiSDThetaX))->Fill(dphiS,dThetaX);
  //
  if (histos->UncheckedAt(kHEtaDPhiS))
    ((TH2F*)histos->UncheckedAt(kHEtaDPhiS))->Fill(eta,dphiS);
  //
  if (histos->UncheckedAt(kHEtaDThetaX))
    ((TH2F*)histos->UncheckedAt(kHEtaDThetaX))->Fill(eta,dThetaX);
  //
  if (histos->UncheckedAt(kHEtaDist))
    ((TH2F*)histos->UncheckedAt(kHEtaDist))->Fill(eta,dist);
  //
  if (histos->UncheckedAt(kHZvDPhiS)) 
    ((TH2F*)histos->UncheckedAt(kHZvDPhiS))->Fill(fESDVtx[2],dphiS);
  //
  if (histos->UncheckedAt(kHZvDThetaX))
    ((TH2F*)histos->UncheckedAt(kHZvDThetaX))->Fill(fESDVtx[2],dThetaX);
  //
  if (histos->UncheckedAt(kHZvDist)) 
    ((TH2F*)histos->UncheckedAt(kHZvDist))->Fill(fESDVtx[2],dist);
  //
  if (dist<=1 && histos->UncheckedAt(kHEtaZvCut)) 
    ((TH2F*)histos->UncheckedAt(kHEtaZvCut))->Fill(eta,fESDVtx[2]);
  //
}
 
//__________________________________________________________________
void AliTrackletTaskUni::FillSpecies(Int_t primsec, Int_t id, double dist)
{
  // fill PDGcode 
  TH1 *hPart=0,*hParent=0;
  if (primsec==0) {
    hPart   = (TH1F*)fHistosCustom->UncheckedAt(kHPrimPDG);
    hParent = (TH1F*)fHistosCustom->UncheckedAt(kHPrimParPDG);
  } 
  else if (primsec==1) {
    hPart   = (TH1F*)fHistosCustom->UncheckedAt(kHSecPDG);
    hParent = (TH1F*)fHistosCustom->UncheckedAt(kHSecParPDG);    
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
  hPart->Fill(pdgBin,dist);
  hParent->Fill(pdgBinPar,dist);
  //
}

//_________________________________________________________________________
Int_t AliTrackletTaskUni::GetPdgBin(Int_t pdgCode)
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
Bool_t AliTrackletTaskUni::HaveCommonParent(const float* clLabs0,const float* clLabs1)
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
void AliTrackletTaskUni::CheckReconstructables()
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
      FillHistosSet(histos,
		    trComp[AliITSMultReconstructor::kTrPhi][ind],trComp[AliITSMultReconstructor::kTrTheta][ind],
		    trComp[AliITSMultReconstructor::kTrDPhi][ind],trComp[AliITSMultReconstructor::kTrDTheta][ind],
		    trComp[5][ind]);
    }
  }
  //
  delete[] clIndL[0];
  delete[] clIndL[1];
}


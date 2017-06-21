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
// Class AliAnalysisTaskdNdEtapp13                                            //
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

#include "TFile.h"
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
#include "TGrid.h"

#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
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
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSMultReconstructor.h"

#include "AliLog.h"

#include "AliPhysicsSelection.h"
#include "AliAnalysisTaskdNdEtapp13.h"
#include "AliITSMultRecBg.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliESDtrackCuts.h"

#include "AliPPVsMultUtils.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskdNdEtapp13)


//                                                     0     1     2     3     4     5    6      7      8        9      10          11        12       13    14
const char* AliAnalysisTaskdNdEtapp13::fgCentSelName[] = {"V0M","V0A","V0C","FMD","TRK","TKL","CL0","SPDClusters1","V0MvsFMD","ZNA","TKLvsV0M","ZEMvsZDC","V0A123","V0A0","V0S", "MB", "RefMult08","SPDTracklets08","SPDTracklets08to15", "V0av"};

const char*  AliAnalysisTaskdNdEtapp13::fgkPDGNames[] = {
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

const int AliAnalysisTaskdNdEtapp13::fgkPDGCodes[] = {
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

TH1 *hMCcalib1 = NULL;

//________________________________________________________________________
/*//Default constructor
AliAnalysisTaskdNdEtapp13::AliAnalysisTaskdNdEtapp13(const char *name)
: AliAnalysisTaskSE(name),
*/
//________________________________________________________________________
AliAnalysisTaskdNdEtapp13::AliAnalysisTaskdNdEtapp13(const char *name)
: AliAnalysisTaskSE(name),
//
fOutput(0),
//
fDoNormalReco(kFALSE),
fDoInjection(kFALSE),
fDoRotation(kFALSE),
//
fUseMC(kFALSE),
fCheckReconstructables(kFALSE),
//
fHistosTrData(0),
fHistosTrInj(0),
fHistosTrRot(0),
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
fEtaBin(0.1),
fZVertexMin(-20),
fZVertexMax( 20),
fZVBin(1.),
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
fMCV0Scale(1.),
//
fTrigSel(AliVEvent::kINT7),
//
fMultReco(0),
fRPTree(0),
fStack(0),
fMCEvent(0),
//
fNPrimMCeta2(0),
fNTreta2(0),
fNPart(0),
fNBColl(0),
fCurrCentBin(-1),
fNCentBins(0),
fCentPerc(),
fUseCentralityVar("V0M"),
fIsSelected(kFALSE),
fVtxOK(kFALSE),
fUseSpecialOutput(kFALSE),
fUseBCMod(kFALSE),
fBCMod4(2),
fCutOnPhi(kFALSE),
fCalibfilePath(""),
fHcalib3(NULL),
//
fWeight(1.),
fPPVsMultUtils(new AliPPVsMultUtils())
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
  fCentPerc.Set(2);
  fCentPerc[0] = 0.0;
  fCentPerc[1] = 100.0;
  //
}

//________________________________________________________________________
AliAnalysisTaskdNdEtapp13::~AliAnalysisTaskdNdEtapp13()
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
  delete fHistosTrRcblPrim;
  delete fHistosTrRcblSec;
  delete fHistosCustom;
  //
  delete fHcalib3;
}

//________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::UserCreateOutputObjects()
{
  //
  if (fUseSpecialOutput) OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  //
  //
  Bool_t needGeom = GetDoNormalReco() || GetDoInjection() || GetDoRotation();
  if (needGeom) {
    AliCDBManager *man = AliCDBManager::Instance();
    if (fUseMC) {
      man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",188359);
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",188359,-1,-1)) AliFatal("Failed to misalign geometry");
    }
    else {
      man->SetDefaultStorage("alien://Folder=/alice/data/2012/OCDB"); //man->SetRun(137366);
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",188359);
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject()); //!!!
      man->UnsetDefaultStorage();
      man->SetDefaultStorage("alien://Folder=/alice/cern.ch/user/s/shahoian/ITS_fixSPD72");

      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",188359,-1,-1)) AliFatal("Failed to misalign geometry");
    }
  }
  //
  // Create histograms
  //---------------------------------------------Standard histos per tracklet type--->>
  UInt_t hPattern = 0xffffffff;
  fHistosTrData                      = BookHistosSet("TrData",hPattern);
  //  hPattern &= ~(BIT(kHEtaZvSPD1));  // fill single clusters for "data" only
  if (GetDoInjection()) fHistosTrInj = BookHistosSet("TrInj",hPattern);
  if (GetDoRotation())  fHistosTrRot = BookHistosSet("TrRot",hPattern);
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
  AliInfo(Form("Centrality type selected: %s\n",fUseCentralityVar.Data()));
  PostData(1, fOutput);
  //
}

//________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::UserExec(Option_t *)
{
  // Main loop
  //
  static Bool_t statOK = kFALSE;
  if (!statOK) {
    RegisterStat();
    statOK = kTRUE;
  }
  //
  //  printf("Do: %d %d %d\n",GetDoNormalReco(),GetDoInjection(),GetDoRotation());
  Bool_t needRecPoints = GetDoNormalReco() || GetDoInjection() || GetDoRotation();
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  fRPTree = 0;
  AliESDInputHandler *handler = (AliESDInputHandler*)anMan->GetInputEventHandler();
  AliESDInputHandlerRP *handRP = 0;
  if (needRecPoints) {
    handRP = (AliESDInputHandlerRP*)handler;
    if (!handRP) { printf("No RP handler\n"); return; }
  }
  AliESDEvent *esd  = (AliESDEvent*)handler->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }
  //
  //
  if (needRecPoints) {
    fRPTree = handRP->GetTreeR("ITS");
    if (!fRPTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  }
  //
  // do we need to initialize the field?
  AliMagF* field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!field && !esd->InitMagneticField()) {printf("Failed to initialize the B field\n");return;}
  //
  AliMCEventHandler* eventHandler = 0;
  fMCEvent = 0;
  fStack = 0;
  //
  // ========================== MC EVENT INFO ====================================>>>
  if (fUseMC) {
    eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
    if (!eventHandler) { printf("ERROR: Could not retrieve MC event handler\n"); return; }
    fMCEvent = eventHandler->MCEvent();
    if (!fMCEvent) { printf("ERROR: Could not retrieve MC event\n"); return; }
    fStack = fMCEvent->Stack();
    if (!fStack) { printf("Stack not available\n"); return; }
  }
  //
  // =============================================================================>>>
  //

  // =============================================================================>>>

  const AliMultiplicity* mult = esd->GetMultiplicity();
  if (!mult) {
    AliDebug(AliLog::kError, "AliMultiplicity not available");
    return;
  }

  /* IMPLEMENT INEL > 0 */
/*
  Float_t totalNch = 0;
  Float_t totalNtr = 0;
  if (fUseMC) {
    for (Int_t i = 0; i < fStack->GetNtrack(); i++) {
      if (!fStack->IsPhysicalPrimary(i))
      continue;
      AliMCParticle* particle = (AliMCParticle*)fMCEvent->GetTrack(i);
      if (!particle)
      continue;
      if (particle->Charge() == 0)
      continue;
      if (TMath::Abs(particle->Eta()) < 1.)
      //printf(" eta =  %f \n", particle->Eta());
      totalNch++;
    }
    if (totalNch < 1) return;
  }

  Bool_t inelgt0 = kFALSE;
  for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i) {
    if (TMath::Abs(mult->GetEta(i)) < 1.) {
      totalNtr++;

      inelgt0 = kTRUE;
    }
  }



  if (fUseMC) {
    // fill matrix
    ((TH2F *)fHistosCustom->UncheckedAt(kHEffMatrix))->Fill(totalNch, totalNtr);
  }

*/


  //
  // MC Generator info
  AliGenEventHeader* mcGenH = 0;
  AliGenHijingEventHeader* hHijing=0;
  AliGenDPMjetEventHeader* hDpmJet=0;
  AliGenPythiaEventHeader* hPythia=0;
  fNPart  = 0;
  fNBColl = 0;
  Int_t processType = 0;
  Float_t processWeight = 1.;
  if (fUseMC) {
    mcGenH = fMCEvent->GenEventHeader();
    //
    if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      TIter next(headers);
      AliGenEventHeader* mcGenH1 = 0;
      while ( (mcGenH1=(AliGenEventHeader*)next()) ) {
        if (mcGenH1->InheritsFrom(AliGenHijingEventHeader::Class())) {
          hHijing = (AliGenHijingEventHeader*)mcGenH1;
          break;
        }
        if (mcGenH1->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
          hDpmJet = (AliGenDPMjetEventHeader*)mcGenH1;
          break;
        }
      }
    }
    //
    if ( hHijing || mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())) {
      if (!hHijing) hHijing = (AliGenHijingEventHeader*)mcGenH;
      fNPart  = (hHijing->ProjectileParticipants()+hHijing->TargetParticipants())/2.;
      fNBColl = hHijing->NN()+hHijing->NNw()+hHijing->NwN()+hHijing->NwNw();
    }
    else if ( hDpmJet || mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
      if (!hDpmJet) hDpmJet = (AliGenDPMjetEventHeader*)mcGenH;
      fNPart  = (hDpmJet->ProjectileParticipants()+hDpmJet->TargetParticipants())/2.;
      fNBColl = hDpmJet->NN()+hDpmJet->NNw()+hDpmJet->NwN()+hDpmJet->NwNw();

      /*
      // Special for SD rejection
      int npProj = hDpmJet->ProjectileParticipants(); // A
      //      int npTarg = hDpmJet->TargetParticipants(); // p
      int nsd1,nsd2,ndd;
      hDpmJet->GetNDiffractive(nsd1,nsd2,ndd);
      if (ndd==0 && (npProj==nsd1 || npProj==nsd2)) {
      hstat->Fill(kNEvSDMC);
      return; // reject SD
    }
    */

    Double_t fsd = 0.13;
    Double_t fdd = 0.08;
    switch (hDpmJet->ProcessType()) {
      case 5: case 6: // single-diffractive
      processType = 1;
      processWeight = 1. / 0.138 * fsd;
      break;
      case 7: // double-diffractive
      processType = 2;
      processWeight = 1. / 0.050 * fdd;
      break;
      default: // the rest
      processType = 0;
      processWeight = 1. / (1. - 0.138 - 0.050) * (1. - fsd - fdd);
      break;
    }

  }
  //
  else if ( hPythia || mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())) {
    if (!hPythia) hPythia = (AliGenPythiaEventHeader *)mcGenH;
    Double_t fsd = 0.193;
    Double_t fdd = 0.130;
    switch (hPythia->ProcessType()) {
      case 92: case 93: // single-diffractive
      processType = 1;
      processWeight = 1. / 0.193 * fsd;
      break;
      case 94: // double-diffractive
      processType = 2;
      processWeight = 1. / 0.130 * fdd;
      break;
      default: // the rest
      processType = 0;
      processWeight = 1. / 0.677 * (1. - fsd - fdd);
      break;
    }
  }
  else {} // unknown generator
  //
  //
  TArrayF vtmc(3);
  mcGenH->PrimaryVertex(vtmc);
  for (int i=3;i--;) fVtxMC[i] = vtmc[i];
}

TH1* hstat = (TH1*)fHistosCustom->UncheckedAt(kHStat);
hstat->Fill(kEvTot0, fWeight); // RS
// =============================================================================>>>
//
Int_t trg = handler->IsEventSelected();
// FIXME: add event selection histo, can probably adapt hstat
Bool_t isPileupfromSPD = esd->IsPileupFromSPD(3, 0.8);
Bool_t pileup = esd->IsPileupFromSPDInMultBins();
fIsSelected = trg & fTrigSel;
if(fIsSelected) hstat->Fill(kEvAfterPhysSel, fWeight);
fIsSelected &= fPPVsMultUtils->IsNotPileupSPDInMultBins(esd);
//fIsSelected &= isPileupfromSPD;
fIsSelected &= fPPVsMultUtils->IsINELgtZERO(esd);
if(fIsSelected) hstat->Fill(kEvAfterPileUp);


// Clusters vs trakclets
Int_t fNofTracklets = esd->GetMultiplicity()->GetNumberOfTracklets();
Int_t fNofITSClusters0 = esd->GetNumberOfITSClusters(0);
Int_t fNofITSClusters1 = esd->GetNumberOfITSClusters(1);
if ( fNofITSClusters0+ fNofITSClusters1 > 65+4*fNofTracklets) fIsSelected = kFALSE;
if(fIsSelected) hstat->Fill(kEvAfterClsVsTrk, fWeight);
// hm tail in V0
AliVVZERO* vzero = fInputEvent->GetVZEROData();
Double_t v0c012 = vzero->GetMRingV0C(0) + vzero->GetMRingV0C(1) + vzero->GetMRingV0C(2);
Double_t v0c3   = vzero->GetMRingV0C(3);

if (fUseBCMod & !fUseMC){
  Int_t  fBC = fInputEvent->GetBunchCrossNumber();
  Int_t  bcmod4 = fBC%4;

  fIsSelected &= (bcmod4 == fBCMod4);
  //if (bcmod4 ==2) fIsSelected &= (bcmod4 == 2);
  //else if (bcmod4 == 0)
}



//fIsSelected &= vzero->GetMTotV0C() < (330. + 100. * TMath::Power(vzero->GetMTotV0A(), .2));
//fIsSelected &= (v0c012 < 160.) || (v0c3 > 12.*TMath::Power(.01*(v0c012 - 160.), 1.7));
if(fIsSelected) hstat->Fill(kEvAfterAsymCut);
/*
if (!fUseMC) {
TString trigStr(esd->GetFiredTriggerClasses());
if (!trigStr.Contains("CINT5-B-")) return;
}
AliVVZERO* esdVZERO = esd->GetVZEROData();
fIsSelected = !((esdVZERO->GetV0ADecision()!=1) || (esdVZERO->GetV0CDecision()!=1));
*/
/*
fIsSelected = kTRUE;
AliVVZERO* esdVZERO = esd->GetVZEROData();
if ((esdVZERO->GetV0ADecision()!=1) || (esdVZERO->GetV0CDecision()==2)) fIsSelected = kFALSE;
if (!fUseMC) {
TString trigStr(esd->GetFiredTriggerClasses());
if (!trigStr.Contains("CINT5-B-")) fIsSelected = kFALSE;
}
*/
//
const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
const AliESDVertex* vtxESD_ = esd->GetPrimaryVertex();
const AliMultiplicity* multESD = esd->GetMultiplicity();

// MF: mult classes using new framework
AliMultSelection *MultSelection = (AliMultSelection*) esd -> FindListObject("MultSelection");
static Bool_t skipCentrality = (fUseCentralityVar == "MB"); // If the centrality var is MB, the centrality selection is skept
Double_t centPercentile = -1;
if(!skipCentrality) {
  if (fUseMC) {
    // TFile *mcCalib1 = NULL;
    //
    // if (hPythia) {
    //   TString str1 = handler->GetTree()->GetCurrentFile()->GetName();
    //
    //   if (str1.Contains("LHC15g3c3")) {
    //     printf("loading PYTHIA6 Perugia-0 V0M calibration (LHC15g3c3)\n");
    //
    //     if (TString(fCalibfilePath).BeginsWith("alien:")) {
    //       TGrid::Connect("alien:");
    //       mcCalib1 = TFile::Open(fCalibfilePath);
    //     }
    //     else    mcCalib1 = TFile::Open(fCalibfilePath);
    //
    //   }
    //   else if (str1.Contains("LHC15g3a3")) {
    //     printf("loading PYTHIA8 Monash V0M calibration (LHC15g3a3)\n");
    //
    //     if (TString(fCalibfilePath).BeginsWith("alien:")) {
    //       TGrid::Connect("alien:");
    //       mcCalib1 = TFile::Open(fCalibfilePath);
    //     }
    //     else   mcCalib1 = TFile::Open(fCalibfilePath);
    //   }
    // }
    //
    // else {
    //   printf("loading EPOS-LHC V0M calibration (LHC16d3)\n");
    //   if (TString(fCalibfilePath).BeginsWith("alien:")) {
    //     TGrid::Connect("alien:");
    //     mcCalib1 = TFile::Open(fCalibfilePath);
    //   }
    //   else    mcCalib1 = TFile::Open(fCalibfilePath);
    // }
    //
    // hMCcalib1 = (TH1 *)mcCalib1->Get("h3");

    Float_t multV01=0;
    AliESDVZERO* esdV01 = esd->GetVZEROData();
    multV01 = esdV01->GetMTotV0A()+esdV01->GetMTotV0C();
    multV01 *= fMCV0Scale;

    //    centPercentile = hMCcalib1->Interpolate(multV01);

    centPercentile = fHcalib3->Interpolate(multV01);
    //  printf(" centPercentile =  %f \n", centPercentile);

    fCurrCentBin = GetCentralityBin(centPercentile);
  }
  else  { centPercentile = MultSelection->GetMultiplicityPercentile(fUseCentralityVar);
    fCurrCentBin = GetCentralityBin(centPercentile);
  }

} else {
  fCurrCentBin = 0;
}

((TH1*)fHistosCustom->UncheckedAt(kHCentDistNoSel))->Fill(centPercentile, fWeight);
if (fIsSelected)
((TH1*)fHistosCustom->UncheckedAt(kHCentDistTrig))->Fill(centPercentile, fWeight);
//
/*
if (centPercentile<1 || fIsSelected) {
printf("Cent%s:%.2f Sel:%d Ntrk:%d SPD1:%d SPD2:%d V0A:%f V0C:%f CentV0A:%.2f CentV0M:%.2f CentCl1:%.2f\n",
fUseCentralityVar.Data(),centPercentile,fIsSelected,multESD->GetNumberOfTracklets(),
multESD->GetNumberOfITSClusters(0),multESD->GetNumberOfITSClusters(1),
esd->GetVZEROData()->GetMTotV0A(),esd->GetVZEROData()->GetMTotV0C(),
centrality->GetCentralityPercentileUnchecked("V0A"),
centrality->GetCentralityPercentileUnchecked("V0M"),
centrality->GetCentralityPercentileUnchecked("CL1"));
}
*/

/*
const double kSafeMargin = 1e-3;
if (centPercentile<-kSafeMargin) return;

if (centPercentile<kSafeMargin)     centPercentile = kSafeMargin;
if (centPercentile>100-kSafeMargin) centPercentile = 100.-kSafeMargin;
*/

//
//
// ==================== STORE SOME GLOBAL INFO FOR ALL EVENTS ==============>>>
Float_t multV0=0;
AliESDVZERO* esdV0 = esd->GetVZEROData();
if (esdV0) {
  multV0 = esdV0->GetMTotV0A()+esdV0->GetMTotV0C();
  if (fUseMC) multV0 *= fMCV0Scale;
}
((TH1*)fHistosCustom->UncheckedAt(kHV0NoSel))->Fill(multV0, fWeight);
if ( trg & fTrigSel)
((TH1*)fHistosCustom->UncheckedAt(kHV0Trig))->Fill(multV0, fWeight);
//
float nSPD2 = multESD->GetNumberOfITSClusters(1);
((TH1*)fHistosCustom->UncheckedAt(kHNClSPD2NoSel))->Fill(nSPD2, fWeight);
//
//------------------------------------------------------
AliESDZDC *esdZDC = esd->GetESDZDC();
float zdcEnergy=0,zemEnergy=0;
if (esdZDC) {
  zdcEnergy = (esdZDC->GetZDCN1Energy() + esdZDC->GetZDCP1Energy() + esdZDC->GetZDCN2Energy()+ esdZDC->GetZDCP2Energy());
  zemEnergy = (esdZDC->GetZDCEMEnergy(0)+ esdZDC->GetZDCEMEnergy(1))/8.;
}
// PutZDCSelection
((TH2*)fHistosCustom->UncheckedAt(kHZDCZEMNoSel))->Fill(zemEnergy,zdcEnergy, fWeight);
//
// ==================== STORE SOME GLOBAL INFO FOR ALL EVENTS ==============<<<
//
//  printf("CentPerc: %f : Bin %d ZDC: %f ZEM: %f\n",centPercentile, fCurrCentBin, zdcEnergy,zemEnergy);
if (fCurrCentBin<0) {
  //printf("Reject: %.1f : V0:%.1f V0Cor:%.1f V0CR:%.1f SPD2c:%.1f\n",mltTst, multV0,multV0Corr,multV0CorrResc,nSPD2Corr);
  return;
}
hstat->Fill(kBinEntries+kEvCentBin + kEntriesPerBin*fCurrCentBin, fWeight);
//
fVtxOK = kFALSE;
for (int i=3;i--;) fESDVtx[i] = 0;
if (vtxESD->GetNContributors()>0) {
  TString vtxTyp = vtxESD->GetTitle();
  /* R+HACK
  if ( !vtxTyp.Contains("vertexer: Z") || (vtxESD->GetDispersion()<0.04 && vtxESD->GetZRes()<0.25)) {
  fVtxOK = kTRUE;
  fESDVtx[0] = vtxESD->GetX();
  fESDVtx[1] = vtxESD->GetY();
  fESDVtx[2] = vtxESD->GetZ();
}
*/

if ( !vtxESD->IsFromVertexerZ() || (vtxESD->GetDispersion()<0.02)) {
  fVtxOK = kTRUE;
  fESDVtx[0] = vtxESD_->GetX(); // R+HACK FIXME: what's this?
  fESDVtx[1] = vtxESD_->GetY(); // R+HACK
  fESDVtx[2] = vtxESD_->GetZ(); // R+HACK
}

}
//
if (fIsSelected) hstat->Fill(kBinEntries+kEvPassPS + kEntriesPerBin*fCurrCentBin, fWeight);
//
if (fVtxOK && fIsSelected) {
  ((TH1*)fHistosCustom->UncheckedAt(kHZVtxNoSel))->Fill(fESDVtx[2], fWeight);
  hstat->Fill(kBinEntries+kEvPassVtx + kEntriesPerBin*fCurrCentBin, fWeight);
}
//
fVtxOK &= (fESDVtx[2] >= fZVertexMin && fESDVtx[2] <= fZVertexMax);
//
//  if (!fVtxOK || !fIsSelected) return;




  Float_t totalNch = 0;
  Float_t totalNtr = 0;
  if (fUseMC && (fVtxMC[2] >= fZVertexMin && fVtxMC[2] <= fZVertexMax)) {
    for (Int_t i = 0; i < fStack->GetNtrack(); i++) {
      if (!fStack->IsPhysicalPrimary(i))
      continue;
      AliMCParticle* particle = (AliMCParticle*)fMCEvent->GetTrack(i);
      if (!particle)
      continue;
      if (particle->Charge() == 0)
      continue;
      if (TMath::Abs(particle->Eta()) < 1.)
      //printf(" eta =  %f \n", particle->Eta());
      totalNch++;
    }
    if (totalNch < 1) return;
  }

  Bool_t inelgt0 = kFALSE;
    if (fVtxOK && fIsSelected){
  for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i) {
    if (TMath::Abs(mult->GetEta(i)) < 1.) {
      totalNtr++;
      inelgt0 = kTRUE;
    }
  }
}


  if (fUseMC) {
    // fill matrix
    ((TH2F *)fHistosCustom->UncheckedAt(kHEffMatrix))->Fill(totalNch, totalNtr);
}




if (fIsSelected) ((TH2F *)fHistosCustom->UncheckedAt(kHCorrMatrixSel+fCurrCentBin))->Fill(totalNch, totalNtr);
((TH2F *)fHistosCustom->UncheckedAt(kHCorrMatrix+fCurrCentBin))->Fill(totalNch, totalNtr);

Float_t totalNtr2 = 0.0;
if (fVtxOK) {
for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i) {
  if (TMath::Abs(mult->GetEta(i)) < 1.) {
    totalNtr2++;
  }
}
}

 ((TH2F *)fHistosCustom->UncheckedAt(kHCorrMatrixSel2+fCurrCentBin))->Fill(totalNch, totalNtr2); // response matrix for trigger efficiency.



if (fUseMC) {
  FillMCPrimaries();
  ((TH2*)fHistosCustom->UncheckedAt(kHProcessMCNoPhSel))->Fill(processType,fCurrCentBin, fWeight);
  ((TH2*)fHistosCustom->UncheckedAt(kHTotalNchNoPhSel))->Fill(totalNch,fCurrCentBin, fWeight);
  ((TH2*)fHistosCustom->UncheckedAt(kHZVtxMCNoPhSel))->Fill(fVtxMC[2],fCurrCentBin, fWeight);
  if (fIsSelected) ((TH2*)fHistosCustom->UncheckedAt(kHZVtxMCNoVtSel))->Fill(fVtxMC[2],fCurrCentBin, fWeight);

  if ((fVtxMC[2] < fZVertexMin || fVtxMC[2] > fZVertexMax)) hstat->Fill(kBinEntries+kEvInMltBin + kEntriesPerBin*fCurrCentBin, fWeight);
  //
}
if (!fVtxOK || !fIsSelected) return;
//
// ===================== Store multiplicity estimators ===============================
double etaRange = TMath::Abs(fEtaMax-fEtaMin)/2.;
double mltE0 = AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTrackletsITSTPC, etaRange);
double mltE1 = AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTrackletsITSSA , etaRange);
double mltE2 = AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTracklets, etaRange);
((TH1F*)fHistosCustom->UncheckedAt(kHMltEstTrITSTPC))->Fill(mltE0, fWeight);
((TH1F*)fHistosCustom->UncheckedAt(kHMltEstTrITSSA))->Fill(mltE1, fWeight);
((TH1F*)fHistosCustom->UncheckedAt(kHMltEstTr))->Fill(mltE2, fWeight);
//
// ===================== Process event passing all selections ===============================
//
((TH1*)fHistosCustom->UncheckedAt(kHCentDist))->Fill(centPercentile, fWeight);
hstat->Fill(kEvTot, fWeight); // RS
((TH1*)fHistosCustom->UncheckedAt(kHStatCentBin))->Fill(fCurrCentBin, fWeight);
((TH1*)fHistosCustom->UncheckedAt(kHStatCent))->Fill(centPercentile, fWeight);
//
// register Ntracklets and ZVertex of the event
if (fUseMC) ((TH2*)fHistosCustom->UncheckedAt(kHZVtxMC))->Fill(fVtxMC[2],fCurrCentBin, fWeight);
((TH2*)fHistosCustom->UncheckedAt(kHZVtx))->Fill(fESDVtx[2],fCurrCentBin, fWeight);
((TH2*)fHistosCustom->UncheckedAt(kHV0))->Fill(multV0,fCurrCentBin, fWeight);
((TH2*)fHistosCustom->UncheckedAt(kHNClSPD2))->Fill(nSPD2,fCurrCentBin, fWeight);
((TH3*)fHistosCustom->UncheckedAt(kHZDCZEM))->Fill(zemEnergy,zdcEnergy,fCurrCentBin, fWeight);
//
// normal reconstruction
hstat->Fill(kBinEntries+kEvProcData + kEntriesPerBin*fCurrCentBin, fWeight);
//
if (GetDoNormalReco() || GetDoInjection()) { // for the injection the normal reco should be done
  InitMultReco();
  fMultReco->Run(fRPTree, fESDVtx);
  printf("Multiplicity Reconstructed:\n");
  AliMultiplicity* mlt = fMultReco->GetMultiplicity();
  if (mlt) mlt->Print();
  if (GetDoNormalReco()) FillHistos(kData,mlt);
  FillClusterInfo();
  //
}
if (!GetDoNormalReco()) {
  FillHistos(kData,multESD); // fill data histos from ESD
  FillClusterInfoFromMult(multESD, fESDVtx[2] );
}
//
// Injection: it must come right after the normal reco since needs its results
if (GetDoInjection()) {
  if (!fMultReco) InitMultReco(); // in principle, not needed, the reco is created above
  fMultReco->SetRecType(AliITSMultRecBg::kBgInj);
  fMultReco->Run(fRPTree, fESDVtx);
  printf("Multiplicity from Injection:\n");
  AliMultiplicity* mlt = fMultReco->GetMultiplicity();
  if (mlt) mlt->Print();
  hstat->Fill(kBinEntries + kEvProcInj + kEntriesPerBin*fCurrCentBin, fWeight);
  FillHistos(kBgInj,mlt);
}
//
// Rotation
if (GetDoRotation()) {
  InitMultReco();
  fMultReco->SetRecType(AliITSMultRecBg::kBgRot);
  fMultReco->SetPhiRotationAngle(fPhiRot);
  fMultReco->Run(fRPTree, fESDVtx);
  printf("Multiplicity from Rotation:\n");
  AliMultiplicity* mlt = fMultReco->GetMultiplicity();
  if (mlt) mlt->Print();
  hstat->Fill(kBinEntries + kEvProcRot + kEntriesPerBin*fCurrCentBin, fWeight);
  FillHistos(kBgRot,mlt);
}
//
// =============================================================================<<<
//
if (fMultReco) delete fMultReco;
fMultReco = 0;
//
}


//________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::RegisterStat()
{
  TH1* hstat;
  TList *lst = dynamic_cast<TList*>(GetOutputData(1));
  if (lst && (hstat=(TH1*)lst->FindObject("hStat"))) {
    Info("Terminate","Registering used settings");
    // fill used settings
    hstat->Fill(kOneUnit,1.);
    hstat->Fill(kCentVar,1);
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
}


//________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::Terminate(Option_t *)
{
  Printf("Terminating...");
  RegisterStat();
}


//_________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::InitMultReco()
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
TObjArray* AliAnalysisTaskdNdEtapp13::BookCustomHistos()
{
  // book custom histos, not related to specific tracklet type
  TObjArray* histos = new TObjArray();
  TH1F* hstat;
  //
  // ------------ job parameters, statistics ------------------------------>>>
  int nbs = kBinEntries + fNCentBins*kEntriesPerBin;
  hstat = new TH1F("hStat","Run statistics",nbs,0.5,nbs+0.5);
  //
  hstat->GetXaxis()->SetBinLabel(kEvTot0, "Ev.Tot0");

  hstat->GetXaxis()->SetBinLabel(kEvAfterPhysSel , "Ev. After PS");
  hstat->GetXaxis()->SetBinLabel(kEvAfterPileUp  , "Ev. After Pup");
  hstat->GetXaxis()->SetBinLabel(kEvAfterClsVsTrk, "Ev. After ClsVsTrk");
  hstat->GetXaxis()->SetBinLabel(kEvAfterAsymCut , "Ev. After V0 Asym");

  hstat->GetXaxis()->SetBinLabel(kEvTot, "Ev.Tot After All Evts. Sel");

  hstat->GetXaxis()->SetBinLabel(kOneUnit,"ScaleMerge");
  hstat->GetXaxis()->SetBinLabel(kNWorkers,"Workers");
  //
  hstat->GetXaxis()->SetBinLabel(kCentVar,  fUseCentralityVar.Data());
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
  hstat->GetXaxis()->SetBinLabel(kNEvSDMC,"MC PureSD");
  //
  for (int i=0;i<fNCentBins;i++) {
    TString bnt = "b"; bnt+= i;
    int offs = kBinEntries + i*kEntriesPerBin;
    hstat->GetXaxis()->SetBinLabel(offs + kEvInMltBin, bnt+" Ev.In.CntBn");
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcData, bnt+" Ev.ProcData");
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcInj,  bnt+" Ev.ProcInj");
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcRot,  bnt+" Ev.ProcRot");
    hstat->GetXaxis()->SetBinLabel(offs + kEvProcMix,  bnt+" Ev.ProcMix");
    hstat->GetXaxis()->SetBinLabel(offs + kEvCentBin,  bnt+" Ev.CentBinNS");
    hstat->GetXaxis()->SetBinLabel(offs + kEvPassPS,   bnt+" Ev.CentBinPSTR");
    hstat->GetXaxis()->SetBinLabel(offs + kEvPassVtx,  bnt+" Ev.CentBinVtx");
    //
  }
  //
  hstat->Fill(kNWorkers);
  //
  AddHisto(histos,hstat,kHStat);
  //
  // ------------------------ events per centrality bin ----------------------
  TH1F* hCentAx = new TH1F("EvCentr","Events per centrality",fNCentBins,fCentPerc.GetArray());
  hCentAx->GetXaxis()->SetTitle("Centrality parameter");
  AddHisto(histos,hCentAx,kHStatCent);
  //
  TH1F* hCentBin = new TH1F("EvCentrBin","Events per centrality bin",fNCentBins,-0.5,fNCentBins-0.5);
  hCentBin->GetXaxis()->SetTitle("Centrality Bin");
  AddHisto(histos,hCentBin,kHStatCentBin);
  //
  TH1F* hCentDistNoSel = new TH1F("EvCentrDistNoSel","Events per centrality Before selection",100,0,100);
  hCentDistNoSel->GetXaxis()->SetTitle("Centrality percentile");
  AddHisto(histos,hCentDistNoSel,kHCentDistNoSel);
  //
  TH1F* hCentDistTrig = new TH1F("EvCentrDistTrig","Events per centrality After trigger",100,0,100);
  hCentDistTrig->GetXaxis()->SetTitle("Centrality percentile");
  AddHisto(histos,hCentDistTrig,kHCentDistTrig);
  //
  TH1F* hCentDist = new TH1F("EvCentrDist","Events per centrality ",100,0,100);
  hCentDist->GetXaxis()->SetTitle("Centrality percentile");
  AddHisto(histos,hCentDist,kHCentDist);
  //

  // ------------ job parameters, statistics ------------------------------<<<
  //
  //
  int nZVBinsS = int((fZVertexMax-fZVertexMin)/fZVBin);
  if (nZVBinsS<1) nZVBinsS = 1;
  fZVBin = (fZVertexMax-fZVertexMin)/nZVBinsS; // readjust
  //
  int nZVBins = nZVBinsS;
  double zMn=fZVertexMin, zMx=fZVertexMax;
  const double kZVWide=30.;
  while (zMn>-kZVWide) {zMn -= fZVBin; nZVBins++;}
  while (zMx<kZVWide)  {zMx += fZVBin; nZVBins++;}
  AliInfoF("Wide ZV histo will have %d bins in %.2f<Zv<%.2f range",nZVBins,zMn,zMx);
  //
  // Z vertex distribution for events before selection
  TH1F* hzvns = new  TH1F("zvNoSel","Z vertex before selection",nZVBins,zMn,zMx);
  hzvns->GetXaxis()->SetTitle("Zvertex");
  AddHisto(histos,hzvns,kHZVtxNoSel);
  //
  // V0 for events before selection
  int nbmltV0 = 250;
  double maxmltV0 = 2500*0.4; //RSmod
  //
  TH1F* hnV0ns = new  TH1F("V0NoSel","V0 signal Before Cent. Selection",nbmltV0,0,maxmltV0);
  hnV0ns->GetXaxis()->SetTitle("V0 signal");
  AddHisto(histos,hnV0ns,kHV0NoSel);
  //
  TH1F* hnV0trg = new  TH1F("V0Trig","V0 signal after trigger",nbmltV0,0,maxmltV0);
  hnV0trg->GetXaxis()->SetTitle("V0 signal");
  AddHisto(histos,hnV0trg,kHV0Trig);
  //
  // N SPD2 clusters
  int nbmltSPD2 = 175;
  double maxmltSPD2 = 7000*0.4; // RSmod
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
  // multiplicity distribution histos
  TH1F* hmlt;
  const int kMaxMlt = 300;
  hmlt = new TH1F("mltTrTPCITS","mltTPCITS",kMaxMlt,0,kMaxMlt);
  AddHisto(histos,hmlt,kHMltEstTrITSTPC);
  //
  hmlt = new TH1F("mltTrITSSA","mltITSSA",kMaxMlt,0,kMaxMlt);
  AddHisto(histos,hmlt,kHMltEstTrITSSA);
  //
  hmlt = new TH1F("mltTracklets","mltTracklets",kMaxMlt,0,kMaxMlt);
  AddHisto(histos,hmlt,kHMltEstTr);
  //
  //
  if (fUseMC) {
    hmlt = new TH1F("mltMCTruth","mltMCTruth",kMaxMlt,0,kMaxMlt);
    AddHisto(histos,hmlt,kHMltMC);

    TH2F *heffMatrix = new TH2F("effMatrix", "effMatrix", kMaxMlt, 0, kMaxMlt, kMaxMlt, 0, kMaxMlt);
    AddHisto(histos, heffMatrix, kHEffMatrix);

    //
    TH2F* hprocessMCNoPS = new  TH2F("processMCNoPS","MC Process  Before Ph.Sel per Cent.Bin",3,0,3, fNCentBins, -0.5,fNCentBins-0.5);
    hprocessMCNoPS->GetXaxis()->SetTitle("MCprocess");
    hprocessMCNoPS->GetYaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hprocessMCNoPS,kHProcessMCNoPhSel);
    //
    //
    TH2F* hTotalNchNoPS = new  TH2F("totalNchNoPS","Total Nch  Before Ph.Sel per Cent.Bin",1000,0,1000, fNCentBins, -0.5,fNCentBins-0.5);
    hprocessMCNoPS->GetXaxis()->SetTitle("MCprocess");
    hprocessMCNoPS->GetYaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hTotalNchNoPS,kHTotalNchNoPhSel);
    //

    //
    TH2F* hzvMCNoPS = new  TH2F("zvMCNoPS","Z MC vertex Before Ph.Sel per Cent.Bin",nZVBins,zMn,zMx, fNCentBins, -0.5,fNCentBins-0.5);
    hzvMCNoPS->GetXaxis()->SetTitle("ZvertexMC");
    hzvMCNoPS->GetYaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hzvMCNoPS,kHZVtxMCNoPhSel);
    //
    TH2F* hzvMCNoVS = new  TH2F("zvMCNoVS","Z MC vertex Before Rec.Vtx.Selection per Cent.Bin",nZVBins,zMn,zMx, fNCentBins, -0.5,fNCentBins-0.5);
    hzvMCNoVS->GetXaxis()->SetTitle("ZvertexMC");
    hzvMCNoVS->GetYaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hzvMCNoVS,kHZVtxMCNoVtSel);
    //
    TH2F* hzvMC = new  TH2F("zvMC","Z vertex MC after Selection per Cent.Bin",nZVBins,zMn,zMx, fNCentBins, -0.5,fNCentBins-0.5);
    hzvMC->GetXaxis()->SetTitle("ZvertexMC");
    hzvMC->GetYaxis()->SetTitle("Cent.Bin ID");
    AddHisto(histos,hzvMC,kHZVtxMC);
    //
    TH1* hNPrimEta2All = new TH1F("nPrimEta2All","NPrim |#eta|<2, allMC",1000,0.,1000.);
    hNPrimEta2All->SetXTitle("NPrim");
    AddHisto(histos,hNPrimEta2All,kHNPrimEta2All);
    //
    TH1* hNPrimEta2Vt = new TH1F("nPrimEta2Vt","NPrim |#eta|<2, VtxMC",1000,0.,1000.);
    hNPrimEta2Vt->SetXTitle("NPrim");
    AddHisto(histos,hNPrimEta2Vt,kHNPrimEta2Vt);
    //
    TH1* hNPrimEta2Sel = new TH1F("nPrimEta2Sel","NPrim |#eta|<2, Sel",1000,0.,1000.);
    hNPrimEta2Sel->SetXTitle("NPrim");
    AddHisto(histos,hNPrimEta2Sel,kHNPrimEta2Sel);
    //
    TH1* hNPrimEta2SelVt = new TH1F("nPrimEta2SelVt","NPrim |#eta|<2, Sel+VtxMC",1000,0.,1000.);
    hNPrimEta2SelVt->SetXTitle("NPrim");
    AddHisto(histos,hNPrimEta2SelVt,kHNPrimEta2SelVt);
    //
  }
  //
  TH2F* hcorrEta2 = new TH2F("nTrNprEta2","Ntr vs NPrim |#eta|<2, Sel+Vtx",100,0.,100.,100,0.,100.);
  hcorrEta2->SetXTitle("NPrim");
  hcorrEta2->SetYTitle("NTracklets");
  AddHisto(histos,hcorrEta2,kHNCorrMCEta2);

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
  TH2* hetaphi = new TH2F("etaphiTracklets","etaphiTracklets",50,-2.5,2.5, 200, 0,2*TMath::Pi());
  hetaphi->SetXTitle("#eta");
  hetaphi->SetYTitle("#phi");
  AddHisto(histos,hetaphi,kHEtaPhi);
  //----------------------------------------------------------------------
  int nEtaBinsS = TMath::Nint((fEtaMax-fEtaMin)/fEtaBin);
  if (nEtaBinsS<1) nEtaBinsS = 1;
  //

  char mybuffn[100],mybufft[500];
  for (int ib1=0;ib1<fNCentBins;ib1++) {
    sprintf(mybuffn,"b%d_corrMatrix",ib1);
    sprintf(mybufft,"bin%d Correlation Matrix",ib1);
    TH2F* hcorr3 = new  TH2F(mybuffn,mybufft, kMaxMlt, 0, kMaxMlt, kMaxMlt, 0, kMaxMlt);
    hcorr3->GetXaxis()->SetTitle("n_{ch}");
    hcorr3->GetYaxis()->SetTitle("n_{tracklets}");
    AddHisto(histos,hcorr3,kHCorrMatrix+ib1);
    //
  }


  char mybuffn2[100],mybufft2[500];
  for (int ib2=0;ib2<fNCentBins;ib2++) {
    sprintf(mybuffn2,"b%d_corrMatrixSel",ib2);
    sprintf(mybufft2,"bin%d Correlation Matrix Selected",ib2);
    TH2F* hcorr4 = new  TH2F(mybuffn2,mybufft2, kMaxMlt, 0, kMaxMlt, kMaxMlt, 0, kMaxMlt);
    hcorr4->GetXaxis()->SetTitle("n_{ch}");
    hcorr4->GetYaxis()->SetTitle("n_{tracklets}");
    AddHisto(histos,hcorr4,kHCorrMatrixSel+ib2);
    //
  }

  char mybuffn3[100],mybufft3[500];
  for (int ib3=0;ib3<fNCentBins;ib3++) {
    sprintf(mybuffn3,"b%d_corrMatrixSel2",ib3);
    sprintf(mybufft3,"bin%d Correlation Matrix2 Selected",ib3);
    TH2F* hcorr5 = new  TH2F(mybuffn3,mybufft3, kMaxMlt, 0, kMaxMlt, kMaxMlt, 0, kMaxMlt);
    hcorr5->GetXaxis()->SetTitle("n_{ch}");
    hcorr5->GetYaxis()->SetTitle("n_{tracklets}");
    AddHisto(histos,hcorr5,kHCorrMatrixSel2+ib3);
    //
  }


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
      //
      sprintf(buffn,"b%d_zvrecEtaPrimMCSel",ib);
      sprintf(bufft,"bin%d ZvertexR vs #eta PrimMC sel evs",ib);
      TH2F* hzvretap = new  TH2F(buffn,bufft, nEtaBinsS,fEtaMin,fEtaMax,nZVBinsS,fZVertexMin,fZVertexMax);
      hzvretap->GetXaxis()->SetTitle("#eta");
      hzvretap->GetYaxis()->SetTitle("Zvertex");
      AddHisto(histos,hzvretap,kHZVrEtaPrimMC+ib);
      //
      sprintf(buffn,"b%d_dZV_ZVGen",ib);
      sprintf(bufft,"bin%d ZvRec-ZvGen vs ZvGen",ib);
      TH2F* hdzv = new  TH2F(buffn,bufft,nZVBinsS,fZVertexMin,fZVertexMax, 200,-1.,1);
      hdzv->GetXaxis()->SetTitle("VZGen");
      hdzv->GetYaxis()->SetTitle("VZRec-VZGen");
      AddHisto(histos,hdzv,kHZVResMC+ib);
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
TObjArray* AliAnalysisTaskdNdEtapp13::BookHistosSet(const char* pref, UInt_t selHistos)
{
  // book standard set of histos attaching the pref in front of the name/title
  //
  const int kNDPhiBins = 100;
  const int kNDThtBins = 100;
  int nDistBins = int(fNStdDev)*5;
  //
  int nEtaBins = TMath::Nint((fEtaMax-fEtaMin)/fEtaBin);
  if (nEtaBins<1) nEtaBins = 1;
  //
  int nZVBins = int((fZVertexMax-fZVertexMin)/fZVBin);
  if (nZVBins<1) nZVBins = 1;
  fZVBin = (fZVertexMax-fZVertexMin)/nZVBins; // readjust
  AliInfoF("%d bin in %.2f<Zv<%.2f (step=%.2f) will be defined",nZVBins,fZVertexMin,fZVertexMax,fZVBin);
  //
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
    /*
    if (selHistos & (0x1<<kHEtaZvSPD1) ) {
    sprintf(buffn,"b%d_%s_ZvEtaSPD1",ib,pref);
    sprintf(bufft,"bin%d (%s) Zv vs Eta SPD1 clusters",ib,pref);
    h2 = new TH2F(buffn,bufft,nEtaBins,fEtaMin,fEtaMax, nZVBins, fZVertexMin,fZVertexMax);
    h2->GetXaxis()->SetTitle("#eta");
    h2->GetYaxis()->SetTitle("Zv");
    AddHisto(histos,h2,offs+kHEtaZvSPD1);
  }
  */
  //
  if (selHistos & (0x1<<kHWDvEta) ) {
    sprintf(buffn,"b%d_%s_WDvsEta",ib,pref);
    sprintf(bufft,"bin%d (%s) Wdist vs Eta",ib,pref);
    h2 = new TH2F(buffn,bufft,nEtaBins,fEtaMin,fEtaMax, 2*nDistBins,0,fNStdDev);
    h2->GetXaxis()->SetTitle("#eta");
    h2->GetYaxis()->SetTitle("Wdist");
    AddHisto(histos,h2,offs+kHWDvEta);
  }
  //
}
//
histos->SetOwner(kFALSE);
return histos;
}

//_________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::AddHisto(TObjArray* histos, TObject* h, Int_t at)
{
  // add single histo to the set
  if (at>=0) histos->AddAtAndExpand(h,at);
  else       histos->Add(h);
  fOutput->Add(h);
}

//_________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::FillHistos(Int_t type, const AliMultiplicity* mlt)
{
  // fill histos of given type
  if (!mlt) return;
  //
  TObjArray* histos = 0;
  if      (type == kData)  histos = fHistosTrData;
  else if (type == kBgInj) histos = fHistosTrInj;
  else if (type == kBgRot) histos = fHistosTrRot;

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
  fNTreta2 = 0;
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
double phi  = mlt->GetPhi(itr);
double eta    = -TMath::Log(TMath::Tan(theta/2));
//
double dtheta = mlt->GetDeltaTheta(itr);
double dThetaX = dtheta;
if (fScaleDTBySin2T) {
  double sint   =  TMath::Sin(theta);
  dThetaX /= (sint*sint);
}

//       if (!(phi > 0 && phi < (2*TMath::Pi()/3))) continue;
//       if (!(phi > (2*TMath::Pi()/3) && phi < (4*TMath::Pi()/3))) continue;
//       if (!(phi > (4*TMath::Pi()/3) && phi < (2*TMath::Pi()))) continue;


if (fCutOnPhi) {
  if ((phi > 0.58 && phi < 0.60 ) || (phi > 1.2 && phi < 1.4 ) || (phi > 1.75 && phi < 2.3 ) || (phi > 4.2 && phi < 4.50 ) || (phi > 4.7 && phi < 5.0 ) || (phi > 5.8 )) continue;  // use only fidutial region by taking out data and mc mismatch regions
}

if (fCutOnDThetaX && TMath::Abs(dThetaX)>fDThetaWindow) continue;
//
//    double phi    = mlt->GetPhi(itr);
double dphi   = mlt->GetDeltaPhi(itr);
double dist   = mlt->CalcDist(itr);
double dphiS  = TMath::Abs(dphi) - fDPhiShift;
if (dphi<0) dphiS = -dphiS;
//
if (dist<fNStdCut && dphiS<fDPhiSCut && TMath::Abs(eta)<2) fNTreta2++;
//
//
if (eta<fEtaMin || eta>fEtaMax) continue;
//
Bool_t isSig = FillHistosSet(histos,eta,/*phi,theta,*/dphi,dtheta,dThetaX,dist);
//
if (type==kData && isSig) {
  ((TH2*)fHistosCustom->UncheckedAt(kHEtaPhi))->Fill(eta,mlt->GetPhi(itr), fWeight);
}
//
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
    if (dphiS<fDPhiSCut) FillSpecies(typeMC, lab0);
  }
  if (fCheckReconstructables) CheckReconstructables();
}
}
//
if (type==kData) {
  TH2* hcorr = ((TH2*)fHistosCustom->UncheckedAt(kHNCorrMCEta2));
  if (hcorr) hcorr->Fill(fNPrimMCeta2,fNTreta2, fWeight);
}
//
//-------------------------------------------------------------TMP RS - singles ------->>>
/*
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
*/
//-------------------------------------------------------------TMP RS - singles -------<<<
//
}

//_________________________________________________________________________
void AliAnalysisTaskdNdEtapp13::FillMCPrimaries()
{

  // fill all MC primaries Zv vs Eta
  if (!fUseMC || !fStack || !fMCEvent) return;
  //
  float zv =  fVtxMC[2]; //fVtxOK ? fESDVtx[2] : fVtxMC[2];
  float zvr = fESDVtx[2];
  /*
  const double kSafeZv = 1e-3
  if (zv<fZVertexMin+kSafeZv) zv = fZVertexMin+kSafeZv;
  if (zv>fZVertexMax-kSafeZv) zv = fZVertexMax-kSafeZv;
  */
  //
  int ntr = fStack->GetNtrack();
  TH2* hprimEtaZ  = (TH2F*)fHistosCustom->UncheckedAt(kHZVEtaPrimMC+fCurrCentBin);
  TH2* hprimEtaZr = (TH2F*)fHistosCustom->UncheckedAt(kHZVrEtaPrimMC+fCurrCentBin);
  if (fVtxOK) {
    TH2* hdzv = (TH2F*)fHistosCustom->UncheckedAt(kHZVResMC+fCurrCentBin);
    hdzv->Fill(fVtxMC[2],fESDVtx[2]-fVtxMC[2], fWeight);
  }
  int nprim = 0;
  fNPrimMCeta2 = 0;
  float nprimSel = 0;
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
    //    if (eta<fEtaMin || eta>fEtaMax) continue;
    if (TMath::Abs(eta)<2) fNPrimMCeta2++;
    hprimEtaZ->Fill(eta, zv, fWeight);
    if (fVtxOK && fIsSelected) {
      hprimEtaZr->Fill(eta, zvr, fWeight);
      if (eta>fEtaMin && eta<fEtaMax) nprimSel++;
    }
    nprim++;
  }
  ((TH1F*)fHistosCustom->UncheckedAt(kHNPrimEta2All))->Fill(fNPrimMCeta2, fWeight);
  if (fVtxOK)                ((TH1F*)fHistosCustom->UncheckedAt(kHNPrimEta2Vt))->Fill(fNPrimMCeta2, fWeight);
  if (fIsSelected)           ((TH1F*)fHistosCustom->UncheckedAt(kHNPrimEta2Sel))->Fill(fNPrimMCeta2, fWeight);
  if (fVtxOK && fIsSelected) ((TH1F*)fHistosCustom->UncheckedAt(kHNPrimEta2SelVt))->Fill(fNPrimMCeta2, fWeight);
  //
  if (fIsSelected)           ((TH1F*)fHistosCustom->UncheckedAt(kHMltMC))->Fill(nprimSel, fWeight);
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
Bool_t AliAnalysisTaskdNdEtapp13::FillHistosSet(TObjArray* histos, double eta,
  //double /*phi*/,double /*theta*/,
  double dphi,double dtheta,double dThetaX,
  double dist)
  {
    // fill standard set of histos
    Bool_t res = kFALSE;
    //
    int offs = fCurrCentBin*kNStandardH;
    //
    if (histos->UncheckedAt(kHWDvEta))
    ((TH2*)histos->UncheckedAt(offs+kHWDvEta))->Fill(eta,dist, fWeight);
    //
    if (dist>fNStdDev) return res;
    //
    double dphiS  = TMath::Abs(dphi) - fDPhiShift;
    if (dphi<0) dphiS = -dphiS;
    //
    if (histos->UncheckedAt(offs+kHDPhiSDThetaX))
    ((TH2*)histos->UncheckedAt(offs+kHDPhiSDThetaX))->Fill(dphiS,dThetaX, fWeight);
    //
    if (histos->UncheckedAt(offs+kHDPhiDTheta))
    ((TH2*)histos->UncheckedAt(offs+kHDPhiDTheta))->Fill(dphi,dtheta, fWeight);
    //
    if (histos->UncheckedAt(kHWDist))
    ((TH2*)histos->UncheckedAt(offs+kHWDist))->Fill(dist, fWeight);
    //
    if (dist<fNStdCut && dphiS<fDPhiSCut && histos->UncheckedAt(offs+kHEtaZvCut)) {
      ((TH2*)histos->UncheckedAt(offs+kHEtaZvCut))->Fill(eta,fESDVtx[2], fWeight);
      res = kTRUE;
    }
    //  else {
    //      printf("--->>> troubles filling histos set for %s\n",histos->UncheckedAt(offs+kHEtaZvCut)->GetName());
    //      printf("dist=%f < fNStdCut=%f | dphiS=%f < fDPhiSCut=%f\n", dist, fNStdCut, dphiS, fDPhiSCut);
    //  }
    //
    return res;
  }

  //__________________________________________________________________
  void AliAnalysisTaskdNdEtapp13::FillSpecies(Int_t primsec, Int_t id)
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
  Int_t AliAnalysisTaskdNdEtapp13::GetCentralityBin(Float_t perc) const
  {
    // calculate centrality bin
    //    for (int i=0;i<fNCentBins;i++) if (perc>=fCentPerc[i] && perc<=fCentPerc[i+1]) return i;
    for (int i=0;i<fNCentBins;i++) if (perc>=fCentPerc[i] && perc<fCentPerc[i+1]) return i; // R+FIX //FIXME: I think this is ok, remove comment?
    return -1;
  }

  //_________________________________________________________________________
  Int_t AliAnalysisTaskdNdEtapp13::GetPdgBin(Int_t pdgCode)
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
  Bool_t AliAnalysisTaskdNdEtapp13::HaveCommonParent(const float* clLabs0,const float* clLabs1)
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
  void AliAnalysisTaskdNdEtapp13::CheckReconstructables()
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
      void AliAnalysisTaskdNdEtapp13::FillClusterInfo()
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
              hclU[il]->Fill( clinf[AliITSMultReconstructor::kClZ], clinf[AliITSMultReconstructor::kClPh], fWeight);
            }
          }
          //
          for (int il=0;il<2;il++) for (int ic=fMultReco->GetNClustersLayer(il);ic--;) {
            Float_t *clinf = fMultReco->GetClusterOfLayer(il,ic);
            hclA[il]->Fill( clinf[AliITSMultReconstructor::kClZ], clinf[AliITSMultReconstructor::kClPh], fWeight);
          }
          //
        }

        //_________________________________________________________________________
        void AliAnalysisTaskdNdEtapp13::FillClusterInfoFromMult(const AliMultiplicity* mlt, double zVertex)
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
            if (goodTracklet) hclU->Fill(z,phi, fWeight);
            hclA->Fill(z,phi, fWeight);
          }
          //
          int ncl = mlt->GetNumberOfSingleClusters();
          for (int icl=ncl;icl--;) {
            double phi   = mlt->GetPhiSingle(icl);
            double z     = kRSPD2/TMath::Tan(mlt->GetThetaSingle(icl)) + zVertex;
            hclA->Fill(z,phi, fWeight);
          }
          //
        }

        //______________________________________________
        void AliAnalysisTaskdNdEtapp13::CheckCentralityVar(const char* var)
        {
          int nv = sizeof(fgCentSelName)/sizeof(char*);
          for (int i=0;i<nv;i++) {
            if (!strcmp(var, fgCentSelName[i])) {
              return;
            }
          }
          AliFatalF("Unknown centrality var: %s",var);
        }

        //______________________________________________
        void AliAnalysisTaskdNdEtapp13::SetCentPercentiles(const Float_t *arr, Int_t nbins)
        {
          // set user defined percentiles
          fCentPerc.Set(nbins+1);
          fNCentBins = nbins;
          AliInfoF("Defining %d centrality bins",fNCentBins);
          for (int i=0;i<=nbins;i++) {
            fCentPerc[i] = arr[i];
            if (i<nbins) AliInfoF("CentBin# %.2f-%.2f",arr[i],arr[i+1]);
          }
          //
        }

        //______________________________________________
        void AliAnalysisTaskdNdEtapp13::SetCentPercentiles(const Double_t *arr, Int_t nbins)
        {
          // set user defined percentiles
          fCentPerc.Set(nbins+1);
          fNCentBins = nbins;
          AliInfoF("Defining %d centrality bins",fNCentBins);
          for (int i=0;i<=nbins;i++) {
            fCentPerc[i] = arr[i];
            if (i<nbins) AliInfoF("CentBin# %.2f-%.2f",arr[i],arr[i+1]);
          }
          //
        }

        //______________________________________________
        void AliAnalysisTaskdNdEtapp13::SetCalibHisto(TH1D *calibhisto)
        {
          // set user calibhisto

          fHcalib3 = (TH1D*)calibhisto->Clone();

          if (!fHcalib3) {
            AliError(Form("No calibration histo found")) ;
          }

        }

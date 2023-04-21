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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "AliVEvent.h"
class TTree;
class TParticle;
class TVector3;
class AliESDtrackCuts;
class AliESDAD; // AD
class AliMCEventHandler;
class AliMCEvent;
class AliStack;
class AliPIDResponse;
class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;
#include "AliESDAD.h" //AD
#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <Riostream.h>

// #include "AliLog.h"

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCascadeVertexer.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliStack.h"
#include "AliV0vertexer.h"

#include "AliAODMCParticle.h"
#include "AliAODcascade.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCFContainer.h"
#include "AliESDUtils.h"
#include "AliESDcascade.h"
#include "AliEventCuts.h"
#include "AliGenEventHeader.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"

using std::cout;
using std::endl;
#include "AliAnalysisTaskFlatenicityLambdaK0sMC.h"

class AliAnalysisTaskFlatenicityLambdaK0sMC;

using namespace std;

ClassImp(AliAnalysisTaskFlatenicityLambdaK0sMC)

    AliAnalysisTaskFlatenicityLambdaK0sMC::
        AliAnalysisTaskFlatenicityLambdaK0sMC()
    : AliAnalysisTaskSE(), fEventCuts(0), fESDtrackCuts(0), fESD(0),
      fOutputList(0), hinvmassK0s(0), fPIDResponse(0),
      hinvmassLambda(0), hinvmassAntiLambda(0), hflat(0), fESDpid(0x0),
      treeK0s(0), treeLambda(0), treeAntiLambda(0), invmK0s(0), invpK0s(0),
      invptK0s(0), invyK0s(0), invmLambda(0), invpLambda(0), invptLambda(0),
      invyLambda(0), invmAntiLambda(0), invpAntiLambda(0), invptAntiLambda(0),
      invyAntiLambda(0), flatenicityK0s(0), flatenicityLambda(0),
      flatenicityAntiLambda(0), hinvmassK0sMC(0), hinvmassLambdaMC(0),
      hinvmassAntiLambdaMC(0), hNeventMC(0), treeK0sMC(0), treeLambdaMC(0),
      treeAntiLambdaMC(0), invmK0sMC(0), invpzK0sMC(0), invptK0sMC(0),
      invyK0sMC(0), invmLambdaMC(0), invpzLambdaMC(0), invptLambdaMC(0),
      invyLambdaMC(0), invmAntiLambdaMC(0), invpzAntiLambdaMC(0),
      invptAntiLambdaMC(0), invyAntiLambdaMC(0), flatenicityK0sMC(0),
      flatenicityLambdaMC(0), flatenicityAntiLambdaMC(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityLambdaK0sMC::AliAnalysisTaskFlatenicityLambdaK0sMC(
    const char *name)
    : AliAnalysisTaskSE(name), fEventCuts(0), fESDtrackCuts(0), fESD(0),
      fOutputList(0), hinvmassK0s(0), fPIDResponse(0),
      hinvmassLambda(0), hinvmassAntiLambda(0), hflat(0), fESDpid(0x0),
      treeK0s(0), treeLambda(0), treeAntiLambda(0), invmK0s(0), invpK0s(0),
      invptK0s(0), invyK0s(0), invmLambda(0), invpLambda(0), invptLambda(0),
      invyLambda(0), invmAntiLambda(0), invpAntiLambda(0), invptAntiLambda(0),
      invyAntiLambda(0), flatenicityK0s(0), flatenicityLambda(0),
      flatenicityAntiLambda(0), hinvmassK0sMC(0), hinvmassLambdaMC(0),
      hinvmassAntiLambdaMC(0), hNeventMC(0), treeK0sMC(0), treeLambdaMC(0),
      treeAntiLambdaMC(0), invmK0sMC(0), invpzK0sMC(0), invptK0sMC(0),
      invyK0sMC(0), invmLambdaMC(0), invpzLambdaMC(0), invptLambdaMC(0),
      invyLambdaMC(0), invmAntiLambdaMC(0), invpzAntiLambdaMC(0),
      invptAntiLambdaMC(0), invyAntiLambdaMC(0), flatenicityK0sMC(0),
      flatenicityLambdaMC(0), flatenicityAntiLambdaMC(0)
{

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityLambdaK0sMC::
    ~AliAnalysisTaskFlatenicityLambdaK0sMC()
{
  // destructor
  if (fOutputList)
  {
    delete fOutputList; // at the end of your task, it is deleted from memory by
                        // calling this function
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityLambdaK0sMC::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  hinvmassK0s = new TH1F("hinvmassK0s", "hinvmassK0s", 100, 0.41, 0.58);
  hinvmassLambda =
      new TH1F("hinvmassLambda", "hinvmassLambda", 100, 1.08, 1.15);
  hinvmassAntiLambda =
      new TH1F("hinvmassAntiLambda", "hinvmassAntiLambda", 100, 1.08, 1.15);
  hflat = new TH1F("hflat", "hflat", 101, -0.005, 1.005);

  treeK0s = new TTree("treeK0s", "treeK0s");
  treeK0s->Branch("invmK0s", &invmK0s, "invmK0s/F");
  treeK0s->Branch("invpK0s", &invpK0s, "invpK0s/F");
  treeK0s->Branch("invptK0s", &invptK0s, "invptK0s/F");
  treeK0s->Branch("invyK0s", &invyK0s, "invyK0s/F");
  treeK0s->Branch("flatenicityK0s", &flatenicityK0s, "flatenicityK0s/F");

  treeLambda = new TTree("treeLambda", "treeLambda");
  treeLambda->Branch("invmLambda", &invmLambda, "invmLambda/F");
  treeLambda->Branch("invpLambda", &invpLambda, "invpLambda/F");
  treeLambda->Branch("invptLambda", &invptLambda, "invptLambda/F");
  treeLambda->Branch("invyLambda", &invyLambda, "invyLambda/F");
  treeLambda->Branch("flatenicityLambda", &flatenicityLambda,
                     "flatenicityLambda/F");

  treeAntiLambda = new TTree("treeAntiLambda", "treeAntiLambda");
  treeAntiLambda->Branch("invmAntiLambda", &invmAntiLambda, "invmAntiLambda/F");
  treeAntiLambda->Branch("invpAntiLambda", &invpAntiLambda, "invpAntiLambda/F");
  treeAntiLambda->Branch("invptAntiLambda", &invptAntiLambda,
                         "invptAntiLambda/F");
  treeAntiLambda->Branch("invyAntiLambda", &invyAntiLambda, "invyAntiLambda/F");
  treeAntiLambda->Branch("flatenicityAntiLambda", &flatenicityAntiLambda,
                         "flatenicityAntiLambda/F");

  fOutputList->Add(hinvmassK0s);
  fOutputList->Add(hinvmassLambda);
  fOutputList->Add(hinvmassAntiLambda);
  fOutputList->Add(hflat);
  fOutputList->Add(treeK0s);
  fOutputList->Add(treeLambda);
  fOutputList->Add(treeAntiLambda);

  // hinvmassK0sMC = new TH1F("hinvmassK0sMC", "hinvmassK0sMC", 100, 0.41, 0.58);
  // hinvmassLambdaMC =
  //     new TH1F("hinvmassLambdaMC", "hinvmassLambdaMC", 100, 1.08, 1.15);
  // hinvmassAntiLambdaMC =
  //     new TH1F("hinvmassAntiLambdaMC", "hinvmassAntiLambdaMC", 100, 1.08, 1.15);
  hNeventMC = new TH1F("hNeventMC", "hNeventMC", 101, -0.005, 1.005);

  treeK0sMC = new TTree("treeK0sMC", "treeK0sMC");
  // treeK0sMC->Branch("invmK0sMC", &invmK0sMC, "invmK0sMC/F");
  treeK0sMC->Branch("invpzK0sMC", &invpzK0sMC, "invpzK0sMC/F");
  treeK0sMC->Branch("invptK0sMC", &invptK0sMC, "invptK0sMC/F");
  treeK0sMC->Branch("invyK0sMC", &invyK0sMC, "invyK0sMC/F");
  // treeK0sMC->Branch("flatenicityK0sMC", &flatenicityK0sMC,
  //                 "flatenicityK0sMC/F");

  treeLambdaMC = new TTree("treeLambdaMC", "treeLambdaMC");
  // treeLambdaMC->Branch("invmLambdaMC", &invmLambdaMC, "invmLambdaMC/F");
  treeLambdaMC->Branch("invpzLambdaMC", &invpzLambdaMC, "invpzLambdaMC/F");
  treeLambdaMC->Branch("invptLambdaMC", &invptLambdaMC, "invptLambdaMC/F");
  treeLambdaMC->Branch("invyLambdaMC", &invyLambdaMC, "invyLambdaMC/F");
  // treeLambdaMC->Branch("flatenicityLambdaMC", &flatenicityLambdaMC,
  //                     "flatenicityLambdaMC/F");

  treeAntiLambdaMC = new TTree("treeAntiLambdaMC", "treeAntiLambdaMC");
  // treeAntiLambdaMC->Branch("invmAntiLambdaMC", &invmAntiLambdaMC,
  //                          "invmAntiLambdaMC/F");
  treeAntiLambdaMC->Branch("invpzAntiLambdaMC", &invpzAntiLambdaMC,
                           "invpzAntiLambdaMC/F");
  treeAntiLambdaMC->Branch("invptAntiLambdaMC", &invptAntiLambdaMC,
                           "invptAntiLambdaMC/F");
  treeAntiLambdaMC->Branch("invyAntiLambdaMC", &invyAntiLambdaMC,
                           "invyAntiLambdaMC/F");
  // treeAntiLambdaMC->Branch("flatenicityAntiLambdaMC", &flatenicityAntiLambdaMC,
  //                        "flatenicityAntiLambdaMC/F");

  // fOutputList->Add(hinvmassK0sMC);
  // fOutputList->Add(hinvmassLambdaMC);
  // fOutputList->Add(hinvmassAntiLambdaMC);
  fOutputList->Add(hNeventMC);
  fOutputList->Add(treeK0sMC);
  fOutputList->Add(treeLambdaMC);
  fOutputList->Add(treeAntiLambdaMC);
  // AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler *inputHandler =
  //     (AliInputEventHandler *)(man->GetInputEventHandler());
  // fPIDResponse = inputHandler->GetPIDResponse();

  fEventCuts.SetManualMode(); // Enable manual mode
  fEventCuts.fRequireTrackVertex = true;
  fEventCuts.fMinVtz = -10.f;
  fEventCuts.fMaxVtz = 10.f;
  fEventCuts.fMaxDeltaSpdTrackAbsolute = 0.5f;
  fEventCuts.fMaxResolutionSPDvertex = 0.25f;
  fEventCuts.fTriggerMask = AliVEvent::kINT7;
  fEventCuts.fRejectDAQincomplete = true;
  fEventCuts.fSPDpileupMinContributors = 3;
  fEventCuts.fSPDpileupMinZdist = 0.8;
  // fEventCuts.fSPDpileupNsigmaZdist = 3.;
  // fEventCuts.fSPDpileupNsigmaDiamXY = 2.;
  // fEventCuts.fSPDpileupNsigmaDiamZ = 5.;
  fEventCuts.fTrackletBGcut = true;
  fEventCuts.AddQAplotsToList(fOutputList);

  // // create track filters
  // fTrackFilter = new AliAnalysisFilter("trackFilter");
  // AliESDtrackCuts *fCuts = new AliESDtrackCuts();
  // fCuts->SetAcceptKinkDaughters(kFALSE);
  // fCuts->SetRequireTPCRefit(kTRUE);
  // fCuts->SetRequireITSRefit(kTRUE);
  // fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  // AliESDtrackCuts::kAny); fCuts->SetDCAToVertex2D(kFALSE);
  // fCuts->SetRequireSigmaToVertex(kFALSE);
  // fCuts->SetEtaRange(-0.8, 0.8);
  // fCuts->SetMinNCrossedRowsTPC(70);
  // fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  // fCuts->SetMaxChi2PerClusterTPC(4);
  // fCuts->SetMaxDCAToVertexZ(2);
  // fCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
  // fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  // fCuts->SetMaxChi2PerClusterTPC(4);
  // fCuts->SetMaxDCAToVertexZ(2);
  // fCuts->SetMaxChi2PerClusterITS(36);
  // fCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  // fCuts->SetMaxChi2PerClusterITS(36);
  // fTrackFilter->AddCuts(fCuts);

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityLambdaK0sMC::UserExec(Option_t *)
{

  AliVEvent *event = InputEvent();
  if (!event)
  {
    Error("UserExec", "Could not retrieve event");
    return;
  }
  AliESDEvent *lESDevent = 0x0;

  lESDevent = dynamic_cast<AliESDEvent *>(event);
  if (!lESDevent)
  {
    AliWarning("ERROR: lESDevent not available \n");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent)
  {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }

  AliStack *lMCstack = NULL;

  lMCstack = mcEvent->Stack();
  if (!lMCstack)
  {
    Printf("ERROR: Could not retrieve stack");
    return;
  }

  TArrayF mcPrimaryVtx;
  AliGenEventHeader *mcHeader = mcEvent->GenEventHeader();
  if (!mcHeader)
    return;
  mcHeader->PrimaryVertex(mcPrimaryVtx);

  Int_t iNumberOfPrimaries = -1;
  iNumberOfPrimaries = lMCstack->GetNprimary();
  if (iNumberOfPrimaries < 1)
    return;

  //////////////////

  //     AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>
  //     (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  //     if (!eventHandler) {
  // Printf("ERROR: Could not retrieve MC event handler");
  // return;
  //     }

  // AliGenEventHeader* mcHeader=mcEvent->GenEventHeader();
  // if(!mcHeader) return;
  // mcHeader->PrimaryVertex(mcPrimaryVtx);

  // if (TMath::Abs(mcPrimaryVtx.At(2)) > cutPrimVertex) return;  /// cut on z
  // of prim. vertex !!!!!!

  ////////////////
  // cout << event->GetFiredTriggerClasses() << endl;

  if (!fEventCuts.AcceptEvent(event))
  {
    // PostData(1, fOutputList);
    return;
  }

  AliVVZERO *esdV0 = lESDevent->GetVZEROData();
  if (!esdV0)
  {
    AliError("AliVVZERO not available");
    return;
  }
  if (lESDevent->IsIncompleteDAQ())
    return;
  AliAnalysisUtils *fUtils = new AliAnalysisUtils();
  if (fUtils->IsSPDClusterVsTrackletBG(lESDevent))
    return;
  if (lESDevent->IsPileupFromSPD(3))
    return;
  //    Bool_t maskIsSelected =
  //    ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  //     Bool_t isSelected = 0;
  //     //isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
  //    // Bool_t maskIsSelected=1;
  //     //pA triggering: CINT7
  //     isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;

  //     //Standard Min-Bias Selection
  //     if ( ! isSelected ) {
  //         return;
  //     }
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man)
  {
    AliInputEventHandler *inputHandler =
        (AliInputEventHandler *)(man->GetInputEventHandler());
    if (inputHandler)
      fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();
  }
  //------------------------------------------------
  // Step 3:  Primary Vertex quality selection
  //------------------------------------------------
  const AliESDVertex *lESDPrimaryTrackingVtx =
      lESDevent->GetPrimaryVertexTracks();
  const AliESDVertex *lESDPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
  const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();

  //------------------------------------------------
  // Step 3.1: reject events if SPDVtx or TrackVtx is not available
  //------------------------------------------------
  if (!(lESDPrimarySPDVtx->GetStatus() &&
        lESDPrimaryTrackingVtx->GetStatus()))
  {
    return;
  }

  //------------------------------------------------
  // Step 3.2: check the spd vertex resolution and reject if not satisfied
  //------------------------------------------------
  if (lESDPrimarySPDVtx->GetStatus() && lESDPrimarySPDVtx->IsFromVertexerZ() &&
      !(lESDPrimarySPDVtx->GetDispersion() < 0.04 &&
        lESDPrimarySPDVtx->GetZRes() < 0.25))
  {
    return;
  }

  //------------------------------------------------
  // Step 3.3: check the proximity between the spd vertex and trak vertex, and
  // reject if not satisfied
  //------------------------------------------------
  if ((TMath::Abs(lESDPrimarySPDVtx->GetZ() - lESDPrimaryTrackingVtx->GetZ()) >
       0.5))
  {
    return;
  }
  const AliESDVertex *lPrimaryBestESDVtx1 = lESDevent->GetPrimaryVertex();
  if (!lPrimaryBestESDVtx1)
  {
    return;
  }
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
  lPrimaryBestESDVtx->GetXYZ(lBestPrimaryVtxPos);

  if (TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0)
  {
    return;
  }
  Double_t lMagneticField = -10;

  Int_t lMultiplicity = -100;
  lMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(
      lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.8);
  if (!(lMultiplicity >= 1))
  {
    return;
  }

  //==========================================MC
  // START=========================================
  Int_t lNbMCPrimary = 0;

  Int_t lPdgcodeCurrentPart = 0;
  Double_t lRapCurrentPart = 0;
  Double_t lPtCurrentPart = 0;

  // Int_t lComeFromSigma      = 0;

  // current mc particle 's mother
  // Int_t iCurrentMother  = 0;
  lNbMCPrimary = lMCstack->GetNprimary();
  Long_t lNumberOfCharged = 0;

  //----- Loop on Stack
  //----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;
       iCurrentLabelStack < (lMCstack->GetNtrack());
       iCurrentLabelStack++)
  { // This is the begining of the loop on tracks
    TParticle *particleOne = lMCstack->Particle(iCurrentLabelStack);
    if (!particleOne)
      continue;
    if (!particleOne->GetPDG())
      continue;
    Double_t lThisCharge = particleOne->GetPDG()->Charge() / 3.;
    if (TMath::Abs(lThisCharge) < 0.001)
      continue;
    if (!(lMCstack->IsPhysicalPrimary(iCurrentLabelStack)))
      continue;

    // Double_t gpt = particleOne -> Pt();
    Double_t geta = particleOne->Eta();

    if (TMath::Abs(geta) < 0.5)
      lNumberOfCharged++;
  } // End of loop on tracks
  //----- End Loop on Stack
  //------------------------------------------------------------

  // Double_t lpArapidityShift = 0.465;
  Bool_t lStackNatural = kTRUE;
  //----- Loop on Lambda, K0Short
  //----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;
       iCurrentLabelStack < (lMCstack->GetNtrack());
       iCurrentLabelStack++)
  { // This is the begining of the loop on tracks

    TParticle *lCurrentParticleForLambdaCheck = 0x0;
    lCurrentParticleForLambdaCheck = lMCstack->Particle(iCurrentLabelStack);
    if (!lCurrentParticleForLambdaCheck)
    {
      Printf("V0s loop %d - MC TParticle pointer to current stack particle = "
             "0x0 ! Skip ...\n",
             iCurrentLabelStack);
      continue;
    }

    //=================================================================================
    // Single-Strange checks
    // Keep only K0s, Lambda and AntiLambda:
    lPdgcodeCurrentPart = lCurrentParticleForLambdaCheck->GetPdgCode();

    if ((lCurrentParticleForLambdaCheck->GetPdgCode() == 310) ||
        (lCurrentParticleForLambdaCheck->GetPdgCode() == 3122) ||
        (lCurrentParticleForLambdaCheck->GetPdgCode() == -3122))
    {
      lRapCurrentPart = MyRapidity(lCurrentParticleForLambdaCheck->Energy(),
                                   lCurrentParticleForLambdaCheck->Pz());
      lPtCurrentPart = lCurrentParticleForLambdaCheck->Pt();

      // Use Close to PV for filling CloseToPV histograms!
      Double_t dx, dy, dz;

      dx = ((mcPrimaryVtx.At(0)) - (lCurrentParticleForLambdaCheck->Vx()));
      dy = ((mcPrimaryVtx.At(1)) - (lCurrentParticleForLambdaCheck->Vy()));
      dz = ((mcPrimaryVtx.At(2)) - (lCurrentParticleForLambdaCheck->Vz()));
      Double_t lDistToPV = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      // if (lDistToPV <= 0.001)
      // {
      //   if (lPdgcodeCurrentPart == 3122)
      //   {
      //     f3dHistPrimCloseToPVPtVsYVsMultLambda->Fill(lPtCurrentPart,
      //     lRapCurrentPart, lMultiplicity);
      //   }
      //   if (lPdgcodeCurrentPart == -3122)
      //   {
      //     f3dHistPrimCloseToPVPtVsYVsMultAntiLambda->Fill(lPtCurrentPart,
      //     lRapCurrentPart, lMultiplicity);
      //   }
      //   if (lPdgcodeCurrentPart == 310)
      //   {
      //     f3dHistPrimCloseToPVPtVsYVsMultK0Short->Fill(lPtCurrentPart,
      //     lRapCurrentPart, lMultiplicity);
      //   }
      // }

      // Use Physical Primaries only for filling PrimRaw Histograms!
      // if (lMCstack->IsPhysicalPrimary(iCurrentLabelStack) != kTRUE)
      //   continue;

      lStackNatural = mcEvent->IsFromBGEvent(iCurrentLabelStack); // Is it?
      if (!lStackNatural)
      {
        if (!(lCurrentParticleForLambdaCheck->GetFirstMother() < 0))
        {
          lStackNatural = kTRUE; // because there are primaries (ALICE
                                 // definition) not produced in the collision
        }
      }

      if (lPdgcodeCurrentPart == 3122)
      {
        if (TMath::Abs(lRapCurrentPart) <= 0.5)
        {
          invpzLambdaMC = lCurrentParticleForLambdaCheck->Pz();
          invptLambdaMC = lPtCurrentPart;
          invyLambdaMC = lRapCurrentPart;
          ((TTree *)(fOutputList->FindObject("treeLambdaMC")))->Fill();
        }
        // f3dHistPrimRawPtVsYVsMultLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicity);
        // f3dHistPrimRawPtVsYCMSVsMultLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicity);

        // f3dHistPrimRawPtVsYVsMultV0ALambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityV0A);
        // f3dHistPrimRawPtVsYCMSVsMultV0ALambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityV0A);
        // f3dHistPrimRawPtVsYVsMultZNALambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityZNA);
        // f3dHistPrimRawPtVsYCMSVsMultZNALambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityZNA);
        // f3dHistPrimRawPtVsYVsMultTRKLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityTRK);
        // f3dHistPrimRawPtVsYCMSVsMultTRKLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityTRK);
        // f3dHistPrimRawPtVsYVsMultSPDLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicitySPD);
        // f3dHistPrimRawPtVsYCMSVsMultSPDLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicitySPD);

        // if (lStackNatural)
        // {
        //   f3dHistPrimRawPtVsYVsMultNonInjLambda->Fill(lPtCurrentPart,
        //   lRapCurrentPart, lMultiplicity);
        // }
        // f3dHistPrimRawPtVsYVsMultMCLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lNumberOfCharged);
        // f3dHistPrimRawPtVsYVsVertexZLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, mcPrimaryVtx.At(2)); if
        // (TMath::Abs(lCurrentParticleForLambdaCheck->Eta()) < 1.2 &&
        // lPtCurrentPart > 2)
        // {
        //   lHasHighPtLambda = kTRUE; // Keep track of events with Lambda
        //   within |eta|<1.2 and pt>2
        // }
      }
      if (lPdgcodeCurrentPart == -3122)
      {
        if (TMath::Abs(lRapCurrentPart) <= 0.5)
        {
          invpzAntiLambdaMC = lCurrentParticleForLambdaCheck->Pz();
          invptAntiLambdaMC = lPtCurrentPart;
          invyAntiLambdaMC = lRapCurrentPart;
          ((TTree *)(fOutputList->FindObject("treeAntiLambdaMC")))->Fill();
        }
        // f3dHistPrimRawPtVsYVsMultAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicity);
        // f3dHistPrimRawPtVsYCMSVsMultAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicity);

        // f3dHistPrimRawPtVsYVsMultV0AAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityV0A);
        // f3dHistPrimRawPtVsYCMSVsMultV0AAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityV0A);
        // f3dHistPrimRawPtVsYVsMultZNAAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityZNA);
        // f3dHistPrimRawPtVsYCMSVsMultZNAAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityZNA);
        // f3dHistPrimRawPtVsYVsMultTRKAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityTRK);
        // f3dHistPrimRawPtVsYCMSVsMultTRKAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityTRK);
        // f3dHistPrimRawPtVsYVsMultSPDAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicitySPD);
        // f3dHistPrimRawPtVsYCMSVsMultSPDAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicitySPD);

        // if (lStackNatural)
        // {
        //   f3dHistPrimRawPtVsYVsMultNonInjAntiLambda->Fill(lPtCurrentPart,
        //   lRapCurrentPart, lMultiplicity);
        // }
        // f3dHistPrimRawPtVsYVsMultMCAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, lNumberOfCharged);
        // f3dHistPrimRawPtVsYVsVertexZAntiLambda->Fill(lPtCurrentPart,
        // lRapCurrentPart, mcPrimaryVtx.At(2));
      }
      if (lPdgcodeCurrentPart == 310)
      {
        if (TMath::Abs(lRapCurrentPart) <= 0.5)
        {
          // invmK0s = v0->GetEffMass();
          invpzK0sMC = lCurrentParticleForLambdaCheck->Pz();
          invptK0sMC = lPtCurrentPart;
          invyK0sMC = lRapCurrentPart;
          ((TTree *)(fOutputList->FindObject("treeK0sMC")))->Fill();
        }
        // f3dHistPrimRawPtVsYVsMultK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicity);
        // f3dHistPrimRawPtVsYCMSVsMultK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicity);

        // f3dHistPrimRawPtVsYVsMultV0AK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityV0A);
        // f3dHistPrimRawPtVsYCMSVsMultV0AK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityV0A);
        // f3dHistPrimRawPtVsYVsMultZNAK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityZNA);
        // f3dHistPrimRawPtVsYCMSVsMultZNAK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityZNA);
        // f3dHistPrimRawPtVsYVsMultTRKK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicityTRK);
        // f3dHistPrimRawPtVsYCMSVsMultTRKK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicityTRK);
        // f3dHistPrimRawPtVsYVsMultSPDK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart, lMultiplicitySPD);
        // f3dHistPrimRawPtVsYCMSVsMultSPDK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart + fpArapidityShift, lMultiplicitySPD);

        // if (lStackNatural)
        // {
        //   f3dHistPrimRawPtVsYVsMultNonInjK0Short->Fill(lPtCurrentPart,
        //   lRapCurrentPart, lMultiplicity);
        // }
        // f3dHistPrimRawPtVsYVsMultMCK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart, lNumberOfCharged);
        // f3dHistPrimRawPtVsYVsVertexZK0Short->Fill(lPtCurrentPart,
        // lRapCurrentPart, mcPrimaryVtx.At(2));
      }
      // Decay Length
      // Acquisition=====================================================
      // Double_t decaylength = -1;
      // Double_t lV0Mass = -1;

      // if (!(lCurrentParticleForLambdaCheck->GetDaughter(0) < 0))
      // {
      //   TParticle *lDght0ofV0 =
      //   lMCstack->Particle(lCurrentParticleForLambdaCheck->GetDaughter(0));
      //   // get first daughter if (lDght0ofV0) { // skip if not defined.
      //     decaylength = TMath::Sqrt(
      //         TMath::Power(lCurrentParticleForLambdaCheck->Vx() -
      //         lDght0ofV0->Vx(), 2) +
      //         TMath::Power(lCurrentParticleForLambdaCheck->Vy() -
      //         lDght0ofV0->Vy(), 2) +
      //         TMath::Power(lCurrentParticleForLambdaCheck->Vz() -
      //         lDght0ofV0->Vz(), 2));
      //     // Need to correct for relativitity! Involves multiplying by mass
      //     and dividing by momentum. if (TMath::Abs(lPdgcodeCurrentPart) ==
      //     3122)
      //     {
      //       lV0Mass = 1.115683;
      //     }
      //     if (TMath::Abs(lPdgcodeCurrentPart) == 310)
      //     {
      //       lV0Mass = 0.497614;
      //     }
      //     if (lCurrentParticleForLambdaCheck->P() + 1e-10 != 0)
      //       decaylength = (lV0Mass * decaylength) /
      //       (lCurrentParticleForLambdaCheck->P() + 1e-10);
      //     if (lCurrentParticleForLambdaCheck->P() + 1e-10 == 0)
      //       decaylength = 1e+5;
      //   }
      // }
      // if (lPdgcodeCurrentPart == 3122)
      //   f3dHistPrimRawPtVsYVsDecayLengthLambda->Fill(lPtCurrentPart,
      //   lRapCurrentPart, decaylength);
      // if (lPdgcodeCurrentPart == -3122)
      //   f3dHistPrimRawPtVsYVsDecayLengthAntiLambda->Fill(lPtCurrentPart,
      //   lRapCurrentPart, decaylength);
      // if (lPdgcodeCurrentPart == 310)
      //   f3dHistPrimRawPtVsYVsDecayLengthK0Short->Fill(lPtCurrentPart,
      //   lRapCurrentPart, decaylength);
    }
  } // End of loop on tracks
  //----- End Loop on Lambda, K0Short
  //------------------------------------------------------------

  //==========================================MC
  // END=========================================

  lMagneticField = lESDevent->GetMagneticField();
  Double_t flat = GetFlatenicityV0();
  ((TH1D *)(fOutputList->FindObject("hflat")))->Fill(1.0 - flat);
  ((TH1D *)(fOutputList->FindObject("hNeventMC")))->Fill(0.5);

  // Run number
  // Int_t fRun = lESDevent->GetRunNumber();
  Int_t lOnFlyStatus = 0; // nv0sOn = 0, nv0sOff = 0;
  Double_t lChi2V0 = 0;
  Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
  Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
  Double_t lV0CosineOfPointingAngle = 0;
  Double_t lV0Radius = 0, lPt = 0;
  Double_t lRapK0Short = 0, lRapLambda = 0;
  Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  Double_t lAlphaV0 = 0, lPtArmV0 = 0;

  Double_t fMinV0Pt = 0;
  Double_t fMaxV0Pt = 100;

  Int_t nv0s = 0;
  nv0s = lESDevent->GetNumberOfV0s();

  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) // extra-crazy test
  {                                      // This is the begining of the V0 loop
    AliESDv0 *v0 = ((AliESDEvent *)lESDevent)->GetV0(iV0);
    if (!v0)
      continue;

    // CheckChargeV0(v0);
    // // Remove like-sign (will not affect offline V0 candidates!)
    // if (v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() > 0)
    // {
    //   continue;
    // }
    // if (v0->GetParamN()->Charge() < 0 && v0->GetParamP()->Charge() < 0)
    // {
    //   continue;
    // }

    lOnFlyStatus = v0->GetOnFlyStatus();
    if (lOnFlyStatus)
      continue;

    Double_t tDecayVertexV0[3];
    v0->GetXYZ(tDecayVertexV0[0], tDecayVertexV0[1], tDecayVertexV0[2]);

    Double_t tV0mom[3];
    v0->GetPxPyPz(tV0mom[0], tV0mom[1], tV0mom[2]);
    Double_t lV0TotalMomentum = TMath::Sqrt(
        tV0mom[0] * tV0mom[0] + tV0mom[1] * tV0mom[1] + tV0mom[2] * tV0mom[2]);

    Double_t lV0Radius = TMath::Sqrt(tDecayVertexV0[0] * tDecayVertexV0[0] +
                                     tDecayVertexV0[1] * tDecayVertexV0[1]);

    Double_t lPt = v0->Pt();
    Double_t lRapK0Short = v0->RapK0Short();
    Double_t lRapLambda = v0->RapLambda();

    UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
    UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

    Double_t lMomPos[3];
    v0->GetPPxPyPz(lMomPos[0], lMomPos[1], lMomPos[2]);
    Double_t lMomNeg[3];
    v0->GetNPxPyPz(lMomNeg[0], lMomNeg[1], lMomNeg[2]);

    // Provisions for cowboy/sailor check
    Double_t lModp1 =
        TMath::Sqrt(lMomPos[0] * lMomPos[0] + lMomPos[1] * lMomPos[1]);
    Double_t lModp2 =
        TMath::Sqrt(lMomNeg[0] * lMomNeg[0] + lMomNeg[1] * lMomNeg[1]);

    // Calculate vec prod with momenta projected to xy plane
    Double_t lVecProd =
        (lMomPos[0] * lMomNeg[1] - lMomPos[1] * lMomNeg[0]) / (lModp1 * lModp2);

    if (lMagneticField < 0)
      lVecProd *= -1; // invert sign

    Bool_t fTreeVariableIsCowboy = kFALSE;
    if (lVecProd < 0)
      fTreeVariableIsCowboy = kTRUE;

    AliESDtrack *pTrack = ((AliESDEvent *)lESDevent)->GetTrack(lKeyPos);
    AliESDtrack *nTrack = ((AliESDEvent *)lESDevent)->GetTrack(lKeyNeg);
    AliExternalTrackParam *fTreeVariablePosTrack; //!
    AliExternalTrackParam *fTreeVariableNegTrack; //!
    fTreeVariablePosTrack = pTrack;
    fTreeVariableNegTrack = nTrack;

    if (!pTrack || !nTrack)
    {
      Printf("ERROR: Could not retreive one of the daughter track");
      continue;
    }
    Int_t fTreeVariablePosPIDForTracking = pTrack->GetPIDForTracking();
    Int_t fTreeVariableNegPIDForTracking = nTrack->GetPIDForTracking();

    const AliExternalTrackParam *innernegv0 = nTrack->GetInnerParam();
    const AliExternalTrackParam *innerposv0 = pTrack->GetInnerParam();
    Float_t lThisPosInnerP = -1;
    Float_t lThisNegInnerP = -1;
    Float_t lThisPosInnerPt = -1;
    Float_t lThisNegInnerPt = -1;
    if (innerposv0)
    {
      lThisPosInnerP = innerposv0->GetP();
    }
    if (innernegv0)
    {
      lThisNegInnerP = innernegv0->GetP();
    }
    if (innerposv0)
    {
      lThisPosInnerPt = innerposv0->Pt();
    }
    if (innernegv0)
    {
      lThisNegInnerPt = innernegv0->Pt();
    }
    Float_t lThisPosdEdx = pTrack->GetTPCsignal();
    Float_t lThisNegdEdx = nTrack->GetTPCsignal();

    // Daughter Eta for Eta selection, afterwards
    Float_t fTreeVariableNegEta = nTrack->Eta();
    Float_t fTreeVariablePosEta = pTrack->Eta();

    if (TMath::Abs(fTreeVariableNegEta) > 0.8 ||
        TMath::Abs(fTreeVariableNegEta) > 0.8)
      continue;

    // Filter like-sign V0 (next: add counter and distribution)
    if (pTrack->GetSign() == nTrack->GetSign())
    {
      continue;
    }

    //________________________________________________________________________
    // Track quality cuts
    Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2, 1);
    Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2, 1);
    Int_t fTreeVariableLeastNbrCrossedRows = (Int_t)lPosTrackCrossedRows;
    if (lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows)
      fTreeVariableLeastNbrCrossedRows = (Int_t)lNegTrackCrossedRows;

    // TPC refit condition (done during reconstruction for Offline but not for
    // On-the-fly)
    if (!(pTrack->GetStatus() & AliESDtrack::kTPCrefit))
      continue;
    if (!(nTrack->GetStatus() & AliESDtrack::kTPCrefit))
      continue;

    // Get status flags
    ULong64_t fTreeVariablePosTrackStatus = pTrack->GetStatus();
    ULong64_t fTreeVariableNegTrackStatus = nTrack->GetStatus();

    Float_t fTreeVariablePosDCAz = GetDCAz(pTrack);
    Float_t fTreeVariableNegDCAz = GetDCAz(nTrack);

    // GetKinkIndex condition
    if (pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0)
      continue;

    // Findable clusters > 0 condition
    if (pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0)
      continue;

    // Compute ratio Crossed Rows / Findable clusters
    // Note: above test avoids division by zero!
    Float_t lPosTrackCrossedRowsOverFindable =
        lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
    Float_t lNegTrackCrossedRowsOverFindable =
        lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));

    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable =
        lPosTrackCrossedRowsOverFindable;
    if (lNegTrackCrossedRowsOverFindable <
        fTreeVariableLeastRatioCrossedRowsOverFindable)
      fTreeVariableLeastRatioCrossedRowsOverFindable =
          lNegTrackCrossedRowsOverFindable;

    // Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
    if (fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8)
      continue;

    // Extra track quality: Chi2/cluster for cross-checks
    Float_t lBiggestChi2PerCluster = -1;

    Float_t lPosChi2PerCluster = 1000;
    Float_t lNegChi2PerCluster = 1000;

    if (pTrack->GetTPCNcls() > 0)
      lPosChi2PerCluster =
          pTrack->GetTPCchi2() / ((Float_t)pTrack->GetTPCNcls());
    if (nTrack->GetTPCNcls() > 0)
      lNegChi2PerCluster =
          nTrack->GetTPCchi2() / ((Float_t)nTrack->GetTPCNcls());

    if (lPosChi2PerCluster > lBiggestChi2PerCluster)
      lBiggestChi2PerCluster = lPosChi2PerCluster;
    if (lNegChi2PerCluster > lBiggestChi2PerCluster)
      lBiggestChi2PerCluster = lNegChi2PerCluster;

    Float_t fTreeVariableMaxChi2PerCluster = lBiggestChi2PerCluster;

    // Extra track quality: min track length
    Float_t lSmallestTrackLength = 1000;
    Float_t lPosTrackLength = -1;
    Float_t lNegTrackLength = -1;

    if (pTrack->GetInnerParam())
      lPosTrackLength = pTrack->GetLengthInActiveZone(
          1, 2.0, 220.0, lESDevent->GetMagneticField());
    if (nTrack->GetInnerParam())
      lNegTrackLength = nTrack->GetLengthInActiveZone(
          1, 2.0, 220.0, lESDevent->GetMagneticField());

    if (lPosTrackLength < lSmallestTrackLength)
      lSmallestTrackLength = lPosTrackLength;
    if (lNegTrackLength < lSmallestTrackLength)
      lSmallestTrackLength = lNegTrackLength;

    if ((((pTrack->GetTPCClusterInfo(2, 1)) < 70) ||
         ((nTrack->GetTPCClusterInfo(2, 1)) < 70)) &&
        lSmallestTrackLength < 80)
      continue;

    // End track Quality Cuts
    //________________________________________________________________________

    lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(
        lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField));

    lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(
        lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField));

    if (lDcaPosToPrimVertex < 0.06 || lDcaNegToPrimVertex < 0.06)
      continue;

    lChi2V0 = v0->GetChi2V0();
    lDcaV0Daughters = v0->GetDcaV0Daughters();
    if (lDcaV0Daughters > 1)
      continue;
    lDcaV0ToPrimVertex = v0->GetD(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1],
                                  lBestPrimaryVtxPos[2]);
    Float_t lV0DecayLength =
        TMath::Sqrt(TMath::Power(tDecayVertexV0[0] - lBestPrimaryVtxPos[0], 2) +
                    TMath::Power(tDecayVertexV0[1] - lBestPrimaryVtxPos[1], 2) +
                    TMath::Power(tDecayVertexV0[2] - lBestPrimaryVtxPos[2], 2));
    if (lV0DecayLength < 0.5)
      continue;
    lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(
        lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2]);
    Float_t fTreeVariableV0CosineOfPointingAngle = lV0CosineOfPointingAngle;

    // Getting invariant mass infos directly from ESD
    // v0->ChangeMassHypothesis(310);
    // lInvMassK0s = v0->GetEffMass();
    // v0->ChangeMassHypothesis(3122);
    // lInvMassLambda = v0->GetEffMass();
    // v0->ChangeMassHypothesis(-3122);
    // lInvMassAntiLambda = v0->GetEffMass();
    // lAlphaV0 = v0->AlphaV0();
    // lPtArmV0 = v0->PtArmV0();

    // Official means of acquiring N-sigmas
    Float_t NSigmasPosProton =
        fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
    Float_t NSigmasPosPion =
        fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
    Float_t NSigmasNegProton =
        fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton);
    Float_t NSigmasNegPion =
        fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);

    // K0Short: Enough to parametrize peak broadening with linear function.
    Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02) * lPt;
    Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02) * lPt;
    // Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
    //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
    Double_t lUpperLimitLambda =
        (1.13688e+00) + (5.27838e-03) * lPt +
        (8.42220e-02) * TMath::Exp(-(3.80595e+00) * lPt);
    Double_t lLowerLimitLambda =
        (1.09501e+00) - (5.23272e-03) * lPt -
        (7.52690e-02) * TMath::Exp(-(3.46339e+00) * lPt);

    v0->ChangeMassHypothesis(310);
    lInvMassK0s = v0->GetEffMass();
    if (lInvMassK0s < lUpperLimitK0Short && lInvMassK0s > lLowerLimitK0Short)
    {
      if (lV0CosineOfPointingAngle > 0.97)
      {
        Float_t lPtK0s = v0->Pt();
        Float_t lPzK0s = v0->Pz();
        if (lPtK0s != 0)
        {
          Float_t lRapK0s = v0->Y(310);
          if (TMath::Abs(lRapK0s) < 0.5)
          {
            Float_t ctauK0s = lV0DecayLength * lInvMassK0s / v0->P();
            if (ctauK0s < 20)
            {
              if (TMath::Abs(NSigmasPosPion) < 5.0 &&
                  TMath::Abs(NSigmasNegPion) < 5.0)
              {
                ((TH1F *)(fOutputList->FindObject("hinvmassK0s")))
                    ->Fill(lInvMassK0s);
                invmK0s = v0->GetEffMass();
                invpK0s = v0->P();
                invptK0s = v0->Pt();
                invyK0s = lRapK0s;
                flatenicityK0s = 1.0 - flat;
                ((TTree *)(fOutputList->FindObject("treeK0s")))->Fill();
              }
            }
          }
        }
      }
    }

    Float_t massLambda = 1.11568;

    v0->ChangeMassHypothesis(3122);
    lInvMassLambda = v0->GetEffMass();
    if (lInvMassLambda < lUpperLimitLambda &&
        lInvMassLambda > lLowerLimitLambda)
    {
      if (lV0CosineOfPointingAngle > 0.995)
      {
        Float_t lPtLambda = v0->Pt();
        Float_t lPzLambda = v0->Pz();
        if (lPtLambda != 0)
        {
          Float_t lRapLambda = v0->Y(3122);
          if (TMath::Abs(lRapLambda) < 0.5)
          {
            Float_t ctauLambda = lV0DecayLength * lInvMassLambda / v0->P();
            if (ctauLambda < 30)
            {
              if (TMath::Abs(NSigmasPosProton) < 5.0 &&
                  TMath::Abs(NSigmasNegPion) < 5.0)
              {
                invmLambda = v0->GetEffMass();
                invpLambda = v0->P();
                invptLambda = v0->Pt();
                invyLambda = lRapLambda;
                flatenicityLambda = 1.0 - flat;

                ((TTree *)(fOutputList->FindObject("treeLambda")))->Fill();
                ((TH1F *)(fOutputList->FindObject("hinvmassLambda")))
                    ->Fill(lInvMassLambda);
              }
            }
          }
        }
      }
    }

    v0->ChangeMassHypothesis(-3122);

    lInvMassAntiLambda = v0->GetEffMass();

    if (lInvMassAntiLambda < lUpperLimitLambda &&
        lInvMassAntiLambda > lLowerLimitLambda)
    {
      if (lV0CosineOfPointingAngle > 0.995)
      {
        Float_t lPtAntiLambda = v0->Pt();
        Float_t lPzAntiLambda = v0->Pz();

        if (lPtAntiLambda != 0)
        {
          Float_t lRapAntiLambda = v0->Y(-3122);
          if (TMath::Abs(lRapAntiLambda) < 0.5)
          {
            Float_t ctauAntiLambda =
                lV0DecayLength * lInvMassAntiLambda / v0->P();
            if (ctauAntiLambda < 30)
            {
              if (TMath::Abs(NSigmasNegProton) < 5.0 &&
                  TMath::Abs(NSigmasPosPion) < 5.0)
              {
                invmAntiLambda = v0->GetEffMass();
                invpAntiLambda = v0->P();
                invptAntiLambda = v0->Pt();
                invyAntiLambda = lRapAntiLambda;
                flatenicityAntiLambda = 1.0 - flat;

                ((TTree *)(fOutputList->FindObject("treeAntiLambda")))->Fill();
                ((TH1F *)(fOutputList->FindObject("hinvmassAntiLambda")))
                    ->Fill(lInvMassAntiLambda);
              }
            }
          }
        }
      }
    }
  }

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityLambdaK0sMC::Terminate(Option_t *)
{
  // terminate
  // called at the END of the anfalysis (when all events are processed)
}
//_____________________________________________________________________________
//________________________________________________________________________

Double_t AliAnalysisTaskFlatenicityLambdaK0sMC::MyRapidity(Double_t rE,
                                                           Double_t rPz) const
{
  // Local calculation for rapidity
  Double_t ReturnValue = -100;
  if ((rE - rPz + 1.e-13) != 0 && (rE + rPz) != 0)
  {
    ReturnValue = 0.5 * TMath::Log((rE + rPz) / (rE - rPz + 1.e-13));
  }
  return ReturnValue;
}

void AliAnalysisTaskFlatenicityLambdaK0sMC::CheckChargeV0(AliESDv0 *v0)
{
  // This function checks charge of negative and positive daughter tracks.
  // If incorrectly defined (onfly vertexer), swaps out.
  if (v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0)
  {
    // V0 daughter track swapping is required! Note: everything is swapped
    // here... P->N, N->P
    Long_t lCorrectNidx = v0->GetPindex();
    Long_t lCorrectPidx = v0->GetNindex();
    Double32_t lCorrectNmom[3];
    Double32_t lCorrectPmom[3];
    v0->GetPPxPyPz(lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2]);
    v0->GetNPxPyPz(lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2]);

    AliExternalTrackParam lCorrectParamN(
        v0->GetParamP()->GetX(), v0->GetParamP()->GetAlpha(),
        v0->GetParamP()->GetParameter(), v0->GetParamP()->GetCovariance());
    AliExternalTrackParam lCorrectParamP(
        v0->GetParamN()->GetX(), v0->GetParamN()->GetAlpha(),
        v0->GetParamN()->GetParameter(), v0->GetParamN()->GetCovariance());
    lCorrectParamN.SetMostProbablePt(v0->GetParamP()->GetMostProbablePt());
    lCorrectParamP.SetMostProbablePt(v0->GetParamN()->GetMostProbablePt());

    // Get Variables___________________________________________________
    Double_t lDcaV0Daughters = v0->GetDcaV0Daughters();
    Double_t lCosPALocal = v0->GetV0CosineOfPointingAngle();
    Bool_t lOnFlyStatusLocal = v0->GetOnFlyStatus();

    // Create Replacement Object_______________________________________
    AliESDv0 *v0correct = new AliESDv0(lCorrectParamN, lCorrectNidx,
                                       lCorrectParamP, lCorrectPidx);
    v0correct->SetDcaV0Daughters(lDcaV0Daughters);
    v0correct->SetV0CosineOfPointingAngle(lCosPALocal);
    v0correct->ChangeMassHypothesis(kK0Short);
    v0correct->SetOnFlyStatus(lOnFlyStatusLocal);

    // Reverse Cluster info..._________________________________________
    v0correct->SetClusters(v0->GetClusters(1), v0->GetClusters(0));

    *v0 = *v0correct;
    // Proper cleanup..._______________________________________________
    v0correct->Delete();
    v0correct = 0x0;

    // Just another cross-check and output_____________________________
    if (v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0)
    {
      AliWarning(
          "Found Swapped Charges, tried to correct but something FAILED!");
    }
    else
    {
      // AliWarning("Found Swapped Charges and fixed.");
    }
    //________________________________________________________________
  }
  else
  {
    // Don't touch it! ---
    // Printf("Ah, nice. Charges are already ordered...");
  }
  return;
}

//________________________________________________________________________
Float_t AliAnalysisTaskFlatenicityLambdaK0sMC::GetDCAz(AliESDtrack *lTrack)
// Encapsulation of DCAz calculation
{
  Float_t b[2];
  Float_t bCov[3];
  lTrack->GetImpactParameters(b, bCov);
  if (bCov[0] <= 0 || bCov[2] <= 0)
  {
    AliDebug(1, "Estimated b resolution lower or equal to zero!");
    bCov[0] = 0;
    bCov[2] = 0;
  }
  // Float_t dcaToVertexXY = b[0];
  Float_t dcaToVertexZ = b[1];

  return dcaToVertexZ;
}

Float_t AliAnalysisTaskFlatenicityLambdaK0sMC::GetFlatenicityV0()
{
  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent)
  {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return -1;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0)
  {
    AliError("AliVVZERO not available");
    return 9999;
  }
  // Flatenicity calculation
  const Int_t nRings = 4;
  const Int_t nSectors = 8;
  Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
  Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
  Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
  Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
  // Grid
  const Int_t nCells = nRings * 2 * nSectors;
  Float_t RhoLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++)
  {
    RhoLattice[iCh] = 0.0;
  }

  Int_t nringA = 0;
  Int_t nringC = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++)
  {
    Float_t detaV0 = -1;
    Float_t mult = lVV0->GetMultiplicity(iCh);
    if (iCh < 32)
    { // V0C
      if (iCh < 8)
      {
        nringC = 0;
      }
      else if (iCh >= 8 && iCh < 16)
      {
        nringC = 1;
      }
      else if (iCh >= 16 && iCh < 24)
      {
        nringC = 2;
      }
      else
      {
        nringC = 3;
      }
      detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
    }
    else
    { // V0A
      if (iCh < 40)
      {
        nringA = 0;
      }
      else if (iCh >= 40 && iCh < 48)
      {
        nringA = 1;
      }
      else if (iCh >= 48 && iCh < 56)
      {
        nringA = 2;
      }
      else
      {
        nringA = 3;
      }
      detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
    }
    RhoLattice[iCh] =
        mult / detaV0; // needed to consider the different eta coverage
  }
  // Filling histos with mult info
  // for (Int_t iCh = 0; iCh < nCells; iCh++) {
  // hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
  //}
  Float_t mRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++)
  {
    mRho += RhoLattice[iCh];
  }
  Float_t multiplicityV0M = mRho;
  // average activity per cell
  mRho /= (1.0 * nCells);
  // get sigma
  Double_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++)
  {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);
  Bool_t fRemoveTrivialScaling = kFALSE;
  if (mRho > 0)
  {
    if (fRemoveTrivialScaling)
    {
      flatenicity = TMath::Sqrt(multiplicityV0M) * sRho / mRho;
    }
    else
    {
      flatenicity = sRho / mRho;
    }
  }
  else
  {
    flatenicity = 9999;
  }
  return flatenicity;
}
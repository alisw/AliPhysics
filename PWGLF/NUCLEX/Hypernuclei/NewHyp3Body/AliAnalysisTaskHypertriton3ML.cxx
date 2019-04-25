/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
// AliAnalysisTaskHypertriton3ML class
// analysis task for the study of the production of hypertriton
// which decays in 3 prongs: d+p+pi^-
// This task is optimized for ESDs.root
//
// Author:
// P. Fecchio, pfecchio@cern.ch
///////////////////////////////////////////////////////////////////////////

/// ROOT includes
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TObjArray.h>
#include <TVector3.h>

/// AliRoot icludes
#include "AliAnalysisTaskSE.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskHypertriton3ML.h"

// #include <string>
#include <Riostream.h>
#include <climits>
#include <vector>

const char *kDataMC[2]                  = {"", "MC"};
const char *kAM                         = "AM";
const char *kSpecNames[4]               = {"pion", "proton", "deuteron", "hyper-triton"};
const AliPID::EParticleType kSpecies[3] = {AliPID::kDeuteron, AliPID::kProton, AliPID::kPion};
const float kDeuteronMass               = 1.87561297416687012f;
const float kProtonMass                 = 9.38271999359130859e-01f;
const float kPionMass                   = 1.39569997787475586e-01f;

ClassImp(AliAnalysisTaskHypertriton3ML)

    using std::cout;
using std::endl;

//________________________________________________________________________
AliAnalysisTaskHypertriton3ML::AliAnalysisTaskHypertriton3ML()
    : AliAnalysisTaskSE(),
      // parameter settings
      fIsMC{0}, fNSigma{0., 0., 0.}, //
      // support objects
      fESDevent{nullptr},      //
      fEventCuts{},            //
      fTrackCuts{nullptr},     //
      fPrimaryVertex{nullptr}, //
      fPIDResponse{nullptr},   //

      fMagneticField{0.0}, //
      // setting parameters
      fCosPoiningAngleLimit{0}, //
      // output objects
      fMInvBackgroundStd{nullptr},  //
      fMInvBackgroundCuts{nullptr}, //
      fHistCosPAngle{nullptr},      //
      fOutputList{nullptr} {
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskHypertriton3ML::AliAnalysisTaskHypertriton3ML(bool isMC, TString taskname)
    : AliAnalysisTaskSE(taskname.Data()),
      // parameter settings
      fIsMC{false}, fNSigma{3., 3., 3.}, //
      // support objects
      fESDevent{nullptr},      //
      fEventCuts{},            //
      fTrackCuts{nullptr},     //
      fPrimaryVertex{nullptr}, //
      fPIDResponse{nullptr},   //
      fMagneticField{0.0},     //
      // setting parameters
      fCosPoiningAngleLimit{0.8}, //
      // output objects
      fMInvBackgroundStd{nullptr},  //
      fMInvBackgroundCuts{nullptr}, //
      fHistCosPAngle{nullptr},      //
      fOutputList{nullptr} {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  // ESD Track cuts
  if (!fTrackCuts) fTrackCuts = new AliESDtrackCuts();
  fTrackCuts->SetMinNClustersTPC(80);
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetMaxChi2PerClusterTPC(4);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetEtaRange(-0.9, 0.9);
}

//________________________________________________________________________
AliAnalysisTaskHypertriton3ML::~AliAnalysisTaskHypertriton3ML() {
  // destructor
  if (fESDevent) delete fESDevent;
  if (fTrackCuts) delete fTrackCuts;
  if (fPrimaryVertex) delete fPrimaryVertex;
  if (fPIDResponse) delete fPIDResponse;
}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3ML::UserCreateOutputObjects() {

  AliAnalysisManager *fMgr = AliAnalysisManager::GetAnalysisManager();
  if (!fMgr) AliFatal("Could not find analysis manager.");
  AliInputEventHandler *fHandl = (AliInputEventHandler *)fMgr->GetInputEventHandler();
  if (!fHandl) AliFatal("No input event handler.");
  fPIDResponse = dynamic_cast<AliPIDResponse *>(fHandl->GetPIDResponse());
  fHandl->SetNeedField();

  // create the output tree
  fHypTree = new TTree("fHypTree", "hyper-triton candidates");
  fHypTree->Branch("MagField", &fMagneticField);

  // Create a TList with Histograms
  fOutputList = new TList();
  fOutputList->SetOwner(true);
  fOutputList->SetName("clistHypertriton");
  fEventCuts.AddQAplotsToList(fOutputList);

  // binning della 3 body?
  fMInvBackgroundStd =
      new TH1D(Form("fMInvBackgroundStd%s", kDataMC[fIsMC]), ";m (dp#pi) [GeV/#it{c}^{2}];Counts", 200, 2.95, 3.35);
  fMInvBackgroundCuts =
      new TH1D(Form("fMInvBackgroundCuts%s", kDataMC[fIsMC]), ";m (dp#pi) [GeV/#it{c}^{2}];Counts", 200, 2.95, 3.35);
  fHistCosPAngle = new TH1D("fCosPointingAngle", ";#it{cos#theta_{pointing}} ;Counts", 200, -1, 1);

  // pt distridution of the daughter particles
  for (int i = 0; i < 3; ++i) {
    fHistPtStd[i]  = new TH1D(Form("fHistPt_%s_std", kSpecNames[i]), "", 100, 0, 10);
    fHistPtCuts[i] = new TH1D(Form("fHistPt_%s_cuts", kSpecNames[i]), "", 100, 0, 10);
    fOutputList->Add(fHistPtStd[i]);
    fOutputList->Add(fHistPtCuts[i]);
  }

  // Adding output to list
  fOutputList->Add(fMInvBackgroundStd);
  fOutputList->Add(fMInvBackgroundCuts);
  fOutputList->Add(fHistCosPAngle);

  PostData(1, fOutputList);
  PostData(2, fHypTree);
}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3ML::UserExec(Option_t *) {
  // main loop called for each analized event

  AliVEvent *fEv = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fEv) {
    ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *fMCEvent = MCEvent();
  if (!fMCEvent && fIsMC) {
    ::Fatal("AliAnalysisTaskStrangenessLifetimes::UserExec", "Could not retrieve MC event");
    return;
  }

  fMagneticField = (float)fEv->GetMagneticField();

  // Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(fEv)) {
    PostData(1, fOutputList);
    PostData(2, fHypTree);
    return;
  }

  fHypTree->Fill();

  PostData(1, fOutputList);
  PostData(2, fHypTree);
}

//________________________________________________________________________
void AliAnalysisTaskHypertriton3ML::Terminate(Option_t *) {
  // Merge output
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList *>(GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: fOutputList not available\n");
    return;
  }

  printf("end of Terminate");
  return;

} // end of Terminate

//________________________________________________________________________
bool AliAnalysisTaskHypertriton3ML::PassPIDSelection(AliESDtrack *track, AliPID::EParticleType specie,
                                                     float nSigmaCut) {
  // Method that return true if the ESD track pass all the PID requirements for the specified specie
  if (abs(fPIDResponse->NumberOfSigmasTPC(track, specie)) < nSigmaCut) {
    return true;
  }
  return false;
}

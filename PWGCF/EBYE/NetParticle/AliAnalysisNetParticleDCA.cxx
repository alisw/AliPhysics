//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliESDtrackCuts.h"

#include "AliAnalysisNetParticleDCA.h"

using namespace std;

// Task for NetParticle checks
// Author: Jochen Thaeder <jochen@thaeder.de>

ClassImp(AliAnalysisNetParticleDCA)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleDCA::AliAnalysisNetParticleDCA() :
  fHelper(NULL),

  fESD(NULL), 
  fESDTrackCuts(NULL),
  fESDTrackCutsBkg(NULL),

  fStack(NULL),
  fMCEvent(NULL),

  fHnDCA(NULL) {
  // Constructor   

  AliLog::SetClassDebugLevel("AliAnalysisNetParticleDCA",10);
}

//________________________________________________________________________
AliAnalysisNetParticleDCA::~AliAnalysisNetParticleDCA() {
  // Destructor

}



/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleDCA::Initialize(AliESDtrackCuts *cuts, AliESDtrackCuts *cutsBkg, AliAnalysisNetParticleHelper* helper) {

  // -- Get Helper
  // ---------------
  fHelper          = helper;

  // -- ESD track cuts
  // -------------------
  fESDTrackCuts    = cuts;
  fESDTrackCutsBkg = cutsBkg;

#if 0
  // -- Get particle species / pdgCode
  // -------------------------
  fPdgCode          = AliPID::ParticleCode(fHelper->GetParticleSpecies());

  // -- MC Labels for efficiency
  // -----------------------------
  fLabelsRec        = new Int_t*[2];
  for (Int_t ii = 0; ii < 2 ; ++ii)
    fLabelsRec[ii] = NULL;

  // -- Create THnSparse Histograms
  // --------------------------------
  CreateHistograms();
#endif  
  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Setup/Reset Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleDCA::SetupEvent(AliESDInputHandler *esdHandler, AliMCEvent *mcEvent) {
  // -- Setup event

  // -- Reset of event
  // -------------------
  ResetEvent();

  // -- Setup of event
  // -------------------
  fESD           = esdHandler->GetEvent();

  fMCEvent       = mcEvent;
  fStack         = fMCEvent->Stack();
#if 0
  fCentralityBin = centBin;

  // -- Create label arrays
  // ------------------------
  fLabelsRec[0] = new Int_t[fNTracks];
  if(!fLabelsRec[0]) {
    AliError("Cannot create fLabelsRec[0]");
    return -1;
  }

  fLabelsRec[1] = new Int_t[fNTracks];
  if(!fLabelsRec[1]) {
    AliError("Cannot create fLabelsRec[1] for TPC");
    return -1;
  }

  for(Int_t ii = 0; ii < fNTracks; ++ii) {
    fLabelsRec[0][ii] = 0;
    fLabelsRec[1][ii] = 0;
  }
#endif
  return 0;
}

//________________________________________________________________________
void AliAnalysisNetParticleDCA::ResetEvent() {
  // -- Reset event
  /*
  for (Int_t ii = 0; ii < 2 ; ++ii) {
    if (fLabelsRec[ii])
      delete[] fLabelsRec[ii];
    fLabelsRec[ii] = NULL;
  }
  */
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDCA::Process() {
  // -- Process event

  // -- Setup (clean, create and fill) MC labels
  // ---------------------------------------------
  //  FillMCLabels();

  // -- Fill  MC histograms for efficiency studies
  // -----------------------------------------------
  //  FillMCEffHist();

  return;
}      

/*
 * ---------------------------------------------------------------------------------
 *                                 Private Methods
 * ---------------------------------------------------------------------------------
 */
#if 0
//________________________________________________________________________
void AliAnalysisNetParticleDCA::CreateHistograms() {
  // -- Create histograms

  Double_t dCent[2] = {-0.5, 8.5};
  Int_t iCent       = 9;
  
  Double_t dEta[2]  = {-0.9, 0.9}; 
  Int_t iEta        = Int_t((dEta[1]-dEta[0]) / 0.075) +1 ; 

  Double_t dRap[2]  = {-0.5, 0.5}; 
  Int_t iRap        = Int_t((dRap[1]-dRap[0]) / 0.075) +1 ; 

  Double_t dPhi[2]  = {0.0, TMath::TwoPi()};
  Int_t iPhi        = 42;

  Double_t dPt[2]   = {0.1, 3.0};
  Int_t iPt         = Int_t((dPt[1]-dPt[0]) / 0.075);

  Double_t dSign[2] = {-1.5, 1.5};
  Int_t iSign       = 3;

  // ------------------------------------------------------------------
  // -- Create THnSparseF - DCA
  // ------------------------------------------------------------------

  //                              cent:   etaMC:     yMC:   phiMC:   ptMC:     sign:contPart:DCArAccepted:DCAr
  Int_t    binHnDCA[9] = {   iCent,    iEta,    iRap,    iPhi,    iPt,    iSign,       4,           2,  50 };     
  Double_t minHnDCA[9] = {dCent[0], dEta[0], dRap[0], dPhi[0], dPt[0], dSign[0],     0.5,        -0.5, -10.};
  Double_t maxHnDCA[9] = {dCent[1], dEta[1], dRap[1], dPhi[1], dPt[1], dSign[1],     4.5,         1.5,  10.};
  fHnDCA = new THnSparseF("fHnDCA", "cent:etaMC:yMC:phiMC:ptMC:sign:contPart:contSign:DCArAccepted:DCAr", 9, binHnDCA, minHnDCA, maxHnDCA);
  
  fHnDCA->Sumw2();
  fHnDCA->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnDCA->GetAxis(1)->SetTitle("#eta_{MC}");                    //  eta  [-0.9,0.9]
  fHnDCA->GetAxis(2)->SetTitle("#it{y}_{MC}");                  //  rapidity  [-0.5, 0.5]
  fHnDCA->GetAxis(3)->SetTitle("#varphi_{MC} (rad)");           //  phi  [ 0. ,2Pi]
  fHnDCA->GetAxis(4)->SetTitle("#it{p}_{T,MC} (GeV/#it{c})");   //  pT   [ 0.1,1.3]
  fHnDCA->GetAxis(5)->SetTitle("sign");                         //  -1 | 0 | +1 
  fHnDCA->GetAxis(6)->SetTitle("contPart");                     //  1  primary | 2 missId | 3 from WeakDecay | 4 p from Material
  fHnDCA->GetAxis(7)->SetTitle("DCArAccepted");                 //  0 not accepted | 1 accepted 
  fHnDCA->GetAxis(8)->SetTitle("DCAr");                         //  DCAr [-10,10]

  // ------------------------------------------------------------------
  
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDCA::FillMCDCA() {
  // Fill MC labels
  // Loop over ESD tracks and fill arrays with MC lables
  //  fLabelsRec[0] : all Tracks
  //  fLabelsRec[1] : all Tracks accepted by PID of TPC
  // Check every accepted track if correctly identified
  //  otherwise check for contamination

  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliESDtrack *track = fESD->GetTrack(idxTrack); 
    
    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted
    if (!fESDTrackCuts->AcceptTrack(track)) 
      continue;

    // -- Check if accepted in rapidity window
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted with thighter DCA cuts
    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;

    Int_t label  = TMath::Abs(track->GetLabel()); 
    
    // -- Fill Label of all reconstructed
    fLabelsRec[0][idxTrack] = label;

    // -- Check if accepted by PID from TPC or TPC+TOF
    Double_t pid[2];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- Fill Label of all reconstructed && recPid_TPC+TOF    
    fLabelsRec[1][idxTrack] = label;    
    
    // -- Check for contamination and fill contamination THnSparse
    CheckContTrack(label, track->GetSign(), idxTrack);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDCA::CheckDCATrack(Int_t label, Float_t sign, Float_t isDCArAccepted, Float_t dcar) {
  // Check if DCA of contamination or correctly identified
  // Fill contamination DCA THnSparse

  TParticle* particle = fStack->Particle(label);
  if (!particle)
    return;

  Bool_t isSecondaryFromWeakDecay = kFALSE;
  Bool_t isSecondaryFromMaterial  = kFALSE;
  
  Int_t contPart = 0;

  // -- Check if correctly identified 
  if (particle->GetPdgCode() == (sign*fPdgCode)) {

    // Check if is physical primary -> all ok 
    if (fStack->IsPhysicalPrimary(label))
      contPart = 1;    
    // -- Check if secondaries from weak decay
    else if(isSecondaryFromWeakDecay = fStack->IsSecondaryFromWeakDecay(label))
      contPart = 3;
    // -- Check if secondaries from material decay
    else if (isSecondaryFromMaterial  = fStack->IsSecondaryFromMaterial(label))
      contPart = 4;
  } 
  // -- MissIdentification
  else
    contPart = 2;
  
  Double_t hnDCA[9] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), sign, contPart, isDCArAccepted, dcar};
  fHnDCA->Fill(hnContDCA);
  
  return;
}


#endif

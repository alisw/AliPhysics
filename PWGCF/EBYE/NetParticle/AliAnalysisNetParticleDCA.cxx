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

/**
 * Class for for NetParticle Distributions
 * -- DCA distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

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
#endif
  // -- Create THnSparse Histograms
  // --------------------------------
  CreateHistograms();

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
  //  DCAr, DCAz, eta phi centrality

  // >  Type : 
  //     Primary MC p + ContMiss p
  //     2nd Material
  //     2nd Weak Decays

  //     IsAcceptedDCA

  //                           cent:   etaMC:     yMC:   phiMC:   ptMC:     sign:contPart:DCArAccepted: DCAr: DCAz
  Int_t    binHnDCA[10] = {   iCent,    iEta,    iRap,    iPhi,    iPt,    iSign,       4,           2,  50 , 50 };     
  Double_t minHnDCA[10] = {dCent[0], dEta[0], dRap[0], dPhi[0], dPt[0], dSign[0],     0.5,        -0.5, -10.,-10.};
  Double_t maxHnDCA[19] = {dCent[1], dEta[1], dRap[1], dPhi[1], dPt[1], dSign[1],     4.5,         1.5,  10., 10.};
  fHnDCA = new THnSparseF("fHnDCA", "cent:etaMC:yMC:phiMC:ptMC:sign:contPart:contSign:DCArAccepted:DCAr:DCAz", 10, binHnDCA, minHnDCA, maxHnDCA);
  
  fHnDCA->Sumw2();
  fHnDCA->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnDCA->GetAxis(1)->SetTitle("#eta_{Rec}");                    //  eta  [-0.9,0.9]
  fHnDCA->GetAxis(2)->SetTitle("#it{y}_{Rec}");                  //  rapidity  [-0.5, 0.5]
  fHnDCA->GetAxis(3)->SetTitle("#varphi_{Rec} (rad)");           //  phi  [ 0. ,2Pi]
  fHnDCA->GetAxis(4)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})");   //  pT   [ 0.1,1.3]
  fHnDCA->GetAxis(5)->SetTitle("sign");                         //  -1 | 0 | +1 

  fHnDCA->GetAxis(6)->SetTitle("contPart");                     //  1  primary | 2 missId | 3 from WeakDecay | 4 p from Material
  fHnDCA->GetAxis(7)->SetTitle("DCArAccepted");                 //  0 not accepted | 1 accepted 
  fHnDCA->GetAxis(8)->SetTitle("DCAr");                         //  DCAr [-10,10]
  fHnDCA->GetAxis(9)->SetTitle("DCAz");                         //  DCAz [-10,10]

  // ------------------------------------------------------------------
  
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDCA::FillDCA() {
  // Fill MC labels
  // Loop over ESD tracks and fill arrays with MC lables
  //  fLabelsRec[0] : all Tracks
  //  fLabelsRec[1] : all Tracks accepted by PID of TPC
  // Check every accepted track if correctly identified
  //  otherwise check for contamination
#if 0
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliESDtrack *track = fESD->GetTrack(idxTrack); 
    
    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted
    if (!fESDTrackCutsBkg->AcceptTrack(track)) 
      continue;

    // -- Check if accepted in rapidity window
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted by PID from TPC or TPC+TOF
    Double_t pid[2];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

  
    Bool_t isDCArAccepted = kTRUE;

    // -- Check if accepted with thighter DCA cuts
    if (!fHelper->IsTrackAcceptedDCA(track))
      isDCArAccepted = kFALSE;

    if (!fESDTrackCuts->AcceptTrack(track)) 
      isDCArAccepted = kFALSE;


    // ?    Int_t label  = ); 
    

    // -- Check for contamination and fill contamination THnSparse
    Int_t contIdx = CheckContTrack(TMath::Abs(track->GetLabel());
    
    Double_t hnDCA[10] = {fCentralityBin, track->Eta(), yP, track->Phi(), track->Pt(), sign, contIdx, isDCArAccepted, dcar, dcaz};
    fHnDCA->Fill(hnContDCA);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
#endif
  return;
}

//________________________________________________________________________
  Int_t AliAnalysisNetParticleDCA::CheckDCATrack(Int_t /*label*/) {
  // Check if DCA of contamination or correctly identified
  // Fill contamination DCA THnSparse
  // return  1  primary | 2 missId | 3 from WeakDecay | 4 p from Material | -1 unknown

  Int_t contIdx = -1;
#if 0
  TParticle* particle = fStack->Particle(label);
  if (!particle)
    return;

  Bool_t isSecondaryFromWeakDecay = kFALSE;
  Bool_t isSecondaryFromMaterial  = kFALSE;
  
  ////  // Check if is physical primary -> all ok 
  if (fStack->IsPhysicalPrimary(label)) {
    
    // -- Check if correctly identified 
    // Check if is physical primary -> all ok 
    if (particle->GetPdgCode() == (sign*fPdgCode))
      contIdx = 1;    
    // -- MissIdentification
    else 
      contIdx = 2;
  }
  // -- Check if secondaries from weak decay
  else if(isSecondaryFromWeakDecay = fStack->IsSecondaryFromWeakDecay(label))
    contIdx = 3;
  // -- Check if secondaries from material decay
  else if (isSecondaryFromMaterial  = fStack->IsSecondaryFromMaterial(label))
    contIdx = 4;
  } 
#endif
  return contIdx;
}




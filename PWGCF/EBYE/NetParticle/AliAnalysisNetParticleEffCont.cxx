//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"
#include "TDatabasePDG.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisNetParticleEffCont.h"

using namespace std;

// Task for NetParticle checks
// Author: Jochen Thaeder <jochen@thaeder.de>

ClassImp(AliAnalysisNetParticleEffCont)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleEffCont::AliAnalysisNetParticleEffCont() :
  fHelper(NULL),
  fPdgCode(2212),

  fESD(NULL), 
  fESDTrackCuts(NULL),

  fAOD(NULL), 
  fArrayMC(NULL),

  fCentralityBin(-1.),
  fNTracks(0),

  fAODtrackCutBit(1024),

  fStack(NULL),
  fMCEvent(NULL),

  fLabelsRec(NULL),

  fHnEff(NULL),
  fHnCont(NULL) {
  // Constructor   

  AliLog::SetClassDebugLevel("AliAnalysisNetParticleEffCont",10);
}

//________________________________________________________________________
AliAnalysisNetParticleEffCont::~AliAnalysisNetParticleEffCont() {
  // Destructor

  if (fLabelsRec[0])
    delete[] (fLabelsRec[0]);
  if (fLabelsRec[1])
    delete[] (fLabelsRec[1]);
  if (fLabelsRec)
    delete[] fLabelsRec;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::Initialize(AliESDtrackCuts *cuts , AliAnalysisNetParticleHelper* helper, Int_t trackCutBit) {

  // -- ESD track cuts
  // -------------------
  fESDTrackCuts     = cuts;

  // -- Get Helper
  // ---------------
  fHelper           = helper;

  // -- Get particle species / pdgCode
  // -------------------------
  fPdgCode          = AliPID::ParticleCode(fHelper->GetParticleSpecies());

  // -- MC Labels for efficiency
  // -----------------------------
  fLabelsRec        = new Int_t*[2];
  for (Int_t ii = 0; ii < 2 ; ++ii)
    fLabelsRec[ii] = NULL;

  // -- AOD track filter bit
  // -------------------------
  fAODtrackCutBit = trackCutBit;


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
Int_t AliAnalysisNetParticleEffCont::SetupEvent(AliESDInputHandler *esdHandler, AliMCEvent *mcEvent) {
  // -- Setup event

  // -- Reset of event
  // -------------------
  ResetEvent();

  // -- Setup of event
  // -------------------
  fESD           = esdHandler->GetEvent();
  fNTracks       = fESD->GetNumberOfTracks();

  fMCEvent       = mcEvent;
  fStack         = fMCEvent->Stack();

  fCentralityBin = fHelper->GetCentralityBin();

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

  return 0;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleEffCont::SetupEvent(AliAODInputHandler *aodHandler) {
  // -- Setup event

  // -- Reset of event
  // -------------------
  ResetEvent();

  // -- Setup of event
  // -------------------
  fAOD           = aodHandler->GetEvent();
  fNTracks       = fAOD->GetNumberOfTracks();

  fArrayMC       = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fArrayMC)
    AliFatal("No array of MC particles found !!!");

  fCentralityBin = fHelper->GetCentralityBin();

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

  return 0;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::ResetEvent() {
  // -- Reset event
  
  for (Int_t ii = 0; ii < 2 ; ++ii) {
    if (fLabelsRec[ii])
      delete[] fLabelsRec[ii];
    fLabelsRec[ii] = NULL;
  }

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::Process() {
  // -- Process event

  // -- Setup (clean, create and fill) MC labels
  // ---------------------------------------------
  if(fESD)
    FillMCLabels();
  else if(fAOD)
    FillMCLabelsAOD();
  
  // -- Fill  MC histograms for efficiency studies
  // -----------------------------------------------
  if(fESD)
    FillMCEffHist();
  else if(fAOD)
    FillMCEffHistAOD();

  return;
}      

/*
 * ---------------------------------------------------------------------------------
 *                                 Private Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::CreateHistograms() {
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
  // -- Create THnSparseF - Eff
  // ------------------------------------------------------------------

  //                           cent:   etaMC:     yMC:  phiMC:   ptMC:     sign:findable:recStatus:pidStatus:   etaRec:  phiRec:  ptRec:deltaEta: deltaPhi: deltaPt
  Int_t    binHnEff[15] = {   iCent,    iEta,    iRap,    iPhi,    iPt,    iSign,       2,      2  ,      2  ,    iEta,    iPhi,    iPt,    iEta, 2*iPhi+1, 2*iPt+1};
  Double_t minHnEff[15] = {dCent[0], dEta[0], dRap[0], dPhi[0], dPt[0], dSign[0],    -0.5,     -0.5,     -0.5, dEta[0], dPhi[0], dPt[0], dEta[0], -dPhi[1], -dPt[1]};
  Double_t maxHnEff[15] = {dCent[1], dEta[1], dRap[0], dPhi[1], dPt[1], dSign[1],     1.5,      1.5,      1.5, dEta[1], dPhi[1], dPt[1], dEta[1],  dPhi[1],  dPt[1]};

  fHnEff = new THnSparseF("fHnEff", "cent:etaMC:yMC:phiMC:ptMC:sign:findable:recStatus:pidStatus:etaRec:phiRec:ptRec:deltaEta:deltaPhi:deltaPt", 
			  15, binHnEff, minHnEff, maxHnEff);

  fHnEff->Sumw2();    

  fHnEff->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnEff->GetAxis(1)->SetTitle("#eta_{MC}");                    //  eta  [-0.9, 0.9]
  fHnEff->GetAxis(2)->SetTitle("#it{y}_{MC}");                  //  rapidity  [-0.5, 0.5]
  fHnEff->GetAxis(3)->SetTitle("#varphi_{MC} (rad)");           //  phi  [ 0. , 2Pi]
  fHnEff->GetAxis(4)->SetTitle("#it{p}_{T,MC} (GeV/#it{c})");   //  pt   [ 0.1, 3.0]
  
  fHnEff->GetAxis(5)->SetTitle("sign");                         //  -1 | 0 | +1 
  fHnEff->GetAxis(6)->SetTitle("findable");                     //  0 not findable      |  1 findable
  fHnEff->GetAxis(7)->SetTitle("recStatus");                    //  0 not reconstructed |  1 reconstructed
  fHnEff->GetAxis(8)->SetTitle("recPid");                       //  0 not accepted      |  1 accepted

  fHnEff->GetAxis(9)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnEff->GetAxis(10)->SetTitle("#varphi_{Rec} (rad)");         //  phi  [ 0. , 2Pi]
  fHnEff->GetAxis(11)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 0.1, 3.0]

  fHnEff->GetAxis(12)->SetTitle("#eta_{MC}-#eta_{Rec}");                      //  eta  [-0.9,0.9]
  fHnEff->GetAxis(13)->SetTitle("#varphi_{MC}-#varphi_{Rec} (rad)");          //  phi  [ 0. ,2Pi]
  fHnEff->GetAxis(14)->SetTitle("#it{p}_{T,MC}-#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 3.0,3.0]

  fHelper->BinLogAxis(fHnEff,  4);
  fHelper->BinLogAxis(fHnEff, 11);

  /* 
     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     
     creation -> findable -> reconstructed -> pid_TPC+TOF
        (1)         (2)          (3)            (4)      
                                                                      ||   findable | recStatus | recPid 
     1) all primary probeParticles_MC                                 ||      -           -         -
     2) all findable primary probeParticles_MC                        ||      x           -         -
     3) all reconstructed primary probeParticles_MC                   ||      x           x         -
     4) all reconstructed primary probeParticles_MC & recPid_TPC+TOF  ||      x           x         x
     
     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  */ 

  // ------------------------------------------------------------------
  // -- Create THnSparseF - Cont
  // ------------------------------------------------------------------

  //                            cent:   etaMC:    yMC:   phiMC:   ptMC:     sign:contPart: contSign:  etaRec:  phiRec:  ptRec:deltaEta: deltaPhi: deltaPt
  Int_t    binHnCont[14] = {   iCent,    iEta,    iRap,    iPhi,    iPt,    iSign,       8,    iSign,    iEta,    iPhi,    iPt,    iEta, 2*iPhi+1, 2*iPt+1};
  Double_t minHnCont[14] = {dCent[0], dEta[0], dRap[0], dPhi[0], dPt[0], dSign[0],     0.5, dSign[0], dEta[0], dPhi[0], dPt[0], dEta[0], -dPhi[1], -dPt[1]};
  Double_t maxHnCont[14] = {dCent[1], dEta[1], dRap[1], dPhi[1], dPt[1], dSign[1],     8.5, dSign[1], dEta[1], dPhi[1], dPt[1], dEta[1],  dPhi[1],  dPt[1]};
  fHnCont = new THnSparseF("fHnCont", "cent:etaMC:yMC:phiMC:ptMC:sign:contPart:contSign:etaRec:phiRec:ptRec:deltaEta:deltaPhi:deltaPt",
			   14, binHnCont, minHnCont, maxHnCont);

  fHnCont->Sumw2();    
  
  fHnCont->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnCont->GetAxis(1)->SetTitle("#eta_{MC}");                    //  eta  [-0.9,0.9]
  fHnCont->GetAxis(2)->SetTitle("#it{y}_{MC}");                  //  rapidity  [-0.5, 0.5]
  fHnCont->GetAxis(3)->SetTitle("#varphi_{MC} (rad)");           //  phi  [ 0. ,2Pi]
  fHnCont->GetAxis(4)->SetTitle("#it{p}_{T,MC} (GeV/#it{c})");   //  pT   [ 0.1,1.3]
  fHnCont->GetAxis(5)->SetTitle("sign");                         //  -1 | 0 | +1 
  
  fHnCont->GetAxis(6)->SetTitle("contPart");                     //  1 pi | 2 K | 3 p | 4 e | 5 mu | 6 other | 7 p from WeakDecay | 8 p from Material
  fHnCont->GetAxis(7)->SetTitle("contSign");                     //  -1 | 0 | +1 
   
  fHnCont->GetAxis(8)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnCont->GetAxis(9)->SetTitle("#varphi_{Rec} (rad)");          //  phi  [ 0. , 2Pi]
  fHnCont->GetAxis(10)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 0.1, 3.0]

  fHnCont->GetAxis(11)->SetTitle("#eta_{MC}-#eta_{Rec}");                      //  eta  [-0.9,0.9]
  fHnCont->GetAxis(12)->SetTitle("#varphi_{MC}-#varphi_{Rec} (rad)");          //  phi  [ 0. ,2Pi]
  fHnCont->GetAxis(13)->SetTitle("#it{p}_{T,MC}-#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 3.0,3.0]

  fHelper->BinLogAxis(fHnCont,  4);
  fHelper->BinLogAxis(fHnCont, 10);

  // ------------------------------------------------------------------
  
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::FillMCLabels() {
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
void AliAnalysisNetParticleEffCont::FillMCLabelsAOD() {
  // Fill MC labels for AOD analysis
  // Loop over ESD tracks and fill arrays with MC lables
  //  fLabelsRec[0] : all Tracks
  //  fLabelsRec[1] : all Tracks accepted by PID of TPC
  // Check every accepted track if correctly identified
  //  otherwise check for contamination

  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliAODTrack *track = fAOD->GetTrack(idxTrack); 
    
    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted
    if (!track->TestFilterBit(fAODtrackCutBit)) 
      continue;

    // -- Check if accepted in rapidity window
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted with thighter DCA cuts (not yet done for AODs XXX)
    //if (!fHelper->IsTrackAcceptedDCA(track))
    //  continue;

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
    CheckContTrackAOD(label, track->Charge(), idxTrack);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::CheckContTrack(Int_t label, Float_t sign, Int_t idxTrack) {
  // Check if particle is contamination or correctly identified for ESDs
  // Fill contamination THnSparse
  
  TParticle* particle = fStack->Particle(label);
  if (!particle)
    return;

  Bool_t isSecondaryFromWeakDecay = kFALSE;
  Bool_t isSecondaryFromMaterial  = kFALSE;
  
  // -- Check if correctly identified 
  if (particle->GetPdgCode() == (sign*fPdgCode)) {

    // Check if is physical primary -> all ok 
    if (fStack->IsPhysicalPrimary(label))
      return;
    
    // -- Check if secondaries from material or weak decay
    isSecondaryFromWeakDecay = fStack->IsSecondaryFromWeakDecay(label);
    isSecondaryFromMaterial  = fStack->IsSecondaryFromMaterial(label);
  } 
  
  // -- Get MC pdg
  Float_t contSign = 0.;
  if      (particle->GetPDG()->Charge() == 0.) contSign =  0.;
  else if (particle->GetPDG()->Charge() <  0.) contSign = -1.;	
  else if (particle->GetPDG()->Charge() >  0.) contSign =  1.;	

  // -- Get contaminating particle
  Float_t contPart = 0;
  if        (isSecondaryFromWeakDecay)                   contPart = 7; // probeParticle from WeakDecay
  else if   (isSecondaryFromMaterial)                    contPart = 8; // probeParticle from Material
  else {
    if      (TMath::Abs(particle->GetPdgCode()) ==  211) contPart = 1; // pion
    else if (TMath::Abs(particle->GetPdgCode()) ==  321) contPart = 2; // kaon
    else if (TMath::Abs(particle->GetPdgCode()) == 2212) contPart = 3; // proton
    else if (TMath::Abs(particle->GetPdgCode()) ==   11) contPart = 4; // electron
    else if (TMath::Abs(particle->GetPdgCode()) ==   13) contPart = 5; // muon
    else                                                 contPart = 6; // other
  }
  
  // -- Get Reconstructed values 
  Float_t etaRec = 0.;
  Float_t phiRec = 0.;
  Float_t ptRec  = 0.;

  AliESDtrack *track = fESD->GetTrack(idxTrack);
  if (track) {
    // if no track present (which should not happen)
    // -> pt = 0. , which is not in the looked at range
    
    // -- Get Reconstructed eta/phi/pt
    etaRec = track->Eta();
    phiRec = track->Phi();	  
    ptRec  = track->Pt();
  } 	

  // -- Fill THnSparse
  Double_t hnCont[14] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), sign, 
			 contPart, contSign, etaRec, phiRec, ptRec, 
			 particle->Eta()-etaRec, particle->Phi()-phiRec, particle->Pt()-ptRec};
  fHnCont->Fill(hnCont);
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::CheckContTrackAOD(Int_t label, Float_t sign, Int_t idxTrack) {
  // Check if particle is contamination or correctly identified for AODs
  // Fill contamination THnSparse
  
  AliAODMCParticle* particle = (AliAODMCParticle*)fArrayMC->At(label);

  if (!particle)
    return;

  Bool_t isSecondaryFromWeakDecay = kFALSE;
  Bool_t isSecondaryFromMaterial  = kFALSE;
  
  // -- Check if correctly identified 
  if (particle->GetPdgCode() == (sign*fPdgCode)) {
    
    // Check if is physical primary -> all ok 
    if(particle->IsPhysicalPrimary())
      return;
    
    // -- Check if secondaries from material or weak decay
    isSecondaryFromWeakDecay = particle->IsSecondaryFromWeakDecay();
    isSecondaryFromMaterial  = particle->IsSecondaryFromMaterial();
    
  } 
  
  // -- Get MC pdg 
  Float_t contSign = 0.;
  if      ((TDatabasePDG::Instance()->GetParticle(particle->PdgCode()))->Charge() == 0.) contSign =  0.;
  else if ((TDatabasePDG::Instance()->GetParticle(particle->PdgCode()))->Charge() <  0.) contSign = -1.;	
  else if ((TDatabasePDG::Instance()->GetParticle(particle->PdgCode()))->Charge() >  0.) contSign =  1.;	

  // -- Get contaminating particle
  Float_t contPart = 0;
  if        (isSecondaryFromWeakDecay)                   contPart = 7; // probeParticle from WeakDecay
  else if   (isSecondaryFromMaterial)                    contPart = 8; // probeParticle from Material
  else {
    if      (TMath::Abs(particle->GetPdgCode()) ==  211) contPart = 1; // pion
    else if (TMath::Abs(particle->GetPdgCode()) ==  321) contPart = 2; // kaon
    else if (TMath::Abs(particle->GetPdgCode()) == 2212) contPart = 3; // proton
    else if (TMath::Abs(particle->GetPdgCode()) ==   11) contPart = 4; // electron
    else if (TMath::Abs(particle->GetPdgCode()) ==   13) contPart = 5; // muon
    else                                                 contPart = 6; // other
  }
  
  // -- Get Reconstructed values 
  Float_t etaRec = 0.;
  Float_t phiRec = 0.;
  Float_t ptRec  = 0.;

  AliAODTrack *track = fAOD->GetTrack(idxTrack);
  
  if (track) {
    // if no track present (which should not happen)
    // -> pt = 0. , which is not in the looked at range
    
    // -- Get Reconstructed eta/phi/pt
    etaRec = track->Eta();
    phiRec = track->Phi();	  
    ptRec  = track->Pt();
  } 	

  // -- Fill THnSparse
  Double_t hnCont[14] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), sign, 
			 contPart, contSign, etaRec, phiRec, ptRec, 
			 particle->Eta()-etaRec, particle->Phi()-phiRec, particle->Pt()-ptRec};
  fHnCont->Fill(hnCont);
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::FillMCEffHist() {
  // Fill efficiency THnSparse for ESDs
  
  Int_t nPart  = fStack->GetNprimary();

  for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
    TParticle* particle = fStack->Particle(idxMC);

    // -- Check basic MC properties -> charged physical primary
    if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
      continue;

    // -- Check rapidity window
    Double_t yMC;
    if (!fHelper->IsParticleAcceptedRapidity(particle, yMC))
      continue;

    // -- Check if probeParticle / anti-probeParticle
    if (TMath::Abs(particle->GetPdgCode()) != fPdgCode)
      continue;
    
    // -- Get sign of particle
    Float_t sign      = (particle->GetPdgCode() < 0) ? -1. : 1.;

    // -- Get if particle is findable 
    Float_t findable  = Float_t(fHelper->IsParticleFindable(idxMC));

    // -- Get recStatus and pidStatus
    Float_t recStatus = 0.;
    Float_t recPid    = 0.;

    // -- Get Reconstructed values 
    Float_t etaRec = 0.;
    Float_t phiRec = 0.;
    Float_t ptRec  = 0.;

    // -- Loop over all labels
    for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {
      if (idxMC == fLabelsRec[0][idxRec]) {
	recStatus = 1.;
	
	if (idxMC == fLabelsRec[1][idxRec])
	  recPid = 1.;

	AliESDtrack *track = fESD->GetTrack(idxRec);
	if (track) {
	  // if no track present (which should not happen)
	  // -> pt = 0. , which is not in the looked at range

	  // -- Get Reconstructed eta/phi/pt
	  etaRec = track->Eta();
	  phiRec = track->Phi();	  
	  ptRec  = track->Pt();
	} 	
	break;
      }
    } // for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {  
    
    // -- Fill THnSparse
    Double_t hnEff[15] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), sign, 
			  findable, recStatus, recPid, etaRec, phiRec, ptRec, 
			  particle->Eta()-etaRec, particle->Phi()-phiRec, particle->Pt()-ptRec};
    fHnEff->Fill(hnEff);

  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::FillMCEffHistAOD() {
  // Fill efficiency THnSparse for AODs
  
  Int_t nPart  = fArrayMC->GetEntriesFast();

  for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
    
    AliAODMCParticle* particle = (AliAODMCParticle*)fArrayMC->At(idxMC);
             
    // -- Check basic MC properties -> charged physical primary
    if (!fHelper->IsParticleAcceptedBasicCharged(particle))
      continue;

    // -- Check rapidity window
    Double_t yMC;
    if (!fHelper->IsParticleAcceptedRapidity((TParticle*)particle, yMC))
      continue;

    // -- Check if probeParticle / anti-probeParticle
    if (TMath::Abs(particle->GetPdgCode()) != fPdgCode)
      continue;
    
    // -- Get sign of particle
    Float_t sign      = (particle->GetPdgCode() < 0) ? -1. : 1.;

    // -- Get if particle is findable 
    Float_t findable  = 1.;// Float_t(fHelper->IsParticleFindable(idxMC)); XXX

    // -- Get recStatus and pidStatus
    Float_t recStatus = 0.;
    Float_t recPid    = 0.;

    // -- Get Reconstructed values 
    Float_t etaRec = 0.;
    Float_t phiRec = 0.;
    Float_t ptRec  = 0.;

    // -- Loop over all labels
    for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {
      if (idxMC == fLabelsRec[0][idxRec]) {
    	recStatus = 1.;
	
    	if (idxMC == fLabelsRec[1][idxRec])
    	  recPid = 1.;

    	AliVTrack *track = NULL;
    	if(fESD)
    	  track = fESD->GetTrack(idxRec);
    	else if(fAOD)
    	  track = fAOD->GetTrack(idxRec);

    	if (track) {
    	  // if no track present (which should not happen)
    	  // -> pt = 0. , which is not in the looked at range

    	  // -- Get Reconstructed eta/phi/pt
    	  etaRec = track->Eta();
    	  phiRec = track->Phi();	  
    	  ptRec  = track->Pt();
    	} 	
    	break;
      }
    } // for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {  
    
    // -- Fill THnSparse
    Double_t hnEff[15] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), sign, 
    			  findable, recStatus, recPid, etaRec, phiRec, ptRec, 
    			  particle->Eta()-etaRec, particle->Phi()-phiRec, particle->Pt()-ptRec};
    fHnEff->Fill(hnEff);

  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  return;
}

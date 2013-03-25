//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"

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

/**
 * Class for for NetParticle Distributions
 * -- Efficiency and contaminations for netParticle distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

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

  if (fLabelsRec) {
    if (fLabelsRec[0])
      delete[] (fLabelsRec[0]);
    if (fLabelsRec[1])
      delete[] (fLabelsRec[1]);
    if (fLabelsRec)
      delete[] fLabelsRec;
  }
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
Int_t AliAnalysisNetParticleEffCont::SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler, AliMCEvent *mcEvent) {
  // -- Setup event

  // -- Reset of event
  // -------------------
  ResetEvent();

  // -- Setup of event
  // -------------------
  
  // -- Get ESD objects
  if(esdHandler) {
    fESD     = esdHandler->GetEvent(); 
    fNTracks = fESD->GetNumberOfTracks();
  }
  
  // -- Get AOD objects
  else if(aodHandler) {
    fAOD     = aodHandler->GetEvent();
    fNTracks = fAOD->GetNumberOfTracks();
    
    fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fArrayMC)
      AliFatal("No array of MC particles found !!!"); // MW  no AliFatal use return values
  }           

  // -- Get MC objects
  if (mcEvent) {
    fMCEvent     = mcEvent;
    if (fMCEvent)
      fStack     = fMCEvent->Stack();
  }

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
    AliError("Cannot create fLabelsRec[1] for PID");
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
  FillMCLabels();
 
  // -- Fill  MC histograms for efficiency studies
  // -----------------------------------------------
  FillMCEffHist();

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

  // ------------------------------------------------------------------
  // -- Create THnSparseF - Eff
  // ------------------------------------------------------------------

  Int_t    binHnEff[13] = {AliAnalysisNetParticleHelper::fgkfHistNBinsCent, AliAnalysisNetParticleHelper::fgkfHistNBinsEta,       //     cent |     etaMC
			   AliAnalysisNetParticleHelper::fgkfHistNBinsRap,  AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,       //      yMC |     phiMC
			   AliAnalysisNetParticleHelper::fgkfHistNBinsPt,   AliAnalysisNetParticleHelper::fgkfHistNBinsSign,      //     ptMC |      sign
			   2,      2  ,      2  ,                                                                                 // findable | recStatus  | pidStatus
			   AliAnalysisNetParticleHelper::fgkfHistNBinsEta,  AliAnalysisNetParticleHelper::fgkfHistNBinsRap,       //   etaRec |      yRec
			   AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,  AliAnalysisNetParticleHelper::fgkfHistNBinsPt};       //   phiRec |     ptRec
  
  Double_t minHnEff[13] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeEta[0], 
			   AliAnalysisNetParticleHelper::fgkfHistRangeRap[0],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[0], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0], 
			   -0.5,     -0.5,     -0.5,
			   AliAnalysisNetParticleHelper::fgkfHistRangeEta[0],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePhi[0],  AliAnalysisNetParticleHelper::fgkfHistRangePt[0]}; 

  Double_t maxHnEff[13] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[1], AliAnalysisNetParticleHelper::fgkfHistRangeEta[1], 
			   AliAnalysisNetParticleHelper::fgkfHistRangeRap[1],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePt[1],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1], 
			   1.5,      1.5,      1.5,
			   AliAnalysisNetParticleHelper::fgkfHistRangeEta[1],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePhi[1],  AliAnalysisNetParticleHelper::fgkfHistRangePt[1]}; 

  fHnEff = new THnSparseF("fHnEff", "cent:etaMC:yMC:phiMC:ptMC:sign:findable:recStatus:pidStatus:etaRec:yRec:phiRec:ptRec", 
			  13, binHnEff, minHnEff, maxHnEff);
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
  fHnEff->GetAxis(10)->SetTitle("#it{y}_{Rec}");                //  rapidity  [-0.5, 0.5]
  fHnEff->GetAxis(11)->SetTitle("#varphi_{Rec} (rad)");         //  phi  [ 0. , 2Pi]
  fHnEff->GetAxis(12)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 0.1, 3.0]

  fHelper->BinLogAxis(fHnEff,  4);
  fHelper->BinLogAxis(fHnEff, 12);

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

  Int_t    binHnCont[13] = {AliAnalysisNetParticleHelper::fgkfHistNBinsCent, AliAnalysisNetParticleHelper::fgkfHistNBinsEta,       //     cent |    etaMC
			    AliAnalysisNetParticleHelper::fgkfHistNBinsRap,  AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,       //      yMC |    phiMC
			    AliAnalysisNetParticleHelper::fgkfHistNBinsPt,   AliAnalysisNetParticleHelper::fgkfHistNBinsSign,      //     ptMC |     sign
			    8,                                               AliAnalysisNetParticleHelper::fgkfHistNBinsSign,      // contPart | contSign
			    AliAnalysisNetParticleHelper::fgkfHistNBinsEta,  AliAnalysisNetParticleHelper::fgkfHistNBinsRap,       //   etaRec |     yRec
			    AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,  AliAnalysisNetParticleHelper::fgkfHistNBinsPt};       //   phiRec |    ptRec
  
  Double_t minHnCont[13] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeEta[0], 
			    AliAnalysisNetParticleHelper::fgkfHistRangeRap[0],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[0], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0], 
			    0.5,                                                AliAnalysisNetParticleHelper::fgkfHistRangeSign[0], 
			    AliAnalysisNetParticleHelper::fgkfHistRangeEta[0],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePhi[0],  AliAnalysisNetParticleHelper::fgkfHistRangePt[0]}; 
  
  Double_t maxHnCont[13] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[1], AliAnalysisNetParticleHelper::fgkfHistRangeEta[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangeRap[1],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePt[1],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1],
			    8.5,                                                AliAnalysisNetParticleHelper::fgkfHistRangeSign[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangeEta[1],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePhi[1],  AliAnalysisNetParticleHelper::fgkfHistRangePt[1]}; 

  fHnCont = new THnSparseF("fHnCont", "cent:etaMC:yMC:phiMC:ptMC:sign:contPart:contSign:etaRec:yRec:phiRec:ptRec",
			   12, binHnCont, minHnCont, maxHnCont);
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
  fHnCont->GetAxis(9)->SetTitle("#it{y}_{Rec}");                 //  rapidity  [-0.5, 0.5]
  fHnCont->GetAxis(10)->SetTitle("#varphi_{Rec} (rad)");         //  phi  [ 0. , 2Pi]
  fHnCont->GetAxis(11)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 0.1, 3.0]

  fHelper->BinLogAxis(fHnCont,  4);
  fHelper->BinLogAxis(fHnCont, 11);

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
    
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

    // -- Check if track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted - ESD
    if (fESD && !fESDTrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))
      continue;
    
    // -- Check if accepted - AOD
    if (fAOD){
      AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
      
      if (!trackAOD) {
	AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
	continue;
      }
      if (!trackAOD->TestFilterBit(fAODtrackCutBit))
	continue;
    }

    // -- Check if accepted in rapidity window -- for identified particles
    //    for charged particles -- eta check is done in the trackcuts
    Double_t yP;
    if (fHelper->GetUsePID() && !fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted with thighter DCA cuts
    // -- returns kTRUE for AODs for now : MW
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
    CheckContTrack(label, track->Charge(), idxTrack);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::CheckContTrack(Int_t label, Float_t sign, Int_t idxTrack) {
  // Check if particle is contamination or correctly identified for ESDs and AODs
  // Check for missidentified primaries and secondaries
  // Fill contamination THnSparse

  AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
  if (!particle)
    return;

  Bool_t isPhysicalPrimary = (fESD) ? fStack->IsPhysicalPrimary(label): (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
  
  // -- Check if correctly identified 
  //    > return if correctly identified -> all ok, no action neededin this method
  //    > if PID required check -> for the correct (signed pdgcode) particle
  //    > no PID just check for primary 
  if (fHelper->GetUsePID()) {
    if (particle->PdgCode() == (sign*fPdgCode))
      if (isPhysicalPrimary)
	return;
  }
  else {
    if (isPhysicalPrimary)
      return;
  }

  // -- Check if secondaries from material or weak decay
  Bool_t isSecondaryFromWeakDecay = (fESD) ? fStack->IsSecondaryFromWeakDecay(label) : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
  Bool_t isSecondaryFromMaterial  = (fESD) ? fStack->IsSecondaryFromMaterial(label)  : (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();

  // -- Get PDG Charge
  Float_t contSign = 0.;
  if      (particle->Charge() == 0.) contSign =  0.;
  else if (particle->Charge() <  0.) contSign = -1.;	
  else if (particle->Charge() >  0.) contSign =  1.;	

  // -- Get contaminating particle
  Float_t contPart = 0;
  if        (isSecondaryFromWeakDecay)                contPart = 7; // probeParticle from WeakDecay
  else if   (isSecondaryFromMaterial)                 contPart = 8; // probeParticle from Material
  else {
    if      (TMath::Abs(particle->PdgCode()) ==  211) contPart = 1; // pion
    else if (TMath::Abs(particle->PdgCode()) ==  321) contPart = 2; // kaon
    else if (TMath::Abs(particle->PdgCode()) == 2212) contPart = 3; // proton
    else if (TMath::Abs(particle->PdgCode()) ==   11) contPart = 4; // electron
    else if (TMath::Abs(particle->PdgCode()) ==   13) contPart = 5; // muon
    else                                              contPart = 6; // other
  }
  
  // -- Get Reconstructed values 
  Float_t etaRec = 0.;
  Float_t phiRec = 0.;
  Float_t ptRec  = 0.;
  Double_t yRec  = 0.;

  AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

  if (track) {
    // if no track present (which should not happen)
    // -> pt = 0. , which is not in the looked at range
    
    // -- Get Reconstructed eta/phi/pt/y
    etaRec = track->Eta();
    phiRec = track->Phi();	  
    ptRec  = track->Pt();
    fHelper->IsTrackAcceptedRapidity(track, yRec); // yRec = y for identified particles | yRec = eta for charged particles
  } 	

  // -- Fill THnSparse
  Double_t hnCont[12] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), sign, 
			 contPart, contSign, etaRec, yRec, phiRec, ptRec};
  fHnCont->Fill(hnCont);
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::FillMCEffHist() {
  // Fill efficiency THnSparse for ESDs

  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Int_t nPart  = (fESD) ? fStack->GetNprimary() : fArrayMC->GetEntriesFast();

  for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
    AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(idxMC) : static_cast<AliVParticle*>(fArrayMC->At(idxMC));

    // -- Check basic MC properties -> charged physical primary
    if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
      continue;

    // -- Check if accepted in rapidity window -- for identified particles
    Double_t yMC;
    if (fHelper->GetUsePID() && !fHelper->IsParticleAcceptedRapidity(particle, yMC))
      continue;

    // -- Check if accepted in eta window -- for charged particles
    if (!fHelper->GetUsePID() && TMath::Abs(particle->Eta()) > etaRange[1])
      continue;

    // -- Check if probeParticle / anti-probeParticle 
    //    > skip check if PID is not required
    if (fHelper->GetUsePID() && TMath::Abs(particle->PdgCode()) != fPdgCode)
      continue;
    
    // -- Get sign of particle
    Float_t sign      = (particle->PdgCode() < 0) ? -1. : 1.;

    // -- Get if particle is findable --- not availible for AODs yet
    Float_t findable  = (fESD) ? Float_t(fHelper->IsParticleFindable(idxMC)) : 1.;

    // -- Get recStatus and pidStatus
    Float_t recStatus = 0.;
    Float_t recPid    = 0.;

    // -- Get Reconstructed values 
    Float_t etaRec = 0.;
    Float_t phiRec = 0.;
    Float_t ptRec  = 0.;
    Double_t yRec  = 0.;

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
	  
          // -- Get Reconstructed eta/phi/pt/y
          etaRec = track->Eta();
          phiRec = track->Phi();         
          ptRec  = track->Pt();
          fHelper->IsTrackAcceptedRapidity(track, yRec); // yRec = y for identified particles | yRec = eta for charged particles
        }      
        break;
      }
    } // for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {  
    
    // -- Fill THnSparse
    Double_t hnEff[15] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), sign, 
			  findable, recStatus, recPid, etaRec, yRec, phiRec, ptRec};
    fHnEff->Fill(hnEff);

  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  return;
}

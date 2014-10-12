//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"

#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisNetParticleEffCont.h"

using namespace std;

/**
 * Class for NetParticle Distributions
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
  AliAnalysisNetParticleBase("EffCont", "EffCont"),
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
void AliAnalysisNetParticleEffCont::Process() {
  // -- Process event

  // -- Setup (clean, create and fill) MC labels
  FillMCLabels();
 
  // -- Fill  MC histograms for efficiency studies
  FillMCEffHist();

  return;
}      

/*
 * ---------------------------------------------------------------------------------
 *                                Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::Init() {
  // -- Init eventwise

  // -- MC Labels for efficiency
  fLabelsRec        = new Int_t*[2];
  for (Int_t ii = 0; ii < 2 ; ++ii)
    fLabelsRec[ii] = NULL;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::CreateHistograms() {
  // -- Create histograms

  // ------------------------------------------------------------------
  // -- Create THnSparse - Eff
  // ------------------------------------------------------------------

  Int_t    binHnEff[19] = {AliAnalysisNetParticleHelper::fgkfHistNBinsCent,    AliAnalysisNetParticleHelper::fgkfHistNBinsEta,       //     cent |     etaMC
			   AliAnalysisNetParticleHelper::fgkfHistNBinsRap,     AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,       //      yMC |     phiMC
			   AliAnalysisNetParticleHelper::fgkfHistNBinsPt,      AliAnalysisNetParticleHelper::fgkfHistNBinsSign,      //     ptMC |    signMC
			   2,      2  ,      2  ,                                                                                    // findable | recStatus  | pidStatus
			   AliAnalysisNetParticleHelper::fgkfHistNBinsEta,     AliAnalysisNetParticleHelper::fgkfHistNBinsRap,       //   etaRec |      yRec
			   AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,     AliAnalysisNetParticleHelper::fgkfHistNBinsPt,        //   phiRec |     ptRec
			   AliAnalysisNetParticleHelper::fgkfHistNBinsSign,                                                          //  signRec
  			   AliAnalysisNetParticleHelper::fgkfHistNBinsEta,     AliAnalysisNetParticleHelper::fgkfHistNBinsRap,       // deltaEta | deltaY
			   2*AliAnalysisNetParticleHelper::fgkfHistNBinsPhi+1, 2*AliAnalysisNetParticleHelper::fgkfHistNBinsPt+1,    // deltaPhi | deltaPt
			   AliAnalysisNetParticleHelper::fgkfHistNBinsSign+2};                                                       // deltaSign

  
  Double_t minHnEff[19] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeEta[0], 
			   AliAnalysisNetParticleHelper::fgkfHistRangeRap[0],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[0], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0], 
			   -0.5,     -0.5,     -0.5,
			   AliAnalysisNetParticleHelper::fgkfHistRangeEta[0],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePhi[0],  AliAnalysisNetParticleHelper::fgkfHistRangePt[0],
			   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0],
			   AliAnalysisNetParticleHelper::fgkfHistRangeEta[0],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0], 
			   -1.*AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], -1.*AliAnalysisNetParticleHelper::fgkfHistRangePt[1],
			   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0]-1.};

  Double_t maxHnEff[19] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[1], AliAnalysisNetParticleHelper::fgkfHistRangeEta[1], 
			   AliAnalysisNetParticleHelper::fgkfHistRangeRap[1],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePt[1],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1], 
			   1.5,      1.5,      1.5,
			   AliAnalysisNetParticleHelper::fgkfHistRangeEta[1],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePhi[1],  AliAnalysisNetParticleHelper::fgkfHistRangePt[1],
			   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1],
			   AliAnalysisNetParticleHelper::fgkfHistRangeEta[1],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1], 
			   AliAnalysisNetParticleHelper::fgkfHistRangePhi[1],  AliAnalysisNetParticleHelper::fgkfHistRangePt[1],
			   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1]+1.};

  fHnEff = new THnSparseF("hnEff", "cent:etaMC:yMC:phiMC:ptMC:signMC:findable:recStatus:pidStatus:etaRec:yRec:phiRec:ptRec:signRec:deltaEta:deltaY:deltaPhi:deltaPt:deltaSign", 
			  19, binHnEff, minHnEff, maxHnEff);
  fHnEff->Sumw2();    

  fHnEff->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnEff->GetAxis(1)->SetTitle("#eta_{MC}");                    //  eta  [-0.9, 0.9]
  fHnEff->GetAxis(2)->SetTitle("#it{y}_{MC}");                  //  rapidity  [-0.5, 0.5]
  fHnEff->GetAxis(3)->SetTitle("#varphi_{MC} (rad)");           //  phi  [ 0. , 2Pi]
  fHnEff->GetAxis(4)->SetTitle("#it{p}_{T,MC} (GeV/#it{c})");   //  pT   [ 0.2, 2.3]
  
  fHnEff->GetAxis(5)->SetTitle("sign");                         //  -1 | 0 | +1 
  fHnEff->GetAxis(6)->SetTitle("findable");                     //  0 not findable      |  1 findable
  fHnEff->GetAxis(7)->SetTitle("recStatus");                    //  0 not reconstructed |  1 reconstructed
  fHnEff->GetAxis(8)->SetTitle("recPid");                       //  0 not accepted      |  1 accepted

  fHnEff->GetAxis(9)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnEff->GetAxis(10)->SetTitle("#it{y}_{Rec}");                //  rapidity  [-0.5, 0.5]
  fHnEff->GetAxis(11)->SetTitle("#varphi_{Rec} (rad)");         //  phi  [ 0. , 2Pi]
  fHnEff->GetAxis(12)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 0.2, 2.3]
  fHnEff->GetAxis(13)->SetTitle("sign");                        //  -1 | 0 | +1 

  fHnEff->GetAxis(14)->SetTitle("#eta_{MC}-#eta_{Rec}");                      //  eta  [-0.9, 0.9]
  fHnEff->GetAxis(15)->SetTitle("#it{y}_{MC}-#it{y}_{Rec}");                  //  rapidity  [-0.5, 0.5]
  fHnEff->GetAxis(16)->SetTitle("#varphi_{MC}-#varphi_{Rec} (rad)");          //  phi  [ -2Pi , 2Pi]
  fHnEff->GetAxis(17)->SetTitle("#it{p}_{T,MC}-#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ -2.3, 2.3]
  fHnEff->GetAxis(18)->SetTitle("sign_{MC}-sign_{Rec}");                      //  -2 | 0 | +2 

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
  // -- Create THnSparse - Cont
  // ------------------------------------------------------------------

  Int_t    binHnCont[17] = {AliAnalysisNetParticleHelper::fgkfHistNBinsCent, AliAnalysisNetParticleHelper::fgkfHistNBinsEta,       //     cent |    etaMC
			    AliAnalysisNetParticleHelper::fgkfHistNBinsRap,  AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,       //      yMC |    phiMC
			    AliAnalysisNetParticleHelper::fgkfHistNBinsPt,   AliAnalysisNetParticleHelper::fgkfHistNBinsSign,      //     ptMC |   signMC=contSign
			    8,                                                                                                     // contPart | 
			    AliAnalysisNetParticleHelper::fgkfHistNBinsEta,  AliAnalysisNetParticleHelper::fgkfHistNBinsRap,       //   etaRec |     yRec
			    AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,  AliAnalysisNetParticleHelper::fgkfHistNBinsPt,        //   phiRec |    ptRec
			    AliAnalysisNetParticleHelper::fgkfHistNBinsSign,                                                       //  signRec  
			    AliAnalysisNetParticleHelper::fgkfHistNBinsEta,     AliAnalysisNetParticleHelper::fgkfHistNBinsRap,    // deltaEta | deltaY
			    2*AliAnalysisNetParticleHelper::fgkfHistNBinsPhi+1, 2*AliAnalysisNetParticleHelper::fgkfHistNBinsPt+1, // deltaPhi | deltaPt
			    AliAnalysisNetParticleHelper::fgkfHistNBinsSign+2};                                                    // deltaSign
  
  Double_t minHnCont[17] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeEta[0], 
			    AliAnalysisNetParticleHelper::fgkfHistRangeRap[0],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[0], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0], 
			    0.5,                                            
			    AliAnalysisNetParticleHelper::fgkfHistRangeEta[0],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePhi[0],  AliAnalysisNetParticleHelper::fgkfHistRangePt[0],
			    AliAnalysisNetParticleHelper::fgkfHistRangeSign[0],
			    AliAnalysisNetParticleHelper::fgkfHistRangeEta[0],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[0], 
			    -1.*AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], -1.*AliAnalysisNetParticleHelper::fgkfHistRangePt[1],
			    AliAnalysisNetParticleHelper::fgkfHistRangeSign[0]-1.};

  
  Double_t maxHnCont[17] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[1], AliAnalysisNetParticleHelper::fgkfHistRangeEta[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangeRap[1],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePt[1],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1],
			    8.5,                        
			    AliAnalysisNetParticleHelper::fgkfHistRangeEta[1],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePhi[1],  AliAnalysisNetParticleHelper::fgkfHistRangePt[1],
			    AliAnalysisNetParticleHelper::fgkfHistRangeSign[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangeEta[1],  AliAnalysisNetParticleHelper::fgkfHistRangeRap[1], 
			    AliAnalysisNetParticleHelper::fgkfHistRangePhi[1],  AliAnalysisNetParticleHelper::fgkfHistRangePt[1],
			    AliAnalysisNetParticleHelper::fgkfHistRangeSign[1]+1.};

  fHnCont = new THnSparseF("hnCont", "cent:etaMC:yMC:phiMC:ptMC:signMC:contPart:etaRec:yRec:phiRec:ptRec:signRec:deltaEta:deltaY:deltaPhi:deltaPt:deltaSign",
			   17, binHnCont, minHnCont, maxHnCont);
  fHnCont->Sumw2();    

  fHnCont->GetAxis(0)->SetTitle("centrality");                   //  0-5|5-10|10-20|20-30|30-40|40-50|50-60|60-70|70-80|80-90 --> 10 bins
  fHnCont->GetAxis(1)->SetTitle("#eta_{MC}");                    //  eta  [-0.9,0.9]
  fHnCont->GetAxis(2)->SetTitle("#it{y}_{MC}");                  //  rapidity  [-0.5, 0.5]
  fHnCont->GetAxis(3)->SetTitle("#varphi_{MC} (rad)");           //  phi  [ 0. ,2Pi]
  fHnCont->GetAxis(4)->SetTitle("#it{p}_{T,MC} (GeV/#it{c})");   //  pT   [ 0.2,2.3]
  fHnCont->GetAxis(5)->SetTitle("sign");                         //  -1 | 0 | +1 
  
  fHnCont->GetAxis(6)->SetTitle("contPart");                     //  1 pi | 2 K | 3 p | 4 e | 5 mu | 6 other | 7 p from WeakDecay | 8 p from Material
 
  fHnCont->GetAxis(7)->SetTitle("#eta_{Rec}");                   //  eta  [-0.9, 0.9]
  fHnCont->GetAxis(8)->SetTitle("#it{y}_{Rec}");                 //  rapidity  [-0.5, 0.5]
  fHnCont->GetAxis(9)->SetTitle("#varphi_{Rec} (rad)");          //  phi  [ 0. , 2Pi]
  fHnCont->GetAxis(10)->SetTitle("#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ 0.2, 2.3]
  fHnCont->GetAxis(11)->SetTitle("sign");                        //  -1 | 0 | +1 

  fHnCont->GetAxis(12)->SetTitle("#eta_{MC}-#eta_{Rec}");                      //  eta  [-0.9, 0.9]
  fHnCont->GetAxis(13)->SetTitle("#it{y}_{MC}-#it{y}_{Rec}");                  //  rapidity  [-0.5, 0.5]
  fHnCont->GetAxis(14)->SetTitle("#varphi_{MC}-#varphi_{Rec} (rad)");          //  phi  [ -2Pi , 2Pi]
  fHnCont->GetAxis(15)->SetTitle("#it{p}_{T,MC}-#it{p}_{T,Rec} (GeV/#it{c})"); //  pt   [ -2.3, 2.3]
  fHnCont->GetAxis(16)->SetTitle("sign_{MC}-sign_{Rec}");                      //  -2 | 0 | +2 

  fHelper->BinLogAxis(fHnCont,  4);
  fHelper->BinLogAxis(fHnCont, 10);

  return;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleEffCont::Setup() {
  // -- Setup eventwise

  // -- Create label arrays
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
void AliAnalysisNetParticleEffCont::Reset() {
  // -- Reset eventwise

  for (Int_t ii = 0; ii < 2 ; ++ii) {
    if (fLabelsRec[ii])
      delete[] fLabelsRec[ii];
    fLabelsRec[ii] = NULL;
  }
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Private Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::FillMCLabels() {
  // Fill MC labels
  // Loop over ESD tracks and fill arrays with MC lables
  //  fLabelsRec[0] : all Tracks
  //  fLabelsRec[1] : all Tracks accepted by PID of TPC
  // Check every accepted track if correctly identified
  //  otherwise check for contamination

  // -- Get ranges for AOD particles
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);

  // -- Track Loop
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

      // -- Check if in pT and eta range (is done in ESDTrackCuts for ESDs)
      if(!(track->Pt() > ptRange[0] && track->Pt() <= ptRange[1] && TMath::Abs(track->Eta()) <= etaRange[1]))
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
    Double_t pid[3];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- Fill Label of all reconstructed && recPid_TPC+TOF    
    fLabelsRec[1][idxTrack] = label;    
    
    // -- Check for contamination and fill contamination THnSparse
    CheckContTrack(track);

  } // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleEffCont::CheckContTrack(AliVTrack *track) {
  // Check if particle is contamination or correctly identified for ESDs and AODs
  // Check for missidentified primaries and secondaries
  // Fill contamination THnSparse

  Int_t label     = TMath::Abs(track->GetLabel()); 
  Float_t signRec = track->Charge();

  
  AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(label) : static_cast<AliVParticle*>(fArrayMC->At(label));
  if (!particle)
    return;

  Bool_t isPhysicalPrimary = (fESD) ? fStack->IsPhysicalPrimary(label): (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
  
  // -- Check if correctly identified 
  //    > return if correctly identified -> all ok, no action neededin this method
  //    > if PID required check -> for the correct (signed pdgcode) particle
  //    > no PID just check for primary 
  if (fHelper->GetUsePID()) {
    if (particle->PdgCode() == (signRec*fPdgCode))
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

  // -- Get PDG Charge of contaminating particle
  Float_t signMC = 0.;
  if      (particle->Charge() == 0.) signMC =  0.;
  else if (particle->Charge() <  0.) signMC = -1.;	
  else if (particle->Charge() >  0.) signMC =  1.;	

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
  
  // -- Get Reconstructed y
  //    yRec = y for identified particles | yRec = eta for charged particles
  Double_t yRec  = 0.;
  fHelper->IsTrackAcceptedRapidity(track, yRec); 

  Double_t deltaPhi = particle->Phi()-track->Phi();
  if (TMath::Abs(deltaPhi) > TMath::TwoPi()) {
    if (deltaPhi < 0)
      deltaPhi += TMath::TwoPi();
    else
      deltaPhi -= TMath::TwoPi();
  }

  // -- Fill THnSparse
  Double_t hnCont[17] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), signMC, 
			 contPart, track->Eta(), yRec, track->Phi(), track->Pt(), signRec,
			 particle->Eta()-track->Eta(), particle->Y()-yRec, deltaPhi, particle->Pt()-track->Pt(), signMC-signRec};
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
    Float_t signMC    = (particle->PdgCode() < 0) ? -1. : 1.;

    // -- Get if particle is findable --- not availible for AODs yet
    Float_t findable  = (fESD) ? Float_t(fHelper->IsParticleFindable(idxMC)) : 1.;

    // -- Get recStatus and pidStatus
    Float_t recStatus = 0.;
    Float_t recPid    = 0.;

    // -- Get Reconstructed values 
    Float_t etaRec  = 0.;
    Float_t phiRec  = 0.;
    Float_t ptRec   = 0.;
    Double_t yRec   = 0.;
    Float_t signRec = 0.;

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
	  
          // -- Get Reconstructed values
          etaRec  = track->Eta();
          phiRec  = track->Phi();         
          ptRec   = track->Pt();
	  signRec = track->Charge();
          fHelper->IsTrackAcceptedRapidity(track, yRec); // yRec = y for identified particles | yRec = eta for charged particles
        }      
        break;
      }
    } // for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {  

    Double_t deltaPhi = particle->Phi()-phiRec;
    if (TMath::Abs(deltaPhi) > TMath::TwoPi()) {
      if (deltaPhi < 0)
	deltaPhi += TMath::TwoPi();
      else
    	deltaPhi -= TMath::TwoPi();
    }
    
    // -- Fill THnSparse
    Double_t hnEff[19] = {fCentralityBin, particle->Eta(), particle->Y(), particle->Phi(), particle->Pt(), signMC, 
			  findable, recStatus, recPid, etaRec, yRec, phiRec, ptRec, signRec,
			  particle->Eta()-etaRec, particle->Y()-yRec, deltaPhi, particle->Pt()-ptRec, signMC-signRec};
    fHnEff->Fill(hnEff);

  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  return;
}

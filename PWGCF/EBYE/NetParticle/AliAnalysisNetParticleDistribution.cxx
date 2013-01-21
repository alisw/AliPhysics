//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"
#include "TSystem.h" 
#include "TProfile.h" 
#include "TH2F.h" 
#include "TH3F.h" 
#include "TFile.h" 
#include "TPRegexp.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliTracker.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"

#include "AliAnalysisNetParticleDistribution.h"

using namespace std;

// Task for NetParticle checks
// Author: Jochen Thaeder <jochen@thaeder.de>

ClassImp(AliAnalysisNetParticleDistribution)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleDistribution::AliAnalysisNetParticleDistribution() :
  fHelper(NULL),

  fOutList(NULL),

  fESDHandler(NULL),
  fPIDResponse(NULL),
  fESD(NULL),
  fAODHandler(NULL),
  fAOD(NULL),
  fIsMC(kFALSE),
  fMCEvent(NULL),
  fStack(NULL),
  fESDTrackCuts(NULL),
  fEtaMax(0.9),
  fPtRange(NULL),
  fAODtrackCutBit(1024),
  fNp(NULL),
  fNCorrNp(1),
  fCorrNp(NULL),
  fNMCNp(5),
  fMCNp(NULL),
  fNControlMCNp(5),
  fControlMCNp(NULL),

  fHnTrackUnCorr(NULL) {
  // Constructor   
  
  AliLog::SetClassDebugLevel("AliAnalysisNetParticleDistribution",10);
}

//________________________________________________________________________
AliAnalysisNetParticleDistribution::~AliAnalysisNetParticleDistribution() {
  // Destructor

  if (fNp) delete[] fNp;  

  for (Int_t ii = 0; ii < fNCorrNp; ++ii) 
    if (fCorrNp[ii]) delete[] fCorrNp[ii];
  if (fCorrNp) delete[] fCorrNp;

  for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    if (fMCNp[ii]) delete[] fMCNp[ii];
  if (fMCNp) delete[] fMCNp;

  for (Int_t ii = 0; ii < fNControlMCNp; ++ii) 
    if (fControlMCNp[ii]) delete[] fControlMCNp[ii];
  if (fControlMCNp) delete[] fControlMCNp;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::Initialize(AliAnalysisNetParticleHelper* helper, AliESDtrackCuts* cuts, 
						     Bool_t isMC, Float_t *ptRange, Float_t etaMax, Int_t trackCutBit, Int_t nCorrNp) {
  // -- Initialize
  
  fHelper = helper;
  fESDTrackCuts = cuts;
  fIsMC = isMC;
  fPtRange = ptRange;
  fEtaMax = etaMax;
  fAODtrackCutBit = trackCutBit;
  fNCorrNp = nCorrNp;

  // ------------------------------------------------------------------
  // -- N particles / N anti-particles
  // ------------------------------------------------------------------
  //  Np          : arr[particle]
  //  CorrNp      : arr[corrSet][particle]
  //  MCNp        : arr[corrSet][particle] - MC
  //  ControlMCNp : arr[corrSet][particle] - Control MC

  fNp     = new Float_t[2];

  fCorrNp = new Float_t*[fNCorrNp];
  for (Int_t ii = 0 ; ii < fNCorrNp; ++ii)
    fCorrNp[ii] = new Float_t[2];

  fMCNp = new Float_t*[fNMCNp];
  for (Int_t ii = 0 ; ii < fNMCNp; ++ii)
    fMCNp[ii] = new Float_t[2];

  fControlMCNp = new Float_t*[fNControlMCNp];
  for (Int_t ii = 0 ; ii < fNControlMCNp; ++ii)
    fControlMCNp[ii] = new Float_t[2];

  ResetEvent();

  return 0;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::CreateHistograms(TList* outList) {
  // -- Add histograms to outlist

  fOutList = outList;
  
  // ------------------------------------------------------------------
  // -- Create net particle histograms
  // ------------------------------------------------------------------

  AddHistSet("fHDist", "Uncorrected");
  AddHistSet("fHDistCorr0", "Corrected [without cross section correction]");
  AddHistSet("fHDistCorr1", "Corrected [with cross section correction]");
  AddHistSet("fHDistCorr2", "Corrected [only cross section correction]");

  if (fIsMC) {
    AddHistSet("fHMCrapidity", "MC primary in |y| < 0.5");
    //    AddHistSet("fHMCptMin",    "MC primary in |y| + #it{p}_{T} > 0.1");
    //    AddHistSet("fHMCpt",       Form("MC primary in |y| < 0.5 + #it{p}_{T} [%.1f,%.1f]", fPtRange[0], fPtRange[1]));
    //    AddHistSet("fHMCeta",      Form("MC primary in |y| < 0.5 + |#eta| < %.1f", fEtaMax));
    AddHistSet("fHMCetapt",    Form("MC primary in |y| < 0.5 + |#eta| < %.1f + #it{p}_{T} [%.1f,%.1f]", fEtaMax, fPtRange[0], fPtRange[1]));
    
    //    AddHistSet("fHControlMCLambdarapidity", "Control MC Lambda primary in |y| < 0.5");
    //    AddHistSet("fHControlMCLambdaptMin",    "Control MC Lambda primary in |y| + #it{p}_{T} > 0.1");
    //    AddHistSet("fHControlMCLambdapt",       Form("Control MC primary in |y| < 0.5 + #it{p}_{T} [%.1f,%.1f]", fPtRange[0], fPtRange[1]));
    //    AddHistSet("fHControlMCLambdaeta",      Form("Control MC primary in |y| < 0.5 + |#eta| < %.1f", fEtaMax));
    //    AddHistSet("fHControlMCLambdaetapt",    Form("Control MC primary in |y| < 0.5 + |#eta| < %.1f + #it{p}_{T} [%.1f,%.1f]", fEtaMax, fPtRange[0], fPtRange[1]));
  }

  // ------------------------------------------------------------------
  // -- Get Probe Particle Container
  // ------------------------------------------------------------------

  Double_t dCent[2] = {-0.5, 8.5};
  Int_t iCent       = 9;
  
  Double_t dEta[2]  = {-0.9, 0.9}; // -> 37 bins
  Int_t iEta        = Int_t((dEta[1]-dEta[0]) / 0.075) +1 ; 

  Double_t dRap[2]  = {-0.5, 0.5}; 
  Int_t iRap        = Int_t((dRap[1]-dRap[0]) / 0.075) +1 ; 

  Double_t dPhi[2]  = {0.0, TMath::TwoPi()};
  Int_t iPhi        = 42;

  Double_t dPt[2]   = {0.1, 3.0}; 
  Int_t iPt         = Int_t((dPt[1]-dPt[0]) / 0.075);

  Double_t dSign[2] = {-1.5, 1.5};
  Int_t iSign       = 3;

  //                          cent:     pt:     sign:     eta:     phi:       y
  Int_t    binShort[6] = {   iCent,    iPt,    iSign,    iEta,    iPhi,    iRap};
  Double_t minShort[6] = {dCent[0], dPt[0], dSign[0], dEta[0], dPhi[0], dRap[0]};
  Double_t maxShort[6] = {dCent[1], dPt[1], dSign[1], dEta[1], dPhi[1], dRap[1]};

  // -- UnCorrected
  fOutList->Add(new THnSparseF("fHnTrackUnCorr", "cent:pt:sign:eta:phi:y", 6, binShort, minShort, maxShort));  
  fHnTrackUnCorr = static_cast<THnSparseF*>(fOutList->Last());
  fHnTrackUnCorr->GetAxis(0)->SetTitle("centrality");
  fHnTrackUnCorr->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHnTrackUnCorr->GetAxis(2)->SetTitle("sign");
  fHnTrackUnCorr->GetAxis(3)->SetTitle("#eta");
  fHnTrackUnCorr->GetAxis(4)->SetTitle("#varphi");
  fHnTrackUnCorr->GetAxis(5)->SetTitle("#it{y}");

  fHelper->BinLogAxis(fHnTrackUnCorr, 1);

  // ------------------------------------------------------------------

  return;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::SetupEvent(AliESDInputHandler *esdHandler, AliAODInputHandler *aodHandler,  AliMCEvent *mcEvent) {
  // -- Setup Event
  
  ResetEvent();

  // -- Get ESD objects
  if(esdHandler){
    fESDHandler  = esdHandler;
    fPIDResponse = esdHandler->GetPIDResponse();
    fESD         = fESDHandler->GetEvent();
  }

  // -- Get AOD objects
  else if(aodHandler){
    fAODHandler  = aodHandler;
    fPIDResponse = aodHandler->GetPIDResponse();
    fAOD         = fAODHandler->GetEvent();
  }

  // -- Get MC objects
  fMCEvent     = mcEvent;
  if (fMCEvent)
    fStack     = fMCEvent->Stack();

  return 0;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::ResetEvent() {
  // -- Reset event
  
  // -- Reset ESD Event
  fESD       = NULL;

  // -- Reset AOD Event
  fAOD       = NULL;

  // -- Reset MC Event
  if (fIsMC)
    fMCEvent = NULL;

  // -- Reset N particles/anti-particles
  for (Int_t jj = 0; jj < 2; ++jj)
    fNp[jj] = 0.;
  
  // -- Reset N corrected particles/anti-particles
  for (Int_t ii = 0; ii < fNCorrNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fCorrNp[ii][jj] = 0.;

  // -- Reset N MC particles/anti-particles
  for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fMCNp[ii][jj] = 0.;

  // -- Reset N control MC particles/anti-particles
  for (Int_t ii = 0; ii < fNControlMCNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fControlMCNp[ii][jj] = 0.;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::Process() {
  // -- Process NetParticle Distributions

  // -- Fill ESD tracks
  if (fESD) 
    ProcessESDTracks();
  
  // -- Fill AOD tracks
  else if (fAOD) 
    ProcessAODTracks();
    
  // -- Fill MC truth particles (missing for AOD XXX)
  if (fIsMC)  {
    ProcessStackParticles();
    //    ProcessStackControlParticles();
  }

  return 0;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::UpdateMinPtForTOFRequired() {
  // -- Update MinPtForTOFRequired, using the pT log-scale

  if (fHelper && fHnTrackUnCorr) {
    Float_t minPtForTOF = fHelper->GetMinPtForTOFRequired();
    TH1D *h =  static_cast<TH1D*>(fHnTrackUnCorr->Projection(1));
      
    for (Int_t ii = 0; ii < h->GetNbinsX(); ++ii)
      if (h->GetBinLowEdge(ii) <= minPtForTOF && h->GetBinLowEdge(ii) + h->GetBinWidth(ii) >  minPtForTOF) {
	minPtForTOF = h->GetBinLowEdge(ii) + h->GetBinWidth(ii);
	fHelper->SetMinPtForTOFRequired(minPtForTOF);
      }
  }

  return ;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessESDTracks() {
  // -- Process ESD tracks and fill histograms

  for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {
    AliESDtrack *track = fESD->GetTrack(idxTrack); 

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Check track cuts
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    // -- Check if charged track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted
    if (!fESDTrackCuts->AcceptTrack(track)) 
      continue;

    // -- Check if accepted in rapidity window
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted bt PID from TPC or TPC+TOF
    Double_t pid[2];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- Check if accepted with thighter DCA cuts
    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Fill Probe Particle
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    Double_t aTrack[6] = {
      Double_t(fHelper->GetCentralityBin()),       //  0 centrality 
      track->Pt(),                    //  1 pt
      track->GetSign(),               //  2 sign
      track->Eta(),                   //  3 eta
      track->Phi(),                   //  4 phi
      yP                              //  5 rapidity
    };
    
    fHnTrackUnCorr->Fill(aTrack);
    
    // -- Count particle / anti-particle 
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle

    Int_t idxPart = (track->GetSign() < 0) ? 0 : 1;
    fNp[idxPart] += 1.;
    
    for (Int_t ii = 0; ii < fNCorrNp; ++ii) 
      fCorrNp[ii][idxPart] += fHelper->GetTrackbyTrackCorrectionFactor(aTrack, ii);      
    
  } // for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  // -- Uncorrected
  FillHistSet("fHDist", fNp);
    
  // -- Corrected 
  for (Int_t ii = 0 ; ii < fNCorrNp; ++ii) 
    FillHistSet(Form("fHDistCorr%d", ii), fCorrNp[ii]); 

  return 0;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessAODTracks() {
  // -- Process AOD tracks and fill histograms

  for (Int_t idxTrack = 0; idxTrack < fAOD->GetNumberOfTracks(); ++idxTrack) {
    AliAODTrack *track = (AliAODTrack*)fAOD->GetTrack(idxTrack); 

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Check track cuts
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    // -- Check if charged track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;

    // -- Check if accepted
    if(!track->TestFilterBit(fAODtrackCutBit)) 
      continue;
    
    // -- Check if in pT and eta range (is done in ESDTrackCuts for ESDs)
    if(!(track->Pt() > fPtRange[0] && track->Pt() < fPtRange[1] && TMath::Abs(track->Eta()) <= fEtaMax))
      continue;

    // -- Check if accepted in rapidity window
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted bt PID from TPC or TPC+TOF
    Double_t pid[2];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- Check if accepted with thighter DCA cuts XXX
    // if (fHelper->IsTrackAcceptedDCA(track))
    //  continue;
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Fill Probe Particle
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    Double_t aTrack[6] = {
      Double_t(fHelper->GetCentralityBin()),       //  0 centrality 
      track->Pt(),                    //  1 pt
      track->Charge(),                //  2 sign
      track->Eta(),                   //  3 eta
      track->Phi(),                   //  4 phi
      yP                              //  5 rapidity
    };
    
    fHnTrackUnCorr->Fill(aTrack);
    
    // -- Count particle / anti-particle 
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle

    Int_t idxPart = (track->Charge() < 0) ? 0 : 1;
    fNp[idxPart] += 1.;
    
    for (Int_t ii = 0; ii < fNCorrNp; ++ii) 
      fCorrNp[ii][idxPart] += fHelper->GetTrackbyTrackCorrectionFactor(aTrack, ii);      
    
  } // for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  // -- Uncorrected
  FillHistSet("fHDist", fNp);
    
  // -- Corrected 
  for (Int_t ii = 0 ; ii < fNCorrNp; ++ii) 
    FillHistSet(Form("fHDistCorr%d", ii), fCorrNp[ii]); 

  return 0;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessStackParticles() {
  // -- Process primary particles from the stack and fill histograms

  Int_t pdgCode    = AliPID::ParticleCode(fHelper->GetParticleSpecies());
  
  for (Int_t idxMC = 0; idxMC < fStack->GetNprimary(); ++idxMC) {
    TParticle* particle = fStack->Particle(idxMC);
    if (!particle) 
      continue;

    // -- Check basic MC properties -> charged physical primary
    if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
      continue;
    
    // -- Check if particle / anti-particle
    if (fHelper->GetUsePID() && TMath::Abs(particle->GetPdgCode()) != pdgCode)
      continue;
    
    // -- Get particle : 0 anti-particle / 1 particle
    Int_t idxPart = (particle->GetPdgCode() < 0) ? 0 : 1;

    // >> NOW only anti-particle / particle 
    // >> With idxPart
    
    // -- Check rapidity window
    Double_t yMC;
    if (!fHelper->GetUsePID()) 
      yMC = fEtaMax;
    if (!fHelper->IsParticleAcceptedRapidity(particle, yMC))
      continue;
    
    fMCNp[0][idxPart] += 1.;         // -> MCrapidity
    
    // -- Check acceptance
    if (particle->Pt() > 0.1 )
      fMCNp[1][idxPart] += 1.;       // -> MCptMin
    
    if (particle->Pt() > fPtRange[0] && particle->Pt() <= fPtRange[1])
      fMCNp[2][idxPart] += 1.;       // -> MCpt
    
    if (TMath::Abs(particle->Eta()) <= fEtaMax)
      fMCNp[3][idxPart] += 1.;       // -> MCeta
    
    if (particle->Pt() > fPtRange[0] && particle->Pt() <= fPtRange[1] && TMath::Abs(particle->Eta()) <= fEtaMax)
      fMCNp[4][idxPart] += 1.;       // -> MCetapt
    
  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSet("fHMCrapidity", fMCNp[0]);
  //  FillHistSet("fHMCptMin",    fMCNp[1]);
  //  FillHistSet("fHMCpt",       fMCNp[2]);
  //  FillHistSet("fHMCeta",      fMCNp[3]);
  FillHistSet("fHMCetapt",    fMCNp[4]);

  return 0;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessStackControlParticles() {
  // -- Process primary particles from the stack and fill histograms

  Int_t pdgCode    = fHelper->GetControlParticleCode();
  Bool_t isNeutral  = fHelper->IsControlParticleNeutral();
  //  const Char_t* name = fHelper->GetControlParticleName().Data();

  for (Int_t idxMC = 0; idxMC < fStack->GetNprimary(); ++idxMC) {
    TParticle* particle = fStack->Particle(idxMC);
    if (!particle) 
      continue;
    
    // -- Check basic MC properties -> neutral or charged physical primary
    if (isNeutral) {
      if (!fHelper->IsParticleAcceptedBasicNeutral(particle, idxMC))
	continue;
    }
    else {
      if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
	continue;
    }
    
    // -- Check if particle / anti-particle
    if (TMath::Abs(particle->GetPdgCode()) != pdgCode)
      continue;

    // -- Get particle : 0 anti-particle / 1 particle
    Int_t idxPart = (particle->GetPdgCode() < 0) ? 0 : 1;

    // >> NOW only anti-particle / particle 
    // >> With idxPart
    
    // -- Check rapidity window
    Double_t yMC;
    if (!fHelper->GetUsePID()) 
      yMC = fEtaMax;
    if (!fHelper->IsParticleAcceptedRapidity(particle, yMC))
      continue;

    fControlMCNp[0][idxPart] += 1.;         // -> ControlMCrapidity

    // -- Check acceptance
    if (particle->Pt() > 0.1 )
      fControlMCNp[1][idxPart] += 1.;       // -> ControlMCptMin
    
    if (particle->Pt() > fPtRange[0] && particle->Pt() <= fPtRange[1])
      fControlMCNp[2][idxPart] += 1.;       // -> ControlMCpt
    
    if (TMath::Abs(particle->Eta()) <= fEtaMax)
      fControlMCNp[3][idxPart] += 1.;       // -> ControlMCeta
    
    if (particle->Pt() > fPtRange[0] && particle->Pt() > fPtRange[1] && TMath::Abs(particle->Eta()) <= fEtaMax)
      fControlMCNp[4][idxPart] += 1.;       // -> ControlMCetapt
  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  //  FillHistSet(Form("fHControlMC%srapidity", name), fControlMCNp[0], 0);
  //  FillHistSet(Form("fHControlMC%sptMin",    name), fControlMCNp[1], 1);
  //  FillHistSet(Form("fHControlMC%spt",       name), fControlMCNp[2], 2);
  //  FillHistSet(Form("fHControlMC%seta",      name), fControlMCNp[3], 3);
  //  FillHistSet(Form("fHControlMC%setapt",    name), fControlMCNp[4], 4);

  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Helper Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void  AliAnalysisNetParticleDistribution::AddHistSet(const Char_t *name, const Char_t *title)  {
  // -- Add histogram sets for particle and anti-particle

  fOutList->Add(new TList);
  TList *list = static_cast<TList*>(fOutList->Last());
  list->SetName(name) ;
  list->SetOwner(kTRUE);

  TString sName(name);

  TString sPtTitle("");
  if (!sName.Contains("fHMC") || !sName.Contains("fHControlMC"))
    sPtTitle += Form("#it{p}_{T} [%.1f,%.1f]", fPtRange[0], fPtRange[1]);
  
  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    list->Add(new TH2F(Form("%s%s", name, fHelper->GetParticleName(idxPart).Data()), 
		       Form("%s %s Dist %s;Centrality;N_{%s}", title, fHelper->GetParticleTitle(idxPart).Data(), sPtTitle.Data(), 
			    fHelper->GetParticleTitleLatex(idxPart).Data()), 
		       24, -0.5, 11.5, 2601, -0.5, 2600.49));
  } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
   
  list->Add(new TH2F(Form("%sNet%s", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s Net %s Dist %s;Centrality;N_{%s} - N_{%s}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			  fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()), 24, -0.5, 11.5, 601, -300.5, 300.49));

  list->Add(new TProfile(Form("%sNet%sM", name, fHelper->GetParticleName(1).Data()), 
			 Form("%s Net %s Dist %s;Centrality;N_{%s} - N_{%s}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			      fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()),  24, -0.5, 11.5));
  list->Add(new TProfile(Form("%sNet%s2M", name, fHelper->GetParticleName(1).Data()), 
			 Form("%s (Net %s)^{2} Dist %s;Centrality;(N_{%s} - N_{%s})^{2}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			      fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()), 24, -0.5, 11.5));
  list->Add(new TProfile(Form("%sNet%s3M", name, fHelper->GetParticleName(1).Data()), 
			 Form("%s (Net %s)^{3} Dist %s;Centrality;(N_{%s} - N_{%s})^{3}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			      fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()), 24, -0.5, 11.5));
  list->Add(new TProfile(Form("%sNet%s4M", name, fHelper->GetParticleName(1).Data()), 
			 Form("%s (Net %s)^{4} Dist %s;Centrality;(N_{%s} - N_{%s})^{4}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			      fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()), 24, -0.5, 11.5));
  list->Add(new TProfile(Form("%sNet%s5M", name, fHelper->GetParticleName(1).Data()), 
			 Form("%s (Net %s)^{5} Dist %s;Centrality;(N_{%s} - N_{%s})^{5}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			      fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()), 24, -0.5, 11.5));
  list->Add(new TProfile(Form("%sNet%s6M", name, fHelper->GetParticleName(1).Data()), 
			 Form("%s (Net %s)^{6} Dist %s;Centrality;(N_{%s} - N_{%s})^{6}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			      fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()), 24, -0.5, 11.5));
  
  list->Add(new TH2F(Form("%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s (Net %s)/ Sum Dist %s;Centrality;(N_{%s} - N_{%s})/(N_{%s} + N_{%s})", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			  fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data(),fHelper->GetParticleTitleLatex(1).Data(),
			  fHelper->GetParticleTitleLatex(0).Data()), 
		     24, -0.5, 11.5, 801, -20.5, 20.49));
  
  if (sName.Contains("fHControlMC")) {
    for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
      list->Add(new TH2F(Form("%s%sOver%s", name, fHelper->GetParticleName(idxPart).Data(), fHelper->GetControlParticleName(idxPart).Data()), 
			 Form("%s %s / %s Dist %s;Centrality;N_{%s}/N_{Tracks}", 
			      title, fHelper->GetParticleTitle(idxPart).Data(), fHelper->GetControlParticleTitle(idxPart).Data(), sPtTitle.Data(), 
			      fHelper->GetParticleTitleLatex(idxPart).Data()), 
			 24, -0.5, 11.5, 101, 0., 1.));
    } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
  }
 
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::FillHistSet(const Char_t *name, Float_t *np, Int_t controlIdx)  {
  // -- Add histogram sets for particle and anti-particle

  TList *list = static_cast<TList*>(fOutList->FindObject(name));

  Float_t centralityBin = fHelper->GetCentralityBin();

  (static_cast<TH2F*>(list->FindObject(Form("%s%s", name, fHelper->GetParticleName(0).Data()))))->Fill(centralityBin, np[0]);
  (static_cast<TH2F*>(list->FindObject(Form("%s%s", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, np[1]);

  Float_t sumNp    = np[1]+np[0];
  Float_t deltaNp  = np[1]-np[0];
  Float_t deltaNp2 = deltaNp * deltaNp;
  Float_t deltaNp3 = deltaNp2 * deltaNp;

  (static_cast<TH2F*>(list->FindObject(Form("%sNet%s",  name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp);

  (static_cast<TProfile*>(list->FindObject(Form("%sNet%sM",  name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp);
  (static_cast<TProfile*>(list->FindObject(Form("%sNet%s2M", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp2);
  (static_cast<TProfile*>(list->FindObject(Form("%sNet%s3M", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp2*deltaNp);
  (static_cast<TProfile*>(list->FindObject(Form("%sNet%s4M", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp2*deltaNp2);
  (static_cast<TProfile*>(list->FindObject(Form("%sNet%s5M", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp3*deltaNp2);
  (static_cast<TProfile*>(list->FindObject(Form("%sNet%s6M", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp3*deltaNp3);

  (static_cast<TH2F*>(list->FindObject(Form("%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp/sumNp);

  TString sName(name);
  if (sName.Contains("fHControlMC") && controlIdx >= 0) {
    (static_cast<TH2F*>(list->FindObject(Form("%s%sOverlambdabar", name, fHelper->GetParticleName(0).Data()))))->Fill(centralityBin, np[0]/fControlMCNp[controlIdx][0]);
    (static_cast<TH2F*>(list->FindObject(Form("%s%sOverlambda",    name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, np[1]/fControlMCNp[controlIdx][1]);
  }

  return;
}


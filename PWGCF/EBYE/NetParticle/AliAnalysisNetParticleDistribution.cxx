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

/**
 * Class for for NetParticle Distributions
 * -- Create input for distributions
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

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
  fOrder(8),
  fAODtrackCutBit(1024),
  fNNp(5),
  fNp(NULL),
  fNMCNp(5),
  fMCNp(NULL),
  fHnTrackUnCorr(NULL) {
  // Constructor   
  
  AliLog::SetClassDebugLevel("AliAnalysisNetParticleDistribution",10);
}

//________________________________________________________________________
AliAnalysisNetParticleDistribution::~AliAnalysisNetParticleDistribution() {
  // Destructor

  for (Int_t ii = 0; ii < fNNp; ++ii) 
    if (fNp[ii]) delete[] fNp[ii];
  if (fNp) delete[] fNp;

  for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    if (fMCNp[ii]) delete[] fMCNp[ii];
  if (fMCNp) delete[] fMCNp;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::Initialize(AliAnalysisNetParticleHelper* helper, AliESDtrackCuts* cuts, 
						     Bool_t isMC, Float_t *ptRange, Float_t etaMax, Int_t trackCutBit) {
  // -- Initialize
  
  fHelper = helper;
  fESDTrackCuts = cuts;
  fIsMC = isMC;
  fPtRange = ptRange;
  fEtaMax = etaMax;
  fAODtrackCutBit = trackCutBit;

  // ------------------------------------------------------------------
  // -- N particles / N anti-particles
  // ------------------------------------------------------------------
  //  Np          : arr[set][particle]
  //  MCNp        : arr[set][particle] - MC

  fNp = new Float_t*[fNNp];
  for (Int_t ii = 0 ; ii < fNNp; ++ii)
    fNp[ii] = new Float_t[2];

  fMCNp = new Float_t*[fNMCNp];
  for (Int_t ii = 0 ; ii < fNMCNp; ++ii)
    fMCNp[ii] = new Float_t[2];

  ResetEvent();

  return 0;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::CreateHistograms(TList* outList) {
  // -- Add histograms to outlist

  fOutList = outList;
  
  // ------------------------------------------------------------------
  // -- Get Probe Particle Container
  // ------------------------------------------------------------------

  Int_t    binHnUnCorr[6] = {AliAnalysisNetParticleHelper::fgkfHistNBinsCent, AliAnalysisNetParticleHelper::fgkfHistNBinsEta,          // cent |  eta
			     AliAnalysisNetParticleHelper::fgkfHistNBinsRap,  AliAnalysisNetParticleHelper::fgkfHistNBinsPhi,          //    y |  phi
			     AliAnalysisNetParticleHelper::fgkfHistNBinsPt,   AliAnalysisNetParticleHelper::fgkfHistNBinsSign};        //   pt | sign
  
  Double_t minHnUnCorr[6] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeEta[0], 
			     AliAnalysisNetParticleHelper::fgkfHistRangeRap[0],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[0], 
			     AliAnalysisNetParticleHelper::fgkfHistRangePt[0],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[0]};
  
  Double_t maxHnUnCorr[6] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[1], AliAnalysisNetParticleHelper::fgkfHistRangeEta[1], 
			     AliAnalysisNetParticleHelper::fgkfHistRangeRap[1],  AliAnalysisNetParticleHelper::fgkfHistRangePhi[1], 
			     AliAnalysisNetParticleHelper::fgkfHistRangePt[1],   AliAnalysisNetParticleHelper::fgkfHistRangeSign[1]};
  
  // -- UnCorrected
  fOutList->Add(new THnSparseF("fHnTrackUnCorr", "cent:eta:y:phi:pt:sign", 6, binHnUnCorr, minHnUnCorr, maxHnUnCorr));  
  fHnTrackUnCorr = static_cast<THnSparseF*>(fOutList->Last());
  fHnTrackUnCorr->Sumw2(); 
  fHnTrackUnCorr->GetAxis(0)->SetTitle("centrality");
  fHnTrackUnCorr->GetAxis(1)->SetTitle("#eta");
  fHnTrackUnCorr->GetAxis(2)->SetTitle("#it{y}");
  fHnTrackUnCorr->GetAxis(3)->SetTitle("#varphi");
  fHnTrackUnCorr->GetAxis(4)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHnTrackUnCorr->GetAxis(5)->SetTitle("sign");

  fHelper->BinLogAxis(fHnTrackUnCorr, 1);

  // ------------------------------------------------------------------
  // -- Create net particle histograms
  // ------------------------------------------------------------------

  TString sTitle("");
  sTitle = (fHelper->GetUsePID()) ? Form("Uncorrected in |y| < %.1f", fHelper->GetRapidityMax()) : Form("Uncorrected in |#eta| < %.1f", fEtaMax);

  AddHistSet("fHDist",       Form("%s + #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), fPtRange[0], fPtRange[1]));
  AddHistSet("fHDistTPC",    Form("%s + #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), fPtRange[0], fHelper->GetMinPtForTOFRequired()));
  AddHistSet("fHDistphi",    Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				  sTitle.Data(), fPtRange[0], fPtRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  AddHistSet("fHDistTPCphi", Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				  sTitle.Data(), fPtRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  
  if (fIsMC) {
    TString sMCTitle("");
    sMCTitle = (fHelper->GetUsePID()) ?  Form("MC primary in |y| < %.1f", fHelper->GetRapidityMax()) : Form("MC primary in |#eta| < %.1f", fEtaMax);

    AddHistSet("fHMC",         Form("%s", sTitle.Data()));
    AddHistSet("fHMCpt",       Form("%s + #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), fPtRange[0], fPtRange[1]));
    AddHistSet("fHMCptTPC",    Form("%s + #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), fPtRange[0], fHelper->GetMinPtForTOFRequired()));
    AddHistSet("fHMCptphi",    Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				    sMCTitle.Data(), fPtRange[0], fPtRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
    AddHistSet("fHMCptTPCphi", Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				    sMCTitle.Data(), fPtRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  }

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
  for (Int_t ii = 0; ii < fNNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fNp[ii][jj] = 0.;
  
  // -- Reset N MC particles/anti-particles
  for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fMCNp[ii][jj] = 0.;
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
  if (fIsMC)
    ProcessStackParticles();

  return 0;
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

    // -- Check if accepted in rapidity window -- for identified particles
    //    for charged - eta check is done in the trackcuts
    Double_t yP;
    if (fHelper->GetUsePID() && !fHelper->IsTrackAcceptedRapidity(track, yP))
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

    if (!fHelper->GetUsePID())  
      yP = track->Eta();
    
    Double_t aTrack[6] = {
      Double_t(fHelper->GetCentralityBin()),  //  0 centrality 
      track->Eta(),                           //  1 eta
      yP,                                     //  2 rapidity
      track->Phi(),                           //  3 phi
      track->Pt(),                            //  4 pt
      track->GetSign()                        //  5 sign
    };
    
    fHnTrackUnCorr->Fill(aTrack);
    
    // -- Count particle / anti-particle 
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle

    Int_t idxPart = (track->Charge() < 0) ? 0 : 1;

    // -- in pt Range
    fNp[0][idxPart] += 1.;

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[1][idxPart] += 1.;
      
    // -- check phi range
    if(!fHelper->IsTrackAcceptedPhi(track))
      continue;

    // -- in pt Range
    fNp[2][idxPart] += 1.;

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[3][idxPart] += 1.;
  
  } // for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSet("fHDist",       fNp[0]);
  FillHistSet("fHDistTPC",    fNp[1]);
  FillHistSet("fHDistphi",    fNp[2]);
  FillHistSet("fHDistTPCphi", fNp[3]);
    
  return 0;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessAODTracks() {
  // -- Process AOD tracks and fill histograms

  for (Int_t idxTrack = 0; idxTrack < fAOD->GetNumberOfTracks(); ++idxTrack) {

    AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(idxTrack)); 

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

    // -- Check if accepted in rapidity window -- for identified particles
    //    for charged particles -- see above
    Double_t yP;
    if (fHelper->GetUsePID() && !fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted bt PID from TPC or TPC+TOF
    Double_t pid[2];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- Check if accepted with thighter DCA cuts -- check MW
    // if (fHelper->IsTrackAcceptedDCA(track))
    //  continue;
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Fill Probe Particle
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    if (!fHelper->GetUsePID())  
      yP = track->Eta();
    
    Double_t aTrack[6] = {
      Double_t(fHelper->GetCentralityBin()), //  0 centrality 
      track->Eta(),                          //  1 eta
      yP,                                    //  2 rapidity
      track->Phi(),                          //  3 phi
      track->Pt(),                           //  4 pt
      track->Charge(),                       //  5 sign
    };
    
    fHnTrackUnCorr->Fill(aTrack);
    
    // -- Count particle / anti-particle 
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle

    Int_t idxPart = (track->Charge() < 0) ? 0 : 1;

    // -- in pt Range
    fNp[0][idxPart] += 1.;

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[1][idxPart] += 1.;
      
    // -- check phi range
    if(!fHelper->IsTrackAcceptedPhi(track))
      continue;

    // -- in pt Range
    fNp[2][idxPart] += 1.;

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[3][idxPart] += 1.;

  } // for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSet("fHDist",       fNp[0]);
  FillHistSet("fHDistTPC",    fNp[1]);
  FillHistSet("fHDistphi",    fNp[2]);
  FillHistSet("fHDistTPCphi", fNp[3]);
    
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
    
    // -- Check rapidity window -- for identfied particles
    Double_t yMC;
    if (fHelper->GetUsePID() && !fHelper->IsParticleAcceptedRapidity(particle, yMC))
      continue;

    // -- Check eta window -- for charged particles
    if (!fHelper->GetUsePID() && TMath::Abs(particle->Eta()) > fEtaMax)
      continue;
    
    // -- Count particle / anti-particle 
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle

    Int_t idxPart = (particle->GetPdgCode() < 0) ? 0 : 1;


    // -- MCrapidity for identfied particles
    //    MCeta for charged particles
    fMCNp[0][idxPart] += 1.;        

    // -- in pt Range
    if (particle->Pt() > fPtRange[0] && particle->Pt() <= fPtRange[1])
      fMCNp[1][idxPart] += 1.;

    // -- in TPC pt Range
    if (particle->Pt() > fPtRange[0] && particle->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[2][idxPart] += 1.;
      
    // -- check phi range
    if(!fHelper->IsParticleAcceptedPhi(particle))
      continue;
    
    // -- in pt Range
    if (particle->Pt() > fPtRange[0] && particle->Pt() <= fPtRange[1])
      fMCNp[3][idxPart] += 1.;

    // -- in TPC pt Range
    if (particle->Pt() > fPtRange[0] && particle->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[4][idxPart] += 1.;

  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSet("fHMC",         fMCNp[0]);
  FillHistSet("fHMCpt",       fMCNp[1]);
  FillHistSet("fHMCptTPC",    fMCNp[2]);
  FillHistSet("fHMCptphi",    fMCNp[3]);
  FillHistSet("fHMCptTPCphi", fMCNp[4]);

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
  list->SetName(name);
  list->SetOwner(kTRUE);

  TString sName(name);

  TString sPtTitle("");
  if (!sName.Contains("fHMC"))
    sPtTitle += Form("#it{p}_{T} [%.1f,%.1f]", fPtRange[0], fPtRange[1]);
  
  // -- Add Particle / Anti-Particle Distributions
  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    list->Add(new TH2F(Form("%s%s", name, fHelper->GetParticleName(idxPart).Data()), 
		       Form("%s %s Dist %s;Centrality;N_{%s}", title, fHelper->GetParticleTitle(idxPart).Data(), sPtTitle.Data(), 
			    fHelper->GetParticleTitleLatex(idxPart).Data()), 
		       24, -0.5, 11.5, 2601, -0.5, 2600.49));
  } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
   
  // -- Add NetParticle Distributions
  list->Add(new TH2F(Form("%sNet%s", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s Net %s Dist %s;Centrality;N_{%s} - N_{%s}", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			  fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data()), 24, -0.5, 11.5, 601, -300.5, 300.49));

  // -- Add NetParticle vs SumParticle
  list->Add(new TH2F(Form("%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s (Net %s)/ Sum Dist %s;Centrality;(N_{%s} - N_{%s})/(N_{%s} + N_{%s})", title, fHelper->GetParticleTitle(1).Data(), sPtTitle.Data(),
			  fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data(),fHelper->GetParticleTitleLatex(1).Data(),
			  fHelper->GetParticleTitleLatex(0).Data()), 
		     24, -0.5, 11.5, 801, -20.5, 20.49));

  // -- Add TProfiles for <NetParticle^k>
  for (Int_t idx = 1; idx <= fOrder; ++idx) {
    list->Add(new TProfile(Form("%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idx), 
			   Form("%s (Net %s)^{%d} Dist %s;Centrality;(N_{%s} - N_{%s})^{%d}", title, fHelper->GetParticleTitle(1).Data(), idx, sPtTitle.Data(),
				fHelper->GetParticleTitleLatex(1).Data(),fHelper->GetParticleTitleLatex(0).Data(), idx), 24, -0.5, 11.5));
  }

  // -- Add TProfiles for <f_ik>
  list->Add(new TList);
  TList *fikList = static_cast<TList*>(list->Last());
  fikList->SetName(Form("%sFik",name));
  fikList->SetOwner(kTRUE);

  for (Int_t ii = 0; ii <= fOrder; ++ii) {
    for (Int_t kk = 0; kk <= fOrder; ++kk) {
      fikList->Add(new TProfile(Form("%sNet%sF%d%d", name, fHelper->GetParticleName(1).Data(), ii, kk), 
				Form("%s (%s) f_{%d%d}Dist %s;Centrality; f_{%d%d} (%s)", title, fHelper->GetParticleTitle(1).Data(), ii, kk, sPtTitle.Data(),
				     ii, kk, fHelper->GetParticleTitleLatex(1).Data()), 24, -0.5, 11.5));
    }
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::FillHistSet(const Char_t *name, Float_t *np)  {
  // -- Add histogram sets for particle and anti-particle

  TList *list = static_cast<TList*>(fOutList->FindObject(name));

  Float_t centralityBin = fHelper->GetCentralityBin();

  (static_cast<TH2F*>(list->FindObject(Form("%s%s", name, fHelper->GetParticleName(0).Data()))))->Fill(centralityBin, np[0]);
  (static_cast<TH2F*>(list->FindObject(Form("%s%s", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, np[1]);

  Float_t sumNp    = np[1]+np[0];
  Float_t deltaNp  = np[1]-np[0];

  (static_cast<TH2F*>(list->FindObject(Form("%sNet%s",  name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp);
  (static_cast<TH2F*>(list->FindObject(Form("%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp/sumNp);

  // -- Fill TProfile for <NetParticle^k>
  Float_t delta = 1.;
  for (Int_t idx = 1; idx <= fOrder; ++idx) {
    delta *= deltaNp;
    (static_cast<TProfile*>(list->FindObject(Form("%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idx))))->Fill(centralityBin, delta);
  }

  // -- Fill TProfiles for <f_ik>
  TList *fikList = static_cast<TList*>(list->FindObject(Form("%sFik",name)));
  for (Int_t ii = 0; ii <= fOrder; ++ii)
    for (Int_t kk = 0; kk <= fOrder; ++kk)
      (static_cast<TProfile*>(fikList->FindObject(Form("%sNet%sF%d%d", name, fHelper->GetParticleName(1).Data(), ii, kk))))->Fill(centralityBin, GetF(ii, kk, np));

  return;
}

// ____________________________________________________________________________________________________
Double_t AliAnalysisNetParticleDistribution::GetF(Int_t ii, Int_t kk, Float_t *np)  {
  // -- Calculate ingredient for faktorial moment
  
  return TMath::Factorial(np[0]) * TMath::Factorial(np[1]) / (TMath::Factorial(np[0] - ii) * TMath::Factorial(np[1] - kk));
}

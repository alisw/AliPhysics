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
#include "AliAODMCParticle.h"

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
  fPdgCode(2212),
  fOutList(NULL),

  fESD(NULL),
  fESDTrackCuts(NULL),
  fPIDResponse(NULL),
  fAOD(NULL),
  fArrayMC(NULL),
  fAODtrackCutBit(1024),
  fIsMC(kFALSE),
  fMCEvent(NULL),
  fStack(NULL),

  fOrder(12),
  fNNp(6),
  fNp(NULL),
  fNMCNp(7),
  fMCNp(NULL),
  fFactp(NULL),

  fCentralityBin(-1.),
  fNTracks(0),

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

  for (Int_t ii = 0; ii <= fOrder; ++ii) 
    if (fFactp[ii]) delete[] fFactp[ii];
  if (fFactp) delete[] fFactp;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::Initialize(AliAnalysisNetParticleHelper* helper, AliESDtrackCuts* cuts, Bool_t isMC, Int_t trackCutBit) {
  // -- Initialize

  // -- Get Helper
  // ---------------
  fHelper           = helper;

  // -- ESD track cuts
  // -------------------
  fESDTrackCuts     = cuts;

  // -- Is MC
  // ----------
  fIsMC             = isMC;

  // -- AOD track filter bit
  // -------------------------
  fAODtrackCutBit   = trackCutBit;

  // -- Get particle species / pdgCode
  // -------------------------
  if (fIsMC)
    fPdgCode        = AliPID::ParticleCode(fHelper->GetParticleSpecies());

  // ------------------------------------------------------------------
  // -- N particles / N anti-particles
  // ------------------------------------------------------------------
  //  Np          : arr[set][particle]
  //  MCNp        : arr[set][particle] - MC
  //  Factorials  : arr[order][particle]
  
  fNp = new Double_t*[fNNp];
  for (Int_t ii = 0 ; ii < fNNp; ++ii)
    fNp[ii] = new Double_t[2];

  fMCNp = new Double_t*[fNMCNp];
  for (Int_t ii = 0 ; ii < fNMCNp; ++ii)
    fMCNp[ii] = new Double_t[2];

  fFactp = new Double_t*[fOrder+1];
  for (Int_t ii = 0 ; ii <= fOrder; ++ii)
    fFactp[ii] = new Double_t[2];

  ResetEvent();

  return 0;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::CreateHistograms(TList* outList) {
  // -- Add histograms to outlist

  fOutList = outList;

  // -- Get ranges for pt and eta 
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);

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
  sTitle = (fHelper->GetUsePID()) ? Form("Uncorrected in |y| < %.1f", fHelper->GetRapidityMax()) : Form("Uncorrected in |#eta| < %.1f", etaRange[1]);

  AddHistSet("fHDist",       Form("%s + #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));
  AddHistSet("fHDistTPC",    Form("%s + #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired()));
  AddHistSet("fHDistTOF",    Form("%s + #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1]));

  AddHistSet("fHDistphi",    Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				  sTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  AddHistSet("fHDistTPCphi", Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				  sTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  AddHistSet("fHDistTOFphi", Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				  sTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  
  if (fIsMC) {
    TString sMCTitle("");
    sMCTitle = (fHelper->GetUsePID()) ?  Form("MC primary in |y| < %.1f", fHelper->GetRapidityMax()) : Form("MC primary in |#eta| < %.1f", etaRange[1]);

    AddHistSet("fHMC",         Form("%s", sTitle.Data()));
    AddHistSet("fHMCpt",       Form("%s + #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], ptRange[1]));
    AddHistSet("fHMCptTPC",    Form("%s + #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired()));
    AddHistSet("fHMCptTOF",    Form("%s + #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1]));

    AddHistSet("fHMCptphi",    Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				    sMCTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
    AddHistSet("fHMCptTPCphi", Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				    sMCTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
    AddHistSet("fHMCptTOFphi", Form("%s + #it{p}_{T} [%.1f,%.1f] + #varphi [%.1f,%.1f]", 
				    sMCTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
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
    fPIDResponse = esdHandler->GetPIDResponse();
    fESD         = esdHandler->GetEvent();
    fNTracks     = fESD->GetNumberOfTracks();
  }

  // -- Get AOD objects
  else if(aodHandler){
    fPIDResponse = aodHandler->GetPIDResponse();
    fAOD         = aodHandler->GetEvent();
    fNTracks     = fAOD->GetNumberOfTracks();
    
    if (fIsMC) {
      fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!"); // MW  no AliFatal use return values
    }
  }

  // -- Get MC objects
  if (fIsMC && mcEvent) {
    fMCEvent     = mcEvent;
    if (fMCEvent)
      fStack     = fMCEvent->Stack();
  }

  // -- Get CentralityBin
  fCentralityBin = fHelper->GetCentralityBin();

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

 // -- Reset factorials for particles/anti-particles
  for (Int_t ii = 0; ii <= fOrder; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fFactp[ii][jj] = 0.;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::Process() {
  // -- Process NetParticle Distributions

  // -- Fill ESD/AOD tracks
  ProcessTracks();
      
  // -- Fill MC truth particles (missing for AOD MW - However AliVParticle already used)
  if (fIsMC)
    ProcessParticles();

  return 0;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessTracks() {
  // -- Process ESD/AOD tracks and fill QA histogram

  // -- Get ranges for AOD particles
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);
  
  // -- Track Loop
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    
    AliVTrack *track = (fESD) ? static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Check track cuts
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

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
    //    for charged - eta check is done in the trackcuts
    Double_t yP;
    if (fHelper->GetUsePID() && !fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted with thighter DCA cuts
    // -- returns kTRUE for AODs for now : MW
    if (!fHelper->IsTrackAcceptedDCA(track))
      continue;

    // -- Check if accepted by PID from TPC or TPC+TOF
    Double_t pid[2];
    if (!fHelper->IsTrackAcceptedPID(track, pid))
      continue;

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Fill Probe Particle
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    if (!fHelper->GetUsePID())  
      yP = track->Eta();
    
    Double_t aTrack[6] = {
      Double_t(fCentralityBin),               //  0 centrality 
      track->Eta(),                           //  1 eta
      yP,                                     //  2 rapidity
      track->Phi(),                           //  3 phi
      track->Pt(),                            //  4 pt
      track->Charge()                         //  5 sign
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

    // -- in TPC+TOF pt Range
    if (track->Pt() > fHelper->GetMinPtForTOFRequired())
      fNp[2][idxPart] += 1.;

    // -- check phi range ----------------------------------------------------------
    if(!fHelper->IsTrackAcceptedPhi(track))
      continue;

    // -- in pt Range
    fNp[3][idxPart] += 1.;

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[4][idxPart] += 1.;

    // -- in TPC+TOF pt Range
    if (track->Pt() > fHelper->GetMinPtForTOFRequired())
      fNp[5][idxPart] += 1.;
  
  } // for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSet("fHDist",       fNp[0]);
  FillHistSet("fHDistTPC",    fNp[1]);
  FillHistSet("fHDistTOF",    fNp[2]);
  FillHistSet("fHDistphi",    fNp[3]);
  FillHistSet("fHDistTPCphi", fNp[4]);
  FillHistSet("fHDistTOFphi", fNp[5]);
    
  return 0;
}

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessParticles() {
  // -- Process primary particles from the stack and fill histograms
  
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);

  for (Int_t idxMC = 0; idxMC < fStack->GetNprimary(); ++idxMC) {
    AliVParticle* particle = (fESD) ? fMCEvent->GetTrack(idxMC) : NULL;

    if (!particle) 
      continue;

    // -- Check basic MC properties -> charged physical primary
    if (!fHelper->IsParticleAcceptedBasicCharged(particle, idxMC))
      continue;
    
    // -- Check if particle / anti-particle
    if (fHelper->GetUsePID() && TMath::Abs(particle->PdgCode()) != fPdgCode)
      continue;
    
    // -- Check rapidity window -- for identfied particles
    Double_t yMC;
    if (fHelper->GetUsePID() && !fHelper->IsParticleAcceptedRapidity(particle, yMC))
      continue;

    // -- Check eta window -- for charged particles
    if (!fHelper->GetUsePID() && TMath::Abs(particle->Eta()) > etaRange[1])
      continue;
    
    // -- Count particle / anti-particle 
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle

    Int_t idxPart = (particle->PdgCode() < 0) ? 0 : 1;

    // -- MCrapidity for identfied particles
    //    MCeta for charged particles
    fMCNp[0][idxPart] += 1.;        

    // -- in pt Range
    if (particle->Pt() > ptRange[0] && particle->Pt() <= ptRange[1])
      fMCNp[1][idxPart] += 1.;

    // -- in TPC pt Range
    if (particle->Pt() > ptRange[0] && particle->Pt() <= fHelper->GetMinPtForTOFRequired())
      fMCNp[2][idxPart] += 1.;

    // -- in TPC+TOF pt Range
    if (particle->Pt() > fHelper->GetMinPtForTOFRequired() && particle->Pt() <= ptRange[1])
      fMCNp[3][idxPart] += 1.;
      
    // -- check phi range ----------------------------------------------------------
    if(!fHelper->IsParticleAcceptedPhi(particle))
      continue;
    
    // -- in pt Range
    if (particle->Pt() > ptRange[0] && particle->Pt() <= ptRange[1])
      fMCNp[4][idxPart] += 1.;

    // -- in TPC pt Range
    if (particle->Pt() > ptRange[0] && particle->Pt() <= fHelper->GetMinPtForTOFRequired())
      fMCNp[5][idxPart] += 1.;

    // -- in TPC+TOF pt Range
    if (particle->Pt() > fHelper->GetMinPtForTOFRequired() && particle->Pt() <= ptRange[1])
      fMCNp[6][idxPart] += 1.;

  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSet("fHMC",         fMCNp[0]);
  FillHistSet("fHMCpt",       fMCNp[1]);
  FillHistSet("fHMCptTPC",    fMCNp[2]);
  FillHistSet("fHMCptTOF",    fMCNp[3]);
  FillHistSet("fHMCptphi",    fMCNp[4]);
  FillHistSet("fHMCptTPCphi", fMCNp[5]);
  FillHistSet("fHMCptTOFphi", fMCNp[6]);

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
  TString sTitle(title);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);

  TString sPtTitle("");
  if (!sName.Contains("fHMC"))
    sPtTitle += Form("#it{p}_{T} [%.1f,%.1f]", ptRange[0], ptRange[1]);
  
  TString sNetTitle(Form("N_{%s} - N_{%s}", fHelper->GetParticleTitleLatex(0).Data(), fHelper->GetParticleTitleLatex(1).Data()));
  TString sSumTitle(Form("N_{%s} + N_{%s}", fHelper->GetParticleTitleLatex(0).Data(), fHelper->GetParticleTitleLatex(1).Data()));

  // -- Add Particle / Anti-Particle Distributions
  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    list->Add(new TH2F(Form("%s%s", name, fHelper->GetParticleName(idxPart).Data()), 
		       Form("N_{%s} : %s;Centrality;N_{%s}", fHelper->GetParticleTitleLatex(idxPart).Data(), sTitle.Data(), 
			    fHelper->GetParticleTitleLatex(idxPart).Data()), 24, -0.5, 11.5, 2601, -0.5, 2600.49));
  } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
   
  // -- Add NetParticle Distributions
  list->Add(new TH2F(Form("%sNet%s", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s : %s;Centrality;%s", sNetTitle.Data(), sTitle.Data(), sNetTitle.Data()), 
		     24, -0.5, 11.5, 601, -300.5, 300.49));

  // -- Add NetParticle vs SumParticle
  list->Add(new TH2F(Form("%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()), 
		     Form("(%s)/(%s) : %s;Centrality;(%s)/(%s)", sNetTitle.Data(), sSumTitle.Data(), sTitle.Data(), 
			  sNetTitle.Data(), sSumTitle.Data()), 24, -0.5, 11.5, 801, -20.5, 20.49));

  // -- Add TProfiles for <NetParticle^k>
  for (Int_t idx = 1; idx <= fOrder; ++idx) {
    list->Add(new TProfile(Form("%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idx), 
			   Form("(%s)^{%d} : %s;Centrality;(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx), 
			   24, -0.5, 11.5));
  }
  
  // -- Add TProfiles for <f_ik>
  list->Add(new TList);
  TList *fikList = static_cast<TList*>(list->Last());
  fikList->SetName(Form("%sFik",name));
  fikList->SetOwner(kTRUE);

  for (Int_t ii = 0; ii <= fOrder; ++ii) {
    for (Int_t kk = 0; kk <= fOrder; ++kk) {
      fikList->Add(new TProfile(Form("%sNet%sF%02d%02d", name, fHelper->GetParticleName(1).Data(), ii, kk),
				Form("f_{%d%d} : %s;Centrality;f_{%d%d}", ii, kk, sTitle.Data(), ii, kk), 24, -0.5, 11.5));
    }
  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::FillHistSet(const Char_t *name, Double_t *np)  {
  // -- Add histogram sets for particle and anti-particle

  TList *list = static_cast<TList*>(fOutList->FindObject(name));

  Float_t centralityBin = fHelper->GetCentralityBin();

  (static_cast<TH2F*>(list->FindObject(Form("%s%s", name, fHelper->GetParticleName(0).Data()))))->Fill(centralityBin, np[0]);
  (static_cast<TH2F*>(list->FindObject(Form("%s%s", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, np[1]);

  Double_t sumNp    = np[1]+np[0];
  Double_t deltaNp  = np[1]-np[0];

  (static_cast<TH2F*>(list->FindObject(Form("%sNet%s",  name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp);
  (static_cast<TH2F*>(list->FindObject(Form("%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp/sumNp);

  // -- Fill TProfile for <NetParticle^k>
  Double_t delta = 1.;
  for (Int_t idx = 1; idx <= fOrder; ++idx) {
    delta *= deltaNp;
    (static_cast<TProfile*>(list->FindObject(Form("%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idx))))->Fill(centralityBin, delta);
  }

  // -- Calculate all factorials only once before filling
  for (Int_t idx = 0; idx <= fOrder; ++ idx) {
    fFactp[idx][0] = TMath::Factorial(np[0] - idx);
    fFactp[idx][1] = TMath::Factorial(np[1] - idx);
  }

  // -- Fill TProfiles for <f_ik> 
  TList *fikList = static_cast<TList*>(list->FindObject(Form("%sFik",name)));
  for (Int_t ii = 0; ii <= fOrder; ++ii) {
    for (Int_t kk = 0; kk <= fOrder; ++kk) {
      Double_t fik = fFactp[0][1]*fFactp[0][0]/(fFactp[ii][1]*fFactp[kk][0]);
      (static_cast<TProfile*>(fikList->FindObject(Form("%sNet%sF%02d%02d", name, fHelper->GetParticleName(1).Data(), ii, kk))))->Fill(centralityBin, fik);
    }
  }

  return;
}


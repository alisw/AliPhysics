//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"
#include "TProfile.h" 
#include "TProfile2D.h" 
#include "TH2D.h" 
#include "TH3D.h" 

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"

#include "AliCentrality.h"
#include "AliTracker.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisNetParticleDistribution.h"

using namespace std;

/**
 * Class for NetParticle Distributions
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
  AliAnalysisNetParticleBase("Dist", "Dist"),

  fOutList(NULL),

  fOrder(8),
  fNNp(6),
  fNp(NULL),
  fNpPt(NULL),
  fNMCNp(7),
  fMCNp(NULL),
  fMCNpPt(NULL),
  fRedFactp(NULL),
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

  for (Int_t ii = 0; ii < 2; ++ii) {
    for (Int_t kk = 0; kk < 2; ++kk)
      if (fNpPt[ii][kk]) delete[] fNpPt[ii][kk];      
    if (fNpPt[ii]) delete[] fNpPt[ii];
  }
  if (fNpPt) delete[] fNpPt;

  // -------------------------------------------
  for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    if (fMCNp[ii]) delete[] fMCNp[ii];
  if (fMCNp) delete[] fMCNp;

  for (Int_t ii = 0; ii < 2; ++ii) {
    for (Int_t kk = 0; kk < 2; ++kk)
      if (fMCNpPt[ii][kk]) delete[] fMCNpPt[ii][kk];      
    if (fMCNpPt[ii]) delete[] fMCNpPt[ii];
  }
  if (fMCNpPt) delete[] fMCNpPt;

  // -------------------------------------------
  for (Int_t ii = 0; ii <= fOrder; ++ii) 
    if (fRedFactp[ii]) delete[] fRedFactp[ii];
  if (fRedFactp) delete[] fRedFactp;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::Process() {
  // -- Process NetParticle Distributions

  // -- Fill ESD/AOD tracks
  ProcessTracks();
      
  // -- Fill MC truth particles (missing for AOD MW - However AliVParticle already used)
  if (fIsMC)
    ProcessParticles();

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::Init() {
  // -- Init eventwise

  // ------------------------------------------------------------------
  // -- N particles / N anti-particles
  // ------------------------------------------------------------------
  //  Np            : arr[set][particle]
  //  MCNp          : arr[set][particle] - MC
  //  Factorials    : arr[order][particle]

  //  NpPt          : arr[phiBin][particle][ptBins]
  //  MCNpPt        : arr[phiBin][particle][ptBins] - MC
  //  FactorialsPt  : arr[order][particle][ptBins]
  
  //   set      ->  fNNp  / fNMCNp
  //   order    ->  fOrder
  //   particle ->  2 {0,1}
  //   phiBin   ->  2 {0,1}

  fNp = new Int_t*[fNNp];
  for (Int_t ii = 0 ; ii < fNNp; ++ii)
    fNp[ii] = new Int_t[2];

  fNpPt = new Int_t**[2];
  for (Int_t ii = 0 ; ii < 2; ++ii) {
    fNpPt[ii] = new Int_t*[2];
    for (Int_t kk = 0 ; kk < 2; ++kk)
      fNpPt[ii][kk] = new Int_t[AliAnalysisNetParticleHelper::fgkfHistNBinsPt];
  }

  // -------------------------------------------
  fMCNp = new Int_t*[fNMCNp];
  for (Int_t ii = 0 ; ii < fNMCNp; ++ii)
    fMCNp[ii] = new Int_t[2];

  fMCNpPt = new Int_t**[2];
  for (Int_t ii = 0 ; ii < 2; ++ii) {
    fMCNpPt[ii] = new Int_t*[2];
    for (Int_t kk = 0 ; kk < 2; ++kk)
      fMCNpPt[ii][kk] = new Int_t[AliAnalysisNetParticleHelper::fgkfHistNBinsPt];
  }

  // -------------------------------------------
  fRedFactp = new Double_t*[fOrder+1];
  for (Int_t ii = 0 ; ii <= fOrder; ++ii)
    fRedFactp[ii] = new Double_t[2];
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::Reset() {
  // -- Reset eventwise

  // -- Reset N particles/anti-particles
  for (Int_t ii = 0; ii < fNNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fNp[ii][jj] = 0;

  for (Int_t ii = 0; ii < 2; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      for (Int_t kk = 0; kk < AliAnalysisNetParticleHelper::fgkfHistNBinsPt; ++kk)
	fNpPt[ii][jj][kk] = 0;

  // -- Reset N MC particles/anti-particles
  for (Int_t ii = 0; ii < fNMCNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fMCNp[ii][jj] = 0;

  for (Int_t ii = 0; ii < 2; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      for (Int_t kk = 0; kk < AliAnalysisNetParticleHelper::fgkfHistNBinsPt; ++kk)
	fMCNpPt[ii][jj][kk] = 0;

  // -- Reset reduced factorials for particles/anti-particles
  for (Int_t ii = 0; ii <= fOrder; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      fRedFactp[ii][jj] = 1.;

}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::CreateHistograms() {
  // -- Add histograms to outlist

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
  fOutList->Add(new THnSparseD("hnTrackUnCorr", "cent:eta:y:phi:pt:sign", 6, binHnUnCorr, minHnUnCorr, maxHnUnCorr));  
  fHnTrackUnCorr = static_cast<THnSparseD*>(fOutList->Last());
  fHnTrackUnCorr->Sumw2(); 
  fHnTrackUnCorr->GetAxis(0)->SetTitle("centrality");
  fHnTrackUnCorr->GetAxis(1)->SetTitle("#eta");
  fHnTrackUnCorr->GetAxis(2)->SetTitle("#it{y}");
  fHnTrackUnCorr->GetAxis(3)->SetTitle("#varphi");
  fHnTrackUnCorr->GetAxis(4)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHnTrackUnCorr->GetAxis(5)->SetTitle("sign");

  fHelper->BinLogAxis(fHnTrackUnCorr, 4, fESDTrackCuts);

  for (Int_t idx = 1 ; idx <= fHnTrackUnCorr->GetAxis(4)->GetNbins(); ++idx)
    printf("%02d |  %f > %f < %f\n", idx, fHnTrackUnCorr->GetAxis(4)->GetBinLowEdge(idx), fHnTrackUnCorr->GetAxis(4)->GetBinCenter(idx), fHnTrackUnCorr->GetAxis(4)->GetBinUpEdge(idx));

  // ------------------------------------------------------------------
  // -- Create net particle histograms
  // ------------------------------------------------------------------

  // -- Get ranges for pt and eta 
  Float_t etaRange[2];
  fESDTrackCuts->GetEtaRange(etaRange[0],etaRange[1]);

  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);

  // ------------------------------------------------------------------

  TString sTitle("");
  sTitle = (fHelper->GetUsePID()) ? Form("|y|<%.1f", fHelper->GetRapidityMax()) : Form("|#eta|<%.1f", etaRange[1]);

  // -- centrality dependent
  AddHistSetCent("Dist",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));
  AddHistSetCent("DistTPC",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired()));
  AddHistSetCent("DistTOF",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1]));

#if USE_PHI
  AddHistSetCent("Distphi",    Form("%s,#it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				     sTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  AddHistSetCent("DistTPCphi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				     sTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
  AddHistSetCent("DistTOFphi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				     sTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
#endif

  // -- centrality + pt dependent
  AddHistSetCentPt("PtDist",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));  
#if USE_PHI
  AddHistSetCentPt("PtDistphi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				      sTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
#endif
  
  if (fIsMC) {
    TString sMCTitle("");
    sMCTitle = (fHelper->GetUsePID()) ?  Form("MC primary in |y| < %.1f", fHelper->GetRapidityMax()) : Form("MC primary in |#eta| < %.1f", etaRange[1]);

    // -- centrality dependent
    AddHistSetCent("MC",         Form("%s", sTitle.Data()));

    AddHistSetCent("MCpt",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], ptRange[1]));
    AddHistSetCent("MCptTPC",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired()));
    AddHistSetCent("MCptTOF",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1]));

#if USE_PHI
    AddHistSetCent("MCptphi",    Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				       sMCTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
    AddHistSetCent("MCptTPCphi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				       sMCTitle.Data(), ptRange[0], fHelper->GetMinPtForTOFRequired(), fHelper->GetPhiMin(), fHelper->GetPhiMax()));
    AddHistSetCent("MCptTOFphi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
				       sMCTitle.Data(), fHelper->GetMinPtForTOFRequired(), ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
#endif

    // -- centrality + pt dependent
    AddHistSetCentPt("PtMCpt",    Form("%s, #it{p}_{T} [%.1f,%.1f]", sMCTitle.Data(), ptRange[0], ptRange[1]));
#if USE_PHI
    AddHistSetCentPt("PtMCptphi", Form("%s, #it{p}_{T} [%.1f,%.1f], #varphi [%.2f,%.2f]", 
					sMCTitle.Data(), ptRange[0], ptRange[1], fHelper->GetPhiMin(), fHelper->GetPhiMax()));
#endif
  }

  // ------------------------------------------------------------------

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleDistribution::ProcessTracks() {
  // -- Process ESD/AOD tracks and fill histograms

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
    Double_t pid[3];
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
      static_cast<Double_t>(track->Charge())                         //  5 sign
    };
    
    fHnTrackUnCorr->Fill(aTrack);
    
    // -- Count particle / anti-particle 
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle
    Int_t idxPart = (track->Charge() < 0) ? 0 : 1;

    //  idxPhi  = 0 -> Full Acceptance
    //  idxPhi  = 1 -> Partial Acceptance
    Int_t idxPhi = 0;

    //  idxPt       -> [0, nBins-1]
    Int_t idxPt = static_cast<TAxis*>(fHnTrackUnCorr->GetAxis(4))->FindBin(track->Pt())-1;

    // -- using pt Bins
    fNpPt[idxPhi][idxPart][idxPt] += 1;

    // -- in pt Range
    fNp[0][idxPart] += 1;

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired())
      fNp[1][idxPart] += 1;

    // -- in TPC+TOF pt Range
    if (track->Pt() > fHelper->GetMinPtForTOFRequired())
      fNp[2][idxPart] += 1;

    // -- check phi range ----------------------------------------------------------
#if USE_PHI
    if(!fHelper->IsTrackAcceptedPhi(track))
      continue;
    idxPhi = 1;

    // -- using pt Bins
    fNpPt[idxPhi][idxPart][idxPt] += 1;

    // -- in pt Range
    fNp[3][idxPart] += 1;

    // -- in TPC pt Range
    if (track->Pt() <= fHelper->GetMinPtForTOFRequired()) 
      fNp[4][idxPart] += 1;
    
    // -- in TPC+TOF pt Range
    if (track->Pt() > fHelper->GetMinPtForTOFRequired())
      fNp[5][idxPart] += 1;
#endif
  } // for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSetCent("Dist",        0, kFALSE);
  FillHistSetCent("DistTPC",     1, kFALSE);
  FillHistSetCent("DistTOF",     2, kFALSE);

  FillHistSetCentPt("PtDist",    0, kFALSE);

#if USE_PHI
  FillHistSetCent("Distphi",     3, kFALSE);
  FillHistSetCent("DistTPCphi",  4, kFALSE);
  FillHistSetCent("DistTOFphi",  5, kFALSE);

  FillHistSetCentPt("PtDistphi", 1, kFALSE);
#endif
 
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

    //  idxPhi  = 0 -> Full Acceptance
    //  idxPhi  = 1 -> Partial Acceptance
    Int_t idxPhi = 0;

    //  idxPt       -> [0, nBins-1]
    Int_t idxPt = static_cast<TAxis*>(fHnTrackUnCorr->GetAxis(4))->FindBin(particle->Pt()) - 1;

    // -- MCrapidity for identfied particles
    //    MCeta for charged particles
    fMCNp[0][idxPart] += 1.;        
     
    // -- Check main pt window
    if (!(particle->Pt() > ptRange[0] && particle->Pt() <= ptRange[1]))
      continue;

    // -- using pt Bin
    fMCNpPt[idxPhi][idxPart][idxPt] += 1;

    // -- in pt Range
    fMCNp[1][idxPart] += 1;
    
    // -- in TPC pt Range
    if (particle->Pt() <= fHelper->GetMinPtForTOFRequired())
      fMCNp[2][idxPart] += 1;

    // -- in TPC+TOF pt Range
    if (particle->Pt() > fHelper->GetMinPtForTOFRequired())
      fMCNp[3][idxPart] += 1;

    // -- check phi range ----------------------------------------------------------
#if USE_PHI
    if(!fHelper->IsParticleAcceptedPhi(particle))
      continue;
    idxPhi = 1;

    // -- using pt Bin
    fMCNpPt[idxPhi][idxPart][idxPt] += 1;
    
    // -- in pt Range
    fMCNp[4][idxPart] += 1;

    // -- in TPC pt Range
    if (particle->Pt() <= fHelper->GetMinPtForTOFRequired())
      fMCNp[5][idxPart] += 1;

    // -- in TPC+TOF pt Range
    if (particle->Pt() > fHelper->GetMinPtForTOFRequired())
      fMCNp[6][idxPart] += 1;
#endif
  } // for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill Particle Fluctuation Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  FillHistSetCent("MC",          0, kTRUE);
  FillHistSetCent("MCpt",        1, kTRUE);
  FillHistSetCent("MCptTPC",     2, kTRUE);
  FillHistSetCent("MCptTOF",     3, kTRUE);

  FillHistSetCentPt("PtMCpt",    0, kTRUE);

#if USE_PHI
  FillHistSetCent("MCptphi",     4, kTRUE);
  FillHistSetCent("MCptTPCphi",  5, kTRUE);
  FillHistSetCent("MCptTOFphi",  6, kTRUE);

  FillHistSetCentPt("PtMCptphi", 1, kTRUE);
#endif

  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Helper Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void  AliAnalysisNetParticleDistribution::AddHistSetCent(const Char_t *name, const Char_t *title)  {
  // -- Add histogram sets for particle and anti-particle
  //    dependence : centrality

  TString sName(name);
  TString sTitle(title);

  // -- Add List
  fOutList->Add(new TList);
  TList *list = static_cast<TList*>(fOutList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);

  // -- Get Bin Ranges
  Int_t nBinsCent         =  AliAnalysisNetParticleHelper::fgkfHistNBinsCent;
  Double_t centBinRange[] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeCent[1]};

  // -- Create Titles
  TString sNetTitle(Form("N_{%s} - N_{%s}", fHelper->GetParticleTitleLatex(1).Data(), fHelper->GetParticleTitleLatex(0).Data()));
  TString sSumTitle(Form("N_{%s} + N_{%s}", fHelper->GetParticleTitleLatex(1).Data(), fHelper->GetParticleTitleLatex(0).Data()));

  // -- Add Particle / Anti-Particle Distributions
  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    list->Add(new TH2D(Form("h%s%s", name, fHelper->GetParticleName(idxPart).Data()), 
		       Form("N_{%s} : %s;Centrality;N_{%s}", fHelper->GetParticleTitleLatex(idxPart).Data(), sTitle.Data(), fHelper->GetParticleTitleLatex(idxPart).Data()),
		       nBinsCent, centBinRange[0], centBinRange[1], 2601, -0.5, 2600.49));
  } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
   
  // -- Add NetParticle Distributions
  list->Add(new TH2D(Form("h%sNet%s", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s : %s;Centrality;%s", sNetTitle.Data(), sTitle.Data(), sNetTitle.Data()), 
		       nBinsCent, centBinRange[0], centBinRange[1], 601, -300.5, 300.49));

  // -- Add NetParticle vs SumParticle
  list->Add(new TH2D(Form("h%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()), 
		     Form("(%s)/(%s) : %s;Centrality;(%s)/(%s)", sNetTitle.Data(), sSumTitle.Data(), sTitle.Data(), sNetTitle.Data(), sSumTitle.Data()), 
		       nBinsCent, centBinRange[0], centBinRange[1], 41, -2.5, 2.49));

  // -----------------------------------------------------------------------------------------------

  // -- Add Particle / Anti-Particle Distributions
  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    list->Add(new TH2D(Form("h%s%sX", name, fHelper->GetParticleName(idxPart).Data()), 
		       Form("N_{%s} : %s;Centrality;N_{%s}", fHelper->GetParticleTitleLatex(idxPart).Data(), sTitle.Data(), fHelper->GetParticleTitleLatex(idxPart).Data()),
		       nBinsCent, centBinRange[0], centBinRange[1], 2601, -0.5, 2600.49));
  } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
   
  // -- Add NetParticle Distributions
  list->Add(new TH2D(Form("h%sNet%sX", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s : %s;Centrality;%s", sNetTitle.Data(), sTitle.Data(), sNetTitle.Data()), 
		       nBinsCent, centBinRange[0], centBinRange[1], 601, -300.5, 300.49));

  // -- Add NetParticle vs SumParticle
  list->Add(new TH2D(Form("h%sNet%sOverSumX", name, fHelper->GetParticleName(1).Data()), 
		     Form("(%s)/(%s) : %s;Centrality;(%s)/(%s)", sNetTitle.Data(), sSumTitle.Data(), sTitle.Data(), sNetTitle.Data(), sSumTitle.Data()), 
		       nBinsCent, centBinRange[0], centBinRange[1], 41, -2.5, 2.49));

  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <NetParticle^k>
  // -----------------------------------------------------------------------------------------------
  for (Int_t idx = 1; idx <= fOrder; ++idx) {
    list->Add(new TProfile(Form("p%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idx), 
			   Form("(%s)^{%d} : %s;Centrality;(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			   nBinsCent, centBinRange[0], centBinRange[1]));
  }

  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <NetParticle^k> for every SubSample
  // -----------------------------------------------------------------------------------------------
  for (Int_t idxSub = 0; idxSub < fHelper->GetNSubSamples(); ++ idxSub) {
    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile(Form("p%sNet%s%dM_%02d", name, fHelper->GetParticleName(1).Data(), idx, idxSub), 
			     Form("(%s)^{%d} : %s;Centrality;(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			     nBinsCent, centBinRange[0], centBinRange[1]));
    }
  }
  
  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <f_ik>
  // -----------------------------------------------------------------------------------------------
  list->Add(new TList);
  TList *fikList = static_cast<TList*>(list->Last());
  fikList->SetName(Form("f%sFik",name));
  fikList->SetOwner(kTRUE);

  for (Int_t ii = 0; ii <= fOrder; ++ii) {
    for (Int_t kk = 0; kk <= fOrder; ++kk) {
      fikList->Add(new TProfile(Form("p%sNet%sF%02d%02d", name, fHelper->GetParticleName(1).Data(), ii, kk),
				Form("f_{%02d%02d} : %s;Centrality;f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
				nBinsCent, centBinRange[0], centBinRange[1]));
    }
  }

  // -- Add counter for number of Non-zero entries
  for (Int_t idxCent = 0; idxCent < nBinsCent; ++idxCent) 
    fikList->Add(new TH2D(Form("p%sNet%sFCounts_%02d", name, fHelper->GetParticleName(1).Data(), idxCent),
			  Form("f_ik counts : %s Cent %d;i;k;", sTitle.Data(), idxCent),
			  fOrder+1, -0.5, Double_t(fOrder)+0.49, fOrder+1, -0.5, Double_t(fOrder)+0.49));
  
  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <f_ik> for every SubSample
  // -----------------------------------------------------------------------------------------------
  for (Int_t idxSub = 0; idxSub < fHelper->GetNSubSamples(); ++ idxSub) {
    list->Add(new TList);
    TList *fikListSub = static_cast<TList*>(list->Last());
    fikListSub->SetName(Form("f%sFik_%02d",name, idxSub));
    fikListSub->SetOwner(kTRUE);

    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	fikListSub->Add(new TProfile(Form("p%sNet%sF%02d%02d_%02d", name, fHelper->GetParticleName(1).Data(), ii, kk, idxSub),
				     Form("f_{%02d%02d} : %s;Centrality;f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
				     nBinsCent, centBinRange[0], centBinRange[1]));
      }
    }

    // -- Add counter for number of Non-zero entries
    for (Int_t idxCent = 0; idxCent < nBinsCent; ++idxCent) 
      fikListSub->Add(new TH2D(Form("p%sNet%sFCounts_%02d_%02d", name, fHelper->GetParticleName(1).Data(), idxCent, idxSub),
			       Form("f_ik counts : %s Cent %d;i;k;", sTitle.Data(), idxCent),
			       fOrder+1, -0.5, Double_t(fOrder)+0.49, fOrder+1, -0.5, Double_t(fOrder)+0.49));
  }
  
  return;
}

//________________________________________________________________________
void  AliAnalysisNetParticleDistribution::AddHistSetCentPt(const Char_t *name, const Char_t *title)  {
  // -- Add histogram sets for particle and anti-particle
  //    dependence : centrality and pt 
  
  TString sName(name);
  TString sTitle(title);

  // -- Add List
  fOutList->Add(new TList);
  TList *list = static_cast<TList*>(fOutList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);

  // -- Get Bin Ranges
  Int_t nBinsPt           =  AliAnalysisNetParticleHelper::fgkfHistNBinsPt;
  Double_t ptBinRange[]   = {-0.5, AliAnalysisNetParticleHelper::fgkfHistNBinsPt - 0.5};

  Int_t nBinsCent         =  AliAnalysisNetParticleHelper::fgkfHistNBinsCent;
  Double_t centBinRange[] = {AliAnalysisNetParticleHelper::fgkfHistRangeCent[0], AliAnalysisNetParticleHelper::fgkfHistRangeCent[1]};

  // -- Create Titles
  TString sNetTitle(Form("N_{%s} - N_{%s}", fHelper->GetParticleTitleLatex(1).Data(), fHelper->GetParticleTitleLatex(0).Data()));
  TString sSumTitle(Form("N_{%s} + N_{%s}", fHelper->GetParticleTitleLatex(1).Data(), fHelper->GetParticleTitleLatex(0).Data()));

  // -- Add Particle / Anti-Particle Distributions
  for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
    list->Add(new TH3D(Form("h%s%s", name, fHelper->GetParticleName(idxPart).Data()), 
		       Form("N_{%s} : %s;Centrality;N_{%s}", fHelper->GetParticleTitleLatex(idxPart).Data(), sTitle.Data(), fHelper->GetParticleTitleLatex(idxPart).Data()),
		       nBinsCent, centBinRange[0], centBinRange[1], nBinsPt, ptBinRange[0], ptBinRange[1], 2601, -0.5, 2600.49));
  } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
   
  // -- Add NetParticle Distributions
  list->Add(new TH3D(Form("h%sNet%s", name, fHelper->GetParticleName(1).Data()), 
		     Form("%s : %s;Centrality;%s", sNetTitle.Data(), sTitle.Data(), sNetTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nBinsPt, ptBinRange[0], ptBinRange[1], 601, -300.5, 300.49));

  // -- Add NetParticle vs SumParticle
  list->Add(new TH3D(Form("h%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()), 
		     Form("(%s)/(%s) : %s;Centrality;(%s)/(%s)", sNetTitle.Data(), sSumTitle.Data(), sTitle.Data(), sNetTitle.Data(), sSumTitle.Data()), 
		     nBinsCent, centBinRange[0], centBinRange[1], nBinsPt, ptBinRange[0], ptBinRange[1], 41, -2.5, 2.49));

  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <NetParticle^k>
  // -----------------------------------------------------------------------------------------------
  for (Int_t idx = 1; idx <= fOrder; ++idx) {
    list->Add(new TProfile2D(Form("p%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idx), 
			     Form("(%s)^{%d} : %s;Centrality;(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			     nBinsCent, centBinRange[0], centBinRange[1], nBinsPt, ptBinRange[0], ptBinRange[1]));
  }

  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <NetParticle^k>  for every SubSample
  // -----------------------------------------------------------------------------------------------
  for (Int_t idxSub = 0; idxSub < fHelper->GetNSubSamples(); ++ idxSub) {
    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile2D(Form("p%sNet%s%dM_%02d", name, fHelper->GetParticleName(1).Data(), idx, idxSub), 
			       Form("(%s)^{%d} : %s;Centrality;(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			       nBinsCent, centBinRange[0], centBinRange[1], nBinsPt, ptBinRange[0], ptBinRange[1]));
    }
  }

  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <f_ik>
  // -----------------------------------------------------------------------------------------------
  list->Add(new TList);
  TList *fikListPt = static_cast<TList*>(list->Last());
  fikListPt->SetName(Form("f%sPtFik",name));
  fikListPt->SetOwner(kTRUE);

  for (Int_t ii = 0; ii <= fOrder; ++ii) {
    for (Int_t kk = 0; kk <= fOrder; ++kk) {
      fikListPt->Add(new TProfile2D(Form("p%sNet%sF%02d%02d", name, fHelper->GetParticleName(1).Data(), ii, kk),
				    Form("f_{%02d%02d} : %s;Centrality;#it{p}_{T} (GeV/#it{c});f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk), 
				    nBinsCent, centBinRange[0], centBinRange[1], nBinsPt, ptBinRange[0], ptBinRange[1]));
    }
  }

  // -- Add counter for number of Non-zero entries
  for (Int_t idxCent = 0; idxCent < nBinsCent; ++idxCent) 
    fikListPt->Add(new TH3D(Form("p%sNet%sFCounts_%02d", name, fHelper->GetParticleName(1).Data(), idxCent),
			    Form("f_ik counts : %s Cent %d;i;k;#it{p}_{T} Bins", sTitle.Data(), idxCent),
			    fOrder+1, -0.5, Double_t(fOrder)+0.49, fOrder+1, -0.5, Double_t(fOrder)+0.49, nBinsPt+1, -0.5, Double_t(nBinsPt)+0.49));
 
  // -----------------------------------------------------------------------------------------------
  // -- Add TProfiles for <f_ik> for every SubSample
  // -----------------------------------------------------------------------------------------------
  for (Int_t idxSub = 0; idxSub < fHelper->GetNSubSamples(); ++ idxSub) {
    list->Add(new TList);
    TList *fikListPtSub = static_cast<TList*>(list->Last());
    fikListPtSub->SetName(Form("f%sPtFik_%02d",name, idxSub));
    fikListPtSub->SetOwner(kTRUE);
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	fikListPtSub->Add(new TProfile2D(Form("p%sNet%sF%02d%02d_%02d", name, fHelper->GetParticleName(1).Data(), ii, kk, idxSub),
					 Form("f_{%02d%02d} : %s;Centrality;#it{p}_{T} (GeV/#it{c});f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk), 
					 nBinsCent, centBinRange[0], centBinRange[1], nBinsPt, ptBinRange[0], ptBinRange[1]));
      }
    }

    // -- Add counter for number of Non-zero entries
    for (Int_t idxCent = 0; idxCent < nBinsCent; ++idxCent) 
      fikListPtSub->Add(new TH3D(Form("p%sNet%sFCounts_%02d_%02d", name, fHelper->GetParticleName(1).Data(), idxCent, idxSub),
				 Form("f_ik counts : %s Cent %d;i;k;#it{p}_{T} Bins", sTitle.Data(), idxCent),
				 fOrder+1, -0.5, Double_t(fOrder)+0.49, fOrder+1, -0.5, Double_t(fOrder)+0.49, nBinsPt+1, -0.5, Double_t(nBinsPt)+0.49));
  }

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::FillHistSetCent(const Char_t *name, Int_t idx, Bool_t isMC)  {
  // -- Fill histogram sets for particle and anti-particle
  //    dependence : centrality 
  
  // -- Get List
  TList *list = static_cast<TList*>(fOutList->FindObject(Form("f%s",name)));
  
  // -- Get Centrality Bin
  Float_t centralityBin = fHelper->GetCentralityBin();

  // -- Select MC or Data
  Int_t **np = (isMC) ? fMCNp : fNp;

  // -----------------------------------------------------------------------------------------------

  Int_t sumNp   = np[idx][1]+np[idx][0];  // p + pbar
  Int_t deltaNp = np[idx][1]-np[idx][0];  // p - pbar

  // -- Fill Particle / Anti-Particle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%s%s", name, fHelper->GetParticleName(0).Data()))))->Fill(centralityBin, np[idx][0]);
  (static_cast<TH2D*>(list->FindObject(Form("h%s%s", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, np[idx][1]);

  // -- Fill NetParticle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%sNet%s",  name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNp);

  // -- Fill NetParticle vs SumParticle
  Double_t deltaNpOverSumNp = (sumNp == 0.) ? 0. : deltaNp/Double_t(sumNp);
  (static_cast<TH2D*>(list->FindObject(Form("h%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNpOverSumNp);

  // -----------------------------------------------------------------------------------------------

  Double_t CENT[7] = {1.157990, 1.153218, 1.147382, 1.142218, 1.139568, 1.138152, 1.134517};

  Double_t sumNpX   = np[idx][1]+(np[idx][0]*CENT[Int_t(centralityBin)]);
  Double_t deltaNpX = np[idx][1]-(np[idx][0]*CENT[Int_t(centralityBin)]);

  // -- Fill Particle / Anti-Particle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%s%sX", name, fHelper->GetParticleName(0).Data()))))->Fill(centralityBin, np[idx][0]*CENT[Int_t(centralityBin)]);
  (static_cast<TH2D*>(list->FindObject(Form("h%s%sX", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, np[idx][1]);

  // -- Fill NetParticle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%sNet%sX",  name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNpX);

  // -- Fill NetParticle vs SumParticle
  Double_t deltaNpXOverSumNpX = (sumNpX == 0.) ? 0. : deltaNpX/sumNpX;
  (static_cast<TH2D*>(list->FindObject(Form("h%sNet%sOverSumX", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, deltaNpXOverSumNpX);

  // -----------------------------------------------------------------------------------------------

  // -- Fill TProfile for <NetParticle^k>
  Double_t delta = 1.;
  for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
    delta *= deltaNp;
    (static_cast<TProfile*>(list->FindObject(Form("p%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idxOrder))))->Fill(centralityBin, delta);
    (static_cast<TProfile*>(list->FindObject(Form("p%sNet%s%dM_%02d", name, fHelper->GetParticleName(1).Data(), idxOrder, fHelper->GetSubSampleIdx()))))->Fill(centralityBin, delta);
  }

  // -- Generate reduced factorials - explictly removing the factorials
  //    - Reset all to 1
  for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
    fRedFactp[idxOrder][0]  = 1.;
    fRedFactp[idxOrder][1]  = 1.;
  }
  
  //    - start at idx 1, as idx 0 = 1 done by reset
  for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
    fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0]  * Double_t(np[idx][0]-(idxOrder-1));
    fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1]  * Double_t(np[idx][1]-(idxOrder-1));
  }

  // -- Fill TProfiles for <f_ik> 
  TList *fikList    = static_cast<TList*>(list->FindObject(Form("f%sFik",name)));
  TList *fikListSub = static_cast<TList*>(list->FindObject(Form("f%sFik_%02d",name, fHelper->GetSubSampleIdx())));
  TH2D  *hCntik     = static_cast<TH2D*>(fikList->FindObject(Form("p%sNet%sFCounts_%02d", name, fHelper->GetParticleName(1).Data(), Int_t(centralityBin))));
  TH2D  *hCntikSub  = static_cast<TH2D*>(fikListSub->FindObject(Form("p%sNet%sFCounts_%02d_%02d", name, fHelper->GetParticleName(1).Data(), Int_t(centralityBin), fHelper->GetSubSampleIdx())));

  for (Int_t ii = 0; ii <= fOrder; ++ii) {   // ii -> p    -> n1
    for (Int_t kk = 0; kk <= fOrder; ++kk) { // kk -> pbar -> n2
      // -- use the reduced factorials only 
      Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   // n1 *n2 -> p * pbar
      (static_cast<TProfile*>(fikList->FindObject(Form("p%sNet%sF%02d%02d", name, fHelper->GetParticleName(1).Data(), ii, kk))))->Fill(centralityBin, fik);
      (static_cast<TProfile*>(fikListSub->FindObject(Form("p%sNet%sF%02d%02d_%02d", 
							  name, fHelper->GetParticleName(1).Data(), ii, kk, fHelper->GetSubSampleIdx()))))->Fill(centralityBin, fik);

      if (fik != 0.) {
	hCntik->Fill(ii, kk);
	hCntikSub->Fill(ii, kk);
      }
    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisNetParticleDistribution::FillHistSetCentPt(const Char_t *name, Int_t idx, Bool_t isMC)  {
  // -- Add histogram sets for particle and anti-particle
  //    dependence : centrality and pt

  // -- Get List
  TList *list = static_cast<TList*>(fOutList->FindObject(Form("f%s",name)));  

  // -- Get Centrality Bin
  Float_t centralityBin = fHelper->GetCentralityBin();

  // -- Select MC or Data
  Int_t ***npPt = (isMC) ? fMCNpPt : fNpPt;

  // -----------------------------------------------------------------------------------------------

  // -- Loop over the pt bins
  for (Int_t idxPt  = 0; idxPt < AliAnalysisNetParticleHelper::fgkfHistNBinsPt; ++idxPt) {
    
    Int_t deltaNp = npPt[idx][1][idxPt]-npPt[idx][0][idxPt]; // p - pbar
    Int_t sumNp   = npPt[idx][1][idxPt]+npPt[idx][0][idxPt]; // p + pbar

    // -- Fill Particle / Anti-Particle Distributions
    (static_cast<TH3D*>(list->FindObject(Form("h%s%s", name, fHelper->GetParticleName(0).Data()))))->Fill(centralityBin, idxPt, npPt[idx][0][idxPt]);
    (static_cast<TH3D*>(list->FindObject(Form("h%s%s", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, idxPt, npPt[idx][1][idxPt]);
    
    // -- Fill NetParticle Distributions
    (static_cast<TH3D*>(list->FindObject(Form("h%sNet%s",  name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, idxPt, deltaNp);
    
    // -- Fill NetParticle vs SumParticle
    Double_t deltaNpOverSumNp = (sumNp == 0.) ? 0. : deltaNp/Double_t(sumNp);
    (static_cast<TH3D*>(list->FindObject(Form("h%sNet%sOverSum", name, fHelper->GetParticleName(1).Data()))))->Fill(centralityBin, idxPt, deltaNpOverSumNp);

    // -----------------------------------------------------------------------------------------------

    // -- Fill TProfile for <NetParticle^k>
    Double_t delta = 1.;
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
      delta *= deltaNp;
      (static_cast<TProfile2D*>(list->FindObject(Form("p%sNet%s%dM", name, fHelper->GetParticleName(1).Data(), idxOrder))))->Fill(centralityBin, idxPt, delta);
      (static_cast<TProfile2D*>(list->FindObject(Form("p%sNet%s%dM_%02d", 
						      name, fHelper->GetParticleName(1).Data(), idxOrder, fHelper->GetSubSampleIdx()))))->Fill(centralityBin, idxPt, delta);
    }
    
    // -- Generate reduced factorials - explictly removing the factorials
    //    - Reset all to 1
    for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0] = 1.;
      fRedFactp[idxOrder][1] = 1.;
    }

    //    - start at idx 1, as idx 0 = 1 done by reset
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
      fRedFactp[idxOrder][0] = fRedFactp[idxOrder-1][0] * Double_t(npPt[idx][0][idxPt]-(idxOrder-1));
      fRedFactp[idxOrder][1] = fRedFactp[idxOrder-1][1] * Double_t(npPt[idx][1][idxPt]-(idxOrder-1));
    }

    // -- Fill TProfiles for <f_ik> 
    TList *fikListPt    = static_cast<TList*>(list->FindObject(Form("f%sPtFik",name)));
    TList *fikListPtSub = static_cast<TList*>(list->FindObject(Form("f%sPtFik_%02d",name, fHelper->GetSubSampleIdx())));
    TH3D  *hCntikPt     = static_cast<TH3D*>(fikListPt->FindObject(Form("p%sNet%sFCounts_%02d", name, fHelper->GetParticleName(1).Data(), Int_t(centralityBin))));
    TH3D  *hCntikPtSub  = static_cast<TH3D*>(fikListPtSub->FindObject(Form("p%sNet%sFCounts_%02d_%02d", name, fHelper->GetParticleName(1).Data(), 
									   Int_t(centralityBin), fHelper->GetSubSampleIdx())));
    for (Int_t ii = 0; ii <= fOrder; ++ii) {   // ii -> p    -> n1
      for (Int_t kk = 0; kk <= fOrder; ++kk) { // kk -> pbar -> n2
	Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   // n1 *n2 -> p * pbar
	(static_cast<TProfile2D*>(fikListPt->FindObject(Form("p%sNet%sF%02d%02d", name, fHelper->GetParticleName(1).Data(), ii, kk))))->Fill(centralityBin, idxPt, fik);
	(static_cast<TProfile2D*>(fikListPtSub->FindObject(Form("p%sNet%sF%02d%02d_%02d", 
								name, fHelper->GetParticleName(1).Data(), ii, kk, fHelper->GetSubSampleIdx()))))->Fill(centralityBin, idxPt, fik);
	
	if (fik != 0.) {
	  hCntikPt->Fill(ii, kk, idxPt);
	  hCntikPtSub->Fill(ii, kk, idxPt);
	}
      }
    }
    
  } // for (Int_t idxPt  = 0; idxPt < AliAnalysisNetParticleHelper::fgkfHistNBinsPt; ++idxPt) {

  return;
}


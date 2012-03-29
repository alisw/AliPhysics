#include "AliAnalysisTaskPIDFluctuation.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliPIDResponse.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include "TFile.h"
#include "TKey.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "AliESDtrackCuts.h"

/* 
 * Event by event PID fluctuation analysis
 * author: Roberto Preghenella (R+)
 * email:  preghenella@bo.infn.it
 *
 */

ClassImp(AliAnalysisTaskPIDFluctuation)

//________________________________________________________________________

const Char_t *AliAnalysisTaskPIDFluctuation::fgkEventCounterName[AliAnalysisTaskPIDFluctuation::kNEventCounters] = {
  "AllEvents",
  "PhysicsSelection",
  "PrimayVertex",
  "PrimayVertexSPD",
  "VertexAccepted",
  "GoodCentrality",
  "AcceptedEvents"
};

const Char_t *AliAnalysisTaskPIDFluctuation::fgkEventCounterTitle[AliAnalysisTaskPIDFluctuation::kNEventCounters] = {
  "all events",
  "physics selection",
  "primary vertex",
  "primary vertex SPD",
  "vertex accepted",
  "good centrality",
  "accepted events"
};

const Char_t *AliAnalysisTaskPIDFluctuation::fgkSparseDataName[AliAnalysisTaskPIDFluctuation::kNSparseData] = {
  "centV0M",
  "centTRK",
  "Nch",
  "Nplus",
  "Nminus",
  "Npi",
  "Npiplus",
  "Npiminus",
  "Nka",
  "Nkaplus",
  "Nkaminus",
  "Npr",
  "Nprplus",
  "Nprminus"
};

const Char_t *AliAnalysisTaskPIDFluctuation::fgkSparseDataTitle[AliAnalysisTaskPIDFluctuation::kNSparseData] = {
  "centrality percentile (V0M)",
  "centrality percentile (TRK)",
  "N_{charged}",
  "N_{+}",
  "N_{-}",
  "N_{#pi}",
  "N_{#pi^{+}}",
  "N_{#pi^{-}}",
  "N_{K}",
  "N_{K^{+}}",
  "N_{K^{-}}",
  "N_{p+#bar{p}}",
  "N_{p}",
  "N_{#bar{p}}"
};

//________________________________________________________________________

AliAnalysisTaskPIDFluctuation::AliAnalysisTaskPIDFluctuation(const Char_t *name) : 
  AliAnalysisTaskSE(name),
  fPIDMethod(kTPCTOF),
  fESDtrackCuts(NULL),
  fAODfilterBit(AliAODTrack::kTrkGlobal),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fPtMin(0.3),
  fPtMax(1.5),
  fPID(NULL),
  fHistoList(NULL),
  fHistoEventCounter(NULL),
  fHistoAcceptedTracks(NULL),
  fHistoTOFMatchedTracks(NULL),
  fHistoTPCdEdx(NULL),
  fHistoTPCdEdx_inclusive(NULL),
  fHistoTOFbeta(NULL),
  fHistoCorrelation(NULL)
{
  
  /*
   * default constructor
   */

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

AliAnalysisTaskPIDFluctuation::~AliAnalysisTaskPIDFluctuation()
{
  
  /*
   * default destructor
   */

  if (fPID) delete fPID;
  if (fHistoList) delete fHistoList;
  
}

//________________________________________________________________________

void AliAnalysisTaskPIDFluctuation::UserCreateOutputObjects() 
{
  
  /*
   * user create output objects
   */

  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);

  fHistoEventCounter = new TH1F("hHistoEventCounter", "", kNEventCounters, 0, kNEventCounters);
  for (Int_t ievc = 0; ievc < kNEventCounters; ievc++)
    fHistoEventCounter->GetXaxis()->SetBinLabel(ievc + 1, fgkEventCounterTitle[ievc]);
  fHistoList->Add(fHistoEventCounter);

  fHistoAcceptedTracks = new TH2F("hHistoAcceptedTracks", ";p_{T} (GeV/c)", 10, 0., 100., 50, 0., 5.);
  fHistoList->Add(fHistoAcceptedTracks);
  
  fHistoTOFMatchedTracks = new TH2F("hHistoTOFMatchedTracks", ";p_{T} (GeV/c)", 10, 0., 100., 50, 0., 5.);
  fHistoList->Add(fHistoTOFMatchedTracks);
  
  fHistoTPCdEdx = new TH3F("hHistoTPCdEdx", ";centrality percentile;p_{TPCin} (GeV/c);dE/dx (au.)", 10, 0., 100., 50, 0., 5., 500, 0., 500.);
  fHistoList->Add(fHistoTPCdEdx);
  
  fHistoTPCdEdx_inclusive = new TH3F("hHistoTPCdEdx_inclusive", ";centrality percentile;p_{TPCin} (GeV/c);dE/dx (au.)", 10, 0., 100., 50, 0., 5., 500, 0., 500.);
  fHistoList->Add(fHistoTPCdEdx_inclusive);
  
  fHistoTOFbeta = new TH3F(Form("hHistoTOFbeta"), ";centrality percentile;p (GeV/c);v/c", 10, 0., 100., 50, 0., 5., 500, 0.1, 1.1);
  fHistoList->Add(fHistoTOFbeta);
  
  /* loop over species */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    
    fHistoTPCdEdx_selected[ipart] = new TH3F(Form("hHistoTPCdEdx_selected_%s", AliPID::ParticleName(ipart)), ";centrality percentile;p_{TPCin} (GeV/c);dE/dx (au.)", 10, 0., 100., 50, 0., 5., 500, 0., 500.);
    fHistoList->Add(fHistoTPCdEdx_selected[ipart]);
    
    fHistoTOFbeta_selected[ipart] = new TH3F(Form("hHistoTOFbeta_selected_%s", AliPID::ParticleName(ipart)), ";centrality percentile;p (GeV/c);v/c", 10, 0., 100., 50, 0., 5., 500, 0.1, 1.1);
    fHistoList->Add(fHistoTOFbeta_selected[ipart]);
    
    fHistoNSigmaTPC[ipart] = new TH3F(Form("hHistoNSigmaTPC_%s", AliPID::ParticleName(ipart)), Form(";centrality percentile;p_{T} (GeV/c); N_{#sigma-TPC}^{%s}", AliPID::ParticleLatexName(ipart)), 10, 0., 100., 50, 0., 5., 200, -10., 10.);
    fHistoList->Add(fHistoNSigmaTPC[ipart]);
    
    fHistoNSigmaTOF[ipart] = new TH3F(Form("hHistoNSigmaTOF_%s", AliPID::ParticleName(ipart)), Form(";centrality percentile;p_{T} (GeV/c); N_{#sigma-TOF}^{%s}", AliPID::ParticleLatexName(ipart)), 10, 0., 100., 50, 0., 5., 200, -10., 10.);
    fHistoList->Add(fHistoNSigmaTOF[ipart]);
    
  } /* end of loop over species */

  Int_t fgSparseDataBins[kNSparseData] = {
    20, 
    20, 
    5000, 
    2500, 
    2500,
    3000,
    1500,
    1500,
    1000,
    500,
    500,
    500,
    250,
    250
  };
  Double_t fgSparseDataMin[kNSparseData] = {
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
  };
  Double_t fgSparseDataMax[kNSparseData] = {
    100., 100., 5000., 2500., 2500., 3000., 1500., 1500., 1000., 500., 500., 500., 250., 250.
  };


  fHistoCorrelation = new THnSparseI("hHistoCorrelation", "", kNSparseData, fgSparseDataBins, fgSparseDataMin, fgSparseDataMax);
  for (Int_t iaxis = 0; iaxis < kNSparseData; iaxis++)
    fHistoCorrelation->GetAxis(iaxis)->SetTitle(fgkSparseDataTitle[iaxis]);
  fHistoList->Add(fHistoCorrelation);
    
  PostData(1, fHistoList);
}

//________________________________________________________________________

void AliAnalysisTaskPIDFluctuation::UserExec(Option_t *) 
{

  /*
   * user exec
   */

  /* get ESD event */
  AliVEvent *event = InputEvent();
  if (!event) return;

  /* accept event */
  if (!AcceptEvent(event)) return;

  /* get centrality object and centrality */
  AliCentrality *centrality = event->GetCentrality();
  Float_t cent_v0m = centrality->GetCentralityPercentileUnchecked("V0M");
  Float_t cent_trk = centrality->GetCentralityPercentileUnchecked("TRK");

  /* init PID */
  InitPID(event);

  Bool_t pidFlag[AliPID::kSPECIES];
  Int_t icharge;
  Int_t chargedCounts[2];
  Int_t pidCounts[AliPID::kSPECIES][2];
  for (Int_t i = 0; i < 2; i++) {
    chargedCounts[i] = 0;
    for (Int_t ii = 0; ii < 5; ii++) {
      pidCounts[ii][i] = 0;
    }
  }

  /* loop over tracks */
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {

    /* get track */
    AliVTrack *track = (AliVTrack *)event->GetTrack(itrk);
    if (!track) continue;

    /* accept track */
    if (!AcceptTrack(track)) continue;

    /* get charge */
    icharge = track->Charge() > 0 ? 0 : 1;
    chargedCounts[icharge]++;

    /* make PID */
    MakePID(track, pidFlag, cent_v0m);
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
      if (pidFlag[ipart])
	pidCounts[ipart][icharge]++;
    
  }

  /* fill histogram */
  Double_t vsparse[kNSparseData];
  vsparse[kCent_V0M] = cent_v0m;
  vsparse[kCent_TRK] = cent_trk;
  vsparse[kNch] = chargedCounts[0] + chargedCounts[1];
  vsparse[kNch_plus] = chargedCounts[0];
  vsparse[kNch_minus] = chargedCounts[1];
  vsparse[kNpi] = pidCounts[AliPID::kPion][0] + pidCounts[AliPID::kPion][1];
  vsparse[kNpi_plus] = pidCounts[AliPID::kPion][0];
  vsparse[kNpi_minus] = pidCounts[AliPID::kPion][1];
  vsparse[kNka] = pidCounts[AliPID::kKaon][0] + pidCounts[AliPID::kKaon][1];
  vsparse[kNka_plus] = pidCounts[AliPID::kKaon][0];
  vsparse[kNka_minus] = pidCounts[AliPID::kKaon][1];
  vsparse[kNpr] = pidCounts[AliPID::kProton][0] + pidCounts[AliPID::kProton][1];
  vsparse[kNpr_plus] = pidCounts[AliPID::kProton][0];
  vsparse[kNpr_minus] = pidCounts[AliPID::kProton][1];
  fHistoCorrelation->Fill(vsparse);

  PostData(1, fHistoList);
}      

//___________________________________________________________

Bool_t
AliAnalysisTaskPIDFluctuation::AcceptEvent(AliVEvent *event) const
{
  /*
   * accept event
   */

  /* fill histo event counter */
  fHistoEventCounter->Fill(kAllEvents);
  
  /* check physics selection */
  if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return kFALSE;
  fHistoEventCounter->Fill(kPhysicsSelection);

  /* check primary vertex */
  const AliVVertex *vertex = event->GetPrimaryVertex();
  if (vertex->GetNContributors() < 1) {
    /* get ESD vertex SPD */
    if (event->InheritsFrom("AliESDEvent")) {
      AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
      vertex = esdevent->GetPrimaryVertexSPD();
    }
    /* get AOD vertex SPD */
    else if (event->InheritsFrom("AliAODEvent")) {
      AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
      vertex = aodevent->GetPrimaryVertexSPD();
    }
    if (vertex->GetNContributors() < 1) return kFALSE;
    fHistoEventCounter->Fill(kPrimaryVertexSPD);
  }
  fHistoEventCounter->Fill(kPrimaryVertex);

  /* check vertex position */
  if (TMath::Abs(vertex->GetZ()) > 10.) return kFALSE;
  fHistoEventCounter->Fill(kVertexAccepted);

  /* get centrality object and check quality */
  AliCentrality *centrality = event->GetCentrality();
  if (centrality->GetQuality() != 0) return kFALSE;
  fHistoEventCounter->Fill(kGoodCentrality);
  
  /* event accepted */
  fHistoEventCounter->Fill(kAcceptedEvents);
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisTaskPIDFluctuation::AcceptTrack(AliVTrack *track) const
{
  /*
   * accept track
   */

  /* check ESD track cuts */
  if (track->InheritsFrom("AliESDtrack")) {
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if (!fESDtrackCuts->AcceptTrack(esdtrack)) return kFALSE;
  }
  /* check AOD filter bit */
  else if (track->InheritsFrom("AliAODTrack")) {
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
    if (!aodtrack->TestFilterBit(fAODfilterBit)) return kFALSE;
  }

  /* check eta range */
  if (track->Eta() < fEtaMin ||
      track->Eta() > fEtaMax) return kFALSE;
  /* check pt range */
  if (track->Pt() < fPtMin ||
      track->Pt() > fPtMax) return kFALSE;

  /* track accepted */
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisTaskPIDFluctuation::HasTPCPID(AliVTrack *track) const
{
  /*
   * has TPC PID
   */
  
  /* check PID signal */
  if (track->GetTPCsignal() <= 0. ||
      track->GetTPCsignalN() == 0) return kFALSE;
  return kTRUE;
  
}
  
//___________________________________________________________

Bool_t
AliAnalysisTaskPIDFluctuation::HasTOFPID(AliVTrack *track) const
{
  /*
   * has TOF PID
   */

  /* check TOF matched track */
  if (!(track->GetStatus() & AliESDtrack::kTOFout)||
      !(track->GetStatus() & AliESDtrack::kTIME)) return kFALSE;
  return kTRUE;

}

//___________________________________________________________

Double_t
AliAnalysisTaskPIDFluctuation::MakeTPCPID(AliVTrack *track, Double_t *nSigma) const
{
  /*
   * make TPC PID
   * returns measured dEdx if PID available, otherwise -1.
   * fills nSigma array with TPC nsigmas for e, mu, pi, K, p
   */
  
  /* check TPC PID */
  if (!HasTPCPID(track)) return -1.;

  /* get TPC info */
  Double_t ptpc = track->GetTPCmomentum();
  Double_t dEdx = track->GetTPCsignal();
  Int_t dEdxN = track->GetTPCsignalN();
    
  /* loop over particles */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    Double_t bethe = fPID->GetTPCResponse().GetExpectedSignal(ptpc, (AliPID::EParticleType)ipart);
    Double_t diff = dEdx - bethe;
    Double_t sigma = fPID->GetTPCResponse().GetExpectedSigma(ptpc, dEdxN, (AliPID::EParticleType)ipart);
    nSigma[ipart] = diff / sigma;
  }
  return dEdx;

}

//___________________________________________________________

Double_t
AliAnalysisTaskPIDFluctuation::MakeTOFPID(AliVTrack *track, Double_t *nSigma) const
{
  /*
   * make TOF PID
   * returns measured beta if PID available, otherwise -1.
   * fills nSigma array with TOF nsigmas for e, mu, pi, K, p
   */
  
  /* check TOF PID */
  if (!HasTOFPID(track)) return -1.;

  /* get TOF info */
  Double_t p = track->P();
  Double_t time = track->GetTOFsignal() - fPID->GetTOFResponse().GetStartTime(p);
  Double_t timei[5];
  track->GetIntegratedTimes(timei);
  
  /* loop over particles */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    Double_t timez = time - timei[ipart];
    Double_t sigma = fPID->GetTOFResponse().GetExpectedSigma(p, timei[ipart], AliPID::ParticleMass(ipart));
    nSigma[ipart] = timez / sigma;
  }
  
  return timei[0] / time;

}

//___________________________________________________________

void
AliAnalysisTaskPIDFluctuation::MakePID(AliVTrack *track, Bool_t *pidFlag, Float_t centrality) const
{
  /*
   * make PID
   * fills PID QA plots
   * fills pidFlag array with PID flags for e, mu, pi, K, p
   */
  
  /* cut definitions
     (better put them as static variables so they can be changed from outside) */
  Double_t fgTPCPIDmomcut[AliPID::kSPECIES] = {0., 0., 0.5, 0.5, 0.7};
  Double_t fgTPCPIDsigmacut[AliPID::kSPECIES] = {0., 0., 2., 2., 2.};
  Double_t fgTPCPIDcompcut[AliPID::kSPECIES] = {0., 0., 3., 3., 3.};
  Double_t fgTOFPIDmomcut[AliPID::kSPECIES] = {0., 0., 1.5, 1.5, 2.0};
  Double_t fgTOFPIDsigmacut[AliPID::kSPECIES] = {0., 0., 2., 2., 2.};

  /* get momentum information */
  Double_t p = track->P();
  Double_t pt = track->Pt();
  Double_t ptpc = track->GetTPCmomentum();

  /* make pid and check if available */
  Double_t nsigmaTPC[AliPID::kSPECIES];
  Double_t nsigmaTOF[AliPID::kSPECIES];
  Double_t dEdx = MakeTPCPID(track, nsigmaTPC);
  Double_t beta = MakeTOFPID(track, nsigmaTOF);
  Bool_t hasTPCPID = dEdx > 0.;
  Bool_t hasTOFPID = beta > 0.;

  /* check PID method */
  if (fPIDMethod == kTPConly) hasTOFPID = kFALSE; // inhibit TOF PID
  if (fPIDMethod == kTOFonly) hasTPCPID = kFALSE; // inhibit TPC PID

  /* fill qa histos */
  fHistoAcceptedTracks->Fill(centrality, pt);
  if (hasTPCPID)
    fHistoTPCdEdx->Fill(centrality, ptpc, dEdx);
  if (hasTOFPID) {
    fHistoTOFMatchedTracks->Fill(centrality, pt);
    fHistoTOFbeta->Fill(centrality, p, beta);
  }
  if (hasTPCPID && !hasTOFPID)
    fHistoTPCdEdx_inclusive->Fill(centrality, ptpc, dEdx);
    
  /* loop over species */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    
    /* reset PID flag */
    pidFlag[ipart] = kFALSE;
    
    /* fill qa histos */
    if (hasTPCPID)
      fHistoNSigmaTPC[ipart]->Fill(centrality, pt, nsigmaTPC[ipart]);
    if (hasTOFPID)
      fHistoNSigmaTOF[ipart]->Fill(centrality, pt, nsigmaTOF[ipart]);
    
    /* combined PID cuts */
    if (hasTPCPID && hasTOFPID) {
      if (pt < fgTOFPIDmomcut[ipart] &&
	  TMath::Abs(nsigmaTOF[ipart]) < fgTOFPIDsigmacut[ipart] &&
	  TMath::Abs(nsigmaTPC[ipart]) < fgTPCPIDcompcut[ipart]) {
	fHistoTOFbeta_selected[ipart]->Fill(centrality, p, beta);
	pidFlag[ipart] = kTRUE;
      }
    }
    /* TPC-only PID cuts */
    else if (hasTPCPID && !hasTOFPID) {
      if (pt < fgTPCPIDmomcut[ipart] &&
	  TMath::Abs(nsigmaTPC[ipart]) < fgTPCPIDsigmacut[ipart]) { 
	fHistoTPCdEdx_selected[ipart]->Fill(centrality, ptpc, dEdx);
	pidFlag[ipart] = kTRUE;
      }
    }
    /* TOF-only PID cuts */
    else if (!hasTPCPID && hasTOFPID) {
      if (pt < fgTOFPIDmomcut[ipart] &&
	  TMath::Abs(nsigmaTOF[ipart]) < fgTOFPIDsigmacut[ipart]) {
	fHistoTOFbeta_selected[ipart]->Fill(centrality, p, beta);
	pidFlag[ipart] = kTRUE;
      }
    }
    
  } /* end of loop over species */

}

//___________________________________________________________

void
AliAnalysisTaskPIDFluctuation::InitPID(AliVEvent *event)
{
  /*
   * init PID
   */

  /* create PID object if not there yet */
  if (!fPID) {

    /* instance object */
    Bool_t mcFlag = kFALSE; /*** WARNING: check whether is MC ***/
    if (event->InheritsFrom("AliESDEvent"))
      fPID = new AliESDpid(mcFlag);
    else if (event->InheritsFrom("AliAODEvent"))
      fPID = new AliAODpidUtil(mcFlag);

    /* set OADB path */
    fPID->SetOADBPath("$ALICE_ROOT/OADB");
  }

  /* init ESD PID */
  Int_t recoPass = 2; /*** WARNING: need to set the recoPass somehow better ***/
  fPID->InitialiseEvent(event, recoPass); /* warning: this apparently sets TOF time        
					      * resolution to some hardcoded value,     
					      * therefore we have to set correct   
					      * resolution value after this call */

  /* set TOF response */
  Double_t tofReso = 85.; /* ps */ /*** WARNING: need to set tofReso somehow better ***/
  fPID->GetTOFResponse().SetTimeResolution(tofReso);
  fPID->GetTOFResponse().SetTrackParameter(0, 0.007);
  fPID->GetTOFResponse().SetTrackParameter(1, 0.007);
  fPID->GetTOFResponse().SetTrackParameter(2, 0.0);
  fPID->GetTOFResponse().SetTrackParameter(3, 30);
  
}

//___________________________________________________________

void
AliAnalysisTaskPIDFluctuation::MeasureNuDyn(const Char_t *filename, Int_t i1, Int_t i2, Int_t centralityEstimator)
{

  /*
   * measure nu
   */
  
  printf("MeasureNuDyn: running for %s vs. %s\n", fgkSparseDataName[i1], fgkSparseDataName[i2]);

  /* get data */
  TFile *filein = TFile::Open(filename);
  /* output */
  TFile *fileout = TFile::Open(Form("MeasureNuDyn_%s_%s.%s", fgkSparseDataName[i1], fgkSparseDataName[i2], filename), "RECREATE");

  /* loop over available containers */
  TList *keylist = filein->GetListOfKeys();
  for (Int_t ikey = 0; ikey < keylist->GetEntries(); ikey++) {
    
    /* get key and check */
    TKey *key = (TKey *)keylist->At(ikey);
    TString contname = key->GetName();
    if (!contname.BeginsWith("PIDFluctuation")) continue;

    /* get data */
    TList *list = (TList *)filein->Get(contname.Data());
    THnSparse *hsparse = (THnSparse *)list->FindObject("hHistoCorrelation");
    
    /* create output directory and cd there */
    fileout->mkdir(contname.Data());
    fileout->cd(contname.Data());

    /* loop over centralities */
    Double_t centBin[kNCentralityBins + 1] = {
      0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.
    };
    TH1D *hNu = new TH1D("hNu", ";centrality percentile;#nu", kNCentralityBins, centBin);
    TH1D *hNuStat = new TH1D("hNuStat", ";centrality percentile;#nu_{stat}", kNCentralityBins, centBin);
    TH1D *hNuDyn = new TH1D("hNuDyn", ";centrality percentile;#nu_{dyn}", kNCentralityBins, centBin);
    for (Int_t icent = 0; icent < kNCentralityBins; icent++) {
      
      /* select centrality range */
      hsparse->GetAxis(centralityEstimator)->SetRangeUser(centBin[icent] + 0.001, centBin[icent + 1] - 0.001);
      /* projections */
      TH1D *hcent = hsparse->Projection(centralityEstimator);
      TH1D *h1 = hsparse->Projection(i1);
      TH1D *h2 = hsparse->Projection(i2);
      TH2D *hcorr = hsparse->Projection(i2, i1);
      TH1D *hnu = new TH1D("hnu", "", 2000, -1., 1.);
      
      Double_t n, a, b, nev, meanx, rmsx, meany, rmsy, meannu, rmsnu;
      Double_t meanx_err, meany_err, meannu_err;
      
      /* compute mean values */
      nev = 0.; meanx = 0.; meany = 0.;
      for (Int_t ibinx = 0; ibinx < hcorr->GetNbinsX(); ibinx++)
	for (Int_t ibiny = 0; ibiny < hcorr->GetNbinsY(); ibiny++) {
	  n = hcorr->GetBinContent(ibinx + 1, ibiny + 1);
	  if (n <= 0.) continue;
	  meanx += n * ibinx;
	  meany += n * ibiny;
	  nev  += n;
	}
      meanx /= nev;
      meany /= nev;
      //    printf("nev = %f, meanx = %f, meany = %f\n", nev, meanx, meany);
      
      /* compute RMS values */
      nev = 0.; rmsx = 0.; rmsy = 0.;
      for (Int_t ibinx = 0; ibinx < hcorr->GetNbinsX(); ibinx++)
	for (Int_t ibiny = 0; ibiny < hcorr->GetNbinsY(); ibiny++) {
	  n = hcorr->GetBinContent(ibinx + 1, ibiny + 1);
	  if (n <= 0.) continue;
	  a = ibinx - meanx;
	  rmsx += n * a * a;
	  a = ibiny - meany;
	  rmsy += n * a * a;
	  nev  += n;
	}
      rmsx /= nev;
      rmsx = TMath::Sqrt(rmsx);
      rmsy /= nev;
      rmsy = TMath::Sqrt(rmsy);
      //    printf("nev = %f, rmsx = %f, rmsy = %f\n", nev, rmsx, rmsy);
      meanx_err = rmsx / TMath::Sqrt(nev);
      meany_err = rmsy / TMath::Sqrt(nev);
      //    printf("nev = %f, meanx_err = %f, meany_err = %f\n", nev, meanx_err, meany_err);
      
      /* compute mean nu */
      nev = 0.; meannu = 0.;
      for (Int_t ibinx = 0; ibinx < hcorr->GetNbinsX(); ibinx++)
	for (Int_t ibiny = 0; ibiny < hcorr->GetNbinsY(); ibiny++) {
	  n = hcorr->GetBinContent(ibinx + 1, ibiny + 1);
	  if (n <= 0.) continue;
	  a = ibinx / meanx - ibiny / meany;
	  meannu += n * a * a;
	  hnu->Fill(a * a, n);
	  nev  += n;
	}
      meannu /= nev;
      //    printf("nev = %f, meannu = %f\n", nev, meannu);
      
      /* compute RMS nu */
      nev = 0.; rmsnu = 0.;
      for (Int_t ibinx = 0; ibinx < hcorr->GetNbinsX(); ibinx++)
	for (Int_t ibiny = 0; ibiny < hcorr->GetNbinsY(); ibiny++) {
	  n = hcorr->GetBinContent(ibinx + 1, ibiny + 1);
	  if (n <= 0.) continue;
	  a = ibinx / meanx - ibiny / meany;
	  b = a * a - meannu;
	  rmsnu += n * b * b;
	  nev  += n;
	}
      rmsnu /= nev;
      rmsnu = TMath::Sqrt(rmsnu);
      //    printf("nev = %f, rmsnu = %f\n", nev, rmsnu);
      meannu_err = rmsnu / TMath::Sqrt(nev);
      //    printf("nev = %f, meannu_err = %f\n", nev, meannu_err);
      
      /* final calculations */
      Double_t nu = meannu;
      Double_t nu_err = meannu_err;
      Double_t nu_stat = 1. / meanx + 1. / meany;
      Double_t meanx4 = meanx * meanx * meanx * meanx;
      Double_t meanx_err2 = meanx_err * meanx_err;
      Double_t meany4 = meany * meany * meany * meany;
      Double_t meany_err2 = meany_err * meany_err;
      Double_t nu_stat_err = TMath::Sqrt(meanx_err2 / meanx4 + meany_err2 / meany4);
      Double_t nu_dyn = nu - nu_stat;
      Double_t nu_dyn_err = TMath::Sqrt(nu_err * nu_err + nu_stat_err * nu_stat_err);
      
      /* setup final plots */
      hNu->SetBinContent(icent + 1, nu);
      hNu->SetBinError(icent + 1, nu_err);
      hNuStat->SetBinContent(icent + 1, nu_stat);
      hNuStat->SetBinError(icent + 1, nu_stat_err);
      hNuDyn->SetBinContent(icent + 1, nu_dyn);
      hNuDyn->SetBinError(icent + 1, nu_dyn_err);
      
      hcent->Write(Form("hcent_cent%d", icent));
      h1->Write(Form("h1_cent%d", icent));
      h2->Write(Form("h2_cent%d", icent));
      hcorr->Write(Form("hcorr_cent%d", icent));
      hnu->Write(Form("hnu_cent%d", icent));
      
      /* clean-up */
      delete hcent;
      delete h1;
      delete h2;
      delete hcorr;
      delete hnu;
    }

    hNu->Write();
    hNuStat->Write();
    hNuDyn->Write();
  }

  fileout->Close();
  filein->Close();

}

//___________________________________________________________

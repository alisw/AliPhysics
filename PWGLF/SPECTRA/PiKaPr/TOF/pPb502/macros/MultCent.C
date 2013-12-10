#include "CommonDefs.C"

MultCent(const Char_t *filename = "data.root", Int_t evMax = kMaxInt)
{
  
  /* include path for ACLic */
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
  /* load libraries */
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  /* build analysis task class */
  gROOT->LoadMacro("AliAnalysisParticle.cxx+g");
  gROOT->LoadMacro("AliAnalysisEvent.cxx+g");
  gROOT->LoadMacro("AliAnalysisTrack.cxx+g");
 
  /* open file, get tree and connect */
  TFile *filein = TFile::Open(filename);
  TTree *treein = (TTree *)filein->Get("aodTree");
  printf("got \"aodTree\": %d entries\n", treein->GetEntries());
  AliAnalysisEvent *analysisEvent = new AliAnalysisEvent();
  TClonesArray *analysisTrackArray = new TClonesArray("AliAnalysisTrack");
  AliAnalysisTrack *analysisTrack = NULL;
  treein->SetBranchAddress("AnalysisEvent", &analysisEvent);
  treein->SetBranchAddress("AnalysisTrack", &analysisTrackArray);

  /* histo */
  TH2F *hMultCent = new TH2F("hMultCent", "", 120, -10., 110., 5000, 0., 5000.);
  TH2F *hRefMultCent = new TH2F("hRefMultCent", "", 120, -10., 110., 5000, 0., 5000.);
  TH2F *hMultTOFCent = new TH2F("hMultTOFCent", "", 100, 0., 100., 5000, 0., 5000.);
  TH2F *hMultCent_T0TOF = new TH2F("hMultCent_T0TOF", "", 100, 0., 100., 5000, 0., 5000.);

  /* loop over events */
  Int_t mult;
  Double_t cent, refmult;
  for (Int_t iev = 0; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 100000 == 0) printf("iev = %d\n", iev);

    /* check accept event */
    if (!analysisEvent->AcceptEvent(acceptEventType)) continue;

    /* count tracks for multiplicity */
    mult = analysisTrackArray->GetEntries();
    cent = analysisEvent->GetCentralityPercentile(centralityEstimator);
    refmult = analysisEvent->GetCentralityPercentile(999);
    /* check centrality percentile (V0M) */
    //    if (cent == 100.) continue;

    /*** ACCEPTED EVENT ***/

    /* fill histo */
    hMultCent->Fill(cent, mult);
    hRefMultCent->Fill(cent, refmult);

    /* check if event has T0-TOF */
    if (analysisEvent->GetTimeZeroTOFSigma(9) > 150.) continue;

    hMultCent_T0TOF->Fill(cent, mult);

    Int_t nTOFtrk = 0, ntrk = 0;
    AliAnalysisTrack *analysisTrack = NULL;
    for (Int_t itrk = 0; itrk < mult; itrk++) {
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack->AcceptTrack()) continue;
      ntrk++;
      if (analysisTrack->HasTOFMatch()) nTOFtrk++;
      
    }

    hMultTOFCent->Fill(cent, nTOFtrk);

  }

  /* output */
  TFile *fileout = new TFile(Form("MultCent.%s", filename), "RECREATE");
  hMultCent->Write();
  hRefMultCent->Write();
  hMultCent_T0TOF->Write();
  hMultTOFCent->Write();
  fileout->Close();

}

MultCent_T0TOFeff(const Char_t *filename)
{

  TFile *filein = TFile::Open(filename);
  TH2F *hMultCent = (TH2F *)filein->Get("hMultCent");
  TH2F *hMultCent_T0TOF = (TH2F *)filein->Get("hMultCent_T0TOF");

  /* projections */
  TH1D *hMultCent_px = hMultCent->ProjectionX();
  TH1D *hMultCent_T0TOF_px = hMultCent_T0TOF->ProjectionX();
  TH1D *hMultCent_py = hMultCent->ProjectionY();
  TH1D *hMultCent_T0TOF_py = hMultCent_T0TOF->ProjectionY();

  /* efficiency centrality */
  TH1D *hEfficiency_centrality = new TH1D(*hMultCent_T0TOF_px);
  hEfficiency_centrality->SetNameTitle("hEfficiency_centrality", ";centrality percentile;T0-TOF efficiency");
  hEfficiency_centrality->Sumw2();
  hEfficiency_centrality->Divide(hEfficiency_centrality, hMultCent_px, 1., 1., "B");

  /* efficiency tracks */
  TH1D *hEfficiency_tracks = new TH1D(*hMultCent_T0TOF_py);
  hEfficiency_tracks->SetNameTitle("hEfficiency_tracks", ";N_{tracks};T0-TOF efficiency");
  hEfficiency_tracks->Sumw2();
  hEfficiency_tracks->Divide(hEfficiency_tracks, hMultCent_py, 1., 1., "B");

  /* draw */
  TCanvas *cCentrality = new TCanvas("cCentrality");
  hEfficiency_centrality->DrawCopy();
  TCanvas *cTracks = new TCanvas("cTracks");
  hEfficiency_tracks->DrawCopy();

  /* output */
  TFile *fileout = new TFile(Form("MultCent_T0TOFeff.%s", filename), "RECREATE");
  hEfficiency_centrality->Write();
  hEfficiency_tracks->Write();
  fileout->Close();


}

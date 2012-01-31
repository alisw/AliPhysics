MultCent(const Char_t *filename, Int_t evMax = kMaxInt)
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
  TH2F *hMultCent = new TH2F("hMultCent", "", 100, 0., 100., 5000, 0., 5000.);
  TH2F *hMultCent_T0TOF = new TH2F("hMultCent_T0TOF", "", 100, 0., 100., 5000, 0., 5000.);

  /* loop over events */
  Int_t mult;
  Double_t cent;
  for (Int_t iev = 0; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 1000 == 0) printf("iev = %d\n", iev);
    /* check vertex */
    if (!analysisEvent->AcceptVertex()) continue;
    /* check collision candidate */
    if (!analysisEvent->IsCollisionCandidate()) continue;
    /* count tracks for multiplicity */
    mult = analysisTrackArray->GetEntries();
    /* check centrality quality */
    if (analysisEvent->GetCentralityQuality() != 0.) continue;
    cent = analysisEvent->GetCentralityPercentile(AliAnalysisEvent::kCentEst_V0M);
    /* check centrality percentile (V0M) */
    if (cent == 100.) continue;

    /*** ACCEPTED EVENT ***/

    /* fill histo */
    hMultCent->Fill(cent, mult);

    /* check if event has T0-TOF */
    if (analysisEvent->GetTimeZeroTOFSigma(9) > 150.) continue;

    hMultCent_T0TOF->Fill(cent, mult);

  }

  /* output */
  TFile *fileout = new TFile(Form("MultCent.%s", filename), "RECREATE");
  hMultCent->Write();
  hMultCent_T0TOF->Write();
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

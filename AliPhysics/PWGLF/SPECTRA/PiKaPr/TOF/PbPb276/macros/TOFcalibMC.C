enum ECharge_t {
  kPositive,
  kNegative,
  kNCharges
};
const Char_t *chargeName[kNCharges] = {
  "positive",
  "negative",
};

TOFcalibMC_centrality(const Char_t *filename, Int_t evMax = kMaxInt)
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

  /* histos */
  TH2F *hTOFcalib_centrality = new TH2F("hTOFcalib_centrality", "", 18, 0., 90., 200, -2440., 2440.);

  /* loop over events */
  Double_t cent, p, time, t0, tof, texp, deltat, timezerocorr, texpcorr;
  for (Int_t iev = 0; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 1000 == 0) printf("iev = %d\n", iev);
    /* check vertex */
    if (!analysisEvent->AcceptVertex()) continue;
    /* check collision candidate */
    if (!analysisEvent->IsCollisionCandidate()) continue;
    /* check centrality quality */
    if (analysisEvent->GetCentralityQuality() != 0.) continue;
    
    /*** ACCEPTED EVENT ***/

    /* apply time-zero TOF correction */
    analysisEvent->ApplyTimeZeroTOFCorrection();
    
    /* get centrality */
    cent = analysisEvent->GetCentralityPercentile(AliAnalysisEvent::kCentEst_V0M);

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get track */
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check TOF PID */
      if (!analysisTrack->HasTOFPID()) continue;

      /*** ACCEPTED TRACK WITH TOF PID ***/

      /* apply expected time correction */
      analysisTrack->ApplyTOFExpectedTimeCorrection();

      p = analysisTrack->GetP();
      time = analysisTrack->GetTOFTime();
      t0 = analysisEvent->GetTimeZeroTOF(p);
      tof = time - t0;
      texp = analysisTrack->GetTOFExpTime(AliPID::kPion);
      deltat = tof - texp;

      hTOFcalib_centrality->Fill(cent, deltat);

      } /* end of loop over particles */
  } /* end of loop over events */

  /* output */
  TFile *fileout = TFile::Open(Form("TOFcalibMC_centrality.%s", filename), "RECREATE");
  hTOFcalib_centrality->Write();
  fileout->Close();

}

TOFcalibMC_texp(const Char_t *filename, Int_t evMax = kMaxInt)
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

  /* histos */
  TH2F *hTOFcalib_texp[AliPID::kSPECIES][kNCharges];
  TH2F *hTOFcalib_texpr[AliPID::kSPECIES][kNCharges];
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hTOFcalib_texp[ipart][icharge] = new TH2F(Form("hTOFcalib_texp_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 100., 0., 5., 400, -4880., 4880.);
      hTOFcalib_texpr[ipart][icharge] = new TH2F(Form("hTOFcalib_texpr_%s_%s", AliPID::ParticleName(ipart), chargeName[icharge]), "", 100., 0., 5., 400, 0.8, 1.2);
    }
  }

  /* loop over events */
  Int_t charge;
  Double_t cent, p, pt, time, t0, tof, texp, deltat, timezerocorr, texpcorr;
  for (Int_t iev = 0; iev < treein->GetEntries() && iev < evMax; iev++) {
    /* get event */
    treein->GetEvent(iev);
    if (iev % 1000 == 0) printf("iev = %d\n", iev);
    /* check vertex */
    if (!analysisEvent->AcceptVertex()) continue;
    /* check collision candidate */
    if (!analysisEvent->IsCollisionCandidate()) continue;
    /* check centrality quality */
    if (analysisEvent->GetCentralityQuality() != 0.) continue;
    
    /*** ACCEPTED EVENT ***/

    /* apply time-zero TOF correction */
    analysisEvent->ApplyTimeZeroTOFCorrection();
    
    /* get centrality */
    cent = analysisEvent->GetCentralityPercentile(AliAnalysisEvent::kCentEst_V0M);
    //    if (cent > 90.) continue;

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < analysisTrackArray->GetEntries(); itrk++) {
      /* get track */
      analysisTrack = (AliAnalysisTrack *)analysisTrackArray->At(itrk);
      if (!analysisTrack) continue;
      /* check eta */
      if (TMath::Abs(analysisTrack->GetEta()) > 0.8) continue;
      /* check TOF PID */
      if (!analysisTrack->HasTOFPID()) continue;
      /* get charge */
      charge = analysisTrack->GetSign() > 0. ? kPositive : kNegative;
      /* check charged primary with defined PID */
      Int_t ipart = analysisTrack->GetMCPID();
      if (!analysisTrack->IsMCPrimary() || ipart < 0 || analysisTrack->GetSign() == 0.) continue;
      /* check mismatch */
      if (analysisTrack->IsMismatchMC()) continue;

      /*** ACCEPTED TRACK WITH TOF PID ***/

      /* apply expected time correction */
      //      analysisTrack->ApplyTOFExpectedTimeCorrection();

      p = analysisTrack->GetP();
      time = analysisTrack->GetTOFTime();
      t0 = analysisEvent->GetTimeZeroTOF(p);
      tof = time - t0;

      texp = analysisTrack->GetTOFExpTime(ipart);
      deltat = tof - texp;
      
      hTOFcalib_texp[ipart][charge]->Fill(p, deltat);
      hTOFcalib_texpr[ipart][charge]->Fill(p, tof / texp);
      
    } /* end of loop over tracks */
  } /* end of loop over events */

  /* output */
  TFile *fileout = TFile::Open(Form("TOFcalibMC_texp.%s", filename), "RECREATE");
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    for (Int_t icharge = 0; icharge < kNCharges; icharge++) {
      hTOFcalib_texp[ipart][icharge]->Write();
      hTOFcalib_texpr[ipart][icharge]->Write();
    }
  }
  fileout->Close();

}


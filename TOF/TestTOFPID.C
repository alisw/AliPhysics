testPID(const Char_t *filename, Bool_t calibrateESD = kTRUE, Bool_t correctTExp = kTRUE, Bool_t useT0TOF = kTRUE, Double_t timeResolution = 100., Bool_t tuneTOFMC = kFALSE)
{
  /* PID analysis */

  /* check MC flag */
  if (tuneTOFMC) calibrateESD = kFALSE;

  /* init ESD */
  TFile *filein = TFile::Open(filename);
  TTree *treein = (TTree *)filein->Get("esdTree");
  AliESDEvent *event = new AliESDEvent();
  event->ReadFromTree(treein);
  /* init OCDB */
  treein->GetEvent(0);
  Int_t run = event->GetRunNumber();
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(run);
  /* init TOF calibration */
  AliTOFcalib *tofCalib = new AliTOFcalib();
  if (correctTExp)
    tofCalib->SetCorrectTExp(kTRUE);
  tofCalib->Init();
  /* init TOF T0-maker */
  AliESDpid *fPIDesd = new AliESDpid();
  AliTOFT0maker *t0maker = new AliTOFT0maker(fPIDesd, tofCalib); 
  t0maker->SetTimeResolution(timeResolution);

  /* pid histos */
  TH2F *hTOFpid[AliPID::kSPECIES];
  TH2F *hTOFpidSigma[AliPID::kSPECIES];
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    hTOFpid[ipart] = new TH2F(Form("hTOFpid_%s", AliPID::ParticleName(ipart)), Form("%s-ID;p (GeV/c);t - t^{%s}_{exp} (ps)", AliPID::ParticleName(ipart), AliPID::ParticleLatexName(ipart)), 200, 0., 10., 2000, -24400., 24400.);
    hTOFpidSigma[ipart] = new TH2F(Form("hTOFpidSigma_%s", AliPID::ParticleName(ipart)), Form("%s-ID;p (GeV/c);(t - t^{%s}_{exp}) / #sigma^{%s}_{exp} (ps)", AliPID::ParticleName(ipart), AliPID::ParticleLatexName(ipart), AliPID::ParticleLatexName(ipart)), 200, 0., 10., 2000, -10., 10.);
  }

  AliESDtrack *track;
  Double_t p, time, timei[AliPID::kSPECIES], sigma[AliPID::kSPECIES];
  /* loop over events */
  for (Int_t iev = 0; iev < treein->GetEntries(); iev++) {
    
    /* get event */
    treein->GetEvent(iev);

    /* calibrate ESD */
    if (calibrateESD)
      tofCalib->CalibrateESD(event);

    /* tune TOF if requested for MC */
    if (tuneTOFMC) {
      t0maker->TuneForMC(event);
    }

    /* compute and apply T0-TOF */
    if (useT0TOF) {
      t0maker->ComputeT0TOF(event);
      t0maker->ApplyT0TOF(event);
      fPIDesd->MakePID(event,kFALSE,0.);
    }

    /* loop over tracks */
    for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
      track = event->GetTrack(itrk);
      /* check TOF match */
      if (!track || 
	  !(track->GetStatus() & AliESDtrack::kTOFout) ||
	  !(track->GetStatus() & AliESDtrack::kTIME)) continue;

      /* get track momentum */
      p = track->P();
      /* get TOF time */
      time = track->GetTOFsignal();
      /* get expected times */
      track->GetIntegratedTimes(timei);

      /* fill PID histos */
      for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++){
	sigma[ipart] = t0maker->GetExpectedSigma(p, timei[ipart], AliPID::ParticleMass(ipart));
	hTOFpid[ipart]->Fill(p, (time - timei[ipart])); 
	hTOFpidSigma[ipart]->Fill(p, (time - timei[ipart]) / sigma[ipart]); 
      }
    }
  }
  
  /* write output */
  TFile *fileout = TFile::Open("testPID.root", "RECREATE");
  for (Int_t ipart = 0; ipart < 5; ipart++){
    hTOFpid[ipart]->Write();
    hTOFpidSigma[ipart]->Write();
  }
  fileout->Close();
}

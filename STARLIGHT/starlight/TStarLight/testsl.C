// -*- C++ -*-
// $Id$

void testsl() {
  gSystem->Load("libStarLight");
  gSystem->Load("libAliStarLight");

  TStarLight* sl = new TStarLight("starlight generator", "title", "");

  sl->SetParameter("baseFileName = slight	#suite of output files will be saved with this base name");
  sl->SetParameter("BEAM_1_Z = 82    #Z of projectile");
  sl->SetParameter("BEAM_1_A = 208   #A of projectile");
  sl->SetParameter("BEAM_2_Z = 82   #Z of target");
  sl->SetParameter("BEAM_2_A = 208   #A of target");
  sl->SetParameter("BEAM_1_GAMMA = 1470.0 #Gamma of the colliding ion 1");
  sl->SetParameter("BEAM_2_GAMMA = 1470.0 #Gamma of the colliding ion 2");
  sl->SetParameter("W_MAX = -1   #Max value of w");
  sl->SetParameter("W_MIN = -1    #Min value of w");
  sl->SetParameter("W_N_BINS = 50    #Bins i w");
  sl->SetParameter("RAP_MAX = 9.    #max y");
  sl->SetParameter("RAP_N_BINS = 200    #Bins i y");
  sl->SetParameter("CUT_PT = 0 #Cut in pT? 0 = (no, 1 = yes)");
  sl->SetParameter("PT_MIN = 1.0 #Minimum pT in GeV");
  sl->SetParameter("PT_MAX = 3.0 #Maximum pT in GeV");
  sl->SetParameter("CUT_ETA = 0 #Cut in pseudorapidity? (0 = no, 1 = yes)");
  sl->SetParameter("ETA_MIN = -10 #Minimum pseudorapidity");
  sl->SetParameter("ETA_MAX = 10 #Maximum pseudorapidity");
  sl->SetParameter("PROD_MODE = 2     #gg or gP switch (1 = 2-photon, 2 = coherent vector meson (narrow), 3 = coherent vector meson (wide), 4 = incoherent vector meson)");
  sl->SetParameter("N_EVENTS = 1000   #Number of events");
  sl->SetParameter("PROD_PID = 443013   #Channel of interest; this is j/psi --> mu+ mu-");
  sl->SetParameter("RND_SEED = 5574533 #Random number seed");
  sl->SetParameter("BREAKUP_MODE = 5     #Controls the nuclear breakup; a 5 here makes no requirement on the breakup of the ions");
  sl->SetParameter("INTERFERENCE = 0     #Interference (0 = off, 1 = on)");
  sl->SetParameter("IF_STRENGTH = 1.    #percent of intefernce (0.0 - 0.1)");
  sl->SetParameter("INT_PT_MAX = 0.24  #Maximum pt considered, when interference is turned on");
  sl->SetParameter("INT_PT_N_BINS =120   #Number of pt bins when interference is turned on");
  sl->SetParameter("XSEC_METHOD = 1 # Set to 0 to use old method for calculating gamma-gamma luminosity");
  sl->SetParameter("PYTHIA_FULL_EVENTRECORD = 0 # Write full pythia information to output (vertex, parents, daughter etc).");

  sl->InitStarLight();
  sl->PrintInputs(std::cout);
  TClonesArray tca("TParticle", 100);

  TLorentzVector v[2], vSum;
  TH1* hM  = new TH1D("hM",  "STARLIGHT;M#(){#pi^{+}#pi^{-}}",    200, 3.0, 3.2);
  TH1* hPt = new TH1D("hPt", "STARLIGHT;P_{T}#(){#pi^{+}#pi^{-}}", 80, 0., 2.);
  TH1* hY  = new TH1D("hY",  "STARLIGHT;Y#(){#pi^{+}#pi^{-}}",    100,-10., 10.);

  std::ofstream ofs("sl.txt");
  TParticle *p;
  for (Int_t counter(0); counter<20000; ) {
    sl->GenerateEvent();
    sl->BoostEvent();
    sl->ImportParticles(&tca, "ALL");
    Bool_t genOK = kTRUE;
    TLorentzVector vSum;
    for (Int_t i=0; i<tca.GetEntries() && genOK; ++i) {
      p = (TParticle*)tca.At(i);
      p->Momentum(v[i]);
      vSum += v[i];
//       genOK = TMath::Abs(v[i].Rapidity()) <= 1.5;
    }
    tca.Clear();
    if (!genOK) continue;
    Printf("%5d %d", counter, genOK);
    ++counter;
    vSum = v[0] + v[1];
    ofs << std::fixed << std::setprecision(4)
	<< vSum.M() << " " << vSum.Perp() << " " << vSum.Rapidity() << " "
	<< v[0].Eta() << " " << v[0].Px() << " " << v[0].Py() << " " << v[0].Pz() << " "
	<< v[1].Eta() << " " << v[1].Px() << " " << v[1].Py() << " " << v[1].Pz()
	<< std::endl;
    hM->Fill(vSum.M());
    hPt->Fill(vSum.Perp());
    hY->Fill(vSum.Rapidity());
  }
  TFile::Open("sl.root", "RECREATE");
  sl->Write();
  gFile->Write();

  hM->Draw();
  c1->SaveAs("SL.pdf(");
  hPt->Draw();
  c1->SaveAs("SL.pdf");
  hY->Draw();
  c1->SaveAs("SL.pdf)");
}

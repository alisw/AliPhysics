// AliROOT DIME MC TDime quick test macro
//
// Mikael.Mieskolainen@cern.ch, 22.9.2015

void test() {

  gSystem->Load("libEVGEN");
  gSystem->Load("libTDime");
  gSystem->Load("libdime");

  // 0. Create TDime object
  TDime* dime = new TDime();

  // 1. Setup random seed by system time
  // (note that only ~1 sec resolution, might be problematic in grid)
  TDatime dt;
  UInt_t seed = dt.Get();
  AliDimeRndm::GetDimeRandom()->SetSeed(seed);
  cout<< "DIME random seed set: "<< AliDimeRndm::GetDimeRandom()->GetSeed() << endl;

  // 3. Setup parameters of DIME MC
  dime->SetEnergyCMS(13000.0);  // sqrt(s)
  dime->SetProcess("pipm");      // Process {"pipm", "pi0", "kpkm", "ks", "rho"}
  dime->SetFormf("orexp");      // Meson-Pomeron form factor {"orexp", "exp", "power"}
  dime->SetFsi("true");         // Exclusive supression {"true", "false"}
  dime->SetYRange(-0.9, 0.9);   // Set rapidity range of mesons
  dime->SetMinPt(0.1);          // Minimum pT of mesons
  
  // 4. Do the initialization, always last
  dime->Initialize();
  
  TClonesArray* particles = new TClonesArray("TParticle", 100);
  
  // Event loop
  Int_t N_events = 10;
  for (Int_t i = 0; i < N_events; ++i) {
    
    dime->GenerateEvent();
    Int_t np = dime->ImportParticles(particles, "All");
    printf("\n Event %d/%d: Imported %3d final state particles \n", i, N_events-1, np);
    
    // Particle loop
    for (Int_t ip = 0; ip < np; ++ip) {
      TParticle* part = (TParticle*) (particles->At(ip));
      part->Print();
    }
  }
}

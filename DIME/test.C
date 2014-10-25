void test()
{
  gSystem->Load("libEG.so");
  gSystem->Load("libdime.so");
  gSystem->Load("libTDime.so");
  
  TDime* dime = new TDime();
  dime->SetEnergyCMS(7000.);
  dime->Initialize();
  
  TClonesArray* particles = new TClonesArray("TParticle", 100);
  
  for (Int_t i = 0; i < 10; i++)
    {    
    dime->GenerateEvent();
    Int_t np = dime->ImportParticles(particles, "All");
    printf("\n Imported %3d particles \n", np);
    for (Int_t ip = 0; ip < np; ip++) {
      TParticle* part = (TParticle*) (particles->At(ip));
      part->Print();
    }
  }

}

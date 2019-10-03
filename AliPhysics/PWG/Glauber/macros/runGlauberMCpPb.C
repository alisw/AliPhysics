void runGlauberMCpPb(Int_t option=0,Int_t N=250000,Bool_t qparts=kFALSE)
{
  //load libraries
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libPWGGlauber");

  //set the random seed from current time
  TTimeStamp time;
  Int_t seed = time.GetSec();
  gRandom->SetSeed(seed);

  Int_t nevents = N; // number of events to simulate 
  // supported systems are e.g. "p", "d", "Si", "Au", "Pb", "U" 
  Option_t *sysA="p"; 
  Option_t *sysB="Pb";
  Double_t signn=70; // inelastic nucleon nucleon cross section
  Double_t mind=0.4;
  const char *fname="glau_ppb_ntuple.root";

  if (qparts) {
    signn/=9;
    mind=0;
  }

  AliGlauberMC mcg(sysA,sysB,signn);
  mcg.SetMinDistance(mind);
  if (qparts) {
    AliGlauberNucleus &na = mcg.GetNucA();
    na.SetN(3*na.GetN());
    AliGlauberNucleus &nb = mcg.GetNucB();
    nb.SetN(3*nb.GetN());
  }
  if (option==1) 
    mcg.SetDoFluc(0.55,78.5,0.82,kTRUE);
  else if (option==2) 
    mcg.SetDoFluc(1.01,72.5,0.74,kTRUE);
  mcg.SetDoPartProduction(kFALSE);
  
  //////////////////
  mcg.SetdNdEtaType(AliGlauberMC::kNBDSV);
  mcg.GetdNdEtaParam()[0] = 2.49;    //npp
  mcg.GetdNdEtaParam()[1] = 1.7;  //ratioSgm2Mu
  mcg.GetdNdEtaParam()[2] = 0.13; //xhard
  //////////////////

  mcg.Run(nevents);

  TNtuple  *nt = mcg.GetNtuple();
  TFile out(fname,"recreate",fname,9);
  if(nt) nt->Write();
  printf("total cross section with a nucleon-nucleon cross section %.4f is %.4f\n\n",signn,mcg.GetTotXSect());
  out.Close();
}

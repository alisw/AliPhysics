void runGlauberMC()
{
  //load libraries
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libPWG2flowTools");

  //set the random seed from current time
  TTimeStamp time;
  Int_t seed = time.GetSec();
  gRandom->SetSeed(seed);

  Int_t nevents = 1000000; // number of events to simulate 
  // supported systems are e.g. "p", "d", "Si", "Au", "Pb", "U" 
  Option_t *sysA="Pb"; 
  Option_t *sysB="Pb";
  Double_t signn=64; // inelastic nucleon nucleon cross section
  //const char *fname="GlauberMC_PbPb_ntuple.root"; // name output file

  // run the code to produce an ntuple:
  //  AliGlauberMC::runAndSaveNucleons(10000,"Pb","Pb",72);
  Double_t mind=0.4;
  //  AliGlauberMC::RunAndSaveNtuple(nevents,sysA,sysB,signn,mind);
  Double_t r=6.62;
  Double_t a=0.546;
  const char *fname="glau_pbpb_ntuple.root";

  AliGlauberMC mcg(sysA,sysB,signn);
  mcg.SetMinDistance(mind);
  mcg.Setr(r);
  mcg.Seta(a);
  mcg.SetDoPartProduction(kTRUE);
  
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

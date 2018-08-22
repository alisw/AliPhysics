/// Author: Eulogio Serradilla
AliGenerator* AddMCGenLightNuclei(const TString& generator="PYTHIA8", 
				  Double_t energyCMS = 7000, 
				  Double_t p0 = 0.100, 
				  Int_t pdg = 1000010020
				  Int_t kProcess=0,
				  Double_t ptHardMin=0,
				  Double_t ptHardMax=1.)
{
// add afterburner to EPOS or PYTHIA8
//
  AliGenCocktail* gener = new AliGenCocktail();

  AliGenerator* genSrc = 0;
  if(generator == "EPOS") {
    AliGenExtExec* genExt = new AliGenExtExec();
    genExt->SetPathScript("$ALICE_PHYSICS/PWG/MCLEGO/CRMC/gen_eposlhc.sh");
    genSrc = dynamic_cast<AliGenerator*>(genExt);
  } else { // Pythia8
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenPythia8.C");
    genSrc = AddMCGenPythia8(energyCMS, kTRUE, 1, kProcess, ptHardMin, ptHardMax);
  }
  
  AliGenLightNuclei* aft = new AliGenLightNuclei();
  aft->SetNucleusPdgCode(pdg);
  aft->SetCoalescenceMomentum(p0);
  gener->AddGenerator(genSrc, "generator", 1);
  gener->AddGenerator(aft, "afterburner", 1);
  return gener;
}

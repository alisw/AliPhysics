AliGenerator* CreatePythia6Gen(Float_t e_cms, Int_t ptHardMin, Int_t ptHardMax, Int_t tune, Int_t cr,Float_t ptWeight) {
    
  if(tune==3)  gSystem->Load("libpythia6_4_28");
   if(tune<3)  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");

  AliGenPythia* genP = new AliGenPythia(1);

  //   vertex position and smearing 
  genP->SetVertexSmear(kPerEvent);

  // structure function
  // use kCTEQ5l for Perugia tunes 
  // except for tunes: Perugia * (325, MRSTLO*), Perugia 6 (326, CTEQ6L),  
  // Perugia 11 M (355, MRST LO**), Perugia 11 C (356, CTEQ6L1) 
  genP->SetStrucFunc(kCTEQ5L); 

  //   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
  if(ptHardMin>0.) {
   
    genP->SetProcess(kPyJets);
    genP->SetPtHard((float)ptHardMin,(float)ptHardMax);
    if(ptWeight>0) genP->SetWeightPower(ptWeight);
  } else
    genP->SetProcess(kPyMb); // Minimum Bias  

  //   Centre of mass energy 
  genP->SetEnergyCMS(e_cms); // in GeV
    
  genP->UseNewMultipleInteractionsScenario(); // for all Pythia versions >= 6.3

  if(tune == 0){ // tune Perugia0
    genP->SetTune(320);
    if(cr==0) genP->SetTune(324);
  }
  if(tune == 1){ // tune Perugia2010
    genP->SetTune(327);
    if(cr==0) genP->SetTune(324);
  } 
  if(tune == 2){ // tune Perugia2011 ('central' Perugia 2011)
    genP->SetTune(350);
    if(cr==0) genP->SetTune(354);
  }
  if(tune == 3){ // tune Perugia2012 ('central' Perugia 2012)
    genP->SetTune(370);
    if(cr==0) genP->SetTune(375);
  }

  genP->Print();
  return genP;
}


AliGenerator* AddMCGenPythia_He3_Eulogio(Float_t e_cms = 2760., Double_t ptHardMin = 0., Double_t ptHardMax = 1., Int_t tune = 2,Int_t cr=1,Float_t ptWeight=0) 
{
  //Add Pythia generator: pt-hard bin or min bias

  AliGenCocktail* gener = new AliGenCocktail();
  
  // pythia
  
  gSystem->Load("liblhapdf");
 
  AliGenerator *genP = NULL;
  genP = CreatePythia6Gen(e_cms, ptHardMin, ptHardMax, tune,cr,ptWeight);

  // nuclei and hyper-nuclei
  
  AliGenLightNuclei* he3 = new AliGenLightNuclei();
  he3->SetNucleusPdgCode(AliGenLightNuclei::kHe3Nucleus);
  he3->SetCoalescenceMomentum(0.127);
  
  AliGenLightNuclei* h3l = new AliGenLightNuclei();
  h3l->SetNucleusPdgCode(AliGenLightNuclei::kHyperTriton);
  h3l->SetCoalescenceMomentum(0.127);
  
  AliGenLightNuclei* d = new AliGenLightNuclei();
  d->SetNucleusPdgCode(AliGenLightNuclei::kDeuteron);
  d->SetCoalescenceMomentum(0.100); // default
  
  AliGenLightNuclei* hd = new AliGenLightNuclei();
  hd->SetNucleusPdgCode(AliGenLightNuclei::kHDibarion);
  hd->SetCoalescenceMomentum(0.100);

  gener->AddGenerator(genP, "PYTHIA", 1);
  gener->AddGenerator(he3, "3He nuclei", 1);
  gener->AddGenerator(h3l, "Hyper-tritons", 1);
  gener->AddGenerator(d, "Deuterons", 1);
  gener->AddGenerator(hd, "H-Dibarions", 1);

  return gener;
}

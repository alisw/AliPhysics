AliGenerator* AddMCGenPythiaLG(Float_t e_cms = 7000., Int_t tune = 2 ,Int_t cr = 0, Int_t npart = 1, Bool_t useNewMultipleInteractionsScenario = kFALSE, Bool_t useFemtoPARJ = kTRUE, Float_t PARJ_1 = 0.1, Float_t PARJ_2 = 0.3, Float_t PARJ_3 = 0.4, Float_t PARJ_4 = 0.05, Float_t PARJ_5 = 0.5, Float_t PARJ_6 = 0.5, Float_t PARJ_7 = 0.5, Float_t PARJ_8 = 0.6, Float_t PARJ_9 = 1.2, Float_t PARJ_10 = 0.6 ) 
{

  //Add Pythia generator: min bias
gSystem->Load("liblhapdf");
  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");

  AliFemtoPARJGenPythia* genP = new AliFemtoPARJGenPythia(npart);


  if(useFemtoPARJ)
    {
      genP->SetUseFemtoPARJ(useFemtoPARJ);
      genP->SetPARJ1(PARJ_1);
      genP->SetPARJ2(PARJ_2);
      genP->SetPARJ3(PARJ_3);
      genP->SetPARJ4(PARJ_4);
      genP->SetPARJ5(PARJ_5);
      genP->SetPARJ6(PARJ_6);
      genP->SetPARJ7(PARJ_7);
      genP->SetPARJ8(PARJ_8);
      genP->SetPARJ9(PARJ_9);
      genP->SetPARJ10(PARJ_10);
    }

  //   vertex position and smearing 
  genP->SetVertexSmear(kPerEvent);  

  // structure function
  // use kCTEQ5l for Perugia tunes 
  // except for tunes: Perugia * (325, MRSTLO*), Perugia 6 (326, CTEQ6L),  
  // Perugia 11 M (355, MRST LO**), Perugia 11 C (356, CTEQ6L1) 
  genP->SetStrucFunc(kCTEQ5L); 

  genP->SetProcess(kPyMb); // Minimum Bias 

  //   Centre of mass energy 
  genP->SetEnergyCMS(e_cms); // in GeV
    
  if(useNewMultipleInteractionsScenario)
    {
      genP->UseNewMultipleInteractionsScenario(); // for all Pythia versions >= 6.3
    }
  
  genP->SetMomentumRange(0,999999);
  genP->SetPhiRange(0., 360.);
  genP->SetThetaRange(0.,180.);
  
  if(tune == 0)
    { // tune Perugia0
      genP->SetTune(320);
      if(cr==0) genP->SetTune(324);
    }
  if(tune == 1)
    { // tune Perugia2010
      genP->SetTune(327);
      if(cr==0) genP->SetTune(324);
    } 
  if(tune == 2)
    { // tune Perugia2011 ('central' Perugia 2011)
      genP->SetTune(350);
      if(cr==0) genP->SetTune(354);
    }
   

  genP->Print();
  return genP;

  
}

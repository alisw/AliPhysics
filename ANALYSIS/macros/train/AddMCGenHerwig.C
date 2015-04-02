AliGenerator* AddMCGenHerwig(Float_t e_beam1=7,Float_t e_beam2=7, Double_t ptHardMin = 0., Double_t ptHardMax = 1.) 
{
  //Add Herwig generator: pt-hard bin or min bias


 
  AliGenerator *genH = NULL;
  genH = CreateHerwigGen(e_beam1,e_beam2, ptHardMin, ptHardMax);
  
  return genH;
}

AliGenerator* CreateHerwigGen(Float_t e_beam1,Float_t e_beam2, Int_t ptHardMin, Int_t ptHardMax) {
    
  gSystem->Load("libTHerwig");
  gSystem->Load("libHERWIG");
  gSystem->Load("liblhapdf");

  AliGenHerwig* genH = new AliGenHerwig(-1);

  //   vertex position and smearing 
  genH->SetVertexSmear(kPerEvent);
  genH->SetProjectile("P");
  genH->SetTarget("P");
   //  Beam momenta 
  genH->SetBeamMomenta(e_beam1,e_beam2); // in GeV
  // structure function
  // this is overwritten by the tune
  genH->SetStrucFunc(kCTEQ5L); 


  if(ptHardMin>0.) {
    genH->SetProcess(1500);
    genH->SetPtHard((float)ptHardMin,(float)ptHardMax);
    //event weight not yet implemented
    //  if(ptWeight>0) genH->SetWeightPower(ptWeight);
  } else
    genH->SetProcess(8000); // Minimum Bias  soft hadron-hadron

      

  genH->Print();
  return genH;
}

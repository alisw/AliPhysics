AliGenerator* AddMCGenHerwig(Float_t e_beam1=3500,Float_t e_beam2=3500, Double_t ptHardMin = 10., Double_t ptHardMax = 100000.,Double_t ptWeight=0) 
{
  //Add Herwig generator: pt-hard bin or min bias


 
  AliGenerator *genH = NULL;
  genH = CreateHerwigGen(e_beam1,e_beam2, ptHardMin, ptHardMax,ptWeight);
  
  return genH;
}

AliGenerator* CreateHerwigGen(Float_t e_beam1,Float_t e_beam2, Int_t ptHardMin, Int_t ptHardMax,Double_t ptWeight) {
    
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
    genH->SetPtHardMin((float)ptHardMin);
    genH->SetPtHardMax((float)ptHardMax);    
    if(ptWeight>0) genH->SetWeightPower(ptWeight);
  } else
    genH->SetProcess(8000); // Minimum Bias  soft hadron-hadron

      

  genH->Print();
  return genH;
}

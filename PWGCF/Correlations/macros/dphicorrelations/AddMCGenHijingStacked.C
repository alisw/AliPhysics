AliGenerator* GetHIJING()
{
  AliGenHijing* gener = new AliGenHijing(-1);
  // centre of mass energy 
  gener->SetEnergyCMS(2760.);
  gener->SetImpactParameterRange(0, 20);  
  // reference frame
  gener->SetReferenceFrame("CMS");
  // projectile
  gener->SetProjectile("A", 208, 82);
  gener->SetTarget    ("A", 208, 82);
  // tell hijing to keep the full parent child chain
  gener->KeepFullEvent();

  // enable shadowing
  gener->SetShadowing(1);
  // Don't track spectators
  gener->SetSpectators(0);
  // kinematic selection
  gener->SetSelectAll(0);   

  gener->SetJetQuenching(0);   
  gener->SetPtHardMin (2.3);
  
  return gener;
}

AliGenerator *AddMCGenHijingStacked(Int_t stackEvents = 1)
{  
  //
  // if stackEvents is > 1, <stackEvents> number of events are put on top of each other
  //
  
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libhijing");       
  gSystem->Load("libTHijing");
  
  if (stackEvents > 1) 
  {
    AliGenCocktail* cocktail = new AliGenCocktail;
    
    for (int i=0; i<stackEvents; i++)
      cocktail->AddGenerator(GetHIJING(), Form("hijing_%d", i), 1.0);
    
    return cocktail;
  }
  else
    return GetHIJING();
}
 
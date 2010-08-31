void SubmitFMDCorrections(const Char_t* filename, Bool_t store, Float_t energy, Int_t trigger, Float_t mag, Int_t collsystem,Bool_t inelGtZero=kFALSE) {
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG2forward");
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  if(energy == 900)
    pars->SetEnergy(AliFMDAnaParameters::k900);
  else if(energy == 7000)
    pars->SetEnergy(AliFMDAnaParameters::k7000);
  else if(energy == 10000)
    pars->SetEnergy(AliFMDAnaParameters::k10000);
  else if(energy == 14000)
    pars->SetEnergy(AliFMDAnaParameters::k14000);
  else if(energy == 5500)
    pars->SetEnergy(AliFMDAnaParameters::k5500);
  
  if(trigger == 0)
    pars->SetTriggerDefinition(AliFMDAnaParameters::kMB1);
  else if(trigger == 1)
    pars->SetTriggerDefinition(AliFMDAnaParameters::kMB2);

  if(mag==0)
    pars->SetMagField(AliFMDAnaParameters::k0G);
  else if(mag==1)
    pars->SetMagField(AliFMDAnaParameters::k5G);
  
  if(collsystem == 0)
    pars->SetCollisionSystem(AliFMDAnaParameters::kPP);
  else if(collsystem == 1)
    pars->SetCollisionSystem(AliFMDAnaParameters::kPbPb);
  
  
  pars->SetInelGtZero(inelGtZero);
  pars->PrintStatus();
  
  std::cout<<"creating background object"<<std::endl;
  AliFMDAnalysisTaskGenerateCorrection t;
  
  t.ReadFromFile(filename,store,0);
  std::cout<<"object created in backgroundFromFile.root "<<std::flush;
  if(store)
    std::cout<<" - and stored!"<<std::endl;
  else
    std::cout<<" - and not stored!"<<std::endl;
}

void SubmitSharingEffCorrection(const Char_t* filename="fmdana.root", Bool_t store, Float_t energy, Int_t trigger, Float_t mag, Int_t collsystem)){

  gSystem->Load("libANALYSIS"); 
  gSystem->Load("libANALYSISalice"); 
  gSystem->Load("libPWGLFforward"); 
  
  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132,"X");
  gStyle->SetLabelFont(132,"Y"); 
  gStyle->SetLabelFont(132,"Z"); 
  gStyle->SetTitleFont(132,"X");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetTitleFont(132,"Z");
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  pars->Init(kTRUE,AliFMDAnaParameters::kBackgroundCorrection);
  if(energy == 900)
    pars->SetEnergy(AliFMDAnaParameters::k900);
  else if(energy == 7000)
    pars->SetEnergy(AliFMDAnaParameters::k7000);
  else if(energy == 10000)
    pars->SetEnergy(AliFMDAnaParameters::k10000);
  else if(energy == 14000)
    pars->SetEnergy(AliFMDAnaParameters::k14000);
  
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
  
  pars->PrintStatus();
  
  std::cout<<"creating sharing efficiency object"<<std::endl;
  AliFMDDndeta t;
  t.SetNbinsToCut(2);
  //    t.SetVtxCut(2);
  t.Init(filename);
  
  t.CreateSharingEfficiency(filename,store);
  if(store)
    std::cout<<" - and stored!"<<std::endl;
  else
    std::cout<<" - and not stored!"<<std::endl;
  
}

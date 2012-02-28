void drawdNdeta(const Char_t* filename="fmdana.root", 
		Int_t sample = 0, 
		Int_t rebin =1,
		Float_t v1 = -10, 
		Float_t v2 = 10, 
		Bool_t realdata=kTRUE, 
		Float_t energy = 900, 
		Int_t magfield=1){
  
  gSystem->Load("libANALYSIS"); 
  gSystem->Load("libANALYSISalice"); 
  gSystem->Load("libPWG0base"); 
  gSystem->Load("libPWG0dep"); 
  gSystem->Load("libPWGLFforward"); 
  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132,"X");
  gStyle->SetLabelFont(132,"Y"); 
  gStyle->SetLabelFont(132,"Z"); 
  gStyle->SetTitleFont(132,"X");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetTitleFont(132,"Z");
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  if(energy == 900)
    pars->SetEnergy(AliFMDAnaParameters::k900);
  else if(energy == 10000)
    pars->SetEnergy(AliFMDAnaParameters::k10000);
  else if(energy == 7000)
    pars->SetEnergy(AliFMDAnaParameters::k7000);
  else {
    std::cout<<"invalid energy - quitting"<<std::endl;
    return; }
  
  if(magfield == 0)
    pars->SetMagField(AliFMDAnaParameters::k0G); 
  else if(magfield == 1)
    pars->SetMagField(AliFMDAnaParameters::k5G); 
  else
    return;
  pars->SetRealData(realdata);
  pars->Init();
  pars->Print();
  AliFMDDndeta t;
  t.SetNbinsToCut(2);
  t.SetVtxCut1(v1);
  t.SetVtxCut2(v2);
  //  t.SetVtxCut(6);
  t.Init(filename);
  if(sample == 0) {
    t.GenerateMult(AliFMDDndeta::kMult);
    t.DrawDndeta(AliFMDDndeta::kMult,rebin,realdata);
  }
  if(sample == 1) {
    t.GenerateMult(AliFMDDndeta::kMultTrVtx);
    t.DrawDndeta(AliFMDDndeta::kMultTrVtx,rebin,realdata);
  }
  if(sample == 2) {
    t.GenerateMult(AliFMDDndeta::kHits);
    t.DrawDndeta(AliFMDDndeta::kHits,rebin);
  }
  if(sample == 3) {
    t.GenerateMult(AliFMDDndeta::kHitsTrVtx);
    t.DrawDndeta(AliFMDDndeta::kHitsTrVtx,rebin);
  }
  if(sample == 4) {
    t.GenerateMult(AliFMDDndeta::kMultNSD);
    t.DrawDndeta(AliFMDDndeta::kMultNSD,rebin,realdata);
  }
  
  
  if(sample == 10) 
    t.CreateSharingEfficiency(filename,kTRUE);
}

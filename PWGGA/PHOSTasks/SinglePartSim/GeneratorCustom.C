AliGenerator* GeneratorCustom(TString opt = ""){
    
  printf("Creating SingleParticle generator: ptMin=%f ptMax=%f \n",ptminConfig,ptmaxConfig) ;  
  AliGenCocktail *cocktail  = GeneratorCocktail("SingleParticle");
  
  AliGenPHOSlib *plib = new AliGenPHOSlib();            
  AliGenParam *genPi0HagPt0 = new AliGenParam(1,new AliGenPHOSlib(),AliGenPHOSlib::kPi0,"");
  genPi0HagPt0->SetPhiRange(240.,330.); ;
  genPi0HagPt0->SetYRange(-0.15,0.15) ;
  genPi0HagPt0->SetPtRange(ptminConfig,ptmaxConfig) ;

  genPi0HagPt0->SetCutOnChild(kTRUE) ;
  genPi0HagPt0->SetChildPhiRange(255.,325.);
  genPi0HagPt0->SetChildYRange(-0.13,0.13) ;
  
//  genPi0HagPt0->SetForceDecay(kGammaEM); // Ensure the decays are photons
  genPi0HagPt0->SetKeepIfOneChildSelected(kTRUE) ;
  genPi0HagPt0->SetKeepParent(kTRUE) ;
  cocktail->AddGenerator(genPi0HagPt0, "Pi0HagPt0", 1);
  
// //Example for eta-meson simulation  
//   AliGenPHOSlib *plib = new AliGenPHOSlib();
//   AliGenParam *genEtaHagPt0 = new AliGenParam(1,new AliGenPHOSlib(),AliGenPHOSlib::kEta,"");
//   genEtaHagPt0->SetPhiRange(240.,330.); ;
//   genEtaHagPt0->SetYRange(-0.15,0.15) ;
//   genEtaHagPt0->SetPtRange(ptminConfig,ptmaxConfig) ;
// 
//   genEtaHagPt0->SetCutOnChild(kTRUE) ;
//   genEtaHagPt0->SetChildPhiRange(255.,325.);
//   genEtaHagPt0->SetChildYRange(-0.13,0.13) ;
// 
//   genEtaHagPt0->SetForceDecay(kGammaEM); // Ensure the decays are photons
//   genEtaHagPt0->SetKeepIfOneChildSelected(kTRUE) ;
//   genEtaHagPt0->SetKeepParent(kTRUE) ;
//   cocktail->AddGenerator(genEtaHagPt0, "EtaHagPt0", 1);
  
  
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  
  
  //return gener;
  return cocktail;

}

AliGenerator * 
GeneratorCustom()
{
  Int_t ipart = atoi(gSystem->Getenv("CONFIG_PARTICLE")) ;
  comment = comment.Append(Form(" | cocktail (%d)",ipart ));

//  if(ipart==221){
//
//    cout << "eta meson is selected. BR eta->gamma+gamma is set to 1." << endl;
//
//    //If Eta, set only eta->2gamma decays
//    double BR_100_gg_pythia = 1.00;
//    int   idc_pythia = -1;
//  
//    AliPythia* py = AliPythia::Instance();
//    Int_t kc =py->Pycomp(221)  ; //Compressed code
//    Int_t idc_pythia = py->GetMDCY(kc, 2); //entry point into branching table
//    py->SetBRAT(idc_pythia,   BR_100_gg_pythia); // new br for gamma gamma
//    py->SetBRAT(idc_pythia+1, 0.); // new br for dalitz
//    py->SetBRAT(idc_pythia+2, 0.); // new br for dalitz
//    py->SetBRAT(idc_pythia+3, 0.); // new br for dalitz
//    py->SetBRAT(idc_pythia+4, 0.); // new br for dalitz
//    py->SetBRAT(idc_pythia+5, 0.); // new br for dalitz
//    py->SetBRAT(idc_pythia+6, 0.); // new br for dalitz
//    py->SetBRAT(idc_pythia+7, 0.); // new br for dalitz
////      KF     KC
////     221    109    eta                                 0    0    0      0.54745     0.00000     0.00000   0.00000E+00    1
////            IDC    
////            590    1    0    0.392300    gamma           gamma                                                           
////            591    1    0    0.321000    pi0             pi0             pi0                                             
////            592    1    0    0.231700    pi+             pi-             pi0                                             
////            593    1    0    0.047800    gamma           pi+             pi-                                             
////            594    1    2    0.004900    gamma           e-              e+                                              
////            595    1    0    0.001300    pi+             pi-             e-              e+                              
////            596    1    0    0.000300    gamma           mu-             mu+                                             
////            597    1    0    0.000700    pi0             gamma           gamma                                           
//
//  }
  
  
  //
  AliGenCocktail *ctl = new AliGenCocktail();
    //=========================//
  // Generator Configuration //
  //=========================//
//  Float_t thmin = (180./TMath::Pi())*2.*atan(exp(-0.15));   // theta min. <---> eta max
//  Float_t thmax = (180./TMath::Pi())*2.*atan(exp(0.15));    // theta max. <---> eta min 
  Double_t ptMin=0.;
  Double_t ptMax=45.;

  if(ipart==111){//pi0
    ptMin=0.5;
    ptMax=45.;
  }
  else if(ipart==221){//eta
    ptMin=1.0;
    ptMax=35.;
  }
  else if(ipart==22){//gamma
    ptMin=0.3;
    ptMax=45.;
  }

  printf("particle:%d will be generated %f < pT < %f GeV/c.\n",ipart,ptMin,ptMax);
  AliGenBox * gPHS = new AliGenBox(1);
  gPHS->SetMomentumRange(ptMin,ptMax);
  gPHS->SetPhiRange(240,330);

//  gPHS->SetThetaRange(thmin,thmax);

  gPHS->SetYRange(-0.15,0.15);
  gPHS->SetPart(ipart) ;

//  AliGenBox * gMod1 = new AliGenBox(1);
//  gMod1->SetMomentumRange(ptMin,ptMax);
//  gMod1->SetPhiRange(240,260);
//  gMod1->SetThetaRange(thmin,thmax);
//  gMod1->SetPart(ipart) ;
//  AliGenBox * gMod2 = new AliGenBox(1);
//  gMod2->SetMomentumRange(ptMin,ptMax);
//  gMod2->SetPhiRange(260,280);
//  gMod2->SetThetaRange(thmin,thmax);
//  gMod2->SetPart(ipart) ;
//  AliGenBox * gMod3 = new AliGenBox(1);
//  gMod3->SetMomentumRange(ptMin,ptMax);
//  gMod3->SetPhiRange(280,300);
//  gMod3->SetThetaRange(thmin,thmax);
//  gMod3->SetPart(ipart) ;
//  AliGenBox * gMod4 = new AliGenBox(1);
//  gMod4->SetMomentumRange(ptMin,ptMax);
//  gMod4->SetPhiRange(300,320);
//  gMod4->SetThetaRange(thmin,thmax);
//  gMod4->SetPart(ipart) ;

  ctl->AddGenerator(gPHS,"PHOSgen",1) ;
//  ctl->AddGenerator(gMod1,"PHOSmod1",1) ;
//  ctl->AddGenerator(gMod2,"PHOSmod2",1) ;
//  ctl->AddGenerator(gMod3,"PHOSmod3",1) ;
//  ctl->AddGenerator(gMod4,"PHOSmod4",1) ;
  return ctl;
}

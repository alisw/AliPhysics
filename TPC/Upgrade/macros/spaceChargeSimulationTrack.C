void CheckCommutator(const char * inputCharge="SC_NeCO2_eps5_50kHz.root"){
  //
  // 1.) create an space charge correction map - using as input 2x2d space charge distribution
  // file 
  /*
    const char * inputCharge="SC_NeCO2_eps5_50kHz.root"
  */
  AliTPCSpaceCharge3D *spaceCharge = new AliTPCSpaceCharge3D;
  spaceCharge->SetSCDataFileName(inputCharge);
  spaceCharge->SetOmegaTauT1T2(0.325,1,1); // Ne CO2
  //spaceCharge->SetOmegaTauT1T2(0.41,1,1.05); // Ar CO2
  spaceCharge->InitSpaceCharge3DDistortion();
  spaceCharge->CreateHistoSCinZR(0.,50,50)->Draw("surf1");
  spaceCharge->CreateHistoDRPhiinZR(0,100,100)->Draw("colz");
  //
  // 2.) Instantiate magnetic field
  //
  ocdb="local://$ALICE_ROOT/OCDB/";
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(0);   
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));   
  //
  // 3.)   OCDB correction  have to be initialized from somewhere
  //
  TFile * f  = TFile::Open("Correction.root");
  TObjArray  * array = AliCDBEntry->GetObject();
  AliTPCComposedCorrection * corrDefault = (AliTPCComposedCorrection *)array->At(0);
  corrDefault->Init();
  //
  // 4.) Create comutators
  //
  TObjArray * array_SC_Default  = new TObjArray(2);
  TObjArray * array_Default_SC  = new TObjArray(2);
  //
  array_SC_Default->AddAt(spaceCharge,0);
  array_SC_Default->AddAt(corrDefault,1);
  //
  array_Default_SC->AddAt(corrDefault,0);
  array_Default_SC->AddAt(spaceCharge,1);
   AliTPCComposedCorrection *corr_SC_Default = new AliTPCComposedCorrection(array_SC_Default,AliTPCComposedCorrection::kQueue);
  AliTPCComposedCorrection *corr_Default_SC = new AliTPCComposedCorrection(array_Default_SC,AliTPCComposedCorrection::kQueue);
  corr_SC_Default->AddVisualCorrection(corr_SC_Default,1);
  corr_Default_SC->AddVisualCorrection(corr_Default_SC,2);
  //
  // 5.) Use TF1, TF2 functionality to visualize comutators
  //
  //  
  TCanvas * canvasNonComute = new TCanvas("SpacechargeDefault"," SpacechargeDefault",1200,700);
  canvasNonComute->Divide(2,1); 
  canvasNonComute->cd(1);

  canvasNonComute->cd(1);gPad->SetRightMargin(0.2);
  TF2 * fdiffXY0 = new TF2("fdiffXY0", "(AliTPCCorrection::GetCorrXYZ(x,y,10,0,1)-AliTPCCorrection::GetCorrXYZ(x,y,10,0,2))*(sqrt(x*x+y*y)>85&&sqrt(x*x+y*y)<245)",-250,250,-250,250);
  fdiffXY0->SetNpx(200);
  fdiffXY0->SetNpy(200);
  fdiffXY0->GetXaxis()->SetTitle("x (cm)");
  fdiffXY0->GetYaxis()->SetTitle("y (cm)");
  fdiffXY0->GetZaxis()->SetTitle("#delta R (cm)");
  fdiffXY0->Draw("colz");
  TPaveText *paveText = new TPaveText(-250,200,0.0,250,"");
  paveText->AddText("[SpaceCharge,Default]");
  paveText->Draw();

  
  canvasNonComute->cd(2);gPad->SetRightMargin(0.2);
  TF2 * fdiffXY1 = new TF2("fdiffXY0", "(AliTPCCorrection::GetCorrXYZ(x,y,10,1,1)-AliTPCCorrection::GetCorrXYZ(x,y,10,1,2))*(sqrt(x*x+y*y)>85&&sqrt(x*x+y*y)<245)",-250,250,-250,250);
  fdiffXY1->SetNpx(200);
  fdiffXY1->SetNpy(200);
  fdiffXY1->GetXaxis()->SetTitle("x (cm)");
  fdiffXY1->GetYaxis()->SetTitle("y (cm)");
  fdiffXY1->GetZaxis()->SetTitle("#delta R#Phi (cm)");
  fdiffXY1->Draw("colz");

}


void spaceChargeSimulationTrack(){
  //
  // 1. Initialzation form space charge maps
  //
  AliTPCSpaceCharge3D *spaceCharge = new AliTPCSpaceCharge3D;
  spaceCharge->SetSCDataFileName("SpaceCharge.root");
  spaceCharge->SetOmegaTauT1T2(0.325,1,1); // Ne CO2
  //spaceCharge->SetOmegaTauT1T2(0.41,1,1.05); // Ar CO2
  spaceCharge->InitSpaceCharge3DDistortion();
  spaceCharge->CreateHistoSCinZR(0.,50,50)->Draw("surf1");
  spaceCharge->CreateHistoDRPhiinZR(0,100,100)->Draw("colz"); 
  ocdb="local://$ALICE_ROOT/OCDB/";
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(0);   
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));   
  //
  // 2. MC generation of the tracks
  //    Distort track and fit distroted track within AliTPCorrection::FitDistortedTrack(*t, refX, dir,  pcstream)
  //
  Double_t etaCuts=1;
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  TF1 fpt("fpt",Form("x*(1+(sqrt(x*x+%f^2)-%f)/([0]*[1]))^(-[0])",mass,mass),0.4,10);
  fpt.SetParameters(7.24,0.120);
  fpt.SetNpx(10000);
  Int_t nTracks=10000; 
  TTreeSRedirector  * pcstream = new TTreeSRedirector("trackDist.root","recreate");
  
  for(Int_t nt=0; nt<nTracks; nt++){
    Double_t phi = gRandom->Uniform(0.0, 2*TMath::Pi());
    Double_t eta = gRandom->Uniform(-etaCuts, etaCuts);
    Double_t pt = 1/(gRandom->Rndm()*5+0.00001); // momentum for f1
    //   printf("phi %lf  eta %lf pt %lf\n",phi,eta,pt);
    Short_t sign=1;
    if(gRandom->Rndm() < 0.5){
      sign =1;
    }else{
      sign=-1;
    }
    
    Double_t theta = 2*TMath::ATan(TMath::Exp(-eta))-TMath::Pi()/2.;
    Double_t pxyz[3];
    pxyz[0]=pt*TMath::Cos(phi);
    pxyz[1]=pt*TMath::Sin(phi);
    pxyz[2]=pt*TMath::Tan(theta);
    Double_t vertex[3]={0,0,0};
    Double_t cv[21]={0};
    AliExternalTrackParam *t= new AliExternalTrackParam(vertex, pxyz, cv, sign);   
    Double_t refX=1.;
    Int_t dir=-1;
    AliExternalTrackParam *td =  spaceCharge->FitDistortedTrack(*t, refX, dir,  pcstream);
  }  
  delete pcstream;
  //
  // Visulalize distortion of of tracks
  //
  TFile * f = TFile::Open("trackDist.root");
  TTree * tree = (TTree*)f->Get("fitDistortSpaceCharge3D");
  tree->SetMarkerStyle(25);
  tree->SetMarkerSize(0.4);
  
  //
  // DCA distortion example
  //
  tree->Draw("track1.GetY()-track0.GetY():track0.GetTgl():abs(track0.fP[4])","track0.fP[4]>0","colz",10000);
  tree->Draw("track1.GetTgl()-track0.GetTgl():track0.GetSigned1Pt():track0.GetTgl()","abs(track0.GetTgl()-0.5)<0.4","colz",10000);

  /*
    problems:
    
   */
}

/*
 */

void DrawDistortions(){
  //
  //
  //
  TFile * f = TFile::Open("trackDist.root");
  TTree * tree = (TTree*)f->Get("fitDistortSpaceCharge3D");
  tree->SetMarkerStyle(25);
  tree->SetMarkerSize(0.4);  
  tree->Draw("track1.GetY()-track0.GetY():abs(track0.GetSigned1Pt()):track0.GetTgl()","track0.GetSigned1Pt()>0&&track0.GetTgl()>0","colz",10000);
tree->Draw("track1.GetY()-track0.GetY():track0.GetSigned1Pt():track0.GetTgl()","track0.GetTgl()>0","colz",10000);

//
 tree->SetMarkerSize(0.25);
 tree->SetMarkerColor(1);
 tree->Draw("point0.fY:point0.fX","point0.fX!=0","",20,0);
 tree->SetMarkerColor(2);
 tree->Draw("point1.fY:point1.fX","point1.fX!=0","same",20,0);

}

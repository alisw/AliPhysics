//
// macro to create a residual OCDB 
//

/*
//
//
//
.x $ALICE_ROOT/TPC/Upgrade/macros/ConfigOCDB.C(1)
.L $ALICE_ROOT/TPC/Upgrade/macros/makeResidualSpaceChargeOCDB.C

Example usage:

ln -sf /hera/alice/wiechula/Upgrade/LUTs_fluctuation/dir1/SpaceChargeFluc10_1.lookup.root current.root
ln -sf /hera/alice/wiechula/Upgrade/LUTs_fluctuation/average/SpaceChargeFluc0_1.lookup.root  mean.root
ln -sf $ALICE_ROOT/OCDB/TPC/Calib/Correction/Run0_999999999_v0_s2.root ocdb.root
aliroot -b -q $ALICE_ROOT/TPC/Upgrade/macros/makeResidualSpaceChargeOCDB.C

*/

void makeResidualSpaceChargeOCDBLookup(const char *ocdbInName="ocdb.root",const char *scCurrentName="current.root", const char *scMeanName="mean.root"){
  //
  // Macro to create a clone original  OCDB entry with space point distortion and add there space point
  // distortions cause by resifual sapce charge
  // Output is stored by default in the current directory.
  // 

  // Parameters to Specify:
  //    ocdbInName      : path to the original OCDB entry
  //    scCurrentName   : path to the fluctuetaed distortion map
  //    scMeanName      : path to the mean distortion map
  //    
  
  /*
    scCurrentName="/hera/alice/wiechula/Upgrade/LUTs_fluctuation/dir1/SpaceChargeFluc10_1.lookup.root";
    scMeanName="/hera/alice/wiechula/Upgrade/LUTs_fluctuation/average/SpaceChargeFluc0_1.lookup.root";
    ocdbInName="$ALICE_ROOT/OCDB/TPC/Calib/Correction/Run0_999999999_v0_s2.root"
  */
  TFile * finCurrent = TFile::Open(scCurrentName);
  TFile * finMean = TFile::Open(scMeanName);
  //
  AliTPCCorrectionLookupTable *spaceChargeCurrent=  (AliTPCCorrectionLookupTable *)finCurrent->Get("map");
  AliTPCCorrectionLookupTable *spaceChargeMean   =  (AliTPCCorrectionLookupTable *)finMean->Get("map");
  AliTPCInverseCorrection * spaceChargeInverseMean = new AliTPCInverseCorrection(spaceChargeMean);
  //
  TObjArray * arraySC = new TObjArray(2);
  arraySC->AddAt(spaceChargeCurrent,0);
  arraySC->AddAt(spaceChargeMean,1);
  AliTPCComposedCorrection *corrSC = new AliTPCComposedCorrection(arraySC,AliTPCComposedCorrection::kParallel);
  AliTPCComposedCorrection *corrSCW = new AliTPCComposedCorrection(arraySC,AliTPCComposedCorrection::kParallel);
  AliTPCComposedCorrection::AddVisualCorrection(spaceChargeCurrent,1);
  AliTPCComposedCorrection::AddVisualCorrection(spaceChargeMean,2);
  AliTPCComposedCorrection::AddVisualCorrection(spaceChargeInverseMean,3);
  AliTPCComposedCorrection::AddVisualCorrection(corrSC,4);
  AliTPCComposedCorrection::AddVisualCorrection(corrSCW,5);
  //
  //
  TLinearFitter fitterR(2,"pol1");
  for (Int_t ipoint =0; ipoint<10000; ipoint++){
    Double_t z0 = -250+gRandom->Rndm()*500;
    Double_t r0   = 85+gRandom->Rndm()*(245-85.);
    Double_t phi0 = gRandom->Rndm()*TMath::TwoPi();
    // some problematic parts to be skipped  - investigated later
    Double_t xvalue[2]={0};
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,1))>50) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,2))>50) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,1))>20) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,2))>20) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,2,1))>50) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,2,2))>50) continue;
    //
    //
    xvalue[0]=AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,2);
    fitterR.AddPoint(xvalue,AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,1));
  }
  fitterR->Eval();
  TVectorD weights(2);
  weights[0]=1;
  weights[1]=-1;
  corrSC->SetWeights( &weights);
  weights[1]=-fitterR->GetParameter(1);
  corrSCW->SetWeights( &weights);
  //
  //
  //
  TF1 *fSC = new TF1("fSC","AliTPCCorrection::GetCorrXYZ(x,0,10,0,4)",85,245);  
  TF1 *fSCW = new TF1("fSCW","AliTPCCorrection::GetCorrXYZ(x,0,10,0,5)",85,245);
  fSC->SetLineColor(2);
  fSC->SetLineColor(4);
  fSC->Draw();
  fSCW->Draw("same");
  gPad->SaveAs("residualMap.pdf");
  
  //
  //
  //
  TFile *fileOCDBIn=TFile::Open(ocdbInName);
  AliCDBEntry   *entry =  ( AliCDBEntry   *)fileOCDBIn->Get("AliCDBEntry");
  TObjArray * array = (TObjArray*)entry->GetObject();
  for (Int_t i=0;  i<3;  i++){
    AliTPCComposedCorrection *corr = ( AliTPCComposedCorrection*)array->At(i);
    if (corr){
      TObjArray* corrArray = corr->GetCorrections();
      corrArray->AddLast(corrSCW);
    }    
  }
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment("Standard+fluctuation");
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/Correction", 0, 9999999999999999);
  //
  TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  pocdbStorage = AliCDBManager::Instance()->GetStorage(localStorage.Data());
  pocdbStorage->Put(array, (*id1), metaData);


}





void makeResidualSpaceChargeOCDB(const char *ocdbInName="ocdb.root",const char *scCurrentName="current.root", const char *scMeanName="mean.root"){
  //
  // Macro to create a clone original  OCDB entry with space point distortion and add there space point
  // distortions cause by resifual sapce charge
  // Output is stored by default in the current directory.
  // 

  // Parameters to Specify:
  //    ocdbInName      : path to the original OCDB entry
  //    scCurrentName   : path to the fluctuetaed distortion map
  //    scMeanName      : path to the mean distortion map
  //    
  
  /*
    scCurrentName="/hera/alice/miranov/SpaceCharge/Fluctuations/PbPbWithGain/dirmergeAll/dEpsilon20/dir0/SpaceChargeFluc10_1.root";
    scMeanName="/hera/alice/miranov/SpaceCharge/Fluctuations/PbPbWithGain/dirmergeAll/dEpsilon20/dir0//SpaceChargeFluc0_1.root";
    ocdbInName="$ALICE_ROOT/OCDB/TPC/Calib/Correction/Run0_999999999_v0_s2.root"
  */
  TFile * finCurrent = TFile::Open(scCurrentName);
  TFile * finMean = TFile::Open(scMeanName);
  //
  AliTPCCorrectionLookupTable *spaceChargeCurrent=  (AliTPCCorrectionLookupTable *)finCurrent->Get("DistRef");
  AliTPCCorrectionLookupTable *spaceChargeMean   =  (AliTPCCorrectionLookupTable *)finMean->Get("DistRef");
  AliTPCInverseCorrection * spaceChargeInverseMean = new AliTPCInverseCorrection(spaceChargeMean);
  //
  TObjArray * arraySC = new TObjArray(2);
  arraySC->AddAt(spaceChargeCurrent,0);
  arraySC->AddAt(spaceChargeMean,1);
  AliTPCComposedCorrection *corrSC = new AliTPCComposedCorrection(arraySC,AliTPCComposedCorrection::kParallel);
  AliTPCComposedCorrection *corrSCW = new AliTPCComposedCorrection(arraySC,AliTPCComposedCorrection::kParallel);
  AliTPCComposedCorrection::AddVisualCorrection(spaceChargeCurrent,1);
  AliTPCComposedCorrection::AddVisualCorrection(spaceChargeMean,2);
  AliTPCComposedCorrection::AddVisualCorrection(spaceChargeInverseMean,3);
  AliTPCComposedCorrection::AddVisualCorrection(corrSC,4);
  AliTPCComposedCorrection::AddVisualCorrection(corrSCW,5);
  //
  //
  TLinearFitter fitterR(2,"pol1");
  for (Int_t ipoint =0; ipoint<10000; ipoint++){
    Double_t z0 = -250+gRandom->Rndm()*500;
    Double_t r0   = 85+gRandom->Rndm()*(245-85.);
    Double_t phi0 = gRandom->Rndm()*TMath::TwoPi();
    // some problematic parts to be skipped  - investigated later
    Double_t xvalue[2]={0};
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,1))>50) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,0,2))>50) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,1))>20) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,2))>20) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,2,1))>50) continue;
    if (TMath::Abs(AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,2,2))>50) continue;
    //
    //
    xvalue[0]=AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,2);
    fitterR.AddPoint(xvalue,AliTPCCorrection::GetCorrXYZ(r0*TMath::Cos(phi0), r0*TMath::Sin(phi0), z0,1,1));
  }
  fitterR->Eval();
  TVectorD weights(2);
  weights[0]=1;
  weights[1]=-1;
  corrSC->SetWeights( &weights);
  weights[1]=-fitterR->GetParameter(1);
  corrSCW->SetWeights( &weights);
  //
  //
  //
  TF1 *fSC = new TF1("fSC","AliTPCCorrection::GetCorrXYZ(x,0,10,1,4)",85,245);  
  TF1 *fSCW = new TF1("fSCW","AliTPCCorrection::GetCorrXYZ(x,0,10,1,5)",85,245);
  fSC->SetLineColor(2);
  fSC->SetLineColor(4);
  fSC->GetXaxis()->SetTitle(" R (cm)");
  fSC->GetYaxis()->SetTitle(" #Delta_{R#phi} (cm)");
  TCanvas * canvasDist = new TCanvas("canvasDist","canvasDist",600,500);
  TF1 *f1sc = new TF1("f1sc","[0]+[1]*x+[2]*x**2",85,245);
  TF1 *f1scw = new TF1("f1scw","[0]+[1]*x+[2]*x**2",85,245);

  {
    fSC->Draw();
    fSCW->Draw();
    fSC->SetMinimum(-1);
    fSC->SetMaximum(2.);
    fSC->GetHistogram()->SetLineColor(2);
    fSCW->GetHistogram()->SetLineColor(4);
    //
    f1sc->SetLineColor(2);
    f1scw->SetLineColor(4);
    f1sc->SetLineWidth(0.5);
    f1scw->SetLineWidth(0.5);
    //
    fSC->GetHistogram()->Fit(f1sc);
    fSCW->GetHistogram()->Fit(f1scw);
    //
    fSC->GetHistogram()->Draw("");
    fSCW->GetHistogram()->Draw("same");
    //
    TF1 *f1scp=new TF1(*f1sc);
    TF1 *f1scm=new TF1(*f1sc);
    f1scp->SetLineStyle(2); f1scm->SetLineStyle(2);
    f1scp->SetParameter(0, f1sc->GetParameter(0)+ 0.3);
    f1scm->SetParameter(0, f1sc->GetParameter(0)- 0.3);
    f1scp->Draw("same");
    f1scm->Draw("same");
    //
    TF1 *f1scwp=new TF1(*f1scw);
    TF1 *f1scwm=new TF1(*f1scw);
    f1scwp->SetLineStyle(2); f1scwm->SetLineStyle(2);
    f1scwp->SetParameter(0, f1scw->GetParameter(0)+ 0.3);
    f1scwm->SetParameter(0, f1scw->GetParameter(0)- 0.3);
    f1scwp->Draw("same");
    f1scwm->Draw("same");
  }
  TLegend * legend = new TLegend(0.3,0.5,0.89,0.89,"Residual #Delta_{R#phi} distortions and helix fit");
  legend->SetBorderSize(0);
  legend->AddEntry(fSC->GetHistogram(),"Stage 0: #Delta_{R#phic}-#Delta_{R#phimean}");
  legend->AddEntry(fSCW->GetHistogram(),"Stage 1: #Delta_{R#phic}-k_{fit}#Delta_{R#phimean}");
  legend->Draw();
  gPad->SaveAs("residualMapTrackFit.pdf");
  
  //
  //
  //
  TFile *fileOCDBIn=TFile::Open(ocdbInName);
  AliCDBEntry   *entry =  ( AliCDBEntry   *)fileOCDBIn->Get("AliCDBEntry");
  TObjArray * array = (TObjArray*)entry->GetObject();
  for (Int_t i=0;  i<3;  i++){
    AliTPCComposedCorrection *corr = ( AliTPCComposedCorrection*)array->At(i);
    if (corr){
      TObjArray* corrArray = corr->GetCorrections();
      corrArray->AddLast(corrSCW);
    }    
  }
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  metaData->SetComment("Standard+fluctuation");
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/Correction", 0, 9999999999999999);
  //
  TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  pocdbStorage = AliCDBManager::Instance()->GetStorage(localStorage.Data());
  pocdbStorage->Put(array, (*id1), metaData);


}



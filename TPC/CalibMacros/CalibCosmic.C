/*
  .x ~/NimStyle.C
  .x ~/UliStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  .L $ALICE_ROOT/TPC/CalibMacros/CalibCosmic.C
  // init
  Init();
  SetDefaultCut();
  
*/  

AliTPCcalibCosmic * cosmic =0;
TObjArray fitArr;

void Init(){
  //
  //
  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  cosmic = ( AliTPCcalibCosmic *)array->FindObject("cosmicTPC");
  TString axisName[9];
  axisName[0]  ="#Delta"; axisName[1]  ="N_{cl}";
  axisName[2]  ="DCA_{r}(cm)";
  axisName[3]  ="z (cm)"; axisName[4]  ="sin(#phi)";
  axisName[5]  ="tan(#theta)"; axisName[6]  ="1/p_{t} (1/GeV)";
  axisName[7]  ="p_{t} (GeV)"; axisName[8]  ="alpha";

  {
  for (Int_t ivar=0;ivar<6;ivar++){
    for (Int_t ivar2=0;ivar2<9;ivar2++){      
      cosmic->fHistoDelta[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      cosmic->fHistoDelta[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
      cosmic->fHistoPull[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      cosmic->fHistoPull[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    }
  }
  }
}

void SetRangeAll(Int_t axis, Float_t xmin, Float_t xmax){

  for (Int_t i=0;i<6;i++){
    //
    cosmic->fHistoDelta[i]->GetAxis(axis)->SetRangeUser(xmin,xmax);
    cosmic->fHistoPull[i]->GetAxis(axis)->SetRangeUser(xmin,xmax);
  }
}


void SetDefaultCut(){
  for (Int_t i=0;i<6;i++){
    //
    cosmic->fHistoDelta[i]->GetAxis(1)->SetRangeUser(120,200);
    cosmic->fHistoPull[i]->GetAxis(1)->SetRangeUser(120,200);
    cosmic->fHistoDelta[i]->GetAxis(3)->SetRangeUser(-0,250);
    cosmic->fHistoPull[i]->GetAxis(3)->SetRangeUser(-0,250);
    cosmic->fHistoDelta[i]->GetAxis(2)->SetRangeUser(0,60);
    cosmic->fHistoPull[i]->GetAxis(2)->SetRangeUser(0,60);
  }
}

TH2 * GetDelta2D(Int_t type, Int_t var){
  TH2 * his = cosmic->fHistoDelta[type]->Projection(0,var);
  his->SetXTitle(cosmic->fHistoDelta[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmic->fHistoDelta[type]->GetAxis(0)->GetName());
  return his;
}


TH1* GetResol2DSigma(Int_t type, Int_t var){
  
  TH2 * his = cosmic->fHistoDelta[type]->Projection(0,var);
  his->SetXTitle(cosmic->fHistoDelta[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmic->fHistoDelta[type]->GetAxis(0)->GetName());
  his->FitSlicesY(0,0,-1,0,"QNR",&fitArr);
  TH1 * hres = (TH1*)(fitArr.At(2)->Clone());
  return hres;
}


TH1 * GetDelta(Int_t type){
  TH1 * his = cosmic->fHistoDelta[type]->Projection(0);
  his->SetXTitle(cosmic->fHistoDelta[type]->GetAxis(0)->GetName());
  return his;
}

TH2 * GetPull2D(Int_t type, Int_t var){
  TH2 * his = cosmic->fHistoPull[type]->Projection(0,var);
  his->SetXTitle(cosmic->fHistoPull[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmic->fHistoPull[type]->GetAxis(0)->GetName());
  return his;
}

TH1* GetPull2DSigma(Int_t type, Int_t var){
  
  TH2 * his = cosmic->fHistoPull[type]->Projection(0,var);
  his->SetXTitle(cosmic->fHistoPull[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmic->fHistoPull[type]->GetAxis(0)->GetName());
  his->FitSlicesY(0,0,-1,0,"QNR",&fitArr);
  TH1 * hres = (TH1*)(fitArr.At(2)->Clone());
  return hres;
}



TH1 * GetPull(Int_t type){
  TH1 * his = cosmic->fHistoPull[type]->Projection(0);
  his->SetXTitle(cosmic->fHistoPull[type]->GetAxis(0)->GetName());
  return his;
}


void DrawResoldEdx(AliTPCcalibCosmic *cosmic){
  //
  //
  //
  Int_t kmicolors[10]={1,2,3,6,7,8,9,10,11,12};
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};
  TH2 *htemp;
  TObjArray arr;
  TH1 * hResolMax[4];
  TH1 * hResolTot[4];
  //  
  for (Int_t ipad=0;ipad<4;ipad++){
    cosmic->fHistodEdxTot[ipad]->GetAxis(4)->SetRangeUser(-0.6,0.6);
    cosmic->fHistodEdxMax[ipad]->GetAxis(4)->SetRangeUser(-0.6,0.6);
  }
  cosmic->fHistodEdxTot[0]->GetAxis(1)->SetRangeUser(30,62);
  cosmic->fHistodEdxTot[1]->GetAxis(1)->SetRangeUser(30,62);
  cosmic->fHistodEdxTot[2]->GetAxis(1)->SetRangeUser(10,35);
  cosmic->fHistodEdxTot[3]->GetAxis(1)->SetRangeUser(10,150);
  cosmic->fHistodEdxMax[0]->GetAxis(1)->SetRangeUser(30,62);
  cosmic->fHistodEdxMax[1]->GetAxis(1)->SetRangeUser(30,62);
  cosmic->fHistodEdxMax[2]->GetAxis(1)->SetRangeUser(10,35);
  cosmic->fHistodEdxMax[3]->GetAxis(1)->SetRangeUser(10,150);
  //

  for (Int_t ipad=0;ipad<4;ipad++){
    htemp = cosmic->fHistodEdxTot[ipad]->Projection(0,1);
    htemp->FitSlicesY(0,0,-1,0,"QNR",&arr);
    hResolTot[ipad] = (TH1*)(arr.At(2)->Clone());
    delete htemp;
    arr.SetOwner(kTRUE);
    arr.Delete();
    hResolTot[ipad]->Scale(1./TMath::Sqrt(2.));
    //
    htemp = cosmic->fHistodEdxMax[ipad]->Projection(0,1);
    htemp->FitSlicesY(0,0,-1,0,"QNR",&arr);
    hResolMax[ipad] = (TH1*)(arr.At(2)->Clone());
    delete htemp;
    arr.SetOwner(kTRUE);
    arr.Delete();
    hResolMax[ipad]->Scale(1./TMath::Sqrt(2.));    
  }
  hResolTot[3]->GetXaxis()->SetRangeUser(0,160);
  hResolTot[3]->SetXTitle("N_{cl}");
  hResolTot[3]->SetYTitle("#sigma(dEdx/dEdx_{d})/#sqrt{2.}");
  hResolTot[3]->SetTitle("Relative dEdx resolution");
  for (Int_t ipad=3;ipad>=0;ipad--){
    hResolTot[ipad]->SetMaximum(0.1);
    hResolTot[ipad]->SetMinimum(0.);
    hResolTot[ipad]->SetMarkerColor(kmicolors[ipad]+0);
    hResolTot[ipad]->SetMarkerStyle(kmimarkers[ipad]+1);
    if (ipad==3)    hResolTot[ipad]->Draw();
    hResolTot[ipad]->Draw("same");
    //
    hResolMax[ipad]->SetMaximum(0.1);
    hResolMax[ipad]->SetMinimum(0.);
    hResolMax[ipad]->SetMarkerColor(kmicolors[ipad]+0);
    hResolMax[ipad]->SetMarkerStyle(kmimarkers[ipad]+4);
    hResolMax[ipad]->Draw("same");
  }
  
}

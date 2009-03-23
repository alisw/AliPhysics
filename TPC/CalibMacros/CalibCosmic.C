/*
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
  axisName[3]  ="z (cm)"; axisName[4]  ="sin(#phi)";
  axisName[5]  ="tan($theta)"; axisName[6]  ="1/p_{t} (1/GeV)";
  axisName[7]  ="p_{t} (1/GeV)"; axisName[8]  ="alpha";

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

void SetDefaultCut(){
  for (Int_t i=0;i<6;i++){
    //
    cosmic->fHistoDelta[i]->GetAxis(1)->SetRangeUser(130,200);
    cosmic->fHistoPull[i]->GetAxis(1)->SetRangeUser(130,200);
    cosmic->fHistoDelta[i]->GetAxis(3)->SetRangeUser(0,250);
    cosmic->fHistoPull[i]->GetAxis(3)->SetRangeUser(0,250);

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

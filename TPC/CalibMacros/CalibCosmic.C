/*
Draw result of perfomance test:

aliroot -b -q  $ALICE_ROOT/TPC/scripts/loadTPCcalib.C $ALICE_ROOT/TPC/CalibMacros/CalibCosmic.C

  //gROOT->Macro("~/NimStyle.C");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  .L $ALICE_ROOT/TPC/CalibMacros/CalibCosmic.C
  // init
  Init();
  SetDefaultCut();  // check defualt cut 
  //
  MakeDefaultPlots();

  gROOT->Macro("~/NimStyle.C");
  TFile f("cosmicPlots.root");
  TBrowser b
  b.Add(CosmicPlots,"CosmicPlot");

*/  

AliTPCcalibCosmic * cosmic =0;
TObjArray fitArr;
Int_t colors[3]={1,3,4};

void CalibCosmic(){
  // init
  Init();
  SetDefaultCut(); 
  //
  MakeDefaultPlots();
}

void Init(){
  //
  //
  TFile fcalib("CalibObjectsTrain1.root");
  cosmic = ( AliTPCcalibCosmic *)fcalib.Get("cosmicTPC");
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
    //cut on number of clusters
    cosmic->fHistoDelta[i]->GetAxis(1)->SetRangeUser(120,200);
    cosmic->fHistoPull[i]->GetAxis(1)->SetRangeUser(120,200);
    //cut on DCA r
    cosmic->fHistoDelta[i]->GetAxis(2)->SetRangeUser(0,15);
    cosmic->fHistoPull[i]->GetAxis(2)->SetRangeUser(0,15);
    //cut on z at
    cosmic->fHistoDelta[i]->GetAxis(3)->SetRangeUser(20,200);
    cosmic->fHistoPull[i]->GetAxis(3)->SetRangeUser(20,200);
  }
  cosmic->fHistoDelta[0]->GetAxis(0)->SetRangeUser(-1,1);
  cosmic->fHistoDelta[1]->GetAxis(0)->SetRangeUser(-1,1);
}

TH2 * GetDelta2D(Int_t type, Int_t var){
  TH2 * his = cosmic->fHistoDelta[type]->Projection(0,var);
  his->SetXTitle(cosmic->fHistoDelta[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmic->fHistoDelta[type]->GetAxis(0)->GetName());
  return his;
}


TH1* GetFit2D(Int_t type, Int_t var, Bool_t resol){
  
  TH2 * his = cosmic->fHistoDelta[type]->Projection(0,var);
  his->SetXTitle(cosmic->fHistoDelta[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmic->fHistoDelta[type]->GetAxis(0)->GetName());
  his->FitSlicesY(0,0,-1,0,"QNR",&fitArr);
  TH1 * hres = 0;
  if (resol) hres = (TH1*)(fitArr.At(2)->Clone());
  if (!resol) hres = (TH1*)(fitArr.At(1)->Clone());
  hres->SetMarkerStyle(20);
  hres->SetMarkerColor(2);
  hres->GetYaxis()->SetTitleOffset(1.8);
  hres->GetYaxis()->SetDecimals(kTRUE);
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

void DrawStat(Int_t coord, TObjArray *array=0){
  //
  //
  //
  TCanvas *cStat = new TCanvas(Form("Cosmic stat%d",coord), Form("CosmicStat%d",coord),800,600);
  Float_t mx0=0.2, mx1=0.05, my0=0.15, my1=0.1;
  //pad->SetMargin(mx0,mx1,my0,my1);
  cStat->Divide(3,3);
  for (Int_t i=0; i<8; i++){
    cStat->cd(i+1);
    cosmic->fHistoDelta[0]->Projection(i)->Draw();
  }
  if (array) array->AddLast(cStat);
}

void SetStylePad(TVirtualPad *pad){
  Float_t mx0=0.2, mx1=0.05, my0=0.15, my1=0.1;
  pad->SetMargin(mx0,mx1,my0,my1);
  pad->SetTicks(1,1);
  pad->SetGrid(1,1); 
}

void MakePlotPt(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCPtResol","TPCPtResol",600,500);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==0) cosmic->fHistoDelta[5]->GetAxis(6)->SetRangeUser(0,2);
    if (i==1) cosmic->fHistoDelta[5]->GetAxis(6)->SetRangeUser(-2,0);
    if (i==2) cosmic->fHistoDelta[5]->GetAxis(6)->SetRangeUser(-2,2);  
    hisRes  = GetFit2D(5,7,kTRUE);
    hisMean = GetFit2D(5,7,kFALSE);
    hisRes->SetMarkerStyle(20);
    hisMean->SetMarkerStyle(20);
    hisRes->SetMarkerColor(colors[i]);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->Scale(100);
    hisMean->Scale(100);
    hisRes->SetMaximum(30);
    hisRes->SetMinimum(0);
    hisMean->SetMaximum(20);
    hisMean->SetMinimum(-20);
    hisRes->SetName("Pt resol");
    hisRes->SetName("p_{t} resolution");
    hisRes->SetYTitle("#sigma_{p_{t}}/p_{t} (%)");
    hisMean->SetYTitle("#Delta_{p_{t}}/p_{t} (%)");
    hisRes->GetXaxis()->SetRangeUser(0,10);
    hisMean->GetXaxis()->SetRangeUser(0,10);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
  }
  if (array) array->AddLast(cptRes);
}


void MakePlotP4(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCP4Resol","TPCP4Resol",600,500);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  //
  TH1 * hisRes  =0;
  TH1 * hisMeanP=0;
  
  for (Int_t i=0; i<3; i++){
    if (i==0) cosmic->fHistoDelta[4]->GetAxis(6)->SetRangeUser(0,2);
    if (i==1) cosmic->fHistoDelta[4]->GetAxis(6)->SetRangeUser(-2,0);
    if (i==2) cosmic->fHistoDelta[4]->GetAxis(6)->SetRangeUser(-2,2);
    hisRes  = GetFit2D(4,7,kTRUE);
    hisMean = GetFit2D(4,7,kFALSE);
    hisRes->SetMarkerStyle(20+i);
    hisMean->SetMarkerStyle(20+i);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->SetMarkerColor(colors[i]);
    hisRes->SetMaximum(0.04);
    hisRes->SetMinimum(-0.0);
    hisMean->SetMaximum(0.02);
    hisMean->SetMinimum(-0.02);
    hisRes->SetName("C resol");
    hisRes->SetName("C resolution");
    hisRes->SetYTitle("#sigma_{C} (1/GeV)");
    hisMean->SetYTitle("#Delta_{C} (1/GeV)");
    hisMean->SetXTitle("p_{t} (GeV)");
    hisRes->SetXTitle("p_{t} (GeV)");
    hisRes->GetXaxis()->SetRangeUser(0,10);
    hisMean->GetXaxis()->SetRangeUser(0,10);      
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
  }
  if (array) array->AddLast(cptRes);
}






void MakePlotPosY(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCPosResolY","TPCPosResolY",600,500);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==0) cosmic->fHistoDelta[0]->GetAxis(6)->SetRangeUser(0,2);
    if (i==1) cosmic->fHistoDelta[0]->GetAxis(6)->SetRangeUser(-2,0);
    if (i==2) cosmic->fHistoDelta[0]->GetAxis(6)->SetRangeUser(-2,2);
    hisRes  = GetFit2D(0,7,kTRUE);
    hisMean = GetFit2D(0,7,kFALSE);
    hisRes->SetMarkerStyle(20+i);
    hisMean->SetMarkerStyle(20+i);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->SetMarkerColor(colors[i]);
    
    //
    hisRes->SetMaximum(0.4);
    hisRes->SetMinimum(0.0);
    hisMean->SetMaximum(0.4);
    hisMean->SetMinimum(-0.4);
    hisRes->SetName("Y resol");
    hisRes->SetName("Y resolution");
    hisRes->SetYTitle("#sigma_{y} (cm)");
    hisMean->SetYTitle("#Delta_{y} (cm)");
    hisRes->GetXaxis()->SetRangeUser(0,10);
    hisMean->GetXaxis()->SetRangeUser(0,10);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
  }
  if (array) array->AddLast(cptRes);
}

void MakePlotPosZ(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCPosResolZ","TPCPosResolZ",600,500);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==0) cosmic->fHistoDelta[1]->GetAxis(6)->SetRangeUser(0,2);
    if (i==1) cosmic->fHistoDelta[1]->GetAxis(6)->SetRangeUser(-2,0);
    if (i==2) cosmic->fHistoDelta[1]->GetAxis(6)->SetRangeUser(-2,2);
    
    hisRes  = GetFit2D(1,7,kTRUE);
    hisMean = GetFit2D(1,7,kFALSE);
    hisRes->SetMaximum(0.4);
    hisRes->SetMinimum(0.0);
    hisMean->SetMaximum(0.2);
    hisMean->SetMinimum(-0.2); 
    hisRes->SetMarkerStyle(20);
    hisMean->SetMarkerStyle(20);
    hisRes->SetMarkerColor(colors[i]);
    hisMean->SetMarkerColor(colors[i]);
 
    hisRes->SetName("Z resol");
    hisRes->SetName("Z resolution");
    hisRes->SetYTitle("#sigma_{z} (cm)");
    hisMean->SetYTitle("#Delta_{z} (cm)");
    hisRes->GetXaxis()->SetRangeUser(0,10);
    hisMean->GetXaxis()->SetRangeUser(0,10);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw();
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw();
  }
  if (array) array->AddLast(cptRes);
}


void  MakeDefaultPlots(){
  //
  //
  //
  TObjArray *picArray = new TObjArray();
  DrawStat(0,picArray);
  MakePlotPt(picArray);
  MakePlotPosY(picArray);
  MakePlotPosZ(picArray);
  MakePlotP4(picArray);
  TFile f("cosmicPlots.root","recreate");
  picArray->Write("CosmicPlots",TObject::kSingleKey);
  f.Close();
  TPostScript *ps=new TPostScript("cosmicPerformance.ps", 112);
  for (Int_t ipad=0;ipad<picArray->GetEntries();ipad++){
    TCanvas *c =dynamic_cast<TCanvas*> picArray->At(ipad);
    if (c) {
      c->SaveAs(Form("pic/cosmic/%s.gif",c->GetName()));
      c->SaveAs(Form("pic/cosmic/%s.eps",c->GetName()));
      ps->NewPage();
      c->Draw();
      c->Update();
    }
  } 
  ps->Close();
  delete ps;
}

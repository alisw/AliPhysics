void ProcessOutput(const char* filename, TCanvas* cLayers, TCanvas* cLayers_prod, Int_t color, Option_t* option= "");
TPaveStats* SetStatPad(TH1* hst,float x1,float x2,float y1,float y2,Int_t stl=-1,Int_t col=-1);
void SetHStyle(TH1* hst,int col=kRed,int mark=20,float mrsize=0.7);
TPaveStats* GetStatPad(TH1* hst);

void MakePlotsTrackRefs(){

  TCanvas* cLayers = new TCanvas("cLayers", "cLayers", 50, 50, 1150, 750);
  cLayers->Divide(3, 2);
  TCanvas* cLayers_prod = new TCanvas("cLayers_prod", "cLayers_prod", 50, 50, 1150, 750);
  cLayers_prod->Divide(3, 2);
  Float_t slopeG3[5];
  Float_t slopeG3_prod[6];
  Float_t slopeG4[5];
  Float_t slopeG4_prod[6];
  ProcessOutput("LHC17l3/TrackRefsOutput_LHC17l3_cent.root", cLayers, cLayers_prod, kRed, "", &slopeG3[0], &slopeG3_prod[0]);
  ProcessOutput("LHC17l4/TrackRefsOutput_LHC17l4_cent.root", cLayers, cLayers_prod, kBlue, "sames", &slopeG4[0], &slopeG4_prod[0]);
  Printf("Comparison layer-to-layer (SPD0 --> SPD1 --> SDD0 --> SDD1 -->  SSD0 --> SSD1)");
  for (Int_t i = 0; i < 5; i++){
    Printf("G3 --> %f, G4 --> %f, G3/G4 = %f", slopeG3[i], slopeG4[i], slopeG3[i]/slopeG4[i]);
  }
  Printf("Comparison with respect to generated value (SPD0, SPD1, SDD0, SDD1, SSD0, SSD1)");
  for (Int_t i = 0; i < 6; i++){
    Printf("G3 --> %f, G4 --> %f, G3/G4 = %f", slopeG3_prod[i], slopeG4_prod[i], slopeG3_prod[i]/slopeG4_prod[i]);
  }

}


//----------------------------------------------

TH2F* hdtgl[6]={0};// = new TH2F*[6];
TH2F* hdtgl_prod[6]={0};// = new TH2F*[6];

void ProcessOutput(const char* filename, TCanvas* cLayers, TCanvas* cLayers_prod, Int_t color, Option_t* opt, Float_t* slope, Float_t* slope_prod){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TFile* file = new TFile(filename);
  TList* listHisto = (TList*)file->Get("OutputHistos");
  TH1D* htemp; 
  TH1D* htemp_prod; 
  for (Int_t i = 0; i < 6; i++){
    if (i < 5){
      TObjArray tmpArr;
      hdtgl[i] = (TH2F*)listHisto->FindObject(Form("hdtgl%d", i));
      hdtgl[i]->FitSlicesY(0, 0, -1, 0, "QNR", &tmpArr);
      tmpArr.SetOwner(kFALSE);
      htemp = (TH1D*)tmpArr[2];//gDirectory->Get(Form("hdtgl%d_2", i));
      cLayers->cd(i+1);
      htemp->SetLineColor(color);
      htemp->Draw(opt);
      gPad->Modified();
      gPad->Update();
      htemp->Fit("pol1", "+");
      TF1* func = htemp->GetFunction("pol1");
      slope[i] = func->GetParameter(1);
      gPad->Modified();
      gPad->Update();
      TString sopt = opt;
      float dy = !sopt.IsNull() ? -0.2 : 0;
      SetStatPad(htemp, 0.1, 0.45, 0.75+dy, 0.9+dy);
      gPad->Modified();
      gPad->Update();

      SetHStyle(htemp,color,20, 1);

      gPad->Modified();
      gPad->Update();
      
      tmpArr.Clear();
      cLayers->Update();
    }
    TObjArray tmpArr_prod;
    hdtgl_prod[i] = (TH2F*)listHisto->FindObject(Form("hdtgl%d_prod", i));
    hdtgl_prod[i]->FitSlicesY(0, 0, -1, 0, "QNR", &tmpArr_prod);
    tmpArr_prod.SetOwner(kFALSE);
    htemp_prod = (TH1D*)tmpArr_prod[2];//gDirectory->Get(Form("hdtgl%d_2", i));
    cLayers_prod->cd(i+1);
    htemp_prod->SetLineColor(color);
    htemp_prod->Draw(opt);
    gPad->Modified();
    gPad->Update();
    htemp_prod->Fit("pol1", "+");
    TF1* func_prod = htemp_prod->GetFunction("pol1");
    slope_prod[i] = func_prod->GetParameter(1);
    gPad->Modified();
    gPad->Update();
    TString sopt = opt;
    float dy = !sopt.IsNull() ? -0.2 : 0;
    SetStatPad(htemp_prod, 0.1, 0.45, 0.75+dy, 0.9+dy);
    gPad->Modified();
    gPad->Update();
    
    SetHStyle(htemp_prod,color,20, 1);
    
    gPad->Modified();
    gPad->Update();
    
    tmpArr_prod.Clear();
    cLayers_prod->Update();    
  }
  return;
}



void SetHStyle(TH1* hst,int col,int mark,float mrsize)
{
  hst->SetLineColor(col);
  hst->SetMarkerColor(col);
  //  hst->SetFillColor(col);
  hst->SetMarkerStyle(mark);
  hst->SetMarkerSize(mrsize);
  TList *lst = hst->GetListOfFunctions();
  if (lst) {
    int nf = lst->GetSize();
    for (int i=0;i<nf;i++) {
      TObject *fnc = lst->At(i);
      if (fnc->InheritsFrom("TF1")) {
	((TF1*)fnc)->SetLineColor(col);
	((TF1*)fnc)->SetLineWidth(1);
	((TF1*)fnc)->ResetBit(TF1::kNotDraw);
      }
      else if (fnc->InheritsFrom("TPaveStats")) {
	((TPaveStats*)fnc)->SetTextColor(col);
      }
    }
  }
}

TPaveStats* SetStatPad(TH1* hst,float x1,float x2,float y1,float y2, Int_t stl, Int_t col)
{
  printf("Set statpad of %s to %f %f %f %f\n",hst->GetName(),x1,x2,y1,y2);
  TPaveStats* pad = GetStatPad(hst);
  if (!pad) return 0;
  pad->SetX1NDC( x1 );
  pad->SetX2NDC( x2 );
  pad->SetY1NDC( y1 );
  pad->SetY2NDC( y2 );
  if (stl>=0) pad->SetFillStyle(stl);
  if (col>=0) pad->SetFillColor(col);
  //
  gPad->Modified();
  return pad;
}

TPaveStats* GetStatPad(TH1* hst)
{
  TList *lst = hst->GetListOfFunctions();
  if (!lst) return 0;
  int nf = lst->GetSize();
  for (int i=0;i<nf;i++) {
    TPaveStats *fnc = (TPaveStats*) lst->At(i);
    if (fnc->InheritsFrom("TPaveStats")) return fnc;
  }
  return 0;
  //
}

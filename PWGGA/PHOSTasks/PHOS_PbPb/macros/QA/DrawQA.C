TFile* file;
const char* prefixToName = "imgs/";
const char* appendToName = ".pdf";

void Draw(const char* name, const char* options = "", double yFrom=0., double yTo=-1.)
{
  TH1* hist = ((TH1*)file->Get(name))->Clone();
  //hist->SetAxisRange(0, 30 );
  hist->GetXaxis()->SetLabelSize(0.011);
  //hist->GetXaxis()->SetTitle("Run");
  hist->GetXaxis()->SetRange(1,84);
  if( yFrom < yTo )
    hist->GetYaxis()->SetRangeUser(yFrom, yTo);


  //if( ! TString(options).Contains("same") )
  TCanvas* canv = new TCanvas;
  canv->SetGrid();
  hist->GetYaxis()->SetNdivisions(16);

  canv->Divide(1,2);

  hist->GetXaxis()->SetLabelSize(0.051);

  canv->cd(1);
  hist->GetXaxis()->SetRange(1,84);
  hist->DrawCopy(options);

  canv->cd(2);
  hist->GetXaxis()->SetRange(85,200);
  hist->DrawCopy(options);


  canv->SaveAs(Form("%s%s%s", prefixToName, hist->GetName(), appendToName ));
  delete hist;
}

void DrawPID()
{
  TH1* grNPhotAll = (TH1*)file->Get("grNPhotAll")->Clone();
  TH1* grNPhotDisp = (TH1*)file->Get("grNPhotDisp")->Clone();
  TH1* grNPhotDisp2 = (TH1*)file->Get("grNPhotDisp2")->Clone();
  TH1* grNPhotCPV = (TH1*)file->Get("grNPhotCPV")->Clone();
  TH1* grNPhotCPV2 = (TH1*)file->Get("grNPhotCPV2")->Clone();

  grNPhotAll->GetXaxis()->SetLabelSize(0.045);
  grNPhotAll->GetYaxis()->SetRangeUser(0, 40);

  grNPhotAll  ->SetMarkerColor(kBlack);
  grNPhotDisp ->SetMarkerColor(kCyan+1);
  grNPhotDisp2->SetMarkerColor(kBlue);
  grNPhotCPV  ->SetMarkerColor(kOrange+1);
  grNPhotCPV2 ->SetMarkerColor(kRed);

  bool scale = false;
  double scaleCPV =  grNPhotAll->Integral() / grNPhotCPV->Integral() *1.01;
  double scaleCPV2 =  grNPhotAll->Integral() / grNPhotCPV2->Integral() *1.02;
  double scaleDisp =  grNPhotAll->Integral() / grNPhotDisp->Integral() *0.99;
  double scaleDisp2 =  grNPhotAll->Integral() / grNPhotDisp2->Integral() *0.98;
  if(scale) {
    grNPhotCPV->Scale( scaleCPV );
    grNPhotCPV2->Scale( scaleCPV2 );
    grNPhotDisp->Scale( scaleDisp );
    grNPhotDisp2->Scale( scaleDisp2 );
    grNPhotAll->GetYaxis()->SetRangeUser(27, 37);
  }
  
  TCanvas* canv = new TCanvas;
  canv->Divide(1,2);

  canv->cd(1);
  grNPhotAll->SetTitle("#LTN_{clusters}^{PID}#GT");
  grNPhotAll->GetXaxis()->SetRange(0, 84);
  grNPhotAll->DrawCopy();
  grNPhotDisp->DrawCopy("same");
  grNPhotDisp2->DrawCopy("same");
  grNPhotCPV->DrawCopy("same");
  grNPhotCPV2->DrawCopy("same");

  canv->cd(2);
  grNPhotAll->SetTitle("");
  grNPhotAll->GetXaxis()->SetRange(85, 200);
  grNPhotAll->DrawCopy();
  grNPhotDisp->DrawCopy("same");
  grNPhotDisp2->DrawCopy("same");
  grNPhotCPV->DrawCopy("same");
  grNPhotCPV2->DrawCopy("same");

  canv->cd(1);
  leg = new TLegend(0.9,0.7,0.99,0.99);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  if( ! scale) {
    leg->AddEntry("grNPhotAll","All","lP");
    leg->AddEntry("grNPhotCPV","CPV","lP");
    leg->AddEntry("grNPhotCPV2","CPV2","lP");
    leg->AddEntry("grNPhotDisp","Disp","lP");
    leg->AddEntry("grNPhotDisp2","Disp2","lP");
  }
  else {
    leg->AddEntry("grNPhotAll",Form("All"),"lP");
    leg->AddEntry("grNPhotCPV",Form("CPV * %f", scaleCPV),"lP");
    leg->AddEntry("grNPhotCPV2",Form("CPV2 * %f", scaleCPV2),"lP");
    leg->AddEntry("grNPhotDisp",Form("Disp * %f", scaleDisp),"lP");
    leg->AddEntry("grNPhotDisp2",Form("Disp2 * %f", scaleDisp2),"lP");
  }
  leg->Draw();
  
  canv->SaveAs(Form("%s%s%s", prefixToName, "nPhotPID", appendToName ));

  
}

void DrawCPVRatio()
{
  TH1* grNPhotAll = (TH1*)file->Get("grNPhotAll")->Clone();
  TH1* grNPhotCPV = (TH1*)file->Get("grNPhotCPV")->Clone();
  TH1* grNPhotCPV2 = (TH1*)file->Get("grNPhotCPV2")->Clone();

  grNPhotCPV->Divide(grNPhotAll);
  grNPhotCPV->SetTitle(Form("%s / %s", grNPhotCPV->GetTitle(), grNPhotAll->GetTitle()));
  grNPhotCPV->GetYaxis()->SetRangeUser(0.7,0.85);

  grNPhotCPV2->Divide(grNPhotAll);

  TCanvas* canv = new TCanvas;
  canv->Divide(1,2);
  canv->cd(1);
  grNPhotCPV->GetXaxis()->SetRange(0, 84);
  grNPhotCPV->DrawCopy();

  canv->cd(2);
  grNPhotCPV->SetTitle("");
  grNPhotCPV->GetXaxis()->SetRange(85, 200);
  grNPhotCPV->DrawCopy();
  
  canv->SaveAs(Form("%s%s%s", prefixToName, "CPVtoAllRatio", appendToName ));
}

const Int_t kNPID = 8+4;
const char* kPIDNames[kNPID] = {"All", "Allwou", "Disp", "Disp2", "Dispwou", "CPV", "CPV2", "Both",
				"Allcore", "Dispcore", "CPVcore", "Bothcore"};
void DrawQA()
{
  gStyle->SetOptStat(0);

  file = TFile::Open("runQA.root", "read");

  Draw("grVtxZ10Cent", "", 0.7, 1.);
  Draw("grNCellsM1", "E");
  Draw("grNCellsM2");
  Draw("grNCellsM3");
  Draw("grECluster", "", 0.5, 0.7);
  Draw("grNCluster", "", 0, 40);
  Draw("grNTracks0", "", 0 , 12000);
  Draw("grNPhotAll", "", 0, 40);
  Draw("grNPhotAllcore", "", 0, 40);
  Draw("grNPhotAllwou", "", 0, 40);
  Draw("grNPhotDisp", "", 0, 40);
  Draw("grNPhotDisp2", "", 0, 40);
  Draw("grNPhotDispwou", "", 0, 40);
  Draw("grNPhotCPV", "", 0, 40);
  Draw("grNPhotCPV2", "", 0, 40);
  Draw("grNPhotBoth", "", 0, 40);
  Draw("grEnAll", "", 0.4, 0.7);
  Draw("grEnAllcore", "", 0.4, 0.7);
  Draw("grEnAllwou", "", 0.4, 0.7);
  Draw("grEnDisp", "", 0.4, 0.7);
  Draw("grEnDisp2", "", 0.4, 0.7);
  Draw("grEnDispcore", "", 0.4, 0.7);
  Draw("grEnDispwou", "", 0.4, 0.7);
  Draw("grEnCPV", "", 0.4, 0.7);
  Draw("grEnCPVcore", "", 0.4, 0.7);
  Draw("grEnCPV2", "", 0.4, 0.7);
  Draw("grEnBoth", "", 0.4, 0.7);
  Draw("grEnBothcore", "", 0.4, 0.7);


  DrawPID();
  DrawCPVRatio();
  file->Close();
}

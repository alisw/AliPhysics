void SetStyles(TH1 *histo,int marker, int color,char *xtitle = "p", char *ytitle = "E");
void PlotMatchedTracksDepositsFromData(TString datafilename="rootFiles/LHC10hPass2/Et.ESD.simPbPb.EMCal.LHC10hPass2.Run139465.root",TString simfilename="rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root", TString det = "EMCal",Bool_t allCells = kTRUE){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString outname = "/tmp/DataVsSimComparisons"+det;
  if(allCells){outname+="AllCells";}
  else{outname+="BkgdOnly";}
  TString textfilename = outname+".txt";
  outname+=".png";
  TString outnameratios = "/tmp/DataVsSimComparisonsRatios"+det;
  if(allCells){outnameratios+="AllCells";}
  else{outnameratios+="BkgdOnly";}
  outnameratios+=".png";
  TObjArray simulation(11);
  TObjArray data(11);
  TObjArray ratios(11);
//   TH1D *simulation[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};
//   TH1D *data[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};
  int colors[] = {TColor::kRed, TColor::kOrange-3,TColor::kGreen+3, TColor::kBlue, TColor::kViolet,
		  TColor::kRed, TColor::kOrange-3,TColor::kGreen+3, TColor::kBlue, TColor::kViolet};
  int markersSim[] = {24,25,26,32,27,24,25,26,32,27};
  int markersData[] = {20,21,22,23,33,20,21,22,23,33};

  float energycut = 0.5;

  TFile *fsim = TFile::Open(simfilename, "READ");
  TList *lsim = (TList*)fsim->Get("out1");

  TH3F  *fHistMatchedTracksEvspTBkgdvsMult = lsim->FindObject("fHistMatchedTracksEvspTBkgdMult");
  TH3F  *fHistMatchedTracksEvspTSignalvsMult = lsim->FindObject("fHistMatchedTracksEvspTSignalMult");
  if(allCells){
    fHistMatchedTracksEvspTBkgdvsMult->Add(fHistMatchedTracksEvspTSignalvsMult);
  }
  int lowbin = fHistMatchedTracksEvspTBkgdvsMult->GetXaxis()->FindBin(energycut);
  int highbin = fHistMatchedTracksEvspTBkgdvsMult->GetXaxis()->GetNbins();
  fHistMatchedTracksEvspTBkgdvsMult->GetXaxis()->SetRange(lowbin,highbin);

  for(int bin = 1; bin<=10;bin++){
    fHistMatchedTracksEvspTBkgdvsMult->GetZaxis()->SetRange(bin,bin);
    simulation[bin] = fHistMatchedTracksEvspTBkgdvsMult->ProjectionY("y");
    ((TH1D*)simulation[bin])->SetName(Form("Sim%i",bin));
    SetStyles((TH1D*)simulation[bin],markersSim[bin-1],colors[bin-1]);
  }
  TFile *f = TFile::Open(datafilename, "READ");
  f->cd();  
  TList *l = (TList*)f->Get("out1");
  TH3F  *fHistMatchedTracksEvspTvsMult = l->FindObject("fHistMatchedTracksEvspTvsMult");
  int lowbin = fHistMatchedTracksEvspTvsMult->GetXaxis()->FindBin(energycut);
  int highbin = fHistMatchedTracksEvspTvsMult->GetXaxis()->GetNbins();
  fHistMatchedTracksEvspTvsMult->GetXaxis()->SetRange(lowbin,highbin);

  for(int bin = 1; bin<=10;bin++){
    fHistMatchedTracksEvspTvsMult->GetZaxis()->SetRange(bin,bin);
    data[bin] = fHistMatchedTracksEvspTvsMult->Project3D("y");
    SetStyles((TH1D*)data[bin],markersData[bin-1],colors[bin-1]);
    ((TH1D*)data[bin])->SetName(Form("Sim%i",bin));
  }
  for(int bin = 1; bin<=10;bin++){
    float scale = ((TH1D*)data[bin])->Integral();
    ((TH1D*)data[bin])->Scale(1.0/scale);
    scale = ((TH1D*)simulation[bin])->Integral();
    ((TH1D*)simulation[bin])->Scale(1.0/scale);
    ratios[bin] = data[bin]->Clone(Form("ratio%i",bin));
    ((TH1D*)ratios[bin])->Divide((TH1D*)simulation[bin]);
  }

  TLegend *leg = new TLegend(0.631661,0.550143,0.731975,0.955587);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.038682);
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLogy();
  ((TH1D *)data[1])->GetXaxis()->SetRange(((TH1D *)data[1])->FindBin(0.3),((TH1D *)data[1])->GetNbinsX());
  ((TH1D *)data[1])->Draw();
  ofstream myfile;
  myfile.open (textfilename.Data());
  myfile<<"data\tsimulation"<<endl;
  for(int bin = 1; bin<=10;bin++){
    float scale = TMath::Power(2,bin);
    ((TH1D *)data[bin])->Scale(1.0/scale);
    ((TH1D *)simulation[bin])->Scale(1.0/scale);
    leg->AddEntry(data[bin],Form("%i<N_{cluster}<%i times %2.1f",(bin-1)*10,bin*10,scale));
    data[bin]->Draw("same");
    simulation[bin]->Draw("same");
    myfile<< ((TH1D *)data[bin])->GetMean() <<"\t"<< ((TH1D *)simulation[bin])->GetMean()<<endl;
  }
  myfile.close();
  leg->Draw();
  c1->SaveAs(outname.Data());

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  ((TH1D *)ratios[1])->GetYaxis()->SetTitle("data/simulation");

  ((TH1D *)ratios[1])->GetXaxis()->SetRange(((TH1D *)ratios[1])->FindBin(0.3),((TH1D *)ratios[1])->GetNbinsX());
  ((TH1D *)ratios[1])->Draw();
  for(int bin = 1; bin<=10;bin++){
    //float scale = 2-0.2*bin;
     //((TH1D *)ratios[bin])->Scale(1.0/scale);
    ratios[bin]->Draw("same");
  }
  ((TH1D *)ratios[1])->SetMaximum(1.5);
  ((TH1D *)ratios[1])->SetMinimum(0.0);
  //leg->Draw();
  c2->SaveAs(outnameratios.Data());
  return;

}

void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle){
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
  histo->Sumw2();
}


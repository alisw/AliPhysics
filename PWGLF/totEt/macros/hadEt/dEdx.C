//Christine Nattrass, University of Tennessee at Knoxville
//This macro is for  plotting the dE/dx from the TPC as used for particle identification for the transverse energy analysis
//Particles identified are colored
//Uses the output of AliAnalysisTaskHadEt
void Rebin(TH2F *histo){
  //cout<<"Nbins x "<<histo->GetNbinsX()<<" y "<<histo->GetNbinsY()<<endl;
  histo->Rebin2D(1,5);
  //histo->SetAxisRange(0.1,1.5);
}

void dEdx(char *prodname = "LHC11b1a PYTHIA 7 TeV p+p", char *shortprodname = "LHC11b1a",bool TPC = true, char *filename="rootFiles/LHC11b1a/Et.ESD.new.sim.LHC11b1a.root",bool data = true){
  TFile *file =  new TFile(filename);
  if(!file){
    cerr<<"Error, no file found"<<endl;
    return;
  }

  char *myname = "ITS";
  if(TPC) myname = "TPCITS";
  cout<<"Using "<<Form("dEdxAll%s",myname)<<endl;
  TString datastring = "";
  if(data) datastring = "Data";
  TH2F *all = out2->FindObject(Form("dEdx%sAll%s",datastring.Data(),myname));
  TH2F *pi = out2->FindObject(Form("dEdx%sPion%s",datastring.Data(),myname));
  TH2F *k = out2->FindObject(Form("dEdx%sKaon%s",datastring.Data(),myname));
  TH2F *p = out2->FindObject(Form("dEdx%sProton%s",datastring.Data(),myname));
  TH2F *e = out2->FindObject(Form("dEdx%sElectron%s",datastring.Data(),myname));
  gStyle->SetPalette(1);
  pi->SetMarkerColor(2);
  k->SetMarkerColor(3);
  p->SetMarkerColor(4);
  e->SetMarkerColor(TColor::kYellow);
  pi->SetLineColor(2);
  k->SetLineColor(3);
  p->SetLineColor(4);
  e->SetLineColor(TColor::kYellow);
  pi->SetFillColor(2);
  k->SetFillColor(3);
  p->SetFillColor(4);
  e->SetFillColor(TColor::kYellow);

  if(!TPC){all->GetXaxis()->SetRange(all->GetXaxis()->FindBin(0.05));}
  else{all->GetXaxis()->SetRange(all->GetXaxis()->FindBin(0.1));}
  //all->GetYaxis()->SetRange(all->GetYaxis()->FindBin(35.0));

  //e->SetMarkerStyle(20);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetLogx();
  all->Draw();
  //e->Draw("same");
  pi->Draw("same");
  k->Draw("same");
  p->Draw("same");
  TLegend *leg = new  TLegend(0.825503,0.768817,0.963087,0.954301);
 leg->AddEntry(pi,"#pi^{#pm}","F");
 leg->AddEntry(k,"K^{#pm}","F");
 leg->AddEntry(p,"p,#bar{p}","F");
 leg->AddEntry(e,"e^{#pm}","F");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
 leg->Draw();

  Rebin(all);
  Rebin(pi);
  Rebin(k);
  Rebin(p);
  Rebin(e);
  //all->SetAxisRange(0.12,1.5);
  if(!TPC){all->GetXaxis()->SetRange(all->GetXaxis()->FindBin(0.05));}
  else{all->GetXaxis()->SetRange(all->GetXaxis()->FindBin(0.1));}

 float y = 141.021;
 if(!TPC) y = 463.693;
 TLatex *tex = new TLatex(0.119617,y,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TString dir = "pics";
  TString detector = "";
  if(!TPC) detector = "ITS";
  gSystem->MakeDirectory(Form("%s/%s",dir.Data(),shortprodname));
  c->Print(Form("%s/%s/dEdx%s.pdf",dir.Data(),shortprodname,detector.Data()));
  c->Print(Form("%s/%s/dEdx%s.png",dir.Data(),shortprodname,detector.Data()));
  c->Print(Form("%s/%s/dEdx%s.eps",dir.Data(),shortprodname,detector.Data()));
  c->Print(Form("%s/%s/dEdx%s.jpg",dir.Data(),shortprodname,detector.Data()));
  return;

  if(TPC){
    c->SaveAs(Form("pics/%s/dEdx.eps",shortprodname));
    c->SaveAs(Form("pics/%s/dEdx.png",shortprodname));
    c->SaveAs(Form("pics/%s/dEdx.C",shortprodname));
  }
  else{
    c->SaveAs(Form("pics/%s/dEdxITS.eps",shortprodname));
    c->SaveAs(Form("pics/%s/dEdxITS.png",shortprodname));
  }
}

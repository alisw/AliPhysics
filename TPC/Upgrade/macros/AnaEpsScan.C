/*

.L $ALICE_ROOT/TPC/Upgrade/macros/AnaEpsScan.C
AnaEpsScan();

*/

void AnaEpsScan(TString dir=".",TString baseFile="0.0_1_2_130_10", Float_t nSigmas=3.)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.025);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  TString files=gSystem->GetFromPipe( Form("ls %s/eps*/*%s*.root", dir.Data(), baseFile.Data() ) );
  TObjArray *arr=files.Tokenize("\n");

  TGraph *grFrac01=new TGraph;
  TGraph *grFrac05=new TGraph;
  TGraph *grFrac10=new TGraph;

  grFrac01->SetNameTitle("grFrac01",";#varepsilon;fraction of tracks");
  grFrac01->SetMarkerStyle(20);
  grFrac01->SetMarkerSize(1);

//   grFrac05->SetNameTitle("grFrac05",";#varepsilon;fraction of tracks");
  grFrac05->SetLineColor(kBlue);
  grFrac05->SetMarkerColor(kBlue);
  grFrac05->SetMarkerStyle(21);
  grFrac05->SetMarkerSize(1);
  
  grFrac10->SetLineColor(kRed);
  grFrac10->SetMarkerColor(kRed);
  grFrac10->SetMarkerStyle(22);
  grFrac10->SetMarkerSize(1);
  
  Int_t colors[7]={kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow};
  Int_t markers[7]={20,21,22,23,24,25,26};
  
  TObjArray arrHists;
  for (Int_t ifile=0; ifile<arr->GetEntriesFast(); ++ifile) {
    TString file=arr->At(ifile)->GetName();
    TString epsilon=gSystem->GetFromPipe(Form("echo %s | sed 's|.*/eps\\([0-9][0-9]\\)/.*|\\1|'",file.Data()));

    printf("%s: %s\n", file.Data(), epsilon.Data());
    TH1F *h=new TH1F(Form("hResY%.0f_%s",10*nSigmas, epsilon.Data()), Form("#varepsilon %d;fraction of clusters more than %.1f#sigma from track;#tracks", nSigmas, epsilon.Atoi()),200,0,1);
    h->SetLineColor(colors[ifile]);
    h->SetMarkerColor(colors[ifile]);
    h->SetMarkerStyle(colors[ifile]);

    arrHists.Add(h);

    TFile f(file);
    gROOT->cd();
    TTree *t=(TTree*)f.Get("Tracks");
//     Float_t clFracY30=0.;

//     t->SetBranchStatus("*",0);
//     t->SetBranchStatus("clFracY30",1);
//     t->SetBranchAddress("clFracY30",&clFracY30);

    t->Draw(Form("clFracY%.0f>>hResY%.0f_%s",10*nSigmas,10*nSigmas,epsilon.Data()),"","goff");

    printf("entries: %d %d %d\n", grFrac05->GetN(), h->GetEntries(), t->GetEntries());

    Double_t frac01 = h->Integral(h->FindBin(.02),h->GetNbinsX())/h->GetEntries();
    Double_t frac05 = h->Integral(h->FindBin(.05),h->GetNbinsX())/h->GetEntries();
    Double_t frac10 = h->Integral(h->FindBin(.10),h->GetNbinsX())/h->GetEntries();

    grFrac01->SetPoint(grFrac01->GetN(), epsilon.Atoi(), frac01);
    grFrac05->SetPoint(grFrac05->GetN(), epsilon.Atoi(), frac05);
    grFrac10->SetPoint(grFrac10->GetN(), epsilon.Atoi(), frac10);

    delete t;
    f.Close();
  }

  TCanvas *c1=new TCanvas("c1");
  c1->cd();
  gPad->SetLogy();

  TLegend *leg = new TLegend(.7,.3,.9,.9);
  leg->SetBorderSize(1);
  leg->SetFillColor(10);
  
  for (Int_t ihist=0; ihist<arrHists.GetEntriesFast();++ihist) {
    TH1F *h=(TH1F*)arrHists.At(ihist);
    h->Draw((ihist==0)?"":"same");
    leg->AddEntry(h,h->GetTitle(),"lp");
  }
  leg->Draw("same");

  c1->SaveAs(Form("~/tmp/epsScan_clFrac_%.0fsigma.png",10*nSigmas));
  c1->SaveAs(Form("~/tmp/epsScan_clFrac_%.0fsigma.eps",10*nSigmas));
  
  TCanvas *c2=new TCanvas("c2");
  c2->cd();
  
  TLegend *leg2 = new TLegend(.1,.75,.6,.95);
  leg2->SetBorderSize(1);
  leg2->SetFillColor(10);

  TH1F *hDummy = new TH1F("hDummy",";#varepsilon;fraction of tracks",100,0,42.5);
  hDummy->SetMinimum(0);
  hDummy->SetMaximum(.21);
  hDummy->GetYaxis()->SetTitleOffset(1.2);
  hDummy->Draw();
  grFrac01->Draw("lp");
  grFrac05->Draw("lp");
  grFrac10->Draw("lp");

  //leg2->AddEntry(grFrac01,Form("%.1f#sigma deviation >2%%",nSigmas),"lp");
  //leg2->AddEntry(grFrac05,Form("%.1f#sigma deviation >5%%",nSigmas),"lp");
  //leg2->AddEntry(grFrac10,Form("%.1f#sigma deviation >10%%",nSigmas),"lp");
  leg2->AddEntry(grFrac01,"fraction of deviating clusters >2%","lp");
  leg2->AddEntry(grFrac05,"fraction of deviating clusters >5%","lp");
  leg2->AddEntry(grFrac10,"fraction of deviating clusters >10%","lp");
  TLatex l;
  l.DrawLatex(3,.14,Form("cluster deviation > %.1f#sigma",nSigmas));
  leg2->Draw("same");
  c2->SaveAs(Form("~/tmp/epsScan_trFrac_eps_%.0fsigma.png",10*nSigmas));
  c2->SaveAs(Form("~/tmp/epsScan_trFrac_eps_%.0fsigma.eps",10*nSigmas));
}


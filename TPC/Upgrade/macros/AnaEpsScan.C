/*

.L $ALICE_ROOT/TPC/Upgrade/macros/AnaEpsScan.C
AnaEpsScan();

*/

void AnaEpsScan(TString dir=".",TString baseFile="0.0_1_2_130_10")
{
  
  TString files=gSystem->GetFromPipe( Form("ls %s/eps*/*%s*.root", dir.Data(), baseFile.Data() ) );
  TObjArray *arr=files.Tokenize("\n");

  TGraph *grFrac05=new TGraph;
  TGraph *grFrac10=new TGraph;

  grFrac05->SetNameTitle("grFrac05",";#varepsilon;fraction of tracks");
  grFrac05->SetMarkerSize(1);
  
  grFrac10->SetLineColor(kRed);
  grFrac10->SetMarkerColor(kRed);
  grFrac10->SetMarkerStyle(21);
  grFrac10->SetMarkerSize(1);
  
  Int_t colors[7]={kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow};
  Int_t markers[7]={20,21,22,23,24,25,26};
  
  TObjArray arrHists;
  for (Int_t ifile=0; ifile<arr->GetEntriesFast(); ++ifile) {
    TString file=arr->At(ifile)->GetName();
    TString epsilon=gSystem->GetFromPipe(Form("echo %s | sed 's|.*/eps\\([0-9][0-9]\\)/.*|\\1|'",file.Data()));

    printf("%s: %s\n", file.Data(), epsilon.Data());
    TH1F *h=new TH1F(Form("hResY30_%s",epsilon.Data()), Form("#varepsilon %d;fraction of clusters more than 3#sigma from track;#tracks", epsilon.Atoi()),100,0,1);
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

    t->Draw(Form("clFracY30>>hResY30_%s",epsilon.Data()),"","goff");

    printf("entries: %d %d %d\n", grFrac05->GetN(), h->GetEntries(), t->GetEntries());

    Double_t frac05 = h->Integral(h->FindBin(.05),h->GetNbinsX())/h->GetEntries();
    Double_t frac10 = h->Integral(h->FindBin(.10),h->GetNbinsX())/h->GetEntries();

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

  TCanvas *c2=new TCanvas("c2");
  c2->cd();
  
  TLegend *leg2 = new TLegend(.1,.7,.5,.9);
  leg2->SetBorderSize(1);
  leg2->SetFillColor(10);
  
  grFrac05->Draw("alp");
  grFrac10->Draw("lp");

  leg2->AddEntry(grFrac05,"3#sigma deviation >5%","lp");
  leg2->AddEntry(grFrac10,"3#sigma deviation >10%","lp");
  leg2->Draw("same");
}


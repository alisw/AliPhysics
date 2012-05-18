void makeTrending(const Char_t *tl, const Char_t *fl, const Char_t *yax="", const Char_t *lab=0)
{
// Make trending of variable list "tl" from trending file list "fl"
// The trending value list should be formated "var1:var2:var3"
// The trending file from the list should be found on a path formated "your path"/runId/TRD.PerformanceTrend.root 

  TKey *key(NULL); 
  Int_t nr(0), run[10000]; Float_t val[10000][10];

  // get list of trending values for this plot
  TString stList(tl);
  TObjArray *ast=stList.Tokenize(":");
  Int_t nt(ast->GetEntries());

  FILE *fp = fopen(fl, "rt");
  TString sfp;
  while(sfp.Gets(fp)){
    TObjArray *afp=sfp.Tokenize("/");
    Int_t idx = afp->GetEntries()-2;
    Int_t rno = ((TObjString*)(*afp)[idx])->GetString().Atoi();
    if(!TFile::Open(sfp.Data())) continue;
    run[nr] = rno;
    for(Int_t it(0); it<nt; it++){
      if(!(key = (TKey*)gFile->Get(Form("TRDresolution_%s", ((TObjString*)(*ast)[it])->GetName() )))) {
        Warning("makeTrending()", "Missing %09d.TRDresolution.%s", rno, ((TObjString*)(*ast)[it])->GetName());
        continue;
      }
      val[nr][it] = atof(key->GetTitle());
    }
    nr++;
  }

//  TH1 *hT = new TH1F("hT", Form("%s;run;<%s>", tn, tn), nr, -0.5, nr-0.5);
  TH1 *hT = new TH1F("hT", Form(";RUN;%s", yax), nr, -0.5, nr-0.5);
  TAxis *ax = hT->GetXaxis(); ax->SetTitleOffset(1.2);ax->CenterTitle();
  TAxis *ay = hT->GetYaxis(); ay->SetTitleOffset(1.2);ay->CenterTitle();
  TGraph *gT[10];
  for(Int_t it(0); it<nt; it++){
    gT[it] = new TGraph(nr);
    gT[it]->SetMarkerStyle(4);gT[it]->SetMarkerColor(it+1);gT[it]->SetLineColor(it+1);
  }
  Float_t min(1.e5), max(-1.e5);
  for(Int_t ir(0); ir<nr; ir++){
    ax->SetBinLabel(ir+1, Form("%09d", run[ir]));
    for(Int_t it(0); it<nt; it++){
      gT[it]->SetPoint(ir, ir, val[ir][it]);
      if(val[ir][it]<min) min = val[ir][it];
      if(val[ir][it]>max) max = val[ir][it];
    }
  }
  hT->Draw("p");
  hT->GetYaxis()->SetRangeUser(0.5*(3*min-max), 0.5*(3*max-min));
  TLegend *leg = new TLegend(.6, .6, .98, .98);
  leg->SetBorderSize(1); leg->SetFillColor(kWhite);
  for(Int_t it(0); it<nt; it++){
    gT[it]->Draw("p"); leg->AddEntry(gT[it], lab?lab:((TObjString*)(*ast)[it])->GetName(), "p");
  }
  leg->Draw();
}
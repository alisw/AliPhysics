
struct Summarizer
{
  TFile* fFile;
  TCanvas* fCanvas;
  TVirtualPad* fTop;
  TVirtualPad* fBody;
  TLatex*      fTitle;
  Summarizer()
    : fFile(0),
      fCanvas(0),
      fTop(0),
      fBody(0),
      fTitle(0)
  {}
  TObject* GetO(TDirectory* d, const char* name, TClass* cls=0)
  {
    if (!d) return 0;
    TObject* o = d->Get(name);
    if (!o) {
      ::Warning("GetO", "No object %s in %s", name, d->GetName());
      return 0;
    }
    if (!cls) return o;

    if (!o->IsA()->InheritsFrom(cls)) {
      Warning("GetO", "Object %s read from %s, is not a %s but a %s",
	      name, d->GetName(), cls->GetName(), o->ClassName());
      return 0;
    }
    return o;
  }
  THStack* GetHS(TDirectory* d, const char* name)
  {
    return static_cast<THStack*>(GetO(d, name, THStack::Class()));
  }
  TH1* GetH1(TDirectory* d, const char* name)
  {
    return static_cast<TH1*>(GetO(d, name, TH1::Class()));
  }
  TH2* GetH2(TDirectory* d, const char* name)
  {
    return static_cast<TH2*>(GetO(d, name, TH2::Class()));
  }
  TDirectory* GetD(TDirectory* d, const char* name)
  {
    if (!d) return;
    TDirectory* sub = d->GetDirectory(name);
    if (!sub) {
      ::Warning("GetO", "No directory %s in %s", name, d->GetName());
      return 0;
    }
    return sub;
  }
  
  void Run(const char* filename)
  {
    fFile = TFile::Open(filename, "READ");
    if (!fFile) return;

    TString  outName(filename); outName.ReplaceAll(".root", ".pdf");
    Int_t    h = 1000;
    Int_t    w = h/TMath::Sqrt(2);

    fCanvas = new TCanvas(outName, outName, w, h);
    fTop    = new TPad("top","top",0,.9,1,1);
    fTop->SetFillColor(kYellow-10);
    fTop->SetFillStyle(1001);
    fCanvas->cd();
    fTop->Draw();
    fBody   = new TPad("body","body",0,0,1,.9);
    fBody->SetFillColor(kWhite);
    fBody->SetFillStyle(1001);
    fBody->SetTopMargin(0.01);
    fBody->SetRightMargin(0.03);
    fCanvas->cd();
    fBody->Draw();

    fTitle = new TLatex(0.5, 0.5, "");
    fTitle->SetTextAlign(22);
    fTitle->SetTextSize(0.3);
    fTop->cd();
    fTitle->Draw();

    fCanvas->Print(Form("%s[", fCanvas->GetName()), "PDF");

    TAxis* centAxis = VisualizeResults();
    for (Int_t i = 1; i <= centAxis->GetNbins(); i++) {
      Double_t c1 = centAxis->GetBinLowEdge(i); 
      Double_t c2 = centAxis->GetBinUpEdge(i);
      
      VisualizeBin(c1, c2);
    }
    
    fCanvas->Print(Form("%s]", fCanvas->GetName()), "PDF");
  }
  TAxis* VisualizeResults()
  {
    THStack* result = GetHS(fFile, "result");
    fBody->SetLogy();
    fBody->SetTicks();
    fBody->cd();
    result->Draw("nostack");
    Print("Results");

    TH1* cent = GetH1(fFile, "cent");
    TH1* mid  = GetH1(fFile, "mid");
    fBody->cd();
    fBody->Divide(1,2,0,0);
    TVirtualPad* p = fBody->cd(1);
    p->SetRightMargin(0.01);
    p->SetTicks();
    cent->Draw("hist");
    p = fBody->cd(2);
    p->SetRightMargin(0.01);
    p->SetTicks();
    mid ->Draw();
    Print("Centralities");
    
    return cent->GetXaxis();
  }
  void VisualizeBin(Double_t c1, Double_t c2)
  {
    TString binName;
    binName.Form("cent%03dd%02d_%03dd%02d",
		 Int_t(c1), (100*Int_t(c1)) % 100, 
		 Int_t(c2), (100*Int_t(c2)) % 100);

    TDirectory* binDir = GetD(fFile, binName);
    if (!binDir) return;

    THStack* summary = GetHS(binDir, "summary");
    fBody->cd();
    fBody->SetRightMargin(0.2);
    summary->Draw("nostack");
    fBody->BuildLegend(1-fBody->GetRightMargin(),
		       fBody->GetBottomMargin(),
		       .99,
		       1-fBody->GetTopMargin());
    Print(Form("%5.2f-%5.2f%% - Calculations", c1, c2));

    TDirectory* detDir = GetD(binDir, "details");
    THStack*    deltas = GetHS(detDir, "deltas");
    fBody->cd();
    fBody->SetLogx();
    fBody->SetLogy();
    fBody->SetRightMargin(0.2);
    deltas->Draw("nostack");
    fBody->BuildLegend(1-fBody->GetRightMargin(),
		       fBody->GetBottomMargin(),
		       .99,
		       1-fBody->GetTopMargin());
    Print(Form("%5.2f-%5.2f%% - #Delta", c1, c2));
    

  }
  void Print(const char* title="")
  {
    fTitle->SetTitle(title);
    fTop->cd();
    fTitle->Draw();
    fCanvas->Modified();
    fCanvas->Update();
    fCanvas->cd();
    
    fCanvas->Print(Form("%s", fCanvas->GetName()), Form("PDF Title=%s", title));

    fCanvas->WaitPrimitive();
    
    fBody->Clear();
    fTop->Clear();
    fBody->SetLogy(false);
    fBody->SetLogx(false);
    fBody->SetTicks();
  }
};

void Summarize(const char* name)
{
  Summarizer* s = new Summarizer;
  s->Run(name);
}

//
// EOF
// 

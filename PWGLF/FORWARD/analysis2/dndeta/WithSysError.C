/** 
 * 
 * 
 * @param system 
 * @param sNN 
 * @param trigger 
 * @param option 
 * @param rebinned 
 * @param empirical 
 */
void
WithSysError(const TString&  system, 
	     UShort_t        sNN, 
	     const TString&  trigger, 
	     const Option_t* option="e5",
	     Bool_t          rebinned=true, 
	     Bool_t          empirical=true)
{
  const char* fwd =
    gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  gROOT->SetMacroPath(Form("%s:%s/dndeta:%s/gse", gROOT->GetMacroPath(),
			   fwd, fwd));
  gSystem->AddIncludePath(Form("-I%s/dndeta -I%s/gse", fwd,fwd));

  if (!gROOT->GetClass("Drawer"))       gROOT->LoadMacro("Drawer.C+g");
  if (!gROOT->GetClass("GraphSysErr"))  gROOT->LoadMacro("GraphSysErr.C+g");
  if (!gROOT->GetClass("SysErrorAdder"))gROOT->LoadMacro("SysErrorAdder.C+g");


  TCanvas* c = new TCanvas("C", "C");
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);

  Double_t    ly1     = 0.7;
  Double_t    ly2     = 0.9;
  const char* trigs[] = { trigger.Data(), 0 };
  const char* exps[]  = { "ALICE", "WIP", 0 };
  TLegend*    tl      = new TLegend(0.15,ly1,0.3,ly2);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  TObjArray   u;
  Double_t*   effs    = 0;
  if (sNN == 8000 && system.EqualTo("pp", TString::kIgnoreCase)) {
    effs = new Double_t[2];
    effs[1] = 0;
    if (trigger.EqualTo("INEL",TString::kIgnoreCase))
      effs[0] = 0.85;
    else if (trigger.EqualTo("NSD",TString::kIgnoreCase) ||
	     trigger.EqualTo("V0AND",TString::kIgnoreCase)) {
      trigger  = "NSD";
      trigs[0] = "NSD";
      effs[0]  = 0.93;
    }
  }

  TLegend* el = new TLegend(tl->GetX2(),ly1,tl->GetX2()+.2,ly2);
  el->SetBorderSize(0);
  el->SetFillStyle(0);
  SysErrorAdder* adder = 0;
  if (trigger.EqualTo("INEL",TString::kIgnoreCase))
    adder = new INELAdder(system, sNN);
  else if (trigger.EqualTo("NSD",TString::kIgnoreCase))
    adder = new NSDAdder(system, sNN);
  else 
    adder = new CENTAdder(system, sNN, trigger);
  
  TPair* dataOther = Drawer::GetDataOther(tl, u, system, sNN, trigs,
					  exps, option, rebinned,
					  empirical, effs);
  
  if (!dataOther || !dataOther->Key()) {
    Error("", "No data found %s", canvas->GetTitle());
    return;
  }
  THStack*     data  = static_cast<THStack*>(dataOther->Key());
  TMultiGraph* other = static_cast<TMultiGraph*>(dataOther->Value());

  TList* outList = new TList;
  TIter next(data->GetHists());
  TH1* h = 0;
  Bool_t first = true;
  const char* opt = "stack stat quad split west west";
  while ((h = static_cast<TH1*>(next()))) {
    TString n(h->GetName());
    if (n.Contains("syserror", TString::kIgnoreCase)) continue;

    GraphSysErr* gse = adder->Make(h, (first ? el : 0));
    gse->SetTitle(""); // Form("%s @ %d (%s)", system.Data(), sNN, trigger.Data()));
    gse->Draw(Form("%s %s", (first ? "axis" : ""), opt));
    outList->Add(gse);
    if (first) { 
      TMultiGraph* axis = gse->GetMulti();
      TH1*         frame = axis->GetHistogram();
      frame->SetMinimum(0);
      frame->SetMaximum(frame->GetMaximum()*1.3);
    }
    // gse->Export(false);
    first = false;
    }
  other->Draw("p");
  
  TParameter<int>* um =
    new TParameter<int>("PWG-LF/GEO - work in progress", 20);
  um->SetUniqueID(kBlack);
  u.Add(um);

  TLegend* ul = Drawer::MakeUniqueLegend(.6,.8,.95,.97, u, 1);
  ul->Draw();
  tl->Draw();
  el->Draw();

  TLatex* ltx = Drawer::MakeTitle(.2, .98, system, sNN, trigger);
  ltx->SetTextAlign(13);

  TString outName(Form("%s_%04d_%s.root", system.Data(), sNN, trigger.Data()));
  TFile* out = TFile::Open(outName.Data(), "RECREATE");
  outList->Write();
  out->Write();
  
}

  

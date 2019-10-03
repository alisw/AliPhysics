/** 
 * 
 * 
 * @param system    System name 
 * @param sNN       Collision energy in GeV 
 * @param trigger   Trigger type 
 * @param option    Options 
 * @param rebinned  Rebinning factor 
 * @param export    Whether to export 
 * @param empirical Whether to apply empirical correction 
 * @param alsoLog   Whether to do log scale 
 */
void
WithSysError(const TString&  system, 
	     UShort_t        sNN, 
	     TString&        trigger, 
	     const Option_t* option="e5",
	     Bool_t          rebinned=true,
	     Bool_t          export=true,
	     Bool_t          empirical=true,
	     Bool_t          alsoLog=true)
{
  // --- Search path -------------------------------------------------
  const char* fwd = 0;
  if (gSystem->Getenv("FWD"))
    fwd = gSystem->Getenv("FWD");
  else if (gSystem->Getenv("ANA_SRC")) 
    fwd = gSystem->Getenv("ANA_SRC");
  else 
    fwd = gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  gROOT->SetMacroPath(Form("%s:%s/dndeta:%s/gse:%s/scripts",
			   gROOT->GetMacroPath(), fwd, fwd, fwd));
  gSystem->AddIncludePath(Form("-I%s/dndeta -I%s/gse", fwd,fwd));

  // --- Load code ---------------------------------------------------
  if (!gROOT->GetClass("Drawer"))       gROOT->LoadMacro("Drawer.C+g");
  if (!gROOT->GetClass("GraphSysErr"))  gROOT->LoadMacro("GraphSysErr.C+g");
  if (!gROOT->GetClass("SysErrorAdder"))gROOT->LoadMacro("SysErrorAdder.C+g");


  // --- Reaction key ------------------------------------------------
  TString reac = system;
  reac.ToUpper();
  reac.Insert(reac.Index("P", 1, 1), " ");
  reac.Append(" --> CHARGED X");
  
  // --- Fix up trigger and efficiency -------------------------------
  const char* trigs[] = { trigger.Data(), 0 };
  const char* exps[]  = { "ALICE", "WIP", 0 };
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

  // --- Select how to add -------------------------------------------
  SysErrorAdder* adder = 0;
  TString trigLegTitle = ""; // trigger;
  if (trigger.EqualTo("INEL",TString::kIgnoreCase))
    adder = new INELAdder(system, sNN);
  if (trigger.EqualTo("INELGt0",TString::kIgnoreCase))
    adder = new INELGt0Adder(system, sNN);
  else if (trigger.EqualTo("NSD",TString::kIgnoreCase))
    adder = new NSDAdder(system, sNN);
  else {
    trigLegTitle = "Centralities:";
    adder = new CENTAdder(system, sNN, trigger);
  }

  // --- Canvas ------------------------------------------------------
  TCanvas* c = new TCanvas("C", "C", 1600, 1000);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.20);
  c->SetLeftMargin(0.12);
  
  // --- Legend parameters -------------------------------------------
  Double_t    ly2     = 0.98;
  Double_t    lx1     = 0.81;
  Double_t    lx2     = 0.98;
  Double_t    ldy     = 0.2;
  Double_t    ly1     = ly2-(trigLegTitle.IsNull() ? ldy : 2*ldy);
  // --- Legend for triggers -----------------------------------------
  // TLegend* tl = new TLegend(0.15,ly1,0.3,ly2);
  TLegend* tl = new TLegend(lx1,ly1,lx2,ly2, trigLegTitle);
  tl->SetBorderSize(0);
  tl->SetFillStyle(0);
  tl->SetTextColor(Drawer::AliceBlue());
  tl->SetTextFont(42);

  // --- Legend for errors -------------------------------------------
  // TLegend* el = new TLegend(tl->GetX2(),ly1,tl->GetX2()+.2,ly2);
  ly2 = ly1;
  ly1 = ly2-ldy;
  TLegend* el = new TLegend(lx1,ly1,lx2,ly2, "Sys. errors:");
  el->SetBorderSize(0);
  el->SetFillStyle(0);
  el->SetTextColor(Drawer::AliceBlue());
  el->SetTextFont(42);

  // --- Legend for unique names -------------------------------------
  // TLegend* ul = new TLegend(.6,.8,.95,.97);
  ly2 = ly1;
  ly1 = ly2-.5*ldy;
  TLegend* ul = new TLegend(lx1, ly1, lx2, ly2);
  ul->SetBorderSize(0);
  ul->SetFillStyle(0);
  ul->SetTextColor(Drawer::AliceBlue());
  ul->SetTextFont(42);

  // --- Get the data ------------------------------------------------
  TObjArray   u;
  TPair* dataOther = Drawer::GetDataOther(tl, u, system, sNN, trigs,
					  exps, option, rebinned,
					  empirical, effs);  
  if (!dataOther || !dataOther->Key()) {
    Error("", "No data found %s", canvas->GetTitle());
    return;
  }
  THStack*     data  = static_cast<THStack*>(dataOther->Key());
  TMultiGraph* other = static_cast<TMultiGraph*>(dataOther->Value());

  // --- Loop over the data ------------------------------------------
  TList* outList = new TList;
  TIter next(data->GetHists());
  TH1* h = 0;
  Bool_t first = true;
  TH1* frame = 0;
  const char* opt = "stack stat quad split west west";
  while ((h = static_cast<TH1*>(next()))) {
    TString n(h->GetName());
    if (n.Contains("syserror", TString::kIgnoreCase)) continue;

    // Color_t col = h->GetMarkerColor();
    // h->SetMarkerColor(col);
    // h->SetLineColor(col);
    // h->SetFillColor(col);
    GraphSysErr* gse = adder->Make(h, (first ? el : 0));
    gse->SetTitle(""); 
    gse->Draw(Form("%s %s", (first ? "axis" : ""), opt));
    gse->SetKey("laboratory", "CERN");
    gse->SetKey("accelerator", "LHC");
    gse->SetKey("detector", Form("FORWARD%s", rebinned ? "" : "_full"));
    gse->SetKey("reackey", reac);
    gse->SetKey("obskey", "DN/DETARAP");
    gse->SetKey("title", "Systematic study of 1/N dNch/deta over widest possible eta region at the LHC");
    gse->SetKey("author", "CHRISTENSEN");
    gse->SetKey("comment", "We present 1/N dNch/deta over widest possible eta region at the LHC");
    gse->SetKey("dscomment", "The pseudo-rapidity density of charged particle");


    if (export) {
      TString trg = adder->GetTriggerString();
      if (trg.EqualTo("INEL>0")) trg = "INELGt0";
      TString dir(Form("out/%s/%05d/%s", system.Data(),
		       sNN, trg.Data()));
      gSystem->Exec(Form("mkdir -p %s", dir.Data()));
      
      TFile* file = TFile::Open(Form("%s/%s.root",
				     dir.Data(),
				     gse->GetKey("detector")),
				"RECREATE");
      Info("", "Writing to %s", file->GetPath());
      gse->Write("data");
      file->Write();
      file->Close();
    }
	
    
    outList->Add(gse);
    if (first) { 
      TMultiGraph* axis = gse->GetMulti();
      frame = axis->GetHistogram();
      frame->SetMinimum(0.1);
      frame->SetMaximum(frame->GetMaximum()*1.1);
      frame->GetYaxis()->SetTitleOffset(1.4);
    }
    // gse->Export(false);
    first = false;
  }

  // --- Draw the others ---------------------------------------------
  if (other) other->Draw("p");

  // --- Build legend of unique names --------------------------------
  TParameter<int>* um =
    new TParameter<int>("PWG-LF/GEO - work in progress", 20);
  um->SetUniqueID(trigLegTitle.IsNull() ? kRed+2 : kBlack);
  u.Add(um);
  Drawer::MakeUniqueLegend(ul, u, 1);
  ul->Draw();
  tl->Draw();
  el->Draw();

  // --- Make a title ------------------------------------------------
  TLatex* ltx = Drawer::MakeTitle(.45, .98, system, sNN, trigger);
  ltx->SetTextAlign(23);
  ltx->SetTextSize(0.04);

  // --- ALICE logo --------------------------------------------------
  if (!gROOT->GetClass("AliceLogo"))
    gROOT->LoadMacro("AliceLogo.C+");
      
  if (gROOT->GetClass("AliceLogo")) {
    ly1 = ly1-ldy;
    c->Range(0,0,1,1);
    gROOT->ProcessLine("AliceLogo* al = new AliceLogo();");
    gROOT->ProcessLine(Form("al->Draw(0,%f,%f,%f, 0, 0xf);", lx1,ly1,ldy));
  }
  
  
  // --- Output to disk ----------------------------------------------
  TString base(Form("%s_%04d_%s", system.Data(), sNN, trigger.Data()));
  TString outName(Form("%s.root", base.Data()));
  TFile* out = TFile::Open(outName.Data(), "RECREATE");
  outList->Write();
  out->Write();
  Info("", "Wrote to %s.root", base.Data());

  c->Modified();
  c->Update();
  c->cd(); 
  c->Print(Form("%s.pdf", base.Data()));
  c->Print(Form("%s.png", base.Data()));
  c->SaveAs(Form("%s_canvas.root", base.Data()));
  /// c->SaveAs(Form("%s_canvas.C", base.Data()));

  if (trigLegTitle.IsNull() || !alsoLog) return;
  
  frame->SetMinimum(5);
  frame->SetMaximum(3*frame->GetMaximum());
  c->SetLogy();
  c->Modified();
  c->Update();
  c->cd(); 
  c->Print(Form("%s_logy.pdf", base.Data()));
  c->Print(Form("%s_logy.png", base.Data()));
 

}

  
//
// EOF
//

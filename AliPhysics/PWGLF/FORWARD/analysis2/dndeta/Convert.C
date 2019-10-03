


TList*
ConvertOne(const TString& system,
	   UShort_t       sNN,
	   const TString& trigger,
	   Bool_t         rebin=true)
{
  TString t(trigger);
  if (t.EqualTo("INEL>0")) t="INELGt0";
  if (t.EqualTo("NSD"))    t="V0AND";
  TString scr(Form("nosec/%s/%04d/%s/%s.C",
		   system.Data(), sNN, t.Data(), (rebin ? "rebin" : "full")));

  if (gSystem->AccessPathName(scr.Data())) {
    Warning("ConvertOne", "Script %s does not exist", scr.Data());
    return;
  }

  THStack* stack = new THStack("", "");
  gROOT->Macro(Form("%s((THStack*)%p,0,20)", scr.Data(), stack));
  if (!stack->GetHists()) {
    Warning("ConvertOne", "Failed to import from %s", scr.Data());
    return;
  }

  SysErrorAdder* adder = 0;
  if      (t.EqualTo("INEL"))    adder = new INELAdder(system,sNN);
  else if (t.EqualTo("INELGt0")) adder = new INELGt0Adder(system,sNN);
  else if (t.EqualTo("V0AND"))   adder = new NSDAdder(system,sNN);
  else if (t.BeginsWith("CENT")) adder = new CENTAdder(system,sNN,t);
  else {
    Error("ConvertOne", "Don't know how to add to %s", t.Data());
    return;
  }

  TList* ret = new TList;
  ret->SetName(Form("%s_%05d_%s", system.Data(), sNN, t.Data()));
  
  TIter  next(stack->GetHists());
  TH1*   hist = 0;
  while ((hist = static_cast<TH1*>(next()))) {
    TString name(hist->GetName());
    if (name.Contains("mirror") || name.Contains("SysError")) {
      // Info("ConvertOne", "Ignoring %s", name.Data());
      continue;
    }

    // Info("ConvertOne", "Converting %s", hist->GetName());
    hist->SetMarkerColor(kBlack);
    hist->SetFillColor(kBlack);
    hist->SetLineColor(kBlack);

    if (system.EqualTo("pp") && sNN == 8000) {
      if (t.EqualTo("INEL"))  hist->Scale(0.85);
      if (t.EqualTo("V0AND")) hist->Scale(0.93);
    }

    GraphSysErr* gse = adder->Make(hist,0);
    gse->SetName("data");
    
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
    // file->Close();
    ret->Add(gse);
  }
  return ret;
}

void
ConvertAndDraw(const TString& system,
	       UShort_t       sNN,
	       const TString& trigger,
	       Bool_t         rebin=true)
{
  TList*   l = ConvertOne(system, sNN, trigger, rebin);
  if (!l) return;
  if (gROOT->IsBatch()) return;
  
  TCanvas* c = new TCanvas(l->GetName(), l->GetName());
  TIter    n(l);
  TObject* o = 0;
  Bool_t   f = true;
  while ((o = n())) {
    if (f) {
      o->Draw("STACK SPLIT QUAD AXIS");
      GraphSysErr* g = static_cast<GraphSysErr*>(o);
      TH1*         h = g->GetMulti()->GetHistogram();
      h->SetMaximum(h->GetMaximum()*1.2);
      h->SetMinimum(0);
    }
    else
      o->Draw("STACK SPLIT QUAD");
    f = false;
  }
}

void
ConvertPbPb()
{
  ConvertAndDraw("PbPb", 2760, "CENT");
}

void
ConvertpPb()
{
 const char*  paCents[] = { "V0X", "ZNX", "V0M", "CL1", 0 };
  const char** ptrCent   = paCents;
  while (*ptrCent) {
    const char*  paSys[]  = { "pPb", "Pbp", 0 };
    const char*  paSide[] = { "A", "C", 0 };
    const char** ptrSys   = paSys;
    const char** ptrSide  = paSide;
    while (*ptrSys) {
      TString trg(Form("CENT%s", *ptrCent));
      trg.ReplaceAll("X", *ptrSide);
      
      ConvertAndDraw(*ptrSys, 5023, trg);

      ptrSys++;
      ptrSide++;
    }
    ptrCent++;
  }
}
void
ConvertPP()
{
  UShort_t  ppSNN[] = { 900, 2760, 7000, 8000, 0 };
  UShort_t* ptrSNN  = ppSNN;
  while (*ptrSNN) {
    const char*  ppTrig[] = { "INEL", "INELGt0", "NSD", 0 };
    const char** ptrTrig  = ppTrig;

    while (*ptrTrig) {
      ConvertAndDraw("pp", *ptrSNN, *ptrTrig);
      ptrTrig++;
    }
    ptrSNN++;
  }
}
  
void
Convert()
{
  const char* fwd = ".";
  gROOT->SetMacroPath(Form("%s:%s/dndeta:%s/gse:%s/scripts",
			   gROOT->GetMacroPath(), fwd, fwd, fwd));
  gSystem->AddIncludePath(Form("-I%s/dndeta -I%s/gse", fwd,fwd));

  if (!gROOT->GetClass("GraphSysErr"))  gROOT->LoadMacro("GraphSysErr.C+g");
  if (!gROOT->GetClass("SysErrorAdder"))gROOT->LoadMacro("SysErrorAdder.C+g");

  gSystem->Exec("rm -rf out");
  ConvertPbPb();
  ConvertpPb();
  ConvertPP();

 
  gSystem->Exec("(cd out && tar -czvf ../fwd.tar.gz .)");
  
  
}


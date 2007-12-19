// void
// RunAnaESD(Bool_t bg=true)
{
  Int_t  v2    = 6;
  Int_t  nBins = 1;
  Bool_t bg    = kTRUE;
  // gROOT->LoadMacro("Compile.C");
  Compile("AliFMDAnaFlowRing.h"); 
  Compile("AliFMDAnaFlow.C"); 
  AliFMDAnaFlow af(nBins,bg); 
  af.Run(); 
  TBrowser b;
  b.Add(&af);
  af.ToFile(Form("flow_v%02d_b%02d_%s.root", v2, nBins, (bg?"bg":"nobg"))); 
}

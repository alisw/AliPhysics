// Draw one result 
void Export(const TString&  system, 
	    UShort_t        sNN, 
	    const TString&  trigger, 
	    Bool_t          rebinned=false, 
	    Bool_t          empirical=true)
{
  if (gSystem->Getenv("FWD"))
    fwd = gSystem->Getenv("FWD");
  else 
    fwd = gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  gROOT->SetMacroPath(Form("%s/dndeta:%s", gROOT->GetMacroPath(),fwd));

  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C+");

  Drawer::Export(system, sNN, trigger, rebinned, empirical);
}
// EOF


  

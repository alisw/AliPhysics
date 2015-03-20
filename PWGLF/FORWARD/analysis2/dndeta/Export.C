// Draw one result 
void Export(const TString&  system, 
	    UShort_t        sNN, 
	    const TString&  trigger, 
	    Bool_t          rebinned=false, 
	    Bool_t          empirical=true)
{
  gROOT->SetMacroPath(Form("%s:%s", gROOT->GetMacroPath(),
			   "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta"));

  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C+");

  Drawer::Export(system, sNN, trigger, rebinned, empirical);
}
// EOF


  

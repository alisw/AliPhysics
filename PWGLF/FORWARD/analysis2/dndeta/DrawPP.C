// Draw one result 
void DrawPP(UShort_t        sNN, 
	    const Option_t* option="e5",
	    Bool_t          rebinned=true)
{
  gROOT->SetMacroPath(Form("%s:%s", gROOT->GetMacroPath(),
			   "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta"));

  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C++g");

  Drawer::DrawPP(sNN, option, rebinned);
}
// EOF

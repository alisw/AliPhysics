// Draw one result 
void DrawPP(UShort_t        sNN, 
	    const Option_t* option="e5",
	    Bool_t          rebinned=true)
{
  if (gSystem->Getenv("FWD"))
    fwd = gSystem->Getenv("FWD");
  else 
    fwd = gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  gROOT->SetMacroPath(Form("%s/dndeta:%s", gROOT->GetMacroPath(),fwd));
  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C++g");

  Drawer::DrawPP(sNN, option, rebinned);
}
// EOF

// Draw one result 
void Draw(const TString&  system, 
	  UShort_t        sNN, 
	  const TString&  trigger, 
	  const Option_t* option="e5",
	  Bool_t          rebinned=true, 
	  Bool_t          empirical=true,
	  Bool_t          symmetrice=false)
{
  gROOT->SetMacroPath(Form("%s:%s", gROOT->GetMacroPath(),
			   "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta"));

  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C+g");

  if (symmetrice && system.EqualTo("pPb")) 
    Drawer::pPbSym(trigger, option, rebinned, empirical);
  else 
    Drawer::Draw(system, sNN, trigger, option,  rebinned, empirical);
}
// EOF


  

// Draw scaled result 
void ScaleBypp(const TString& system, 
	       UShort_t       sNN, 
	       UShort_t       ppsNN, 
	       const TString& trigger, 
	       const TString& ppTrigger, 
	       Double_t       etaShift=0,
	       Bool_t         rebinned=true, 
	       Bool_t         empirical=true,
	       Bool_t         symmetrice=false,
	       Bool_t         write=false)
{
  if (gSystem->Getenv("FWD"))
    fwd = gSystem->Getenv("FWD");
  else 
    fwd = gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  gROOT->SetMacroPath(Form("%s/dndeta:%s", gROOT->GetMacroPath(),fwd));

  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C+");

  if (symmetrice && system.EqualTo("pPb"))
    Drawer::SymScaleBypp(ppsNN, trigger, ppTrigger, "",
			 rebinned, empirical, write);
  else
    Drawer::ScaleBypp(system, sNN, ppsNN, trigger, ppTrigger, 
		      etaShift, rebinned, empirical, write);
}
void ScaleBypp()
{
  ScaleBypp("pPb", 5023, 7000, "CENTZNX", "NSD",  0, true, true, true, true);
  ScaleBypp("pPb", 5023, 7000, "CENTV0X", "NSD",  0, true, true, true, true);
  ScaleBypp("pPb", 5023, 7000, "CENTZNX", "INEL", 0, true, true, true, true);
  ScaleBypp("pPb", 5023, 7000, "CENTV0X", "INEL", 0, true, true, true, true);
}
  
// EOF


  
  

void
DrawMultCuts(ULong_t       runNo=999, 
	     UShort_t      sys=0, 
	     UShort_t      sNN=0, 
	     Short_t       field=999, 
	     Bool_t        mc=false, 
	     const Char_t* local=0)
{
  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  if (gSystem->Getenv("ANA_SRC"))
    fwd = gSystem->Getenv("ANA_SRC");
  gROOT->SetMacroPath(Form("%s/scripts:%s/corrs:%s", 
			   gROOT->GetMacroPath(), fwd, fwd));
  if (!gROOT->GetClass("AliForwardCorrectionManager"))
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  gSystem->AddIncludePath(Form("-I$ALICE_ROOT/include "
			       "-I$ALICE_PHYSICS/include "
			       "-I%s/scripts -I%s/corrs -I%s",
			       fwd, fwd, fwd));
  gROOT->LoadMacro(Form("%s/scripts/SummaryDrawer.C++", fwd));
  gROOT->LoadMacro(Form("%s/corrs/MultCutDrawer.C++", fwd));

  MultCutDrawer* mcd = new MultCutDrawer();
  mcd->Run(runNo, sys, sNN, field, mc, local);
}

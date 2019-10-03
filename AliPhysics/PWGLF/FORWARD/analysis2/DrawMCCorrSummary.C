/**
 * @file   DrawMCCorrSummary.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 30 09:47:30 2012
 * 
 * @brief  Script to draw summary of AOD pass into a PDF 
 * 
 * 
 */
//____________________________________________________________________
void DrawMCCorrSummary(const char* fname="forward_mccorr.root", 
		       UShort_t what=0x20F)
{
  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  gROOT->SetMacroPath(Form("%s:%s/scripts",
			   gROOT->GetMacroPath(), fwd));
  gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));
  gROOT->LoadMacro(Form("%s/scripts/SummaryMCCorrDrawer.C+g",fwd));
  
  SummaryMCCorrDrawer d;
  d.Run(fname, what);
}
//
// EOF
//

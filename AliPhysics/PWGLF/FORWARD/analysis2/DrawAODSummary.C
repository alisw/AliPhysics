/**
 * @file   DrawAODSummary.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 30 09:47:30 2012
 * 
 * @brief  Script to draw summary of AOD pass into a PDF 
 * 
 * 
 */
//____________________________________________________________________
void DrawAODSummary(const char* fname="forward.root", UShort_t what=0x27F)
{
  gROOT->SetMacroPath(Form("%s:$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts",
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("SummaryAODDrawer.C+g");
  
  SummaryAODDrawer d;
  d.Run(fname, what);
}
//
// EOF
//

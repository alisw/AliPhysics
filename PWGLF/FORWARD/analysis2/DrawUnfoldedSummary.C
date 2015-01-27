/**
 * @file   DrawUnfoldedSummary.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 30 09:47:30 2012
 * 
 * @brief  Script to draw summary of Unfolding pass into a PDF 
 * 
 * 
 */
//____________________________________________________________________
void DrawUnfoldedSummary(const char* fname="forward_unfolded.root")
{
  gROOT->SetMacroPath(Form("%s:$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts",
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("SummaryUnfoldedDrawer.C++g");
  
  SummaryUnfoldedDrawer d;
  d.Run(fname);
}

//
// EOF
//

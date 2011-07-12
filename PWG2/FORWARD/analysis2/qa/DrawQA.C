/**
 * @file   DrawQA.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 13:45:26 2011
 * 
 * @brief  Script to draw most QA stuff 
 * 
 * @ingroup pwg2_forward_scripts_qa
 */
/** 
 * Draw most QA stuff 
 * 
 * @param file File to read 
 * @param full If true, assume output of analysis 
 */
void
DrawQA(const char* file, bool full=false)
{

  gROOT->SetMacroPath(Form(".:$(ALICE_ROOT)/PWG2/FORWARD/analysis2/qa:"
			   "$(ALICE_ROOT)/PWG2/FORWARD/analysis2/corrs:%s",
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("DrawBeforeAfter.C");
  gROOT->LoadMacro("DrawELossPoisson.C");
  gROOT->LoadMacro("DrawNeighbors.C");
  gROOT->LoadMacro("DrawOccupancy.C");
  gROOT->LoadMacro("DrawRecAnaEloss.C");

  DrawBeforeAfter(file);
  DrawELossPoisson(file);
  DrawNeighbors(file);
  DrawOccupancy(file);
  DrawRecAnaEloss(file);

  if (!full) { 
    gROOT->LoadMacro("DrawAnaELoss.C");
    DrawAnaELoss(file);
  }
  else { 
    gROOT->LoadMacro("DrawSteps.C");
    DrawSteps(file);
  }
}
//
// EOF
//

  

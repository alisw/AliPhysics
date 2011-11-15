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
  gROOT->LoadMacro("DrawBeforeAfter.C+g");
  gROOT->LoadMacro("DrawELossPoisson.C+g");
  gROOT->LoadMacro("DrawNeighbors.C+g");
  gROOT->LoadMacro("DrawOccupancy.C+g");
  gROOT->LoadMacro("DrawRecAnaEloss.C+g");
  gROOT->LoadMacro("Draw123.C+g");

  Info("DrawQA", "Drawing before-after");
  DrawBeforeAfter(file);
  Info("DrawQA", "Drawing singles, doubles, tripples");
  Draw123(file);
  Info("DrawQA", "Drawing Neighbors");
  DrawNeighbors(file);
  Info("DrawQA", "Drawing raw and analysed energy loss");
  DrawRecAnaEloss(file);
  Info("DrawQA", "Drawing poisson vs energy loss");
  DrawELossPoisson(file);
  Info("DrawQA", "Drawing Occupancies");
  DrawOccupancy(file);

  if (!full) { 
    Info("DrawQA", "Drawing fit results");
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

  

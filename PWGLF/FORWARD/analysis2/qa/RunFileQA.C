/**
 * @file   RunFileQA.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jan  6 11:43:37 2012
 * 
 * @brief  Script to run a run QA
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
/** 
 * script to run a run QA.  Note, on errors, the ROOT session is terminated
 * 
 * @param input   Input file 
 * @param output  Output tree file (optional)
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
void
RunFileQA(const char* input, const char* output=0)
{
  int ret = 0;
  gROOT->SetMacroPath(Form(".:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2/qa:"
			   "$(ALICE_ROOT)/PWGLF/FORWARD/analysis2/corrs:%s",
			   gROOT->GetMacroPath()));
  gSystem->AddIncludePath("-I\${ALICE_ROOT}/PWGLF/FORWARD/analysis2/qa");
  gSystem->Load("libGpad");
  gSystem->Load("libTree");
  
  gROOT->LoadMacro("QABase.h+g");
  gROOT->LoadMacro("QATrender.C+g");
  
  QATrender t(true, false);
  t.AddFile(input);
  if (output && output[0] != '\0')
    t.SetOutputName(output);
  if (!t.Run()) exit(1);
}
//
// EOF
//

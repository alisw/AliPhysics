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
 * @param prodYear Production year 
 * @param prodLetter Production letter
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
void
RunFileQA(const char* input, const char* output=0,
	  Int_t prodYear=0, const char* prodLetter="")
{
  int ret = 0;
  gROOT->SetMacroPath(Form(".:%s", gROOT->GetMacroPath()));
  gSystem->Load("libGpad");
  gSystem->Load("libTree");
  
  gROOT->LoadMacro("QABase.h+g");
  gROOT->LoadMacro("QATrender.C+g");
  
  QATrender t(true, false, prodYear, prodLetter[0]);
  t.AddFile(input);
  if (output && output[0] != '\0')
    t.SetOutputName(output);
  if (!t.Run()) exit(1);
}
//
// EOF
//

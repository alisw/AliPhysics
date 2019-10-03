/**
 * @file   CheckQAFile.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jan  6 11:39:22 2012
 * 
 * @brief  Script to check a QA file
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
/** 
 * Script to check a QA file.  Note that this terminates the ROOT session. 
 * 
 * @param filename File to read 
 * @param type What type 
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
void CheckQAFile(const char* filename, const char* type="")
{
  int    ret  = 0;
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("CheckQAFile", "No such file %s", filename);
    exit(1);
  }
  TString check    = Form("Forward%s",type);
  TString check2;
  TString check3;
  TObject* forward1 = file->Get(check);
  if (!forward1) {
    check2    = Form("Forward%sSums",type);
    forward1 = file->Get(Form("Forward%sSums",type));
    if (!forward1) {
      check3 = "ForwardSums";
      forward1 = file->Get(check3);
      if (!forward1) { 
	Error("CheckQAFile", "No %s, %s, or %s object found in %s", 
	      check.Data(),check2.Data(), 
	      check3.Data(), filename);
	ret |= 2;
      }
    }
  } 
  check = Form("Forward%sResults", type);
  TObject* forward2 = file->Get(check);
  if (!forward2) {
    check2   = Form("forward%sResults",type);
    forward2 = file->Get(check2);
    if (!forward2) { 
      check3 = "ForwardResults";
      forward2 = file->Get(check3);
      if (!forward2) { 
	Error("CheckQAFile", "No %s, %s, or %s object found in %s", 
	      check.Data(), check2.Data(), 
	      check3.Data(), filename);
	file->ls();
	ret |= 4;
      }
    }
  } 
  file->Close();
  exit(ret);
}
//
// EOF
//

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
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
void CheckQAFile(const char* filename)
{
  int    ret  = 0;
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("CheckQAFile", "No such file %s", filename);
    exit(1);
  }
  TObject* forward1 = file->Get("Forward");
  if (!forward1) {
    Error("CheckQAFile", "No Forward object found in %s", filename);
    ret |= 2;
  } 
  TObject* forward2 = file->Get("ForwardResults");
  if (!forward2) {
    Error("CheckQAFile", "No ForwardResults object found in %s", filename);
    ret |= 4;
  } 
  file->Close();
  exit(ret);
}
//
// EOF
//

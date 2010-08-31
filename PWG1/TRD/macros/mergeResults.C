void mergeResults(Char_t *files, Char_t *file="QAresults.root")
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWG1.so");
  
  AliTRDpwg1Helper::MergeProd(file, files, 10);
}
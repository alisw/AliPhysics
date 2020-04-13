void treeinfo(const char* fileName)
{
  TFile::Open(fileName);
  ((TTree*) gFile->Get("aodTree"))->Print();
}

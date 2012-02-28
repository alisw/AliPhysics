void mergeResults(Char_t *files, Char_t *file="QAresults.root")
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libPWGmuon.so");

  TDatime dt; gRandom->SetSeed(dt.Get());
  gSystem->Exec("mkdir -p merge; rm -rf merge/*");
  AliTRDpwgppHelper::MergeProd(file, files, 5);
  gSystem->Exec("rm -rfv merge");
}

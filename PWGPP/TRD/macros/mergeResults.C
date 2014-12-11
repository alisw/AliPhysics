void mergeResults(Char_t *files, Char_t *file="QAresults.root")
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGmuon");

  TDatime dt; gRandom->SetSeed(dt.Get());
  gSystem->Exec("mkdir -p merge; rm -rf merge/*");
  AliTRDpwgppHelper::MergeProd(file, files, 5);
  gSystem->Exec("rm -rfv merge");
}

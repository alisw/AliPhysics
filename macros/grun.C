Int_t ievent;
void grun (Int_t nevent=1, const char *config="Config.C")
{
  //
  // Simple macro to run aliroot in a batch mode
  //
  ievent=nevent;
  gAlice->Init(config);
  gAlice->Run(ievent);
}

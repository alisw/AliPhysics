void runkalman()
{

  AliL3Kalman *k = new AliL3Kalman("./",0,0);
  k->Init();
  k->LoadTracks(0,kTRUE);
  k->WriteFiles();
  k->ProcessTracks();
  delete k;
}

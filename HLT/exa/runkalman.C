void runkalman(Char_t *path = "/tmp/tvik")
{ 
  Bool_t binary=kFALSE; //Assume input is RLE binary files, or rootfile.
  Bool_t pileup=kFALSE; //Assume input is pileup event = non RLE binary files.
  Int_t npatches = 1;   //Options; 1, 2 and 6.
  Char_t trackparams[] = "SetTrackingParameters_1000bf04.C"; //Set this to correspond 
                                                             //with mult. and BField
  //for aliroot the path should point to a file 
  //containing the tpc geometry called alirunfile.root
  Bool_t isinit=AliHLTTransform::Init("./",!binary);

  AliHLTKalman *k = new AliHLTKalman(path,0,0);
  k->Init();
  //k->DoMakeSeed();
  k->LoadTracks(0,kTRUE);
  k->WriteFiles(path);
  k->ProcessTracks();
  delete k;
}

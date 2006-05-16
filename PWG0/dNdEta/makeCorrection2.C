//
// Script to make correction maps for dndeta measurements using the
// dNdEtaCorrection class.
//
// implementation with TSelector
//

#include "../CreateESDChain.C"

void makeCorrection2(Char_t* dataDir, Int_t nRuns=20)
{
  gSystem->Load("libPWG0base");

  /*gSystem->Load("../esdTrackCuts/libESDtrackQuality.so");
  gSystem->Load("libdNdEta.so");
  gSystem->Load("../AliSelector_cxx.so");*/
  //gSystem->Load("AlidNdEtaEffSelector_cxx.so");

  TChain* chain = CreateESDChain(dataDir, nRuns);
  chain->Process("AlidNdEtaEffSelector.cxx+");
}

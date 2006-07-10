/* $Id$ */

//
// Script to test the class AlidNdEtaVertexRecEffSelector that creates
// a vertex reconstruction efficiency plot
//
// implementation with TSelector
//

#include "../CreateESDChain.C"

runVertexRecEff(Char_t* dataDir, Int_t nRuns=20, Int_t offset=0)
{
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");

  TChain* chain = CreateESDChain(dataDir, nRuns, offset);

  TString selectorName = "AlidNdEtaVertexRecEffSelector";

  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  TStopwatch timer;
  timer.Start();

  chain->Process(selectorName + ".cxx+");

  timer.Stop();
  timer.Print();
}

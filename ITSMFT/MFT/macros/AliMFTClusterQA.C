//================================================================================================================================

void AliMFTClusterQA(const Char_t *readDir= ".",
		     const Char_t *outDir = ".",
		     Int_t nEventsToAnalyze = -1) {
  
  AliMFTClusterQA *myAnalysis = new AliMFTClusterQA();
  myAnalysis -> Init(readDir, outDir, nEventsToAnalyze);

  while (myAnalysis->LoadNextEvent()) continue;

  myAnalysis->Terminate();

}

//================================================================================================================================


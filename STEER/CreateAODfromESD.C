#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TFile.h"
#include "AliReconstruction.h"

#endif

// THIS MACRO IS NOT SUPPORTED ANYMORE. Goodbye. MDO 26/03/2008
void CreateAODfromESD(const char *inFileName = "AliESDs.root",
		      const char *outFileName = "AliAOD.root") {
  
  // open input and output files
  TFile *esdFile = TFile::Open(inFileName, "READ");
  TFile *aodFile = TFile::Open(outFileName, "RECREATE");
  
  // call conversion method
  AliReconstruction reco;
  reco.ESDFile2AODFile(esdFile, aodFile);

  // close files
  esdFile->Close();
  aodFile->Close();
}

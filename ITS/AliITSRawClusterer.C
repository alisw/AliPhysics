#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSclustererV2.h"
#include "AliRawReaderRoot.h"
#include "AliRunLoader.h"
#include <TFile.h>
#include <TTree.h>
#endif

void AliITSRawClusterer(const char* fileNameRawData = "event.root",
			Int_t iEvent = 0,
			const char* fileNameGalice = "galice.root")
{
// To run the cluster finder on raw data a galice.root file is needed.
// This file has to contain a run loader, a loader for ITS and
// the geometry of the ITS

  AliRawReaderRoot rawReader(fileNameRawData, iEvent);

  AliRunLoader* runLoader = AliRunLoader::Open(fileNameGalice);
  runLoader->CdGAFile();
  AliITSgeom* geom = (AliITSgeom*) gFile->Get("AliITSgeom");
  AliITSclustererV2 clusterer(geom);

  runLoader->GetLoader("ITSLoader")->LoadRecPoints("recreate");
  runLoader->SetEventNumber(iEvent);

  clusterer.Digits2Clusters(&rawReader);
}

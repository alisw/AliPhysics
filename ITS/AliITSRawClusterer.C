#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSclustererV2.h"
#include <TFile.h>
#include <TTree.h>
#endif

void AliITSRawClusterer(const char* fileNameParam = "its.digits.root",
			const char* fileNameClusters = "its.clusters.root")
{
  delete gAlice;
  TFile* file = TFile::Open(fileNameParam);
  AliRun* gAlice = (AliRun*) file->Get("gAlice");
  AliITS* its = (AliITS*) gAlice->GetModule("ITS");
  AliITSgeom* geom = (AliITSgeom*) its->GetITSgeom();
  AliITSclustererV2 clusterer(geom);

  TFile* out = TFile::Open(fileNameClusters, "recreate");
  geom->Write();

  clusterer.Digits2Clusters(out);

  out->Close();
  delete out;
  file->Close();
  delete file;
}

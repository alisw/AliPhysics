#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTPC.h"
#include "AliTPCclustererMI.h"
#include <TFile.h>
#include <TTree.h>
#endif

void AliTPCRawClusterer(const char* fileNameParam,
			const char* fileNameClusters = "tpc.clusters.root",
			const char* fileNameRawData = "AltroFormatDDL.dat")
{
  TFile* fileParam = TFile::Open(fileNameParam);
  AliTPCclustererMI clusterer(AliTPC::LoadTPCParam(fileParam));
  TFile* fileClusters = TFile::Open(fileNameClusters, "recreate");
  TTree* output = new TTree("TreeC_TPC_0", "TreeC_TPC_0"); 

  clusterer.SetOutput(output);
  clusterer.Digits2Clusters(fileNameRawData);

  fileClusters->Close();
  delete fileClusters;
  fileParam->Close();
  delete fileParam;
}

// ----------------------------------------------------//
//                                                     //
//    This macro does Digits to Reconstructed Points   //
//                                                     //
// ----------------------------------------------------//

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include <stdlib.h>

void AliPMDDigits2Recpoints(Int_t nevt=1) 
{
  TStopwatch timer;
  timer.Start();

  // Open the AliRoot file
  AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root");
  if (!fRunLoader)
    {
      cerr<<"Can't load RunLoader"<<endl;
      return 1;
    }
  fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();

  printf(" Do you want reconstruction from Digits file or RAW data \n");
  printf(" If RAW,    type 0 \n");
  printf(" If Digits, type 1 \n");
  Int_t itype;
  cin >> itype;

  // Create the PMD Cluster Finder 
  AliPMDClusterFinder *clus = new AliPMDClusterFinder(fRunLoader);
  clus->SetDebug(1);
  if (itype == 0)
    {
      clus->Load();
    }
  else if (itype == 1)
    {
      clus->LoadClusters();
    }  


  for (Int_t ievt = 0; ievt < nevt; ievt++)
    {
      if (itype == 0)
	{
	  // from digits data
	  clus->Digits2RecPoints(ievt);
	}
      else if (itype == 1)
	{
	  // from raw data
	  AliRawReaderFile reader(ievt);
	  clus->Digits2RecPoints(ievt, &reader);
	}
    }
  if (itype == 0)
    {
      clus->UnLoad();
    }
  else if (itype == 1)
    {
      clus->UnLoadClusters();
    }
  timer.Stop();
  timer.Print();
}


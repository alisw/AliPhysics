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


  // Create the PMD Cluster Finder 
  AliPMDClusterFinder *clus = new AliPMDClusterFinder(fRunLoader);
  
  clus->SetDebug(1);
  clus->Load();



  for (Int_t ievt = 0; ievt < nevt; ievt++)
    {
      // from digits data
      //      clus->Digits2RecPoints(ievt);

      // from raw data
      AliRawReaderFile reader(ievt);
      clus->Digits2RecPoints(ievt, &reader);
    }
  clus->UnLoad();

  timer.Stop();
  timer.Print();
}


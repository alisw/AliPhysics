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

void AliPMDReconstruction() 
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

  if (itype == 0)
    {
      AliRawReaderFile reader;
      AliPMDReconstructor pmdreco;
      pmdreco.Reconstruct(fRunLoader, &reader);
    }
  else if (itype == 1)
    {
      AliPMDReconstructor pmdreco;
      pmdreco.Reconstruct(fRunLoader);
    }

  timer.Stop();
  timer.Print();
}


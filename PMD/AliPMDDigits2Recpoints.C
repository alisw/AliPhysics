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

  // Input file name
  Char_t *alifile = "galice.root"; 

  // Create the PMD Cluster Finder 
  AliPMDClusterFinder *clus = new AliPMDClusterFinder();
  
  // Open the AliRoot file
  clus->OpengAliceFile(alifile,"DR");
  
  Int_t ievt;
  
  for (ievt = 0; ievt < nevt; ievt++)
    {
      clus->Digits2RecPoints(ievt);
    }
  clus->UnLoad("R");
  timer.Stop();
  timer.Print();
}


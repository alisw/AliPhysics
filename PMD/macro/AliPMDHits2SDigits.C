// ----------------------------------------------------//
//                                                     //
//       This macro does Hits to SDigits               //
//                                                     //
// ----------------------------------------------------//

#include <stdlib.h>
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"

void AliPMDHits2SDigits(Int_t nevt=1) 
{
  TStopwatch timer;
  timer.Start();
  Float_t zpos = 361.5;

  // Input (and output) file name
  Char_t *alifile = "galice.root"; 

  // Create the PMD digitzer 
  AliPMDDigitizer *digitizer = new AliPMDDigitizer();
                                                  
  // Open the AliRoot file
  digitizer->OpengAliceFile(alifile,"HS");

  digitizer->SetZPosition(zpos);
  
  // Create the sdigits

  for (Int_t ievt=0; ievt<nevt; ievt++)
    {
      digitizer->Hits2SDigits(ievt);
    }

  digitizer->UnLoad("S");
  timer.Stop();
  timer.Print();
  
}


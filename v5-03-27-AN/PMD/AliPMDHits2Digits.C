// ----------------------------------------------------//
//                                                     //
//       This macro does Hits to Digits                //
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

void AliPMDHits2Digits(Int_t nevt=1) 
{
  // Input file name
  Char_t *alifile = "galice.root"; 

  // Create the PMD digitzer 
  AliPMDDigitizer *digitizer = new AliPMDDigitizer();
                                                  
  // Open the AliRoot file
  digitizer->OpengAliceFile(alifile,"HS");

  // Create the digits
  for (Int_t ievt=0; ievt<nevt; ievt++)
    {
      digitizer->Hits2Digits(ievt);
    }
  digitizer->UnLoad("S");
}


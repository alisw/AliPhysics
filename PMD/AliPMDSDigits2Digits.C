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


void AliPMDSDigits2Digits(Int_t nevt=1) 
{
  // Input (and output) file name
  Char_t *alifile = "galice.root"; 

  // Create the PMD digitzer 
  AliPMDDigitizer *digitizer = new AliPMDDigitizer();
                                                  
  // Open the AliRoot file
  digitizer->OpengAliceFile(alifile,"SD");

  Int_t ievt;

  for (ievt = 0; ievt < nevt; ievt++)
    {
      digitizer->SDigits2Digits(ievt);
    }
  digitizer->UnLoad("S");

}


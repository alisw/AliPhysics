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
#include "AliVZERODigitizer.h"

void VZEROHits2Digits() 
{
   
  // Input file name
  Char_t *alifile = "galice.root";   

  // Create the VZERO digitizer 
  AliVZERODigitizer *manager = new AliVZERODigitizer();

  // Open the AliRoot file
  manager->OpengAliceFile(alifile);
                                                   
  // Create the digits

  manager->Exec();

}

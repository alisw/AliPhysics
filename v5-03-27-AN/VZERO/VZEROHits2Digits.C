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
#include "AliRunDigitizer.h"
#include "AliVZERODigitizer.h"

void VZEROHits2Digits() 
{
   
  // Input file name
  Char_t *alifile = "galice.root";   

  // Create the run digitizer 
  AliRunDigitizer* manager = new AliRunDigitizer(1, 1);
  manager->SetInputStream(0, alifile);
  manager->SetOutputFile("H2Dfile");

  // Create the VZERO digitizer 
  new AliVZERODigitizer(manager);

  // Create the digits
  manager->Exec("");

}

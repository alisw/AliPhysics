#include "AliEmptyPreprocessor.h"

// This preprocessor is used as a placeholder for non-existing preprocessors
// during the FDR. Its task is just to fail, so that the run does not stay
// in processing state forever.

ClassImp(AliEmptyPreprocessor)

//______________________________________________________________________________________________
AliEmptyPreprocessor::AliEmptyPreprocessor(AliShuttleInterface* shuttle, const char* detector) :
  AliPreprocessor(detector, shuttle)
{
  // constructor
}

//______________________________________________________________________________________________
AliEmptyPreprocessor::~AliEmptyPreprocessor()
{
  // destructor
}


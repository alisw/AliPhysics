// $Id$
//
// Emcal particle trigger class, which can contain either
//
// Author: J.Kral

#include "AliEmcalTriggerSetupInfo.h"
#include "AliLog.h"

//_________________________________________________________________________________________________
AliEmcalTriggerSetupInfo::AliEmcalTriggerSetupInfo() :
  TNamed(),
  fThresholds(),
  fThresholdsSimple()
{
  // Default constructor.
  for( int i = 0; i < 4; i++ )
    fThresholds[i] = -1;
  for( int i = 0; i < 4; i++ )
    fThresholdsSimple[i] = -1;
}

  
//_________________________________________________________________________________________________
AliEmcalTriggerSetupInfo::AliEmcalTriggerSetupInfo(const AliEmcalTriggerSetupInfo &p) :
  TNamed(p)
{
  // Copy constructor.
  for( int i = 0; i < 4; i++ )
    fThresholds[i] = p.fThresholds[i];
  for( int i = 0; i < 4; i++ )
    fThresholdsSimple[i] = p.fThresholdsSimple[i];
}

//_________________________________________________________________________________________________
AliEmcalTriggerSetupInfo::~AliEmcalTriggerSetupInfo()
{
  // Destructor.
}

//_________________________________________________________________________________________________
AliEmcalTriggerSetupInfo &AliEmcalTriggerSetupInfo::operator=(const AliEmcalTriggerSetupInfo &p)
{
  // Assignment operator.

  if (this != &p) {
    for( int i = 0; i < 4; i++ )
      fThresholds[i] = p.fThresholds[i];
    for( int i = 0; i < 4; i++ )
      fThresholdsSimple[i] = p.fThresholdsSimple[i];
  }

  return *this;
}

//_________________________________________________________________________________________________
void AliEmcalTriggerSetupInfo::Clean(){
  // cleaner
  for( int i = 0; i < 4; i++ )
    fThresholds[i] = -1;
  for( int i = 0; i < 4; i++ )
    fThresholdsSimple[i] = -1;
}


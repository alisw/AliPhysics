/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliEmcalTriggerSetupInfo.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerSetupInfo)
/// \endcond

/**
 * Default constructor
 */
AliEmcalTriggerSetupInfo::AliEmcalTriggerSetupInfo() :
  TNamed(),
  fThresholds(),
  fThresholdsSimple()
{
  for( int i = 0; i < 4; i++ )
    fThresholds[i] = -1;
  for( int i = 0; i < 4; i++ )
    fThresholdsSimple[i] = -1;
}

/**
 * Copy constructor
 * @param p Reference for the copy
 */
AliEmcalTriggerSetupInfo::AliEmcalTriggerSetupInfo(const AliEmcalTriggerSetupInfo &p) :
  TNamed(p)
{
  for( int i = 0; i < 4; i++ )
    fThresholds[i] = p.fThresholds[i];
  for( int i = 0; i < 4; i++ )
    fThresholdsSimple[i] = p.fThresholdsSimple[i];
}

/**
 * Destructor
 */
AliEmcalTriggerSetupInfo::~AliEmcalTriggerSetupInfo()
{
}

/**
 * Assignment operator
 * @param p Reference for the assignment
 * @return This object after assignment
 */
AliEmcalTriggerSetupInfo &AliEmcalTriggerSetupInfo::operator=(const AliEmcalTriggerSetupInfo &p)
{
  if (this != &p) {
    for( int i = 0; i < 4; i++ )
      fThresholds[i] = p.fThresholds[i];
    for( int i = 0; i < 4; i++ )
      fThresholdsSimple[i] = p.fThresholdsSimple[i];
  }

  return *this;
}

/**
 * Cleaning function, set all thresholds to 0
 */
void AliEmcalTriggerSetupInfo::Clean(){
  for( int i = 0; i < 4; i++ )
    fThresholds[i] = -1;
  for( int i = 0; i < 4; i++ )
    fThresholdsSimple[i] = -1;
}


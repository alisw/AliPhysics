// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2002

//=============================================================================
// particle (track) pair
// Allows to calculate the kinematic pair variables typically used in 
// two-particle correlation analyses. 
//=============================================================================

#include "AliDPair.h"

ClassImp(AliDPair)

//=============================================================================
AliDPair::AliDPair() : fP0(), fP1(), fP(), fQ(), fBeta(), fBetat(), fBetaz(), fUbeta(), 
		 fUbetat(), fUbetaz(), fCMp(), fCMq(), fBuf() 
{
  // constructor

  printf("AliDPair object created\n");
}
//=============================================================================

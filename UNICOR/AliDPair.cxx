// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2002

//=============================================================================
// particle (track) pair
// Allows to calculate the kinematic pair variables typically used in 
// two-particle correlation analyses. 
//=============================================================================

#include "AliDPair.h"

ClassImp(AliDPair)

//=============================================================================
AliDPair::AliDPair() : p0(), p1(), p(), q(), beta(), betat(), betaz(), ubeta(), 
		 ubetat(), ubetaz(), CMp(), CMq(), buf() 
{
  // constructor

  printf("AliDPair object created\n");
}
//=============================================================================

/***************************************************************************
 *
 * $Id$
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.8  2005/07/06 18:49:56  fisyak
 * Replace AliFmHelixD, AliFmLorentzVectorD,AliFmLorentzVectorF,AliFmMatrixD,AliFmMatrixF,AliFmPhysicalHelixD,AliFmThreeVectorD,AliFmThreeVectorF by templated version
 *

****************************************************************************/
#ifndef ST_HELIX_D_HH
#define ST_HELIX_D_HH
#include "AliFmThreeVectorD.h"
#include "AliFmHelix.h"
#include <utility>
typedef AliFmHelix AliFmHelixD;
typedef pair<double,double> pairD;
#endif

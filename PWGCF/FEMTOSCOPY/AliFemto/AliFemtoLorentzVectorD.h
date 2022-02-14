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
 * Revision 1.6  2005/07/06 18:49:56  fisyak
 * Replace AliFemtoHelixD, AliFemtoLorentzVectorD,AliFemtoLorentzVectorF,AliFemtoMatrixD,AliFemtoMatrixF,AliFemtoPhysicalHelixD,AliFemtoThreeVectorD,AliFemtoThreeVectorF by templated version
 *

****************************************************************************/
#ifndef ALIFEMTO_LORENTZ_VECTOR_D_HH
#define ALIFEMTO_LORENTZ_VECTOR_D_HH
//#include "AliFemtoThreeVectorF.hh"
#include "AliFmThreeVectorD.h"
#include "AliFmLorentzVector.h"
typedef AliFmLorentzVector<double> AliFmLorentzVectorD;
#endif

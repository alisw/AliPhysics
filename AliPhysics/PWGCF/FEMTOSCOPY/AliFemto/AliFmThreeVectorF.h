/***************************************************************************
 *
 * $Id$
 *
 * Author: Thomas Ullrich, Jan 1999
 ***************************************************************************
 *
 * Description:  
 *
 * Remarks:   This is a 'handmade' specialisation of AliFmThreeVector<T>
 *            for T=float. This code contains no templates.
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.14  2005/07/06 18:49:57  fisyak
 * Replace AliFmHelixD, AliFmLorentzVectorD,AliFmLorentzVectorF,AliFmMatrixD,AliFmMatrixF,AliFmPhysicalHelixD,AliFmThreeVectorD,AliFmThreeVectorF by templated version
 *
 ****************************************************************************/
#ifndef ST_THREE_VECTOR_F_HH
#define ST_THREE_VECTOR_F_HH
#include "AliFmThreeVector.h"
typedef  AliFmThreeVector<float> AliFmThreeVectorF;
#endif

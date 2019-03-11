/***************************************************************************
 *
 * $Id$
 *
 * Author: Brian Lasiuk, Thomas Ullrich, April 1998
 ***************************************************************************
 *
 * Description:
 *
 * Remarks:   Since not all compilers support member templates
 *            we have to specialize the templated member on these
 *            platforms. If member templates are not supported the
 *            ST_NO_MEMBER_TEMPLATES flag has to be set. tu.
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
 * Revision 1.15  2005/09/22 20:09:20  fisyak
 * Make AliFmLorentzVector persistent
 *
 * Revision 1.14  2005/07/19 22:27:11  perev
 * Cleanup
 *
 * Revision 1.13  2005/07/06 18:49:57  fisyak
 * Replace AliFmHelixD, AliFmLorentzVectorD,AliFmLorentzVectorF,AliFmMatrixD,AliFmMatrixF,AliFmPhysicalHelixD,AliFmThreeVectorD,AliFmThreeVectorF by templated version
 *
 * Revision 1.12  2005/03/28 06:03:41  perev
 * Defence FPE added
 *
 * Revision 1.11  2004/12/02 20:07:32  fine
 * define the valid method for both flavor of AliFmThreeVector
 *
 * Revision 1.10  2003/10/30 20:06:46  perev
 * Check of quality added
 *
 * Revision 1.9  2003/09/02 17:59:35  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.8  2002/06/21 17:47:37  genevb
 * Added pseudoProduct
 *
 * Revision 1.7  2000/01/04 19:56:05  ullrich
 * Added cpp macro for CINT.
 *
 * Revision 1.6  1999/12/21 15:14:31  ullrich
 * Modified to cope with new compiler version on Sun (CC5.0).
 *
 * Revision 1.5  1999/10/15 15:46:54  ullrich
 * Changed output format in operator<<
 *
 * Revision 1.4  1999/06/04 18:00:05  ullrich
 * Added new constructor which takes C-style array as argument.
 * New operators operator() and operator[] which can be used
 * as lvalues.
 *
 * Revision 1.3  1999/02/17 11:42:19  ullrich
 * Removed specialization for 'long double'.
 *
 * Revision 1.2  1999/02/14 23:11:48  fisyak
 * Fixes for Rootcint
 *
 * Revision 1.1  1999/01/30 03:59:05  fisyak
 * Root Version of AliFmarClassLibrary
 *
 * Revision 1.1  1999/01/23 00:28:04  ullrich
 * Initial Revision
 *
 **************************************************************************/
/*//
//// General class for a three-vector
///*/

#include "AliFmThreeVector.h"


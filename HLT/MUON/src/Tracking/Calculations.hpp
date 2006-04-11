////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "BasicTypes.hpp"


/* Computes the Pt (transverse mementum) based on the equations given in the
   ALICE dimuon spectrometer Technical Design Report (TDR-5): trigger section.
   
   Reference: 
     "CERN/LHCC 2000-046
      Addendum 1 to ALICE TDR 5
      15 Dec 2000"
     
     Section 3.1.2 pages 144 and 145.
     
   Input can be in meters, cm or mm. 
   Output is in GeV.
 */
Float AliHLTMUONCoreCalculatePt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2
	);

/* Performs the same calculation as above however alows the zf and qBL
   parameters to be specified.
 */
Float AliHLTMUONCoreCalculatePt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		register Float zf, register Float qBL
	);

/* The same Pt calculation as above however the sign of the result indicates
   the sign of the particle.
 */
Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2
	);

/* Performs the same calculation as above however alows the zf and qBL
   parameters to be specified.
 */
Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		register Float zf, register Float qBL
	);

/* The same Pt calculation as above however the sign of the result indicates
   the sign of the particle. The momentum is also computed and returned.
 */
Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		Float& p
	);

/* Performs the same calculation as above however alows the zf and qBL
   parameters to be specified.
 */
Float AliHLTMUONCoreCalculateSignedPt(
		register Float x1,
		register Float y1, register Float y2,
		register Float z1, register Float z2,
		register Float zf, register Float qBL,
		Float& p
	);


/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers This class
// contains the base procedures for the Forward Multiplicity detector
// Detector consists of 5 Si volumes covered pseudorapidity interval
// from 1.7 to 5.1.
// 
// This contains the coarse version of the FMD - that is, the
// simulation produces no hits in the FMD volumes. 
//                                                                          
// The actual code is done by various separate classes.   Below is
// diagram showing the relationship between the various FMD classes
// that handles the geometry 
//
//
//       +----------+   +----------+   
//       | AliFMDv1 |	| AliFMDv1 |   
//       +----------+   +----------+   
//            |              |
//       +----+--------------+
//       |
//       |           +------------+ 1  +---------------+
//       |        +- | AliFMDRing |<>--| AliFMDPolygon | 
//       V     2  |  +------------+    +---------------+   
//  +--------+<>--+        |
//  | AliFMD |             ^                       
//  +--------+<>--+        V 1..2                     
//	       3  | +-------------------+ 
//	          +-| AliFMDSubDetector | 
//	  	    +-------------------+
//                           ^              
//                           |
//             +-------------+-------------+
//             |             |             |	      
//        +---------+   +---------+   +---------+
//        | AliFMD1 |   | AliFMD2 |   | AliFMD3 |
//        +---------+   +---------+   +---------+
//      
//
// See also the AliFMD class for a more detailed description of the
// various components. 
//

#include "AliFMDv0.h"		// ALIFMDV0_H

//____________________________________________________________________
ClassImp(AliFMDv0);

//___________________________________________________________________
//
// EOF
//

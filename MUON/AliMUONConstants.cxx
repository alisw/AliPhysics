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

/*
$Log$
Revision 1.7  2001/03/14 17:22:15  pcrochet
Geometry of the trigger chambers : a vertical gap of has been introduced around x=0 according fig.3.27 of the TDR (P.Dupieux)

Revision 1.6  2001/01/30 12:19:39  morsch
Update chamber positions (AdTDR version update 4/12/2000).

Revision 1.5  2000/10/18 13:26:10  morsch
New z-positions of chambers after Erice

Revision 1.4  2000/10/06 09:09:56  morsch
Outer radius of chambers adjusted to accomodate slat chambers (to be checked and updated).

Revision 1.3  2000/10/02 16:58:29  egangler
Cleaning of the code :
-> coding conventions
-> void Streamers
-> some useless includes removed or replaced by "class" statement

Revision 1.2  2000/06/27 09:46:57  morsch
kMAXZOOM global constant now in AliMUONConstants

Revision 1.1  2000/06/26 14:02:38  morsch
Add class AliMUONConstants with MUON specific constants using static memeber data and access methods.

*/

#include "AliMUONConstants.h"


ClassImp(AliMUONConstants)

Int_t   AliMUONConstants::fgNCh = 14;
Int_t   AliMUONConstants::fgNTrackingCh = 10;
Int_t   AliMUONConstants::fgNTriggerCh = 4;
Int_t   AliMUONConstants::fgNTriggerCircuit = 234;
Float_t AliMUONConstants::fgDefaultChamberZ[14] =
{533.5, 546.5, 678.5, 693.5, 964.0, 986.0, 1251.5, 1278.5, 1416.5, 1443.5,
		   1610, 1625., 1710., 1725.}; 

Float_t  AliMUONConstants::fgDmin[7] = {  36.4,  46.2,  66.0,   80.,   80., 100., 100.};    
Float_t  AliMUONConstants::fgDmax[7]  = {183., 245., 395.,  560.,  563., 850., 900.};  
Int_t   AliMUONConstants::fgMaxZoom = 20;


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
Revision 1.8  2004/03/30 14:09:22  kowal2
Changes due to the coding conventions

Revision 1.7  2003/11/24 09:48:28  kowal2
Changes to obey the coding conventions

Revision 1.6  2003/11/24 09:43:03  kowal2
Obsolete - removed

Revision 1.5  2003/09/29 11:27:39  kowal2
new classes added

Revision 1.3  2002/11/15 14:27:45  hristov
First version of the parallel TPC tracking (M.Ivanov)

Revision 1.2  2001/02/05 14:43:13  hristov
Compare() declared const

Revision 1.1  2000/10/05 16:17:27  kowal2
New class replacing AliCluster


*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber clusters objects                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTPCCluster.gif">
*/
//End_Html
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
// **************
//aaaaaaaaa
//aaaaaaaaa
//aaaaaaaaa
//aaaaaaaa
//aaaaaaaaaa

#include "AliComplexCluster.h"


ClassImp(AliComplexCluster)
  //__________________________________________________________________
  AliComplexCluster::AliComplexCluster()
                    :TObject(),
		     fX(0.),
		     fY(0.),
                     fQ(0.),
		     fSigmaX2(0.),
		     fSigmaY2(0.),
		     fSigmaXY(0.),
		     fArea(0.),
		     fMax(0.)
{
  //
  // default constructor
  //
   fTracks[0]=fTracks[1]=fTracks[2]=0;
}
//_____________________________________________________________________________
Int_t AliComplexCluster::Compare(const TObject * o) const
{
  //
  // compare two clusters according y coordinata
  AliComplexCluster *cl= (AliComplexCluster *)o;
  if (fY<cl->fY) return -1;
  if (fY==cl->fY) return 0;
  return 1;  
}

Bool_t AliComplexCluster::IsSortable() const
{
  //
  //make AliComplexCluster sortabale
  return kTRUE; 
}

ClassImp(AliTPCExactPoint)
ClassImp(AliTPCClusterPoint)
ClassImp(AliTPCTrackerPoint)
ClassImp(AliTPCTrackPoint)
ClassImp(AliTPCTrackPoint2)
ClassImp(AliTPCTrackPointRef)
  //_______________________________________________________________


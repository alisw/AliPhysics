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
Revision 1.1.4.2  2000/04/10 11:34:02  kowal2

Clusters handling in a new data structure

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

#include "AliCluster.h"

ClassImp(AliCluster)
//_____________________________________________________________________________
Int_t AliCluster::Compare(TObject * o)
{
  //
  // compare two clusters according y coordinata
  AliCluster *cl= (AliCluster *)o;
  if (fY<cl->fY) return -1;
  if (fY==cl->fY) return 0;
  return 1;  
}

Bool_t AliCluster::IsSortable() const
{
  //
  //make AliCluster sortabale
  return kTRUE; 
}

ClassImp(AliDigitCluster)
ClassImp(AliDifCluster)

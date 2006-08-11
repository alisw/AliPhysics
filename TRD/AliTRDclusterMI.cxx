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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster, alternative version                                         //
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 

#include "AliTRDclusterMI.h"
#include "AliTRDrecPoint.h"

ClassImp(AliTRDclusterMI)

//___________________________________________________________________________
AliTRDclusterMI::AliTRDclusterMI() 
  :AliTRDcluster() 
  ,fRmsY(0)
  ,fNPads(0)
  ,fRelPos(0)
{ 
  //
  // Default constructor
  //

}

//___________________________________________________________________________
AliTRDclusterMI::AliTRDclusterMI(const AliTRDcluster &c)
  :AliTRDcluster(c)
  ,fRmsY(0)
  ,fNPads(0)
  ,fRelPos(0)
{
  //
  // Copy constructor
  //

}

//_____________________________________________________________________________
AliTRDclusterMI::AliTRDclusterMI(const AliTRDrecPoint &p)
  :AliTRDcluster(p)
  ,fRmsY(0)
  ,fNPads(0)
  ,fRelPos(0)
{
  //
  // Constructor from AliTRDrecPoint
  //
  
}


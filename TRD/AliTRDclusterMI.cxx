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
//  TRD cluster                                                              //
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 

#include "AliTRDclusterMI.h"
#include "AliTRDrecPoint.h"

ClassImp(AliTRDclusterMI)

  //___________________________________________________________________________

AliTRDclusterMI::AliTRDclusterMI() : AliTRDcluster() 
{ 
  fNPads  = 0;
  fRelPos = 0;
  fRmsY   = 0;
}

AliTRDclusterMI::AliTRDclusterMI(AliTRDcluster & cl): AliTRDcluster(cl)
{
  fNPads  = 0;
  fRelPos = 0;
  fRmsY   = 0;
}
//_____________________________________________________________________________
  AliTRDclusterMI::AliTRDclusterMI(const AliTRDrecPoint &p):AliTRDcluster(p)
{
  //
  // Constructor from AliTRDrecPoint
  //
  fNPads = 0;
  fRelPos =0;
  
}


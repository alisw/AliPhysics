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
//     Base class for a detector segment                                     // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDsegmentID.h"

ClassImp(AliTRDsegmentID)

//_____________________________________________________________________________
AliTRDsegmentID::AliTRDsegmentID()
                :fSegmentID(0)
{
  //
  // AliTRDsegmentID default constructor
  //

}

//_____________________________________________________________________________
AliTRDsegmentID::AliTRDsegmentID(Int_t index)
                :fSegmentID(index)
{
  //
  // Defines a detector segment
  //

}

//_____________________________________________________________________________
AliTRDsegmentID::~AliTRDsegmentID()
{
  //
  // AliTRDsegmentID destructor
  //

}

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
Revision 1.1  2000/02/28 19:03:35  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//     Base class for a detector segment                                     // 
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDsegmentID.h"

ClassImp(AliTRDsegmentID)

//_____________________________________________________________________________
AliTRDsegmentID::AliTRDsegmentID():TObject()
{
  //
  // AliTRDsegmentID default constructor
  //

  fSegmentID = 0;

}

//_____________________________________________________________________________
AliTRDsegmentID::AliTRDsegmentID(Int_t index):TObject()
{
  //
  // Defines a detector segment
  //

  fSegmentID = index;

}

//_____________________________________________________________________________
AliTRDsegmentID::~AliTRDsegmentID()
{
  //
  // AliTRDsegmentID destructor
  //

}

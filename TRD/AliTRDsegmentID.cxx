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
Revision 1.3.12.2  2002/07/24 10:09:31  alibrary
Updating VirtualMC

Revision 1.4  2002/03/28 14:59:07  cblume
Coding conventions

Revision 1.3  2000/11/01 14:53:21  cblume
Merge with TRD-develop

Revision 1.1.4.1  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.2  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.1  2000/02/28 19:03:35  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     Base class for a detector segment                                     // 
//                                                                           //
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

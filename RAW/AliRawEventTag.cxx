/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// The AliRawEventTag class handles the raw-data event-oriented tag         //
// information. One object for each raw-data event is stored in a ROOT      //
// tree inside the file controled by AliTagDB class.                        //
// For the moment the tag information includes the raw-data event header +  //
// the raw-data file GUID and the event index.                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "AliRawEventTag.h"
#include "AliRawEventHeaderBase.h"

ClassImp(AliRawEventTag)

//______________________________________________________________________________
AliRawEventTag::AliRawEventTag() :
  fHeader(NULL),
  fGUID(""),
  fEvent(-1)
{
  // Default constructor
}

//___________________________________________________________________________
AliRawEventTag::AliRawEventTag(const AliRawEventTag & tag) :
  TObject(tag),
  fHeader(tag.fHeader),
  fGUID(tag.fGUID),
  fEvent(tag.fEvent)
{
  // copy constructor
}

//___________________________________________________________________________
AliRawEventTag & AliRawEventTag::operator=(const AliRawEventTag &tag) {
  // assignment operator
  if (this != &tag) {
    TObject::operator=(tag);

    SetHeader(tag.GetHeader());
    SetGUID(tag.GetGUID());
    SetEventNumber(tag.GetEventNumber());
  }
  return *this;
}

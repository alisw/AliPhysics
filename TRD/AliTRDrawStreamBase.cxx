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

/* $Id: AliTRDrawStreamBase.cxx 23387 2008-01-17 17:25:16Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class defines access to TRD digits in raw data.                      //
//                                                                           //
// It loops over all TRD digits in the raw data given by the AliRawReader.   //
// The Next method goes to the next digit. If there are no digits left       //
// it returns kFALSE.                                                        //
// Several getters provide information about the current digit.              //
//                                                                           //
// Author: M. Ploskon (ploskon@ikf.uni-frankfurt.de)                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"

#include "AliTRDrawOldStream.h"
#include "AliTRDrawStream.h"

#include "AliTRDrawStreamBase.h"

//--------------------------------------------------------
ClassImp(AliTRDrawStreamBase)

Int_t AliTRDrawStreamBase::fgRawStreamVersion = AliTRDrawStreamBase::kTRDrealStream;

//_____________________________________________________________________________
AliTRDrawStreamBase::AliTRDrawStreamBase()
  : TObject()
{
  //
  // this is just for API
  //
  ;
}

//_____________________________________________________________________________
AliTRDrawStreamBase::AliTRDrawStreamBase(AliRawReader */*rawReader*/)
  : TObject()
{
  //
  // this is just for API
  //
  ;
}

//_____________________________________________________________________________
AliTRDrawStreamBase::AliTRDrawStreamBase(const AliTRDrawStreamBase& /*st*/)
  : TObject()
{
  //
  // copy
  //
  TRD_NOIMP();
  ;
}

//_____________________________________________________________________________
AliTRDrawStreamBase::~AliTRDrawStreamBase()
{
  //
  // destructor
  //
  ;
}

//_____________________________________________________________________________
AliTRDrawStreamBase &
AliTRDrawStreamBase::operator=(const AliTRDrawStreamBase &)
{
  //
  // we are not using this functionality
  //
  TRD_NOIMP();
  return *this;
}

//_____________________________________________________________________________
AliTRDrawStreamBase *AliTRDrawStreamBase::GetRawStream()
{
  //
  // Returns the selected raw stream implementation
  //

  if (fgRawStreamVersion == kTRDoldStream)
    return new AliTRDrawOldStream();

  if (fgRawStreamVersion == kTRDrealStream)
    return new AliTRDrawStream();

  if (fgRawStreamVersion == kTRDsimStream)
    return new AliTRDrawStream();
  
  return new AliTRDrawStreamBase;
}

//_____________________________________________________________________________
AliTRDrawStreamBase *AliTRDrawStreamBase::GetRawStream(AliRawReader *reader)
{
  //
  // Returns the selected raw stream implementation
  //

  if (fgRawStreamVersion == kTRDoldStream)
    return new AliTRDrawOldStream(reader);

  if (fgRawStreamVersion == kTRDrealStream)
    return new AliTRDrawStream(reader);

  if (fgRawStreamVersion == kTRDsimStream)
    return new AliTRDrawStream(reader);

  return new AliTRDrawStreamBase;
}

//_____________________________________________________________________________
void AliTRDrawStreamBase::SetRawStreamVersion(const char *opt) 
{ 
  //
  // Sets the raw stream version
  //

  fgRawStreamVersion = 0; 

  if (strstr(opt, "sim" ) != 0 || strstr(opt, "SIM") != 0) 
    fgRawStreamVersion = kTRDsimStream; 

  if (strstr(opt, "tb" ) != 0 || strstr(opt, "TB") != 0) 
    fgRawStreamVersion = kTRDrealStream; 

  if (strstr(opt, "real" ) != 0 || strstr(opt, "REAL") != 0) 
    fgRawStreamVersion = kTRDrealStream; 

  if (strstr(opt, "old" ) != 0 || strstr(opt, "OLD") != 0) 
    fgRawStreamVersion = kTRDoldStream; 
      
}

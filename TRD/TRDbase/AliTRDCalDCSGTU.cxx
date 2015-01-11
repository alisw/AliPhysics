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

/* $Id: AliTRDCalDCSGTU.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS GTU parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSGTU.h"
#include "AliTRDCalDCSGTUTgu.h"
#include <TObjArray.h>

ClassImp(AliTRDCalDCSGTU)

//_____________________________________________________________________________
AliTRDCalDCSGTU::AliTRDCalDCSGTU()
  :TNamed()
    ,fRunNumber(0)
    ,fSORFlag(0)
    ,fSerial(0)
    ,fDNR(-1)
    ,fSegmentsArr(new TObjArray())
    ,fTgu(new AliTRDCalDCSGTUTgu())
{
  //
  // AliTRDCalDCSGTU default constructor
  //
  fSegmentsArr->SetOwner();
}

//_____________________________________________________________________________
AliTRDCalDCSGTU::AliTRDCalDCSGTU(const char *name, const char *title)
  :TNamed(name,title)
    ,fRunNumber(0)
    ,fSORFlag(0)
    ,fSerial(0)
    ,fDNR(-1)
    ,fSegmentsArr(new TObjArray())
    ,fTgu(new AliTRDCalDCSGTUTgu())
{
  //
  // AliTRDCalDCSGTU constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSGTU::AliTRDCalDCSGTU(const AliTRDCalDCSGTU&)
  :TNamed("","")
    ,fRunNumber(0)
    ,fSORFlag(0)
    ,fSerial(0)
    ,fDNR(-1)
    ,fSegmentsArr(new TObjArray())
    ,fTgu(new AliTRDCalDCSGTUTgu())
{
  //
  // AliTRDCalDCSGTU constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSGTU::~AliTRDCalDCSGTU()
{
  //
  // AliTRDCalDCSGTU destructor
  //

  if (fSegmentsArr) {
    fSegmentsArr->Delete();
    delete fSegmentsArr;
    fSegmentsArr = 0x0;
  }

  if (fTgu) {
    delete fTgu;
    fTgu = 0x0;
  }

}

//_____________________________________________________________________________
AliTRDCalDCSGTU& AliTRDCalDCSGTU::operator=(const AliTRDCalDCSGTU& sh)
{
  //
  // AliTRDCalDCSGTU constructor
  //
  if (&sh == this) return *this;
  
  new (this) AliTRDCalDCSGTU(sh);
  return *this;
}




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

/* $Id: AliTRDCalDCSGTUTgu.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS GTU parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSGTUTgu.h"
#include <TObjArray.h>
#include "AliTRDCalDCSGTUBoardInfo.h"

ClassImp(AliTRDCalDCSGTUTgu)

//_____________________________________________________________________________
AliTRDCalDCSGTUTgu::AliTRDCalDCSGTUTgu()
  :TNamed()
    ,fFromRunNum(0)
    ,fFromSORFlag(0)
    ,fFromChild(0)
    ,fSegmentMask("")
    ,fBusyMask("")
    ,fContribMask("")
    ,fBoardInfo(new AliTRDCalDCSGTUBoardInfo())
    ,fCtpOpcArr(new TObjArray())
{
  //
  // AliTRDCalDCSGTU default constructor
  //
  fCtpOpcArr->SetOwner();
}

//_____________________________________________________________________________
AliTRDCalDCSGTUTgu::AliTRDCalDCSGTUTgu(const char *name, const char *title)
  :TNamed(name,title)
    ,fFromRunNum(0)
    ,fFromSORFlag(0)
    ,fFromChild(0)
    ,fSegmentMask("")
    ,fBusyMask("")
    ,fContribMask("")
    ,fBoardInfo(new AliTRDCalDCSGTUBoardInfo())
    ,fCtpOpcArr(new TObjArray())
{
  //
  // AliTRDCalDCSGTU constructor
  //
}

//_____________________________________________________________________________
AliTRDCalDCSGTUTgu::AliTRDCalDCSGTUTgu(const AliTRDCalDCSGTUTgu&)
  :TNamed("","")
    ,fFromRunNum(0)
    ,fFromSORFlag(0)
    ,fFromChild(0)
    ,fSegmentMask("")
    ,fBusyMask("")
    ,fContribMask("")
    ,fBoardInfo(0)
    ,fCtpOpcArr(0)
{
  //
  // AliTRDCalDCSGTU constructor
  //

}

//_____________________________________________________________________________
AliTRDCalDCSGTUTgu::~AliTRDCalDCSGTUTgu()
{
  //
  // AliTRDCalDCSGTU destructor
  //

  if (fBoardInfo) {
    delete fBoardInfo;
    fBoardInfo = 0x0;
  }

  if (fCtpOpcArr) {
    fCtpOpcArr->Delete();
    delete fCtpOpcArr;
    fCtpOpcArr = 0x0;
  }

}

//_____________________________________________________________________________
AliTRDCalDCSGTUTgu& AliTRDCalDCSGTUTgu::operator=(const AliTRDCalDCSGTUTgu& sh)
{
  //
  // AliTRDCalDCSGTU constructor
  //
  if (&sh == this) return *this;

  new (this) AliTRDCalDCSGTUTgu(sh);
  return *this;
}




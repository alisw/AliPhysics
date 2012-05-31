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

/* $Id:  $ */

///////////////////////////////
//                           //
// AliTOFDCSmaps class       //
// container for HV&&LV maps //
// as found during a run     //
//                           //
// Author: A. De Caro        //
// Email: decaro@sa.innf.it  //
//                           //
///////////////////////////////

#include "AliTOFDCSmaps.h"

ClassImp(AliTOFDCSmaps)

////////////////////////////////////////////////////////////////////////
AliTOFDCSmaps::AliTOFDCSmaps() :
  TObject(),
  fTime(0)
{
  //
  // default ctor
  //
  for (Int_t ii=0; ii<91*96*18; ii++) fArray[ii]=-1;

}

////////////////////////////////////////////////////////////////////////
AliTOFDCSmaps::AliTOFDCSmaps(Int_t time, Short_t array[]) :
  TObject(),
  fTime(time)
{
  //
  // ctor
  //
  for (Int_t ii=0; ii<91*96*18; ii++) fArray[ii]=array[ii];
}

////////////////////////////////////////////////////////////////////////
AliTOFDCSmaps::AliTOFDCSmaps(const AliTOFDCSmaps & digit):
  TObject(digit),
  fTime(digit.fTime)
{
  // 
  // copy ctor
  //
  for (Int_t ii=0; ii<91*96*18; ii++) fArray[ii]=digit.fArray[ii];

}

////////////////////////////////////////////////////////////////////////
AliTOFDCSmaps& AliTOFDCSmaps::operator=(const AliTOFDCSmaps & digit)
{
  // 
  // operator =
  //

  if (this == &digit)
    return *this;

  TObject::operator=(digit);

  fTime = digit.fTime;
  for (Int_t ii=0; ii<91*96*18; ii++) fArray[ii]=digit.fArray[ii];
  return *this;

}

////////////////////////////////////////////////////////////////////////
AliTOFDCSmaps::~AliTOFDCSmaps()
{
  //
  // dtor
  //
}

////////////////////////////////////////////////////////////////////////
void AliTOFDCSmaps::Update(AliTOFDCSmaps *object)
{
  //
  // Update aready exsisting AliTOFDCSmap
  // with the value of the passed object
  //

  Short_t value = -1;

  for (Int_t ii=0; ii<91*96*18; ii++) {
    value = object->GetCellValue(ii);
    if (fArray[ii]==-1 && value!=-1)
      fArray[ii] = value;
  }

}


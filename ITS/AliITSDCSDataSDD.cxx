/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class containing SDD DCS data           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSDCSDataSDD.h"
#include "AliLog.h"

ClassImp(AliITSDCSDataSDD)

//---------------------------------------------------------------
AliITSDCSDataSDD::AliITSDCSDataSDD():
TObject(),
fNpts(0),
fSetPoints(0),
fTimeStamp(),
fDriftField(),
fTemperatureLeft(),
fTemperatureRight()
{
  // default constructor
}
//---------------------------------------------------------------
AliITSDCSDataSDD::AliITSDCSDataSDD(Int_t npts):
TObject(),
fNpts(npts),
fSetPoints(0),
fTimeStamp(npts),
fDriftField(npts),
fTemperatureLeft(npts),
fTemperatureRight(npts)
{
  // standard constructor
}
//---------------------------------------------------------------
void AliITSDCSDataSDD::SetNPoints(Int_t npts){
  // Redimension and resets arrays
  fNpts=npts;
  fTimeStamp.Set(npts);
  fDriftField.Set(npts);
  fTemperatureLeft.Set(npts);
  fTemperatureRight.Set(npts);
}

//---------------------------------------------------------------
void AliITSDCSDataSDD::SetValues(Int_t time, Float_t field, Float_t templ, Float_t tempr){
  // Set values of next array elements
  if(fSetPoints>=fNpts){
    AliWarning("Try to write outside array range");
    return;
  }
  fTimeStamp.AddAt(time,fSetPoints);
  fDriftField.AddAt(field,fSetPoints);
  fTemperatureLeft.AddAt(templ,fSetPoints);
  fTemperatureRight.AddAt(tempr,fSetPoints);
  fSetPoints++;
}
//---------------------------------------------------------------
void AliITSDCSDataSDD::Compress(){
  // Redimension arrays removing the empty elements at the end
  SetNPoints(fSetPoints);
}
//---------------------------------------------------------------
void AliITSDCSDataSDD::PrintValues() const {
  // Printout array contents
  for(Int_t i=0;i<fNpts;i++){
    printf("TimeStamp=%d   Drift Field=%f  Temperatures: Left Hybr=%f  Right Hybr=%f\n",fTimeStamp.At(i),fDriftField.At(i),fTemperatureLeft.At(i),fTemperatureRight.At(i));
  }
}

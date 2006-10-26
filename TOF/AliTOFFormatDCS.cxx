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
*/  

#include "AliTOFFormatDCS.h"

// AliTOFFormatDCS class
// contining the format of a DCS calibration objects
// consisting in 10 floats
// (5 data points + 5 corresponding timestamps) 
// and a short_integer

ClassImp(AliTOFFormatDCS)

//---------------------------------------------------------------
AliTOFFormatDCS::AliTOFFormatDCS(): TObject(),fShort(0){
  // main constructor
  for (Int_t i=0;i<3;i++){
    fFloats[i]=0;
    fTimeStampsFloat[i]=0;
  }
  for (Int_t i=0;i<2;i++){
    fDeltas[i]=0;
    fTimeStampsDelta[i]=0;
  }
}
//---------------------------------------------------------------

AliTOFFormatDCS::AliTOFFormatDCS(const AliTOFFormatDCS & format):
  TObject(),
  fShort(format.fShort)
{ 
  // copy constructor
  for (Int_t i=0;i<3;i++){
    this->fFloats[i]=format.GetFloat(i);
    this->fTimeStampsFloat[i]=format.GetTimeStampFloat(i);
  }
  for (Int_t i=0;i<2;i++){
    this->fDeltas[i]=format.GetDelta(i);
    this->fTimeStampsDelta[i]=format.GetTimeStampFloat(i);
  }
}
//---------------------------------------------------------------

AliTOFFormatDCS& AliTOFFormatDCS:: operator=(const AliTOFFormatDCS & format) { 

  // assignment operator
  for (Int_t i=0;i<3;i++){
    this->fFloats[i]=format.GetFloat(i);
    this->fTimeStampsFloat[i]=format.GetTimeStampFloat(i);
  }
  for (Int_t i=0;i<2;i++){
    this->fDeltas[i]=format.GetDelta(i);
    this->fTimeStampsDelta[i]=format.GetTimeStampFloat(i);
  }
  this->fShort=format.GetShort();
  return *this;
}
//---------------------------------------------------------------

AliTOFFormatDCS::~AliTOFFormatDCS(){

}


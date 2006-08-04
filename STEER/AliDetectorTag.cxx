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

//-----------------------------------------------------------------
//           Implementation of the DetectorTag class
//   This is the class to deal with the tags in the detector level
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "AliDetectorTag.h"

ClassImp(AliDetectorTag)

//___________________________________________________________________________
  AliDetectorTag::AliDetectorTag() :
    TObject(),
    fITS(kFALSE),
    fTPC(kFALSE),
    fTRD(kFALSE),
    fTOF(kFALSE),
    fHMPID(kFALSE),
    fPHOS(kFALSE),
    fZDC(kFALSE),
    fMUON(kFALSE),
    fPMD(kFALSE),
    fEMCAL(kFALSE),
    fVZERO(kFALSE),
    fTZERO(kFALSE)
{
  // Default constructor
}

//___________________________________________________________________________
AliDetectorTag::AliDetectorTag(const AliDetectorTag & detTag) : TObject(detTag) {
  // DetectorTag copy constructor
  SetITS();
  SetTPC();
  SetTRD();
  SetTOF();
  SetHMPID();
  SetPHOS();
  SetZDC();
  SetMUON();
  SetPMD();
  SetEMCAL();
  SetVZERO();
  SetTZERO();
}

//___________________________________________________________________________
AliDetectorTag & AliDetectorTag::operator=(const AliDetectorTag &detTag) {
  //DetectorTag assignment operator
  if (this != &detTag) {
    TObject::operator=(detTag);
    
    SetITS();
    SetTPC();
    SetTRD();
    SetTOF();
    SetHMPID();
    SetPHOS();
    SetZDC();
    SetMUON();
    SetPMD();
    SetEMCAL();
    SetVZERO();
    SetTZERO();
  }
  return *this;
}


//___________________________________________________________________________
AliDetectorTag::~AliDetectorTag() {
  // Destructor
}

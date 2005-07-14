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

//______________________________________________________________________________
AliDetectorTag::AliDetectorTag()
{
  fITS = 0;
  fTPC = 0;
  fTRD = 0;
  fTOF = 0;
  fHMPID = 0;
  fPHOS = 0;
  fZDC = 0;
  fMUON = 0;
  fABSORBER = 0;
  fPMD = 0;
  fRICH = 0;
  fEMCAL = 0;
  fVZERO = 0;
  fTZERO = 0;
}

//______________________________________________________________________________
AliDetectorTag::AliDetectorTag(AliDetectorTag *DetTag)
{
  // DetectorTag copy constructor
  SetITS(DetTag->GetITS());
  SetTPC(DetTag->GetTPC());
  SetTRD(DetTag->GetTRD());
  SetTOF(DetTag->GetTOF());
  SetHMPID(DetTag->GetHMPID());
  SetPHOS(DetTag->GetPHOS());
  SetZDC(DetTag->GetZDC());
  SetMUON(DetTag->GetMUON());
  SetABSORBER(DetTag->GetABSORBER());
  SetPMD(DetTag->GetPMD());
  SetRICH(DetTag->GetRICH());
  SetEMCAL(DetTag->GetEMCAL());
  SetVZERO(DetTag->GetVZERO());
  SetTZERO(DetTag->GetTZERO());
}

//______________________________________________________________________________
AliDetectorTag::~AliDetectorTag()
{
}

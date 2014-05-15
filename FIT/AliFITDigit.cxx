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

/////////////////////////////////////////////////////////////////////////
//  Class AliFITDigits for FIT digits
//  fTimeCFD, fTimeLED   for PMT fNPMT
// amplitude start fTimeQT0 and stop fTimeQT1
///////////////////////////////////////////////////////////////////////

#include "AliFITDigit.h"
#include "AliFIT.h"

ClassImp(AliFITDigit)

//-----------------------------------------------
  AliFITDigit::AliFITDigit() :AliDigit(),
  fTimeCFD(-1),   
  fTimeLED(-1),   
  fTimeQT0(-1),   
  fTimeQT1(-1),
  fNPMT(-1)
{
  // default contructor    
}

//_____________________________________________________________________________
AliFITDigit::AliFITDigit(Int_t npmt,  
			 Int_t timeCFD, Int_t timeLED, Int_t timeQT0, 
			 Int_t timeQT1 , Int_t *labels):
  AliDigit(),
  fTimeCFD(timeCFD),   
  fTimeLED(timeLED),   
  fTimeQT0(timeQT0),   
  fTimeQT1(timeQT1),
  fNPMT(npmt)
{
  // Constructor
  // Used in the digitizer
 if (labels)
    for(Int_t iTrack = 0; iTrack < 3; ++iTrack) fTracks[iTrack] = labels[iTrack];
}

//_____________________________________________________________________________
AliFITDigit::~AliFITDigit() {
  // destructor
 
}

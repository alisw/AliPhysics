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
//           Implementation of the ITS PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include <TMath.h>

#include "AliITSpidESD.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"

ClassImp(AliITSpidESD)


//______________________________________________________________________
AliITSpidESD::AliITSpidESD():TObject(){
  //Default constructor
 
}

Double_t AliITSpidESD::Bethe(Double_t p,Double_t mass) {
  // returns AliExternalTrackParam::BetheBloch normalized to 1 at the minimum
  Double_t density=2.33; // g/cm3
  Double_t thickness=0.03; // cm
  Double_t meanMIPSi=116.24; // keV in 300 microns of Si
  Double_t conv=density*1E6*thickness/meanMIPSi;
  Float_t betagamma=p/mass;
  return conv*AliExternalTrackParam::BetheBlochSolid(betagamma);

}

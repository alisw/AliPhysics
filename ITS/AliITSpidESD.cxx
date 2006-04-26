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

//-----------------------------------------------------------------
//           Implementation of the ITS PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include <TMath.h>

#include "AliITSpidESD.h"
#include "AliESDtrack.h"

ClassImp(AliITSpidESD)


//______________________________________________________________________
AliITSpidESD::AliITSpidESD():TObject(){
  //Default constructor
 
}
Double_t AliITSpidESD::Bethe(Double_t bg) {
  //
  // This is the Bethe-Bloch function normalised to 1 at the minimum
  //
  Double_t bg2=bg*bg;
  Double_t bethe;
  /*if (bg<3.5) 
     bethe=(1.+ bg2)/bg2*(log(5940*bg2) - bg2/(1.+ bg2));
  else*/  // not 100% clear why...
  bethe=(1.+ bg2)/bg2*(TMath::Log(3.5*5940*bg) - bg2/(1.+ bg2));
  return bethe/11.091;
}


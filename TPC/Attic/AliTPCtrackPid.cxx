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

//
// Class to determine the particle ID for TPC tracks
// More comments and a description of the class will be added here
// missing description
// as above
// as above

#include "AliTPCtrackPid.h"

ClassImp(AliTPCtrackPid)

AliTPCtrackPid::AliTPCtrackPid()
  :TObject(),
     fWpi(0.),     
     fWk(0.),      
     fWp(0.),      
     fSignal(0.),  
     fMom(0.),     
     fPhi(0.),     
     fLam(0.),     
     fPcode(0),   
     fLabel(0),
     fGSignal(0.), 
     fGMom(0.),    
     fGpx(0.),     
     fGpy(0.),     
     fGpz(0.),     
     fGx(0.),      
     fGy(0.),      
     fGz(0.),      
     fGcode(0),   
     fGlab(0)       
{
  //
  // default costructor
  //
}

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
// Event header for GeVSim event generator
// support event plane and elliptic flow
// in next release will suport full differential 
// directed and elliptic flow
//
// Sylwester Radomski, GSI
// mail: S.Radomski@gsi
// 31 Oct, 2002
//
//

#include "AliGenGeVSimEventHeader.h"


ClassImp(AliGenGeVSimEventHeader)


////////////////////////////////////////////////////////////////////////////////

AliGenGeVSimEventHeader::AliGenGeVSimEventHeader()
  :AliGenEventHeader() {
  //
  // Default Constructor 
  //

  fEventPlane = 0;
  fEllipticFlow = 0;
}


////////////////////////////////////////////////////////////////////////////////

AliGenGeVSimEventHeader::AliGenGeVSimEventHeader(const char *name)
  :AliGenEventHeader(name) {
  //
  // Standard constructor
  //
  
  fEventPlane = 0;
  fEllipticFlow = 0;
}

////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSimEventHeader::SetEventPlane(Float_t psi) {
  //
  // Psi in rad.
  //

  fEventPlane = psi;
}

////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSimEventHeader::SetEllipticFlow(Float_t v2) {
  //
  // Set elliptic flow
  //

  fEllipticFlow = v2;
}

////////////////////////////////////////////////////////////////////////////////

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


ClassImp(AliGenGeVSimEventHeader);


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

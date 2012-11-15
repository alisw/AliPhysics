#ifndef ALIITSUAUX
#define ALIITSUAUX

#include <TObject.h>
#include <TMath.h>

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Namespace AliITSUAux                                             //
//  Set of utilities for the ITSU classes                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////


class AliITSUGeomTGeo;
class AliITSsegmentation;
using namespace TMath;


namespace AliITSUAux {
  void   BringTo02Pi(double &phi);
  Bool_t OKforPhiMin(double phiMin,double phi);
  Bool_t OKforPhiMax(double phiMax,double phi);
  //
  const Double_t kNominalBz = 5.01;           // nominal field
  const Double_t kPionMass = 1.3957e-01;
}

//_________________________________________________________________________________
inline void AliITSUAux::BringTo02Pi(double &phi) {   
  // bring phi to 0-2pi range
  if (phi<0) phi+=TwoPi(); else if (phi>TwoPi()) phi-=TwoPi();
}

//_________________________________________________________________________________
inline Bool_t AliITSUAux::OKforPhiMin(double phiMin,double phi) {
  // check if phi is above the phiMin, phi's must be in 0-2pi range
  double dphi = phi-phiMin;
  return ((dphi>0 && dphi<Pi()) || dphi<-Pi()) ? kTRUE:kFALSE;
}

//_________________________________________________________________________________
inline Bool_t AliITSUAux::OKforPhiMax(double phiMax,double phi) {
  // check if phi is below the phiMax, phi's must be in 0-2pi range
  double dphi = phi-phiMax;
  return ((dphi<0 && dphi>-Pi()) || dphi>Pi()) ? kTRUE:kFALSE;
}


#endif

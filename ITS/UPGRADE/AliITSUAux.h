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
  UInt_t PackCluster(Int_t lr, Int_t clID);
  Int_t  UnPackCluster(UInt_t p, Int_t &lr);
  Bool_t IsCluster(UInt_t p);
  //
  const Double_t kNominalBz = 5.01;           // nominal field
  const Double_t kPionMass  = 1.3957e-01;
  const UInt_t   kLrBitMax  = 5;                            // layer mask highest bit
  const UInt_t   kMaxLayers = UInt_t(Power(2.,int(kLrBitMax)-1));  // max number of active layers
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

//_________________________________________________________________________________
inline UInt_t AliITSUAux::PackCluster(Int_t lr, Int_t clID) {
  // pack layer/cluster into single uint
  UInt_t p = clID<0 ? 0 : clID+1;
  p<<=kLrBitMax;
  return p + lr;
}

//_________________________________________________________________________________
inline Int_t AliITSUAux::UnPackCluster(UInt_t p, Int_t &lr) {
  // unpack layer/cluster
  lr = p&kMaxLayers;
  p>>=kLrBitMax;
  return int(p)-1;
}

//_________________________________________________________________________________
inline Bool_t AliITSUAux::IsCluster(UInt_t p) {
  // does it correspond to cluster?
  return p>kMaxLayers;
}


#endif

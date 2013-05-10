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

#define _ITSU_TUNING_MODE_

class AliITSUGeomTGeo;
class AliITSsegmentation;
using namespace TMath;



namespace AliITSUAux {
  void   BringTo02Pi(double &phi);
  Bool_t OKforPhiMin(double phiMin,double phi);
  Bool_t OKforPhiMax(double phiMax,double phi);
  UInt_t PackCluster(Int_t lr, Int_t clID);
  Int_t  UnpackCluster(UInt_t p, Int_t &lr);
  Int_t  UnpackLayer(UInt_t p);
  Int_t  UnpackCluster(UInt_t p);
  Bool_t IsCluster(UInt_t p);
  Int_t  NumberOfBitsSet(UInt_t x);
  void   PrintBits(ULong64_t patt, Int_t maxBits);
  //
  const Double_t kNominalBz = 5.01;           // nominal field
  const Double_t kPionMass  = 1.3957e-01;
  const UInt_t   kLrBitLow  = 28;             // layer mask lowest bit
  const UInt_t   kLrMask    = 0xf0000000;     // layer mask
  const UInt_t   kClMask    = 0x0fffffff;     // cluster mask
  const UInt_t   kMaxLayers = 15;             // max number of active layers
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
  UInt_t p = (clID<0 ? 0 : clID+1) + (lr<<=kLrBitLow);
  return p;
}

//_________________________________________________________________________________
inline Int_t AliITSUAux::UnpackCluster(UInt_t p, Int_t &lr) {
  // unpack layer/cluster
  lr = (p&kLrMask)>>kLrBitLow;
  p &= kClMask;
  return int(p)-1;
}

//_________________________________________________________________________________
inline Int_t AliITSUAux::UnpackLayer(UInt_t p) {
  // unpack layer
  return (p&kLrMask)>>kLrBitLow;
}

//_________________________________________________________________________________
inline Int_t AliITSUAux::UnpackCluster(UInt_t p) {
  // unpack cluster
  return int(p&kClMask)-1;
}

//_________________________________________________________________________________
inline Bool_t AliITSUAux::IsCluster(UInt_t p) {
  // does it correspond to cluster?
  return (p&kClMask);
}

//_________________________________________________________________________________
inline Int_t AliITSUAux::NumberOfBitsSet(UInt_t x) {
  // count number of non-0 bits in 32bit word
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}


#endif

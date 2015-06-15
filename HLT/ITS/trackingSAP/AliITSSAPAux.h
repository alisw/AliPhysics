#ifndef ALIITSSAPAUX_H
#define ALIITSSAPAUX_H

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Namespace AliITSSAPAux                                             //
//  Set of utilities for the XXX classes                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////


namespace AliITSSAPAux {
  void   BringTo02Pi(double &phi);
  void   BringTo02Pi(float  &phi);
  bool   OKforPhiMin(double phiMin,double phi);
  bool   OKforPhiMax(double phiMax,double phi);
  double MeanPhiSmall(double phi0, double phi1);
  double DeltaPhiSmall(double phi0, double phi1);
  //
  bool   OKforPhiMin(float phiMin,float phi);
  bool   OKforPhiMax(float phiMax,float phi);
  float  MeanPhiSmall(float phi0, float phi1);
  float  DeltaPhiSmall(float phi0, float phi1);

  unsigned int PackCluster(int lr, int clID);
  int    UnpackCluster(unsigned int p, int &lr);
  int    UnpackLayer(unsigned int p);
  int    UnpackCluster(unsigned int p);
  bool   IsCluster(unsigned int p);
  int    NumberOfBitsSet(unsigned int x);
  void   PrintBits(unsigned long long patt, int maxBits);
  //
  const double   kNominalBz = 5.01;           // nominal field
  const double   kPionMass  = 1.3957e-01;
  const double   kPi  = 3.14159265358979312e+00;
  const double   k2Pi = 2*kPi;
  const unsigned int   kLrBitLow  = 28;             // layer mask lowest bit
  const unsigned int   kLrMask    = 0xf0000000;     // layer mask
  const unsigned int   kClMask    = 0x0fffffff;     // cluster mask
  const unsigned int   kMaxLayers = 15;             // max number of active layers
  const unsigned int   kMaxLrMask = 0x7fff;         // bitmask for allowed layers
}

//_________________________________________________________________________________
inline void AliITSSAPAux::BringTo02Pi(double &phi) {   
  // bring phi to 0-2pi range
  while (phi<0) phi+=k2Pi; 
  while (phi>k2Pi) phi-=k2Pi;
}

//_________________________________________________________________________________
inline void AliITSSAPAux::BringTo02Pi(float &phi) {   
  // bring phi to 0-2pi range
  while (phi<0) phi+=k2Pi; 
  while (phi>k2Pi) phi-=k2Pi;
}

//_________________________________________________________________________________
inline bool AliITSSAPAux::OKforPhiMin(double phiMin,double phi) {
  // check if phi is above the phiMin, phi's must be in 0-2pi range
  double dphi = phi-phiMin;
  return ((dphi>0 && dphi<kPi) || dphi<-kPi) ? true:false;
}

//_________________________________________________________________________________
inline bool AliITSSAPAux::OKforPhiMin(float phiMin,float phi) {
  // check if phi is above the phiMin, phi's must be in 0-2pi range
  float dphi = phi-phiMin;
  return ((dphi>0 && dphi<kPi) || dphi<-kPi) ? true:false;
}

//_________________________________________________________________________________
inline bool AliITSSAPAux::OKforPhiMax(double phiMax,double phi) {
  // check if phi is below the phiMax, phi's must be in 0-2pi range
  double dphi = phi-phiMax;
  return ((dphi<0 && dphi>-kPi) || dphi>kPi) ? true:false;
}

//_________________________________________________________________________________
inline bool AliITSSAPAux::OKforPhiMax(float phiMax,float phi) {
  // check if phi is below the phiMax, phi's must be in 0-2pi range
  float dphi = phi-phiMax;
  return ((dphi<0 && dphi>-kPi) || dphi>kPi) ? true:false;
}

//_________________________________________________________________________________
inline unsigned int AliITSSAPAux::PackCluster(int lr, int clID) {
  // pack layer/cluster into single uint
  unsigned int p = (clID<0 ? 0 : clID+1) + (lr<<=kLrBitLow);
  return p;
}

//_________________________________________________________________________________
inline int AliITSSAPAux::UnpackCluster(unsigned int p, int &lr) {
  // unpack layer/cluster
  lr = (p&kLrMask)>>kLrBitLow;
  p &= kClMask;
  return int(p)-1;
}

//_________________________________________________________________________________
inline int AliITSSAPAux::UnpackLayer(unsigned int p) {
  // unpack layer
  return (p&kLrMask)>>kLrBitLow;
}

//_________________________________________________________________________________
inline int AliITSSAPAux::UnpackCluster(unsigned int p) {
  // unpack cluster
  return int(p&kClMask)-1;
}

//_________________________________________________________________________________
inline bool AliITSSAPAux::IsCluster(unsigned int p) {
  // does it correspond to cluster?
  return (p&kClMask);
}

//_________________________________________________________________________________
inline int AliITSSAPAux::NumberOfBitsSet(unsigned int x) {
  // count number of non-0 bits in 32bit word
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//_________________________________________________________________________________
inline double AliITSSAPAux::MeanPhiSmall(double phi0, double phi1) {
  // return mean phi, assume phis in 0:2pi
  double phi;
  if (!OKforPhiMin(phi0,phi1)) {phi=phi0; phi0=phi1; phi1=phi;}
  if (phi0>phi1) phi = (phi1 - (k2Pi-phi0))/2; // wrap
  else           phi = (phi0+phi1)/2;
  BringTo02Pi(phi);
  return phi;
}

//_________________________________________________________________________________
inline float AliITSSAPAux::MeanPhiSmall(float phi0, float phi1) {
  // return mean phi, assume phis in 0:2pi
  float phi;
  if (!OKforPhiMin(phi0,phi1)) {phi=phi0; phi0=phi1; phi1=phi;}
  if (phi0>phi1) phi = (phi1 - (k2Pi-phi0))/2; // wrap
  else           phi = (phi0+phi1)/2;
  BringTo02Pi(phi);
  return phi;
}

//_________________________________________________________________________________
inline double AliITSSAPAux::DeltaPhiSmall(double phi0, double phi1) {
  // return delta phi, assume phis in 0:2pi
  double del;
  if (!OKforPhiMin(phi0,phi1)) {del=phi0; phi0=phi1; phi1=del;}
  del = phi1 - phi0;
  if (del<0) del += k2Pi;
  return del;
}

//_________________________________________________________________________________
inline float AliITSSAPAux::DeltaPhiSmall(float phi0, float phi1) {
  // return delta phi, assume phis in 0:2pi
  float del;
  if (!OKforPhiMin(phi0,phi1)) {del=phi0; phi0=phi1; phi1=del;}
  del = phi1 - phi0;
  if (del<0) del += k2Pi;
  return del;
}


#endif

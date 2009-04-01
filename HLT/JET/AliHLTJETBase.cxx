//-*- Mode: C++ -*-
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTJETBase.cxx
    @author Jochen Thaeder
    @date   
    @brief  Base functionality for HLT JET package
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliHLTJETBase.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETBase)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
  AliHLTJETBase::AliHLTJETBase() {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


//##################################################################################
AliHLTJETBase::~AliHLTJETBase() {
  // see header file for class documentation
}

#if 0


//##################################################################################
void AliHLTJETBase::XYZtoRPhiEta( const Double_t *xyz, Double_t *rpe ) {
  // see header file for class documentation

  // Spherical [r,phi,eta] - rpe
  // Cartesian [x,y,z] - xyz
  
  rpe[0] = TMath::Sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
  rpe[1] = AliHLTJETBase::GetPhiFromXYZ( xyz );
  rpe[2] = 0.5 * TMath::Log( (rpe[0]+xyz[2])/(rpe[0]-xyz[2]) );
}

//##################################################################################
void AliHLTJETBase::XYZEtoRPhiEtaPt( const Double_t *xyze, Double_t *rpep ) {
  // see header file for class documentation

  // Spherical [r,phi,eta,pt] - rpep
  // Cartesian [x,y,z,E] - xyze
  
  rpep[0] = TMath::Sqrt( xyze[0]*xyze[0] + xyze[1]*xyze[1] + xyze[2]*xyze[2] );
  rpep[1] = AliHLTJETBase::GetPhiFromXYZ( xyze );
  rpep[2] = 0.5 * TMath::Log( (rpep[0]+xyze[2])/(rpep[0]-xyze[2]) );
  rpep[3] = AliHLTJETBase::GetPtFromXYZ( xyze );
}

//##################################################################################
void AliHLTJETBase::XYZEtoRPhiEtaPt( const Float_t *xyze, Double_t *rpep ) {
  // see header file for class documentation

  // Spherical [r,phi,eta,pt] - rpep
  // Cartesian [x,y,z,E] - xyze
  Double_t x = (Double_t) xyze[0];
  Double_t y = (Double_t) xyze[1];
  Double_t z = (Double_t) xyze[2];
  
  rpep[0] = TMath::Sqrt( x*x + y*y + z*z );
  rpep[1] = AliHLTJETBase::GetPhiFromXYZ( xyze );
  rpep[2] = 0.5 * TMath::Log( (rpep[0]+z)/(rpep[0]-z) );
  rpep[3] = AliHLTJETBase::GetPtFromXYZ( xyze );
}

//##################################################################################
void AliHLTJETBase::XYZEtoEPhiEtaPt( const Double_t *xyze, Double_t *epep ) {
  // see header file for class documentation

  // Spherical [r,phi,eta,pt] - rpep
  // Cartesian [x,y,z,E] - xyze

  Double_t r = TMath::Sqrt( xyze[0]*xyze[0] + xyze[1]*xyze[1] + xyze[2]*xyze[2] );
  
  epep[0] = xyze[3];
  epep[1] = AliHLTJETBase::GetPhiFromXYZ( xyze );
  epep[2] = 0.5 * TMath::Log( (r+xyze[2])/(r-xyze[2]) );
  epep[3] = AliHLTJETBase::GetPtFromXYZ( xyze );
}

//##################################################################################
void AliHLTJETBase::XYZEtoEPhiEtaPt( const Float_t *xyze, Double_t *epep ) {
  // see header file for class documentation

  // Spherical [r,phi,eta,pt] - rpep
  // Cartesian [x,y,z,E] - xyze
  Double_t x = (Double_t) xyze[0];
  Double_t y = (Double_t) xyze[1];
  Double_t z = (Double_t) xyze[2];
  Double_t e = (Double_t) xyze[3];
  Double_t r = TMath::Sqrt( x*x + y*y + z*z );
  
  epep[0] = e;
  epep[1] = AliHLTJETBase::GetPhiFromXYZ( xyze );
  epep[2] = 0.5 * TMath::Log( (r+z)/(r-z) );
  epep[3] = AliHLTJETBase::GetPtFromXYZ( xyze );
}

//##################################################################################
Double_t AliHLTJETBase::GetPtFromXYZ( const Double_t *pxpypz ) {
  // see header file for class documentation
  // Cartesian [px,py,pz] 
  
  return TMath::Sqrt( pxpypz[0]*pxpypz[0] + pxpypz[1]*pxpypz[1] );
}

//##################################################################################
Double_t AliHLTJETBase::GetPtFromXYZ( const Float_t *pxpypz ) {
  // see header file for class documentation
  // Cartesian [px,py,pz] 
  
  Double_t px = (Double_t) pxpypz[0];
  Double_t py = (Double_t) pxpypz[1];

  return TMath::Sqrt( px*px + py*py );
}

//##################################################################################
Double_t AliHLTJETBase::GetPhiFromXYZ( const Double_t *xyz ) {
  // see header file for class documentation
  // Cartesian [x,y,z,E] - xyz
 
  Double_t phi = TMath::ATan2( xyz[1], xyz[0] );  
  if ( phi < 0.0 ) phi += TMath::TwoPi();

  return phi;
}

//##################################################################################
Double_t AliHLTJETBase::GetPhiFromXYZ( const Float_t *xyz ) {
  // see header file for class documentation
  // Cartesian [x,y,z,E] - xyz
 
  Double_t x = (Double_t) xyz[0];
  Double_t y = (Double_t) xyz[1];

  Double_t phi = TMath::ATan2( y, x );  
  if ( phi < 0.0 ) phi += TMath::TwoPi();

  return phi;
}

//##################################################################################
Double_t AliHLTJETBase::GetEtaFromXYZ( const Double_t *xyz ) {
  // see header file for class documentation
  // Cartesian [x,y,z,E] - xyz
  
  Double_t r = TMath::Sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
  return 0.5 * TMath::Log( (r+xyz[2])/(r-xyz[2]) );
}
//##################################################################################
Double_t AliHLTJETBase::GetEtaFromXYZ( const Float_t *xyz ) {
  // see header file for class documentation
  // Cartesian [x,y,z,E] - xyz
  
  Double_t x = (Double_t) xyz[0];
  Double_t y = (Double_t) xyz[1];
  Double_t z = (Double_t) xyz[2];
  
  Double_t r = TMath::Sqrt( x*x + y*y + z*z );
  return 0.5 * TMath::Log( (r+z)/(r-z) );
}

//################################################################################## 
Double_t AliHLTJETBase::GetDistance2( const Double_t eta1, const Double_t phi1, 
					     const Double_t eta2, const Double_t phi2) {
  // see header file for class documentation
  
  return ( (eta1-eta2)*(eta1-eta2) ) + ( (phi1-phi2)*(phi1-phi2) );
}


#endif

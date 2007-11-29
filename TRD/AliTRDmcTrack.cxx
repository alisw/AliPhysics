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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD MC track                                                          //
//  Used for efficiency estimates and matching of reconstructed tracks    //
//  to MC particles                                                       //                    
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDmcTrack.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDmcTrack)

//_____________________________________________________________________________
AliTRDmcTrack::AliTRDmcTrack() 
  :TObject()
  ,fLab(-1)
  ,fSeedLab(-1)
  ,fPrimary(kFALSE)
  ,fMass(0)
  ,fCharge(0)
  ,fPDG(0)
  ,fN(0) 
{ 
  //
  // Default constructor 
  //

  for (Int_t ltb = 0; ltb < kMAXTB; ltb++) {
    for (Int_t plane = 0; plane < 6; plane++) {
      fIndex[ltb][plane][0] = -1;
      fIndex[ltb][plane][1] = -1;
    }
  }

  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 3; j++) { 
      fPin[i][j]    = 0.0; 
      fPout[i][j]   = 0.0;
      fXYZin[i][j]  = 0.0; 
      fXYZout[i][j] = 0.0;
    }
  }

}

//_____________________________________________________________________________
AliTRDmcTrack::AliTRDmcTrack(Int_t label, Int_t seedLabel, Bool_t primary
			   , Float_t mass, Int_t charge, Int_t pdg) 
  :TObject()
  ,fLab(label)
  ,fSeedLab(seedLabel)
  ,fPrimary(primary)
  ,fMass(mass)
  ,fCharge(charge)
  ,fPDG(pdg)
  ,fN(0) 
{ 
  //
  // Main constructor 
  //
  
  for (Int_t ltb = 0; ltb < kMAXTB; ltb++) {
    for (Int_t plane = 0; plane < 6; plane++) {
      fIndex[ltb][plane][0] = -1;
      fIndex[ltb][plane][1] = -1;
    }
  }
  
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 3; j++) { 
      fPin[i][j]    = 0.0; 
      fPout[i][j]   = 0.0;
      fXYZin[i][j]  = 0.0; 
      fXYZout[i][j] = 0.0;
    }
  }

}

//_____________________________________________________________________________
void AliTRDmcTrack::GetPxPyPzXYZ(Double_t& px, Double_t& py, Double_t& pz
			       , Double_t&  x, Double_t&  y, Double_t&  z 
			       , Int_t opt) const 
{
  //
  // Returns track momentum components and coordinates at the entrance 
  // (opt >= 0), or exit (opt < 0) of TRD. 
  //

  Int_t i = 0;

  if (opt >= 0) {

    for (i = 0; i < AliTRDgeometry::Nplan(); i++) {
      if ((fPin[i][0]*fPin[i][0]
         + fPin[i][1]*fPin[i][1]
         + fPin[i][2]*fPin[i][2]) > 0.0005) {
        break;
      }
    }

    px = fPin[i][0];   
    py = fPin[i][1];
    pz = fPin[i][2];
    x  = fXYZin[i][0];   
    y  = fXYZin[i][1];
    z  = fXYZin[i][2];

  }
  else {

    for (i = AliTRDgeometry::Nplan() - 1; i >= 0; i--) {
      if ((fPout[i][0]*fPout[i][0]
         + fPout[i][1]*fPout[i][1]
         + fPout[i][2]*fPout[i][2]) > 0.0005) {
        break;
      }
    }

    px = fPout[i][0];
    py = fPout[i][1];
    pz = fPout[i][2];
    x  = fXYZout[i][0];
    y  = fXYZout[i][1];
    z  = fXYZout[i][2];

  }

  return;
}

//_____________________________________________________________________________
void AliTRDmcTrack::GetPlanePxPyPz(Double_t& px, Double_t& py, Double_t& pz
                                 , Int_t plane, Int_t opt) const 
{
  //
  // Returns momentum components at the entrance (opt >= 0), or
  // exit (opt < 0) of TRD plane <plane>. 
  //

  if (opt >= 0) {
    px = fPin[plane][0];
    py = fPin[plane][1];
    pz = fPin[plane][2];
  }
  else {
    px = fPout[plane][0];
    py = fPout[plane][1];
    pz = fPout[plane][2];
  }

  return;

}

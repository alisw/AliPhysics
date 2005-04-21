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
 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD clusterCorrection                                                    //
//  marian.ivanov@cern.ch                                                    //
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 

#include "AliTRDclusterCorrection.h"
#include "TFile.h"

ClassImp(AliTRDclusterCorrection)

AliTRDclusterCorrection * gAliTRDclusterCorrection=0;

//_____________________________________________________________________________
AliTRDclusterCorrection::AliTRDclusterCorrection()
{
  //
  //default constructor for correction
  //

  for (Int_t iplane=0;iplane<6;iplane++) {
    for (Int_t itime=0;itime<30;itime++) {
      for(Int_t iangle=0;iangle<20;iangle++){	
	fCorrections[iplane][itime][iangle][0]=0;
	fCorrections[iplane][itime][iangle][1]=0;
      }
    }
  }

  fOffsetAngle =0;

}

//_____________________________________________________________________________
void AliTRDclusterCorrection::SetCorrection(Int_t plane,Int_t timebin, Float_t angle, 
					    Float_t value, Float_t sigma)
{
  //
  // Set the correction factors
  //

  Int_t iangle = int( (angle-fOffsetAngle+1.)*10.+0.5);
  if (iangle<0) return;
  if (iangle>=20) return;
  fCorrections[plane][timebin][iangle][0] = value;
  fCorrections[plane][timebin][iangle][1] = sigma;

}

//_____________________________________________________________________________
Float_t AliTRDclusterCorrection::GetCorrection(Int_t plane, Int_t timebin, Float_t angle) const
{
  //
  // Get the correction factors
  //

  Int_t iangle = int( (angle-fOffsetAngle+1.)*10.+0.5);
  if (iangle<0) return 0.;
  if (iangle>=20) return 0.;
  return fCorrections[plane][timebin][iangle][0];

}

//_____________________________________________________________________________
Float_t AliTRDclusterCorrection::GetSigma(Int_t plane, Int_t timebin, Float_t angle) const
{
  //
  // Returns the sigma
  //

  Int_t iangle = int( (angle-fOffsetAngle+1.)*10.+0.5);
  if (iangle<0) return 1.;
  if (iangle>=20) return 1.;
  return fCorrections[plane][timebin][iangle][1];

}

//_____________________________________________________________________________
AliTRDclusterCorrection *  AliTRDclusterCorrection::GetCorrection()
{
  //
  // Return an instance of AliTRDclusterCorrection
  //

  if (gAliTRDclusterCorrection!=0) return gAliTRDclusterCorrection;
  //
  TFile * f  = new TFile("$ALICE_ROOT/TRD/TRDcorrection.root");
  if (!f){
    ////
    gAliTRDclusterCorrection = new AliTRDclusterCorrection();
    return gAliTRDclusterCorrection;
  }
  gAliTRDclusterCorrection = (AliTRDclusterCorrection*) f->Get("TRDcorrection");
  if (gAliTRDclusterCorrection==0)  
    gAliTRDclusterCorrection = new AliTRDclusterCorrection();

  return gAliTRDclusterCorrection;
  
}

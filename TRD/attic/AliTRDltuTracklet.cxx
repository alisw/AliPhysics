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
//                                                                           //
//  TRD chamber local track (LTU, tracklet)                                  //
//                                                                           //
//  Author:                                                                  //
//    Bogdan Vulpescu                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>

#include "AliTRDltuTracklet.h"

ClassImp(AliTRDltuTracklet)

//_____________________________________________________________________________
AliTRDltuTracklet::AliTRDltuTracklet()
  :TObject()
  ,fX(0)
  ,fY(0)
  ,fSlope(0)
  ,fRowz(0)
  ,fDetector(0)
  ,fRow(0)
  ,fNclusters(0)
  ,fLabel(0)
  ,fQ(0)
{
  //
  // AliTRDltuTracklet default constructor
  //

}

//_____________________________________________________________________________
AliTRDltuTracklet::AliTRDltuTracklet(Int_t det, Int_t row, Float_t rowz
				   , Float_t slope, Float_t offset, Float_t time
				   , Int_t ncl, Int_t label, Float_t q) 
  :TObject()
  ,fX(time)
  ,fY(offset)
  ,fSlope(slope)
  ,fRowz(rowz)
  ,fDetector(det)
  ,fRow(row)
  ,fNclusters(ncl)
  ,fLabel(label)
  ,fQ(q)
{
  //
  // AliTRDltuTracklet constructor
  //

}

//_____________________________________________________________________________
AliTRDltuTracklet::~AliTRDltuTracklet()
{
  //
  // Destructor
  //

}

//_____________________________________________________________________________
Float_t AliTRDltuTracklet::GetPt(Float_t field) const
{
  //
  // Transverse momentum calculation
  // Curvature R = (fX*fX + fY*fY) / (2 * sin(alpha))
  // alpha = angle deviation from "infinite momentum"
  //
  // Consistent with AliTRDmcmTracklet::GetPt(...)
  //

  Float_t infSlope = TMath::ATan(fY/fX) / TMath::Pi()*180.0;    
  Float_t alpha    = fSlope - infSlope;
  Float_t r        = TMath::Sqrt(fX*fX + fY*fY)
                   / (2.0 * TMath::Sin(alpha/180.0*TMath::Pi()));
  
  Float_t pt       = 0.3 * field * 0.01 * r;
 
  return pt;
 
}

//_____________________________________________________________________________
Int_t AliTRDltuTracklet::Compare(const TObject *o) const
{
  //
  // compare two LTU tracklets according to the intercept point Y1
  //

  AliTRDltuTracklet *ltutrk = (AliTRDltuTracklet *) o;

  if (fRow      != ltutrk->fRow) {
    return +1;
  }
  if (fDetector != ltutrk->fDetector) {
    return +1;
  }

  if (fY <  ltutrk->fY) {
    return -1;
  }
  if (fY == ltutrk->fY) {
    return  0;
  }

  return 1;

}

//_____________________________________________________________________________
Float_t AliTRDltuTracklet::GetYproj(Float_t xpl) const
{
  //
  // y-projection (bending plane) onto the median plane
  //

  Float_t yproj = fY + TMath::Tan(fSlope/180.0*TMath::Pi()) * (xpl - fX);

  return yproj;

}

//_____________________________________________________________________________
Float_t AliTRDltuTracklet::GetZproj(Float_t xpl) const
{
  //
  // z-projection (using pad row center) onto the median plane
  //

  Float_t zproj = fRowz * xpl / fX;

  return zproj;

}


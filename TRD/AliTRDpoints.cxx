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

/*
$Log$
Revision 1.2  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.2.1  2000/09/18 13:44:21  cblume
New class AliTRDpoints to display the TR photon hits

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the TRD points for the ALICE event display           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TPad.h>
#include <TView.h>

#include "AliRun.h"
#include "AliDetector.h"

#include "AliTRDpoints.h"
 
ClassImp(AliTRDpoints)

//_____________________________________________________________________________
AliTRDpoints::AliTRDpoints():AliPoints()
{
  //
  // Default constructor
  //

  fNTRpoints    = 0;
  fTRpolyMarker = 0;

}

//_____________________________________________________________________________
AliTRDpoints::AliTRDpoints(Int_t nhitsE, Int_t nhitsT):AliPoints(nhitsE)
{
  //
  // Standard constructor
  //

  fNTRpoints    = nhitsT;
  fTRpolyMarker = 0;

}
	 
//_____________________________________________________________________________
AliTRDpoints::AliTRDpoints(const AliTRDpoints &p)
{
  //
  // Copy contructor
  //
 
  ((AliTRDpoints &) p).Copy(*this);

}

//_____________________________________________________________________________
AliTRDpoints::~AliTRDpoints()
{
  //
  // Default destructor
  //

  if (fTRpolyMarker) delete fTRpolyMarker;

}

//_____________________________________________________________________________
AliTRDpoints &AliTRDpoints::operator=(const AliTRDpoints &p)
{
  //
  // Assignment operator 
  //

  if (this != &p) ((AliTRDpoints &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDpoints::Copy(TObject &p)
{
  //
  // Copy function
  //

  ((AliTRDpoints &) p).fNTRpoints = fNTRpoints; 
  for (Int_t i = 0; i < 3*fNTRpoints; i++) {
    ((AliTRDpoints &) p).fTRpoints[i] = fTRpoints[i];
  }

}

//_____________________________________________________________________________
void AliTRDpoints::Draw(Option_t *option) 
{
  //
  // Draws a TRD point
  //

  AliPoints::Draw(option);

  //if (fTRpolyMarker) delete fTRpolyMarker;
  if (fNTRpoints) {
    fTRpolyMarker = new TPolyMarker3D(fNTRpoints,fTRpoints,29);
    fTRpolyMarker->SetMarkerColor(2); 
    fTRpolyMarker->SetMarkerSize(0.8);
    fTRpolyMarker->Draw(option);
  }

}

//_____________________________________________________________________________
void AliTRDpoints::SetTRpoints(Int_t n, Float_t *coor) 
{
  //
  // Sets the number and the coordinates of the photon hits
  //

  if (kNTRpoints >= 3 * n) {
    fNTRpoints = n;
    for (Int_t i = 0; i < 3*n; i++) {
      fTRpoints[i] = coor[i];
    } 
  }
  else {
    printf("AliTRDpoints::SetTRpoints -- ");
    printf("Boundary error: %d/%d\n",3*n,kNTRpoints);
  }

}

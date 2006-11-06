/**************************************************************************
 * Copyright(c) 2001-2002, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////
// Class used by TPC tracking parameterization to handle to tracking
// parameters (efficiencies, etc...) on a kinematic grid [pt,eta].
// User has to provide the grid steps and the values of the parameter
// in the points of the grid. The function GetValueAt(pt,eta) returns
// the result of a linear interpolation at the point [pt,eta].
//
//  Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it
////////////////////////////////////////////////////////////////////////

//-- standard headers -----
#include "Riostream.h"
//--- Root headers --------
#include <TMatrixD.h>
#include <TArrayD.h>
//-- AliRoot headers ------
#include "AliTPCkineGrid.h"
//-------------------------

ClassImp(AliTPCkineGrid)

//------------------------------------------------------------------------
  AliTPCkineGrid::AliTPCkineGrid()
    :TNamed(),
     fNpt(0),
     fNeta(0),
     fPt(0),
     fEta(0),
     fParams(0) 
{
//------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------

}
//------------------------------------------------------------------------
AliTPCkineGrid::AliTPCkineGrid(Int_t npt,Int_t neta,
			       Double_t* pt,Double_t* eta)
  :TNamed(),
    fNpt(npt),
    fNeta(neta),
    fPt(0),
    fEta(0),
    fParams(0) 
{
//------------------------------------------------------------------------
// Standard constructor
//------------------------------------------------------------------------


  fPt   = new TArrayD(fNpt);
  fEta  = new TArrayD(fNeta);
 
  for(Int_t i=0; i<npt; i++)  (*fPt)[i]  = pt[i];
  for(Int_t i=0; i<neta; i++) (*fEta)[i] = eta[i];

  fParams = new TMatrixD(fNpt,fNeta);
}
//-------------------------------------------------------------------------
AliTPCkineGrid::AliTPCkineGrid(const AliTPCkineGrid& grid):TNamed(grid),
     fNpt(0),
     fNeta(0),
     fPt(0),
     fEta(0),
     fParams(0) 
{
//-------------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------------
  fNpt = grid.fNpt;
  fNeta = grid.fNeta;
  fPt = new TArrayD(fNpt);
  for(Int_t i=0; i<fNpt; i++) (*fPt)[i] = grid.fPt->At(i);
  fEta = new TArrayD(fNeta);
  for(Int_t i=0; i<fNeta; i++) (*fEta)[i] = grid.fEta->At(i);

  fParams = new TMatrixD(fNpt,fNeta);
  for(Int_t i=0; i<fNpt; i++) {
    for(Int_t j=0; j<fNeta; j++) (*fParams)(i,j)=(*grid.fParams)(i,j);
  }
}
//--------------------------------------------------------------------------
AliTPCkineGrid::~AliTPCkineGrid() {
//--------------------------------------------------------------------------
// Destructor
//--------------------------------------------------------------------------
  delete fPt;
  delete fEta;
  delete fParams;
}
//__________________________________________________________________________
AliTPCkineGrid & AliTPCkineGrid::operator =(const AliTPCkineGrid & param)
{
  //
  // assignment operator - dummy
  //
  fNpt=param.fNpt;
  return (*this);
}
//--------------------------------------------------------------------------
void AliTPCkineGrid::GetArrayEta(Double_t* eta) const {
//--------------------------------------------------------------------------
// This functions returns an array with the eta points
//--------------------------------------------------------------------------
  for(Int_t i=0;i<fNeta;i++) eta[i] = fEta->At(i);
  return;
}
//--------------------------------------------------------------------------
void AliTPCkineGrid::GetArrayPt(Double_t* pt) const {
//--------------------------------------------------------------------------
// This functions returns an array with the pt points
//--------------------------------------------------------------------------
  for(Int_t i=0;i<fNpt;i++) pt[i] = fPt->At(i);
  return;
}
//--------------------------------------------------------------------------
Int_t AliTPCkineGrid::GetBin(Double_t pt,Double_t eta) const {
//--------------------------------------------------------------------------
// This functions tells in which bin of the grid a certain point falls
//--------------------------------------------------------------------------

  Int_t etaBin=0,ptBin=0,bin=0;
  eta = TMath::Abs(eta);  

  /*   this is how bins are numbered

             ...  ... .
             ---+---+---
       ^      6 | 7 | 8
       |     ---+---+---
              3 | 4 | 5
       Pt    ---+---+---
              0 | 1 | 2

                eta ->
   */

  if(eta < fEta->At(0)) {
    etaBin = 0;
  } else if(eta > fEta->At(fNeta-1)) {
    etaBin = fNeta;
  } else {
    for(Int_t i=0; i<fNeta; i++) {
      if(eta < fEta->At(i)) {
	etaBin = i;
	break;
      } 
    }
  }
  if(pt < fPt->At(0)) {
    ptBin = 0;
  } else if(pt > fPt->At(fNpt-1)) {
    ptBin = fNpt;
  } else {
    for(Int_t i=0; i<fNpt; i++) {
      if(pt < fPt->At(i)) {
	ptBin = i;
	break;
      } 
    }
  }

  bin = ptBin*(fNeta+1) + etaBin;

  return bin;
}
//--------------------------------------------------------------------------
Double_t AliTPCkineGrid::GetParam(Int_t i) const {
//--------------------------------------------------------------------------
// This functions allows to get parameters using only one index
//--------------------------------------------------------------------------
  Int_t ipt = (Int_t)i/fNeta;
  Int_t ieta = i-ipt*fNeta;
  return GetParam(ipt,ieta);
}
//--------------------------------------------------------------------------
Double_t AliTPCkineGrid::GetValueAt(Double_t pt,Double_t eta) const {
//--------------------------------------------------------------------------
// This functions makes a linear interpolation at the point [pt,eta]
//--------------------------------------------------------------------------

  // determine the points to be used in the interpolation: 
  //
  // eta
  eta = TMath::Abs(eta);
  Int_t etaLow=0,etaUp=0;
  if(eta < fEta->At(0)) {
    etaLow = 0;
    etaUp = 1;
  } else if(eta >= fEta->At(fNeta-1)) {
    etaLow = fNeta-2;
    etaUp = fNeta-1;
  } else {
    for(Int_t i=0; i<fNeta; i++) {
      if(eta < fEta->At(i)) {
	etaLow = i-1;
	etaUp = i;
	break;
      } 
    }
  }
  //
  // pt
  Int_t ptLow=0,ptUp=0;
  if(pt < fPt->At(0)) {
    ptLow = 0;
    ptUp = 1;
  } else if(pt >= fPt->At(fNpt-1)) {
    ptLow = fNpt-2;
    ptUp = fNpt-1;
  } else {
    for(Int_t i=0; i<fNpt; i++) {
      if(pt < fPt->At(i)) {
	ptLow = i-1;
	ptUp = i;
	break;
      } 
    }
  }

  //cerr<<" Pt = ("<<ptLow<<","<<ptUp<<")   Eta = ("<<etaLow<<","<<etaUp<<")\n";

  Double_t intValue=0,intValueEtaLow=0,intValueEtaUp=0;
  // interpolate, at etaLow, between ptLow and ptUp
  intValueEtaLow = (*fParams)(ptLow,etaLow)+
                  ((*fParams)(ptUp,etaLow)-(*fParams)(ptLow,etaLow))/
                  (fPt->At(ptUp)-fPt->At(ptLow))*(pt-fPt->At(ptLow));
  // interpolate, at etaUp, between ptLow and ptUp
  intValueEtaUp  = (*fParams)(ptLow,etaUp)+
                  ((*fParams)(ptUp,etaUp)-(*fParams)(ptLow,etaUp))/
                  (fPt->At(ptUp)-fPt->At(ptLow))*(pt-fPt->At(ptLow));
  // interpolate, at pt,  between etaLow and etaUp
  intValue       = intValueEtaLow+
                  (intValueEtaUp-intValueEtaLow)/
                  (fEta->At(etaUp)-fEta->At(etaLow))*(eta-fEta->At(etaLow));

  if(intValue<0.) intValue=0.;
  return intValue;
}
//--------------------------------------------------------------------------
void AliTPCkineGrid::SetParam(Int_t i,Double_t par) {
//--------------------------------------------------------------------------
// This functions allows to set parameters using only one index
//--------------------------------------------------------------------------
  Int_t ipt = (Int_t)i/fNeta;
  Int_t ieta = i-ipt*fNeta;
  SetParam(ipt,ieta,par);

  return;
}















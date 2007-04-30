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
Revision 0.02 2005/07/27 A. De Caro
         Implement new ctor AliTOFcluster(Double_t *, Int_t *)
	 for cluster construction from raw data

Revision 0.01 2005/07/25 A. De Caro
         Implement new class statring from
	 class definition inside AliTOFtracker class
	 (originally implemented by S. Arcelli and C. Zampolli)
*/

////////////////////////////////////////////////////////
//                                                    //
//    AliTOFcluster Class                             //
//    Description: class for TOF cluster definition   //
//                                                    //
////////////////////////////////////////////////////////

#include "AliTOFcluster.h"

ClassImp(AliTOFcluster)

AliTOFcluster::AliTOFcluster():
  fIdx(-1),
  fR(0),
  fPhi(0),
  fZ(0),
  fQuality(-100), 
  fTDC(0),
  fToT(0),
  fADC(0),
  fTdcND(0),
  fTdcRAW(0),
  fStatus(kTRUE) 
 {
  //
  // default ctor
  //

  Int_t ii;
  for (ii=0; ii<3; ii++) fLab[ii]      = -1;
  for (ii=0; ii<5; ii++) fdetIndex[ii] = -1;
}
//-------------------------------------------------------------------------

AliTOFcluster::AliTOFcluster(Double_t *h,  Int_t *ind, Int_t *par, Bool_t status,Int_t *l, Int_t idx):
  TObject(),
  fIdx(idx),
  fR(h[0]),
  fPhi(h[1]),
  fZ(h[2]),
  fQuality(-100), 
  fTDC(par[0]),
  fToT(par[1]),
  fADC(par[2]),
  fTdcND(par[3]),
  fTdcRAW(par[4]),
  fStatus(status) 
 {
  //
  // constructor
  //

  Int_t ii;
  for (ii=0; ii<3; ii++) fLab[ii]      = l[ii];
  for (ii=0; ii<5; ii++) fdetIndex[ii] = ind[ii];
}
//-------------------------------------------------------------------------

AliTOFcluster::AliTOFcluster(Double_t *h, Int_t *ind, Int_t *par):
  TObject(),
  fIdx(-1),
  fR(h[0]),
  fPhi(h[1]),
  fZ(h[2]),
  fQuality(-100), 
  fTDC(par[0]),
  fToT(par[1]),
  fADC(par[2]),
  fTdcND(par[3]),
  fTdcRAW(par[4]),
  fStatus(kTRUE) 
 {
  //
  // constructor
  //

  Int_t ii;
  for (ii=0; ii<3; ii++) fLab[ii]      = -1;
  for (ii=0; ii<5; ii++) fdetIndex[ii] = ind[ii];
}
//-------------------------------------------------------------------------

AliTOFcluster::AliTOFcluster(const AliTOFcluster & cluster):
  TObject(),
  fIdx(-1),
  fR(0),
  fPhi(0),
  fZ(0),
  fQuality(-100), 
  fTDC(0),
  fToT(0),
  fADC(0),
  fTdcND(0),
  fTdcRAW(0),
  fStatus(kTRUE) 
 {
  //
  // copy ctor for AliTOFcluster object
  //

  Int_t ii;
  fR        = cluster.fR;
  fPhi      = cluster.fPhi;
  fZ        = cluster.fZ;
  fQuality    = cluster.fQuality; 
  for (ii=0; ii<3; ii++) fLab[ii]      = cluster.fLab[ii];
  fIdx      = cluster.fIdx;
  for (ii=0; ii<5; ii++) fdetIndex[ii] = cluster.fdetIndex[ii];
  fTDC      = cluster.fTDC;
  fToT    = cluster.fToT; 
  fADC      = cluster.fADC;
  fTdcND    = cluster.fTdcND; 
  fTdcRAW    = cluster.fTdcRAW; 
  fStatus    = cluster.fStatus; 
}
//-------------------------------------------------------------------------

AliTOFcluster::~AliTOFcluster() {
  //
  // dtor
  //

  //delete fLab;
  //delete fdetIndex;

}

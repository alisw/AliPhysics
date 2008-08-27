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
//#include "AliLog.h"
//#include "AliGeomManager.h"
//#include "TGeoMatrix.h"

ClassImp(AliTOFcluster)

AliTOFcluster::AliTOFcluster():
  AliCluster3D(),
  fIdx(-1),
  fQuality(-100), 
  fR(0),
  fPhi(0),
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
  for (ii=0; ii<5; ii++) fdetIndex[ii] = -1;
}
//-------------------------------------------------------------------------

AliTOFcluster::AliTOFcluster(UShort_t volId, 
   Float_t x,   Float_t y,   Float_t z,
   Float_t sx2, Float_t sxy, Float_t sxz,
                Float_t sy2, Float_t syz,
                             Float_t sz2, Int_t *lab, Int_t *ind, Int_t *par, Bool_t status, Int_t idx):
  AliCluster3D(volId,x,y,z,sx2,sxy,sxz,sy2,syz,sz2,lab),
  fIdx(idx),
  fQuality(-100), 
  fR(0),
  fPhi(0),
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
   for (ii=0; ii<5; ii++) fdetIndex[ii] = ind[ii];
   
   Float_t xyz[3];
   GetGlobalXYZ(xyz);
   fR   = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);   
   fPhi = TMath::ATan2(xyz[1], xyz[0]);

}
//-------------------------------------------------------------------------

AliTOFcluster::AliTOFcluster(const AliTOFcluster & cluster):
  AliCluster3D(cluster),
  fIdx(cluster.fIdx),
  fQuality(cluster.fQuality), 
  fR(cluster.fR),
  fPhi(cluster.fPhi),
  fTDC(cluster.fTDC),
  fToT(cluster.fToT),
  fADC(cluster.fADC),
  fTdcND(cluster.fTdcND),
  fTdcRAW(cluster.fTdcRAW),
  fStatus(cluster.fStatus) 
 {
  //
  // copy ctor for AliTOFcluster object
  //

  Int_t ii;
  for (ii=0; ii<5; ii++) fdetIndex[ii] = cluster.fdetIndex[ii];
}
//-------------------------------------------------------------------------

AliTOFcluster::~AliTOFcluster() {

  //
  // dtor
  //

}



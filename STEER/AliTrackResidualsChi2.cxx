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

//-----------------------------------------------------------------
//   Implementation of the derived class for track residuals
//   based on Chi2 minimization
//
//-----------------------------------------------------------------

#include <TMinuit.h>
#include <TGeoMatrix.h>

#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliTrackResidualsChi2.h"

ClassImp(AliTrackResidualsChi2)

//______________________________________________________________________________
void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // This function is called by minuit
  // The corresponding member method is called
  // using SetObjectFit/GetObjectFit methods of TMinuit
  AliTrackResidualsChi2* dummy = (AliTrackResidualsChi2 *)gMinuit->GetObjectFit();
  dummy->Chi2(npar, gin, f, par, iflag);
}

//______________________________________________________________________________
Bool_t AliTrackResidualsChi2::Minimize()
{
  // Implementation of Chi2 based minimization
  // of track residuala sum
  TMinuit *gMinuit = new TMinuit(6);  //initialize TMinuit
  gMinuit->SetObjectFit(this);
  gMinuit->SetFCN(Fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // Set starting values and step sizes for parameters
  Double_t tr[3],rot[3];
  fAlignObj->GetPars(tr,rot);
  static Double_t step[6] = {0,0,0,0,0,0};
  gMinuit->mnparm(0, "x", tr[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "y", tr[1], step[1], 0,0,ierflg);
  gMinuit->mnparm(2, "z", tr[2], step[2], 0,0,ierflg);
  gMinuit->mnparm(3, "psi", rot[0], step[3], 0,0,ierflg);
  gMinuit->mnparm(4, "theta", rot[1], step[4], 0,0,ierflg);
  gMinuit->mnparm(5, "phi", rot[2], step[5], 0,0,ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);

  return kTRUE;
}

//______________________________________________________________________________
void AliTrackResidualsChi2::Chi2(Int_t & /* npar */, Double_t * /* gin */, Double_t &f, Double_t *par, Int_t /* iflag */)
{
  // Chi2 function to be minimized
  // Sums all the track residuals
  Double_t chi2 = 0;

  fAlignObj->SetPars(par[0],par[1],par[2],par[3],par[4],par[5]);
  TGeoHMatrix m;
  fAlignObj->GetMatrix(m);
  Double_t *rot = m.GetRotationMatrix();
  Double_t *tr  = m.GetTranslation();

  AliTrackPoint p1,p2;
  Double_t dR[3];
  Float_t xyzvol[3],xyztr[3];
  Float_t covvol[6],covtr[6];

  for (Int_t itrack = 0; itrack < fLast; itrack++) {
    if (!fVolArray[itrack] || !fTrackArray[itrack]) continue;
    for (Int_t ipoint = 0; ipoint < fVolArray[itrack]->GetNPoints(); ipoint++) {
      fVolArray[itrack]->GetPoint(p1,ipoint);
      fTrackArray[itrack]->GetPoint(p2,ipoint);
      p1.GetXYZ(xyzvol,covvol);
      p2.GetXYZ(xyztr,covtr);

      for (Int_t i = 0; i < 3; i++)
	dR[i] = tr[i] 
               + xyzvol[0]*rot[3*i]
               + xyzvol[1]*rot[3*i+1]
               + xyzvol[2]*rot[3*i+2]
	       - xyztr[i];
      chi2 += dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];
    }
  }
  f = chi2;
}

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

//-----------------------------------------------------------------
//           Implementation of the TPC cluster class
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "AliTPCcluster.h"

ClassImp(AliTPCcluster)
//
AliTPCcluster:: AliTPCcluster(Int_t *lab, Float_t *hit)
              :AliCluster(lab,hit),
               fQ(0.)
{
  //
  // constructor
  //
  fQ=hit[4];
}
 
//____________________________________________________________________________
Double_t AliTPCcluster::SigmaY2(Double_t r, Double_t tgl, Double_t pt)
{
  //
  // Parametrised error of the cluster reconstruction (pad direction)
  //
  // Sigma rphi
  const Float_t kArphi=0.41818e-2;
  const Float_t kBrphi=0.17460e-4;
  const Float_t kCrphi=0.30993e-2;
  const Float_t kDrphi=0.41061e-3;

  pt=TMath::Abs(pt)*1000.;
  Double_t x=r/pt;
  tgl=TMath::Abs(tgl);
  Double_t s=kArphi - kBrphi*r*tgl + kCrphi*x*x + kDrphi*x;
  if (s<0.4e-3) s=0.4e-3;
  s*=1.3; //Iouri Belikov

  return s;
}


//____________________________________________________________________________
Double_t AliTPCcluster::SigmaZ2(Double_t r, Double_t tgl)
{
  //
  // Parametrised error of the cluster reconstruction (drift direction)
  //
  // Sigma z
  const Float_t kAz=0.39614e-2;
  const Float_t kBz=0.22443e-4;
  const Float_t kCz=0.51504e-1;


  tgl=TMath::Abs(tgl);
  Double_t s=kAz - kBz*r*tgl + kCz*tgl*tgl;
  if (s<0.4e-3) s=0.4e-3;
  s*=1.3; //Iouri Belikov

  return s;
}


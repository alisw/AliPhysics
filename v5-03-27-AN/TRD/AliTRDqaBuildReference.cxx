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

/* $Id: AliTRDqaBuildReference.cxx 26344 2008-06-03 10:28:50Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Build the reference histograms                                        //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TH1D.h"
#include "TFile.h"

#include "AliTRDqaBuildReference.h"

//////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
AliTRDqaBuildReference &AliTRDqaBuildReference::operator=(const AliTRDqaBuildReference &qadm)
{
  //
  // Assignment operator
  //

  if (this != &qadm) {
    ((AliTRDqaBuildReference &) qadm).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDqaBuildReference::Copy(TObject &/*qadm*/) const
{
  //
  // Copy function
  //

}

void  AliTRDqaBuildReference::BuildRefHistos(TFile *file) const
{
  //
  // Build the reference histograms
  //

  // check
  if (!file) return;
  if (file->IsZombie()) return;
  
  // recpoints
  Int_t cd = file->cd("TRD/RecPoints");
  if (cd) {
    
    // MPV distribution

    TH1D *r = new TH1D("qaTRD_recPoints_ampMPV_ref", "", 150, 0, 150);
    for(Int_t i=0; i<r->GetNbinsX(); i++) {
      Double_t x = r->GetBinCenter(i+1);
      if (x < 10 || x > 70) r->Fill(x, 1);
      else if (x < 25 || x > 55) r->Fill(x, 0.5);
    }
    
    r->Fill(-1e3, 1);
    r->Fill(1e3, 1);
    
    // number of clusters
    r = new TH1D("qaTRD_recPoints_nCls_ref", "", 500, -0.5, 499.5);
    for(Int_t i=0; i<r->GetNbinsX(); i++) {
      Double_t x = r->GetBinCenter(i+1);
      if ( (i+4)%22 > 9 || i < 10) r->Fill(x, 0.5);
      if (x > 350) r->Fill(x, 1);
    }

    file->Write();
  }

  // esds

}


//////////////////////////////////////////////////////////////////////////////////////

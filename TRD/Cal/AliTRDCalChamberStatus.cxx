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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for status of chambers                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TH2D.h"
#include "AliTRDCalChamberStatus.h"

ClassImp(AliTRDCalChamberStatus)

//_____________________________________________________________________________
AliTRDCalChamberStatus::AliTRDCalChamberStatus()
  :TNamed()
{
  //
  // AliTRDCalChamberStatus default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fStatus[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalChamberStatus::AliTRDCalChamberStatus(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
{
  //
  // AliTRDCalChamberStatus constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fStatus[idet] = 0;
  }

}
//_____________________________________________________________________________
TH2D* AliTRDCalChamberStatus::Plot(Int_t sm, Int_t rphi) 
{
  //
  // Plot chamber status for supermodule and halfchamberside 
  // as a function of layer and stack
  //

  TH2D *h2 = new TH2D(Form("sm_%d_rphi_%d",sm,rphi),Form("sm_%d_rphi_%d",sm,rphi),5,0.0,5.0,6,0.0,6.0);
  
  h2->SetXTitle("stack");
  h2->SetYTitle("layer");

  Int_t start = sm*30;
  Int_t end = (sm+1)*30;
  
  for(Int_t i=start; i<end; i++) {
    Int_t layer  = i%6;
    Int_t stack  = static_cast<int>((i-start)/6.);
    Int_t status = GetStatus(i);
    if(rphi == 0) {    
      if(status!=4) h2->Fill(stack,layer,status);
    }
    else if(rphi == 1) {
      if(status!=3) h2->Fill(stack,layer,status);
    }
    
  }

  return h2;

}
//_____________________________________________________________________________
TH2D* AliTRDCalChamberStatus::Plot(Int_t sm) 
{
  //
  // Plot chamber status for supermodule and halfchamberside 
  // as a function of layer and stack
  //

  TH2D *h2 = new TH2D(Form("sm_%d",sm),Form("sm_%d",sm),5,0.0,5.0,6,0.0,6.0);
  
  h2->SetXTitle("stack");
  h2->SetYTitle("layer");

  Int_t start = sm*30;
  Int_t end = (sm+1)*30;
  
  for(Int_t i=start; i<end; i++) {
    Int_t layer  = i%6;
    Int_t stack  = static_cast<int>((i-start)/6.);
    Int_t status = GetStatus(i);
    h2->Fill(stack,layer,status);
  }

  return h2;

}

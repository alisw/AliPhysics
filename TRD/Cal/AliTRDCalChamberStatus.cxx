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
void AliTRDCalChamberStatus::SetStatus(Int_t det, Char_t status) 
{

  //
  // set the chamber status
  //
  //

  switch(status)
    {
    case AliTRDCalChamberStatus::kGood:
      SETBIT(fStatus[det], kGood);
      CLRBIT(fStatus[det], kNoData);
      CLRBIT(fStatus[det], kBadCalibrated);
      break;
    case AliTRDCalChamberStatus::kNoData:
      CLRBIT(fStatus[det], kGood);
      SETBIT(fStatus[det], kNoData);
      SETBIT(fStatus[det], kNoDataHalfChamberSideA);
      SETBIT(fStatus[det], kNoDataHalfChamberSideB);
      //      SETBIT(fStatus[det], kBadCalibrated);
      break;
    case AliTRDCalChamberStatus::kNoDataHalfChamberSideA:
      SETBIT(fStatus[det], kNoDataHalfChamberSideA);
      if(TESTBIT(fStatus[det], kNoDataHalfChamberSideB)){
	SETBIT(fStatus[det], kNoData);
	CLRBIT(fStatus[det], kGood);
      }
      break;
    case AliTRDCalChamberStatus::kNoDataHalfChamberSideB:
      SETBIT(fStatus[det], kNoDataHalfChamberSideB);
      if(TESTBIT(fStatus[det], kNoDataHalfChamberSideA)) {
	CLRBIT(fStatus[det], kGood);
	SETBIT(fStatus[det], kNoData);
      }
      break;
    case AliTRDCalChamberStatus::kBadCalibrated:
      CLRBIT(fStatus[det], kGood);
      SETBIT(fStatus[det], kBadCalibrated);
      break;
    default:
      CLRBIT(fStatus[det], kGood);
      CLRBIT(fStatus[det], kNoData);
      CLRBIT(fStatus[det], kNoDataHalfChamberSideA);
      CLRBIT(fStatus[det], kNoDataHalfChamberSideB);
      CLRBIT(fStatus[det], kBadCalibrated);
    }

}
//_____________________________________________________________________________
void AliTRDCalChamberStatus::UnsetStatusBit(Int_t det, Char_t status) 
{

  //
  // unset the chamber status bit
  //
  //

  switch(status)
    {
    case AliTRDCalChamberStatus::kGood:
      CLRBIT(fStatus[det], kGood);
      break;
    case AliTRDCalChamberStatus::kNoData:
      CLRBIT(fStatus[det], kNoData);
      break;
    case AliTRDCalChamberStatus::kNoDataHalfChamberSideA:
      CLRBIT(fStatus[det], kNoDataHalfChamberSideA);
      break;
    case AliTRDCalChamberStatus::kNoDataHalfChamberSideB:
      CLRBIT(fStatus[det], kNoDataHalfChamberSideB);
      break;
    case AliTRDCalChamberStatus::kBadCalibrated:
      CLRBIT(fStatus[det], kBadCalibrated);
      break;
    default:
      CLRBIT(fStatus[det], kGood);
      CLRBIT(fStatus[det], kNoData);
      CLRBIT(fStatus[det], kNoDataHalfChamberSideA);
      CLRBIT(fStatus[det], kNoDataHalfChamberSideB);
      CLRBIT(fStatus[det], kBadCalibrated);
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
    Int_t stackn  = static_cast<int>((i-start)/6.);
    Int_t status = GetStatus(i);
    h2->Fill(stackn,layer,status);
    if(rphi == 0) {    
      if(!TESTBIT(fStatus[i], kNoDataHalfChamberSideB)) h2->Fill(stackn,layer,status);
    }
    else if(rphi == 1) {
      if(!TESTBIT(fStatus[i], kNoDataHalfChamberSideA)) h2->Fill(stackn,layer,status);
    }
  }

  return h2;

}
//_____________________________________________________________________________
TH2D* AliTRDCalChamberStatus::PlotNoData(Int_t sm, Int_t rphi) 
{
  //
  // Plot chamber data status for supermodule and halfchamberside 
  // as a function of layer and stack
  //

  TH2D *h2 = new TH2D(Form("sm_%d_rphi_%d_data",sm,rphi),Form("sm_%d_rphi_%d_data",sm,rphi),5,0.0,5.0,6,0.0,6.0);
  
  h2->SetXTitle("stack");
  h2->SetYTitle("layer");

  Int_t start = sm*30;
  Int_t end = (sm+1)*30;
  
  for(Int_t i=start; i<end; i++) {
    Int_t layer  = i%6;
    Int_t stackn  = static_cast<int>((i-start)/6.);
    if(rphi == 0) {    
      if(TESTBIT(fStatus[i], kNoDataHalfChamberSideB)) h2->Fill(stackn,layer,1);
      if(TESTBIT(fStatus[i], kNoData)) h2->Fill(stackn,layer,1);
    }
    else if(rphi == 1) {
      if(!TESTBIT(fStatus[i], kNoDataHalfChamberSideA)) h2->Fill(stackn,layer,1);
      if(!TESTBIT(fStatus[i], kNoData)) h2->Fill(stackn,layer,1);
    }
  }

  return h2;

}
//_____________________________________________________________________________
TH2D* AliTRDCalChamberStatus::PlotBadCalibrated(Int_t sm, Int_t rphi) 
{
  //
  // Plot chamber calibration status for supermodule and halfchamberside 
  // as a function of layer and stack
  //

  TH2D *h2 = new TH2D(Form("sm_%d_rphi_%d_calib",sm,rphi),Form("sm_%d_rphi_%d_calib",sm,rphi),5,0.0,5.0,6,0.0,6.0);
  
  h2->SetXTitle("stack");
  h2->SetYTitle("layer");

  Int_t start = sm*30;
  Int_t end = (sm+1)*30;
  
  for(Int_t i=start; i<end; i++) {
    Int_t layer  = i%6;
    Int_t stackn  = static_cast<int>((i-start)/6.);
    if(rphi == 0) {    
      if(TESTBIT(fStatus[i], kBadCalibrated)) h2->Fill(stackn,layer,1);
    }
    else if(rphi == 1) {
      if(TESTBIT(fStatus[i], kBadCalibrated)) h2->Fill(stackn,layer,1);
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
    Int_t stackn  = static_cast<int>((i-start)/6.);
    Int_t status = GetStatus(i);
    h2->Fill(stackn,layer,status);
  }

  return h2;

}

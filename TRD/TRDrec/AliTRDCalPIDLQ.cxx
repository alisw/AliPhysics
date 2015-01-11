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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Container for the distributions of dE/dx and the time bin of the     //
// max. cluster for electrons and pions                                 //
//                                                                      //
// Authors:                                                             //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>               //
//   Alex Bercuci <a.bercuci@gsi.de>                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliPID.h"

#include "TKDPDF.h"

#include "AliTRDCalPIDLQ.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDCalPIDLQ)

Float_t AliTRDCalPIDLQ::fgTrackSegLength[kNLength] = { 3.7, 3.9, 4.2, 5.0 };

//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ()
  :AliTRDCalPID("pid", "LQ PID references for TRD")
{
  //
  //  The Default constructor
  //

  //Init();

}

//_________________________________________________________________________
AliTRDCalPIDLQ::AliTRDCalPIDLQ(const Text_t *name, const Text_t *title)
  :AliTRDCalPID(name,title)
{
  //
  //  The main constructor
  //
  
  Init();

}

//_________________________________________________________________________
Bool_t AliTRDCalPIDLQ::LoadReferences(Char_t *refFile)
{
  //
  // Read the TRD dEdx histograms.
  //

  if(gSystem->AccessPathName(refFile, kReadPermission)){
    AliError(Form("File %s.root not readable", refFile));
    return kFALSE;
  }
  if(!TFile::Open(refFile)){
    AliError(Form("File %s corrupted", refFile));
    return kFALSE;
  }
  TObjArray *pdf(NULL);
  if (!( pdf = dynamic_cast<TObjArray*>(gFile->Get("PDF_LQ")))) {
    AliError("PID data not available");
    return kFALSE;
  }
  fModel=(TObjArray*)pdf->Clone("PDF");
  gFile->Close();
  return kTRUE;
}

//_________________________________________________________________________
TObject* AliTRDCalPIDLQ::GetModel(Int_t ip, Int_t is, Int_t slices) const
{
  //
  // Returns one selected dEdx histogram
  //

  if (is < 0 || is >= AliPID::kSPECIES) return NULL;
  if(ip<0 || ip>= kNMom ) return NULL;

  AliDebug(2, Form("Retrive dEdx distribution for %s @ p=%5.2f [GeV/c].", AliPID::ParticleShortName(is), fgTrackMomentum[ip]));
  return fModel->At(GetModelID(ip, is, slices));
}

//_________________________________________________________________________
Int_t AliTRDCalPIDLQ::GetNrefs()
{
// Returns the number of PDF distribution 
  return AliTRDCalPID::kNMom*AliPID::kSPECIES*2;
}

//_________________________________________________________________________
Double_t AliTRDCalPIDLQ::GetProbability(Int_t spec, Float_t mom
                                      , const Float_t * const dedx
                                      , Float_t length
                                      , Int_t slices) const
{
//
// Core function of AliTRDCalPID class for calculating the
// likelihood for species "spec" (see AliTRDtrack::kNspecie) of a
// given momentum "mom" and a given dE/dx (slice "dedx") yield per
// layer
//

  if (spec < 0 || spec >= AliPID::kSPECIES) return 0.;

  Bool_t k2D(slices==2);
  Double_t x[]={0., 0.};
  if(!CookdEdx(dedx, x, k2D)) return 0.;
    
  // find the interval in momentum and track segment length which applies for this data
/*  Int_t ilength = 1;
  while(ilength<kNLength-1 && length>fgTrackSegLength[ilength]){
    ilength++;
  }*/
  Int_t imom = 1;
  while(imom<kNMom-1 && mom>fgTrackMomentum[imom]) imom++;

  Double_t p[2], e[2], r;
  TKDPDF *pdf(NULL);

  AliDebug(2, Form("Looking %cD for %s@p=%6.4f[GeV/c] log(dEdx)={%7.2f %7.2f}[a.u.] l=%4.2f[cm].", k2D?'2':'1', AliPID::ParticleShortName(spec), mom, x[0], x[1], length));
  if(!(pdf = dynamic_cast<TKDPDF*>(fModel->At(GetModelID(imom-1, spec, slices))))) {
    AliError(Form("%cD Ref data @ idx[%d] not found in DB.", k2D?'2':'1', GetModelID(imom-1, spec, slices)));
    return 0.;
  }
  if(!pdf->GetSize()) pdf->Bootstrap();
  pdf->Eval(x, r, e[0], kFALSE);
  p[0]=TMath::Abs(r); // conversion from interpolation to PDF
  AliDebug(2, Form("LQ=%6.3f+-%5.2f%% @ %4.1f[GeV/c]", p[0], 1.E2*e[0]/p[0], fgTrackMomentum[imom-1]));

  if(!(pdf = dynamic_cast<TKDPDF*>(fModel->At(GetModelID(imom, spec, slices))))){
    AliError(Form("%cD Ref data @ idx[%d] not found in DB.", k2D?'2':'1', GetModelID(imom, spec, slices)));
    return p[0];
  }
  if(!pdf->GetSize()) pdf->Bootstrap();
  pdf->Eval(x, r, e[1], kFALSE);
  p[1]=TMath::Abs(r); // conversion from interpolation to PDF
  AliDebug(2, Form("LQ=%6.3f+-%5.2f%% @ %4.1f[GeV/c]", p[1], 1.E2*e[1]/p[1], fgTrackMomentum[imom]));
  
  // return interpolation over momentum binning
  if(mom < fgTrackMomentum[0]) return p[0];
  else if(mom > fgTrackMomentum[kNMom-1]) return p[1];
  else{ 
    Double_t lmom[2] = {fgTrackMomentum[imom-1],  fgTrackMomentum[imom]};
    return p[0] + (p[1] - p[0])*(mom - lmom[0])/(lmom[1] - lmom[0]);
  }
}

//_________________________________________________________________________
void AliTRDCalPIDLQ::Init()
{
//
// PID LQ list initialization
//
  fModel = new TObjArray(AliPID::kSPECIES  * kNMom * 2);
  fModel->SetOwner();
}

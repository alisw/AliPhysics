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

////////////////////////////////////////////////
//                                            //
//  Manager class for TRD   hits              //
//                                            //
////////////////////////////////////////////////

#include <Riostream.h>

#include <TClonesArray.h>

#include "AliTRDtrackHits.h"
#include "AliTRDhit.h"    

ClassImp(AliTRDtrackHits)
 
//_____________________________________________________________________________
void  AliTRDtrackHits::AddHitTRD(Int_t volumeID, Int_t trackID, Double_t x
		               , Double_t y, Double_t z, Int_t q, Bool_t inDrift)
{
  //
  // Add one TRD hit
  //

  if (inDrift) {
    q = 2 * q + 1;
  }
  else {
    q = 2 * q;
  }

  AddHitKartez(volumeID,trackID,x,y,z,q,0);

}

//_____________________________________________________________________________
Bool_t AliTRDtrackHits::First()
{
  //
  // Set Current hit for the first hit
  //

  if (fArray->GetSize() <= 0) {
    fCurrentHit->fStatus = kFALSE;
    return kFALSE;
  }

  AliTrackHitsParamV2 *param = (AliTrackHitsParamV2 *) fArray->At(0);
  if (!fHit) {
    fHit = new AliTRDhit;
  }
  if (!(param)) {
    fCurrentHit->fStatus = kFALSE;
    return kFALSE;
  }

  fCurrentHit->fParamIndex = 0;
  fCurrentHit->fStackIndex = 0;

  ((AliTRDhit *) fHit)->SetDetector(param->fVolumeID);
  ((AliTRDhit *) fHit)->SetTrack(param->fTrackID);
  ((AliTRDhit *) fHit)->SetX(param->fR * TMath::Cos(param->fFi));
  ((AliTRDhit *) fHit)->SetY(param->fR * TMath::Sin(param->fFi));
  ((AliTRDhit *) fHit)->SetZ(param->fZ); 
  ((AliTRDhit *) fHit)->SetQ(param->fCharge[0] / 2);  
  if ((param->fCharge[0] % 2) == 0) {
    ((AliTRDhit *) fHit)->SetAmplification(); 
  }
  else {
    ((AliTRDhit *) fHit)->SetDrift();
  }
  fCurrentHit->fR = param->fR;
  
  return fCurrentHit->fStatus = kTRUE;

}
//set current hit to next
//_____________________________________________________________________________
Bool_t AliTRDtrackHits::Next()
{
  //
  // Set current hit to next
  //

  if (!(fCurrentHit->fStatus)) { 
    return kFALSE;
  }
  fCurrentHit->fStackIndex++;

  AliTrackHitsParamV2 *param = (AliTrackHitsParamV2 *) 
                                 fArray->At(fCurrentHit->fParamIndex);

  if (fCurrentHit->fStackIndex >= param->fNHits) {
    fCurrentHit->fParamIndex++;
    if (fCurrentHit->fParamIndex >= fArray->GetEntriesFast()) {
      fCurrentHit->fStatus = kFALSE;
      return kFALSE;
    }
    param = (AliTrackHitsParamV2 *) fArray->At(fCurrentHit->fParamIndex);
    fCurrentHit->fStackIndex = 0; 
    fCurrentHit->fR          = param->fR;
  }

  Double_t ratio;
  Double_t dfi2 = param->fAn;
  dfi2 *= dfi2 * fCurrentHit->fR * fCurrentHit->fR;
  Double_t ddz2 = param->fTheta;
  ddz2 *= ddz2;
  ratio = TMath::Sqrt(1.0 + dfi2 + ddz2);  

  fCurrentHit->fR += fStep * param->fHitDistance[fCurrentHit->fStackIndex] / ratio;

  Double_t dR = fCurrentHit->fR - param->fR;
  Double_t fi = param->fFi + (param->fAn    * dR + param->fAd     * dR*dR);
  Double_t z  = param->fZ  + (param->fTheta * dR + param->fThetaD * dR*dR);

  ((AliTRDhit *) fHit)->SetQ(param->fCharge[fCurrentHit->fStackIndex] / 2);   
  if ((param->fCharge[fCurrentHit->fStackIndex] % 2) ==0) {
    ((AliTRDhit *) fHit)->SetAmplification();
  }
  else {
    ((AliTRDhit *) fHit)->SetDrift();
  }
  ((AliTRDhit *) fHit)->SetX(fCurrentHit->fR * TMath::Cos(fi));
  ((AliTRDhit *) fHit)->SetY(fCurrentHit->fR * TMath::Sin(fi));
  ((AliTRDhit *) fHit)->SetZ(z);   
  ((AliTRDhit *) fHit)->SetDetector(param->fVolumeID);
  ((AliTRDhit *) fHit)->SetTrack(param->fTrackID);

  return kTRUE;

}
  

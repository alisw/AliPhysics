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

#include <TMath.h>
#include <TMatrixF.h>
#include <TStopwatch.h>

#include <AliLog.h>
#include <AliTPCROC.h>

#include "AliTPCCorrection.h"

#include "AliTPCCorrectionLookupTable.h"

ClassImp(AliTPCCorrectionLookupTable)

//_________________________________________________________________________________________
AliTPCCorrectionLookupTable::AliTPCCorrectionLookupTable()
: AliTPCCorrection()
, fNR(90)
, fNPhi(144)
, fNZ(130)
, fLimitsR()
, fLimitsPhi()
, fLimitsZ()
, fLookUpDxDist(0x0)
, fLookUpDyDist(0x0)
, fLookUpDzDist(0x0)
, fLookUpDxCorr(0x0)
, fLookUpDyCorr(0x0)
, fLookUpDzCorr(0x0)
{
  //
  //
  //
}
//_________________________________________________________________________________________
AliTPCCorrectionLookupTable::~AliTPCCorrectionLookupTable()
{
  //
  // dtor
  //

  ResetTables();
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Get interpolated Distortion
  //

  GetInterpolation(x,roc,dx,fLookUpDxDist,fLookUpDyDist,fLookUpDzDist);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Get interplolated correction
  //
  GetInterpolation(x,roc,dx,fLookUpDxCorr,fLookUpDyCorr,fLookUpDzCorr);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetInterpolation(const Float_t x[],const Short_t roc,Float_t dx[],
                                                   TMatrixF **mDx, TMatrixF **mDy, TMatrixF **mDz)
{
  //
  // Calculates the correction/distotring from a lookup table
  // type: 0 = correction
  //       1 = distortion
  //

//   Float_t typeSign=-1;
//   if (type==1) typeSign=1;
  
  Int_t   order     = 1 ;    // FIXME: hardcoded? Linear interpolation = 1, Quadratic = 2
  
  Double_t r, phi, z ;
  Int_t    sign;
  
  r      =  TMath::Sqrt( x[0]*x[0] + x[1]*x[1] ) ;
  phi    =  TMath::ATan2(x[1],x[0]) ;
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  z      =  x[2] ;                                         // Create temporary copy of x[2]
  
  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }
  
  if ( sign==1  && z <  fgkZOffSet ) z =  fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -fgkZOffSet ) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");
  
  // Get the Er and Ephi field integrals plus the integral over Z
    dx[0] = Interpolate3DTable(order, r, z, phi,
                               fNR, fNZ, fNPhi,
                               fLimitsR.GetMatrixArray(),
                               fLimitsZ.GetMatrixArray(),
                               fLimitsPhi.GetMatrixArray(),
                               mDx  );
    dx[1] = Interpolate3DTable(order, r, z, phi,
                               fNR, fNZ, fNPhi,
                               fLimitsR.GetMatrixArray(),
                               fLimitsZ.GetMatrixArray(),
                               fLimitsPhi.GetMatrixArray(),
                               mDy);
    dx[2] = Interpolate3DTable(order, r, z, phi,
                               fNR, fNZ, fNPhi,
                               fLimitsR.GetMatrixArray(),
                               fLimitsZ.GetMatrixArray(),
                               fLimitsPhi.GetMatrixArray(),
                               mDz   );
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::CreateLookupTable(AliTPCCorrection &tpcCorr, Float_t stepSize/*=5.*/)
{
  //
  //
  //

  // create distortion lookup table

  TStopwatch s;
  
  ResetTables();
  InitTables();
  
  const Float_t delta=stepSize; // 5cm
  Float_t x[3]={0.,0.,0.};
  Float_t dx[3]={0.,0.,0.};
  
  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    Double_t phi=fLimitsPhi(iPhi);
    //
    TMatrixF &mDxDist   = *fLookUpDxDist[iPhi];
    TMatrixF &mDyDist   = *fLookUpDyDist[iPhi];
    TMatrixF &mDzDist   = *fLookUpDzDist[iPhi];
    //
    TMatrixF &mDxCorr   = *fLookUpDxCorr[iPhi];
    TMatrixF &mDyCorr   = *fLookUpDyCorr[iPhi];
    TMatrixF &mDzCorr   = *fLookUpDzCorr[iPhi];
    
    for (Int_t ir=0; ir<fNR; ++ir){
      Double_t r=fLimitsR(ir);
      x[0]=r * TMath::Cos(phi);
      x[1]=r * TMath::Sin(phi);
      
      for (Int_t iz=0; iz<fNZ; ++iz){
        Double_t z=fLimitsZ(iz);
        x[2]=z;
        //TODO: change hardcoded value for r>133.?
        Int_t roc=TMath::Nint(phi*TMath::RadToDeg()/20.)%18;
        if (r>133.) roc+=36;
        if (z<0)    roc+=18;

        if (delta>0)
          tpcCorr.GetDistortionIntegralDz(x,roc,dx,delta);
        else
          tpcCorr.GetDistortion(x,roc,dx);
        mDxDist(ir,iz)=dx[0];
        mDyDist(ir,iz)=dx[1];
        mDzDist(ir,iz)=dx[2];

        if (delta>0)
          tpcCorr.GetCorrectionIntegralDz(x,roc,dx,delta);
        else
          tpcCorr.GetCorrection(x,roc,dx);
        mDxCorr(ir,iz)=dx[0];
        mDyCorr(ir,iz)=dx[1];
        mDzCorr(ir,iz)=dx[2];
        
      }
    }
  }

  s.Stop();
  AliInfo("Required time for lookup table creation: %.2f (%.2f) sec. real (cpu)");
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::InitTables()
{
  //
  //
  //
  fLookUpDxCorr = new TMatrixF*[fNPhi];
  fLookUpDyCorr = new TMatrixF*[fNPhi];
  fLookUpDzCorr = new TMatrixF*[fNPhi];
  
  fLookUpDxDist = new TMatrixF*[fNPhi];
  fLookUpDyDist = new TMatrixF*[fNPhi];
  fLookUpDzDist = new TMatrixF*[fNPhi];

  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    fLookUpDxCorr[iPhi] = new TMatrixF(fNR,fNZ);
    fLookUpDyCorr[iPhi] = new TMatrixF(fNR,fNZ);
    fLookUpDzCorr[iPhi] = new TMatrixF(fNR,fNZ);
    
    fLookUpDxDist[iPhi] = new TMatrixF(fNR,fNZ);
    fLookUpDyDist[iPhi] = new TMatrixF(fNR,fNZ);
    fLookUpDzDist[iPhi] = new TMatrixF(fNR,fNZ);
  }

  SetupDefaultLimits();
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::ResetTables()
{
  //
  // Reset the lookup tables
  //

  if (!fLookUpDxCorr) return;
  
  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    delete fLookUpDxCorr[iPhi];
    delete fLookUpDyCorr[iPhi];
    delete fLookUpDzCorr[iPhi];

    delete fLookUpDxDist[iPhi];
    delete fLookUpDyDist[iPhi];
    delete fLookUpDzDist[iPhi];
  }

  delete [] fLookUpDxCorr;
  delete [] fLookUpDyCorr;
  delete [] fLookUpDzCorr;
  
  delete [] fLookUpDxDist;
  delete [] fLookUpDyDist;
  delete [] fLookUpDzDist;
  
  fLookUpDxCorr   = 0x0;
  fLookUpDyCorr = 0x0;
  fLookUpDzCorr    = 0x0;
  
  fLookUpDxDist   = 0x0;
  fLookUpDyDist = 0x0;
  fLookUpDzDist    = 0x0;

  fLimitsR.ResizeTo(1);
  fLimitsPhi.ResizeTo(1);
  fLimitsZ.ResizeTo(1);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::SetupDefaultLimits()
{
  //
  // Set default limits for tables
  //

  fNR   = kNR;
  fNPhi = kNPhi;
  fNZ   = kNZ;
  fLimitsR.  ResizeTo(fNR);
  fLimitsPhi.ResizeTo(fNPhi);
  fLimitsZ.  ResizeTo(fNZ);
  fLimitsR.  SetElements(fgkRList);
  fLimitsPhi.SetElements(fgkPhiList);
  fLimitsZ.  SetElements(fgkZList);
}

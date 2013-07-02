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
, fLimitsRows()
, fLimitsPhiSlices()
, fLimitsColumns()
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

  GetInterpolation(x,roc,dx,fLookUpDxDist,fLookUpDyDist,fLookUpDzDist,1);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Get interplolated correction
  //
  GetInterpolation(x,roc,dx,fLookUpDxCorr,fLookUpDyCorr,fLookUpDzCorr,0);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetInterpolation(const Float_t x[],const Short_t roc,Float_t dx[],
                                                   TMatrixF **mDx, TMatrixF **mDy, TMatrixF **mDz,
                                                   Int_t type)
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
                               fLimitsRows.GetMatrixArray(),
                               fLimitsColumns.GetMatrixArray(),
                               fLimitsPhiSlices.GetMatrixArray(),
                               mDx  );
    dx[1] = Interpolate3DTable(order, r, z, phi,
                               fNR, fNZ, fNPhi,
                               fLimitsRows.GetMatrixArray(),
                               fLimitsColumns.GetMatrixArray(),
                               fLimitsPhiSlices.GetMatrixArray(),
                               mDy);
    dx[2] = Interpolate3DTable(order, r, z, phi,
                               fNR, fNZ, fNPhi,
                               fLimitsRows.GetMatrixArray(),
                               fLimitsColumns.GetMatrixArray(),
                               fLimitsPhiSlices.GetMatrixArray(),
                               mDz   );
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::CreateLookupTable(AliTPCCorrection &tpcCorr, Int_t /*rows*//*=90*/, Int_t /*phiSlices*//*=144*/, Int_t /*columns*//*=130*/ )
{
  //
  //
  //

  // create distortion lookup table

  const Float_t delta=5.; // 5cm
  Float_t x[3]={0.,0.,0.};
  Float_t dx[3]={0.,0.,0.};
  
  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    Double_t phi=fLimitsPhiSlices(iPhi);
    //
    TMatrixF &mDxDist   = *fLookUpDxDist[iPhi];
    TMatrixF &mDyDist   = *fLookUpDyDist[iPhi];
    TMatrixF &mDzDist   = *fLookUpDzDist[iPhi];
    //
    TMatrixF &mDxCorr   = *fLookUpDxCorr[iPhi];
    TMatrixF &mDyCorr   = *fLookUpDyCorr[iPhi];
    TMatrixF &mDzCorr   = *fLookUpDzCorr[iPhi];
    
    for (Int_t ir=0; ir<fNR; ++ir){
      Double_t r=fLimitsRows(ir);
      x[0]=r * TMath::Cos(phi);
      x[1]=r * TMath::Sin(phi);
      
      for (Int_t iz=0; iz<fNZ; ++iz){
        Double_t z=fLimitsColumns(iz);
        x[2]=z;
        //TODO: change hardcoded value for r>133.?
        Int_t roc=TMath::Nint(phi*TMath::RadToDeg()/20.)%18;
        if (r>133.) roc+=36;
        if (z<0)    roc+=18;

        tpcCorr.GetDistortionIntegralDz(x,roc,dx,delta);
        mDxDist(ir,iz)=dx[0];
        mDyDist(ir,iz)=dx[1];
        mDzDist(ir,iz)=dx[2];
        
        tpcCorr.GetCorrectionIntegralDz(x,roc,dx,delta);
        mDxCorr(ir,iz)=dx[0];
        mDyCorr(ir,iz)=dx[1];
        mDzCorr(ir,iz)=dx[2];
        
      }
    }
  }

  // create correction lookup table
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

  fLimitsRows.ResizeTo(1);
  fLimitsPhiSlices.ResizeTo(1);
  fLimitsColumns.ResizeTo(1);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::SetupDefaultLimits()
{
  //
  // Set default limits for tables
  //

  fLimitsRows.ResizeTo(fNR);
  fLimitsPhiSlices.ResizeTo(fNPhi);
  fLimitsColumns.ResizeTo(fNZ);

  Double_t *limRList   = fLimitsRows.GetMatrixArray();
  Double_t *limPhiList = fLimitsPhiSlices.GetMatrixArray();
  Double_t *limZList   = fLimitsColumns.GetMatrixArray();
  
  AliTPCROC * roc = AliTPCROC::Instance();
  const Double_t rLow =  TMath::Floor(roc->GetPadRowRadii(0,0))-1; // first padRow plus some margin
  
  // fulcrums in R
  limRList[0] = rLow;
  for (Int_t i = 1; i<fNR; i++) {
    limRList[i] = limRList[i-1] + 3.5;     // 3.5 cm spacing
    if (limRList[i]<90 ||limRList[i]>245){
      limRList[i] = limRList[i-1] + 0.5; // 0.5 cm spacing
    } else if (limRList[i]<100 || limRList[i]>235){
      limRList[i] = limRList[i-1] + 1.5;  // 1.5 cm spacing
    } else if (limRList[i]<120 || limRList[i]>225){
      limRList[i] = limRList[i-1] + 2.5;  // 2.5 cm spacing
    }
  }

  // fulcrums in Z
  limZList[0] = -249.5;
  limZList[fNZ-1] = 249.5;
  for (Int_t j = 1; j<fNZ/2; j++) {
    limZList[j] = limZList[j-1];
    if      (TMath::Abs(limZList[j])< 0.15){
      limZList[j] = limZList[j-1] + 0.09; // 0.09 cm spacing
    } else if(TMath::Abs(limZList[j])< 0.6){
      limZList[j] = limZList[j-1] + 0.4; // 0.4 cm spacing
    } else if      (TMath::Abs(limZList[j])< 2.5 || TMath::Abs(limZList[j])>248){
      limZList[j] = limZList[j-1] + 0.5; // 0.5 cm spacing
    } else if (TMath::Abs(limZList[j])<10 || TMath::Abs(limZList[j])>235){
      limZList[j] = limZList[j-1] + 1.5;  // 1.5 cm spacing
    } else if (TMath::Abs(limZList[j])<25 || TMath::Abs(limZList[j])>225){
      limZList[j] = limZList[j-1] + 2.5;  // 2.5 cm spacing
    } else{
      limZList[j] = limZList[j-1] + 4;  // 4 cm spacing
    }

    limZList[fNZ-j-1] = -limZList[j];
  }
  
  // fulcrums in phi
  for (Int_t k = 0; k<fNPhi; k++)
    limPhiList[k] = TMath::TwoPi()*k/(fNPhi-1);
  
}

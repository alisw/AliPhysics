/// \class AliTPCCorrectionLookupTable

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
#include <TString.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <THashList.h>
#include <TVector2.h>
#include <TLinearFitter.h>

#include <AliLog.h>
#include <AliTPCROC.h>

#include "AliTPCCorrection.h"

#include "AliTPCCorrectionLookupTable.h"

/// \cond CLASSIMP
ClassImp(AliTPCCorrectionLookupTable)
/// \endcond

//_________________________________________________________________________________________
AliTPCCorrectionLookupTable::AliTPCCorrectionLookupTable()
: AliTPCCorrection()
, fNR(0)
, fNPhi(0)
, fNZ(0)
, fCorrScaleFactor(-1)
, fFillCorrection(kTRUE)
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
  /// dtor

  ResetTables();
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]) {
  /// Get interpolated Distortion

  GetInterpolation(x,roc,dx,fLookUpDxDist,fLookUpDyDist,fLookUpDzDist);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  /// Get interplolated correction

  GetInterpolation(x,roc,dx,fLookUpDxCorr,fLookUpDyCorr,fLookUpDzCorr);

  if (fCorrScaleFactor>0) {
    dx[0]*=fCorrScaleFactor;
    dx[1]*=fCorrScaleFactor;
    dx[2]*=fCorrScaleFactor;
  }
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::GetInterpolation(const Float_t x[],const Short_t roc,Float_t dx[],
                                                   TMatrixF **mDx, TMatrixF **mDy, TMatrixF **mDz)
{
  /// Calculates the correction/distotring from a lookup table
  /// type: 0 = correction
  ///       1 = distortion

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
  /// create lookup table for all phi,r,z bins

  if (fNR==0) {
    AliError("Limits are not set yet. Please use one of the Set..Limits functions first");
    return;
  }

  TStopwatch s;

  ResetTables();
  InitTables();

  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    CreateLookupTablePhiBin(tpcCorr,iPhi,stepSize);
  }

  s.Stop();
  AliInfo(Form("Required time for lookup table creation: %.2f (%.2f) sec. real (cpu)",s.RealTime(),s.CpuTime()));
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::CreateLookupTableSinglePhi(AliTPCCorrection &tpcCorr, Int_t iPhi, Float_t stepSize)
{
  /// Lookup table for only one phi bin. Can be used for parallel processing

  if (fNR==0) {
    AliError("Limits are not set yet. Please use one of the Set..Limits functions first");
    return;
  }

  TStopwatch s;

  ResetTables();
  InitTableArrays();
  InitTablesPhiBin(iPhi);
  CreateLookupTablePhiBin(tpcCorr,iPhi,stepSize);

  s.Stop();
  AliInfo(Form("Required time for lookup table creation: %.2f (%.2f) sec. real (cpu)",s.RealTime(),s.CpuTime()));
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::CreateLookupTablePhiBin(AliTPCCorrection &tpcCorr, Int_t iPhi, Float_t stepSize)
{
  ///

  if (iPhi<0||iPhi>=fNPhi) return;

  const Float_t delta=stepSize; // 5cm
  Float_t x[3]={0.,0.,0.};
  Float_t dx[3]={0.,0.,0.};

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

      if (fFillCorrection) {
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

}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::InitTables()
{
  /// Init all tables

  InitTableArrays();
  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    InitTablesPhiBin(iPhi);
  }
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::CreateLookupTableFromResidualDistortion(THn &resDist)
{
  /// create lookup table from residual distortions stored in a 3d histogram
  /// assume dimensions are r, phi, z

  if (fNR==0) {
    AliError("Limits are not set yet. Please use one of the Set..Limits functions first");
    return;
  }

  ResetTables();
  InitTables();

  Double_t x[3]={0.,0.,0.};

  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    const Double_t phi=fLimitsPhi(iPhi);
    x[1]=phi;
    //
    TMatrixF &mDxDist   = *fLookUpDxDist[iPhi];
    TMatrixF &mDyDist   = *fLookUpDyDist[iPhi];
    TMatrixF &mDzDist   = *fLookUpDzDist[iPhi];
    //
    TMatrixF &mDxCorr   = *fLookUpDxCorr[iPhi];
    TMatrixF &mDyCorr   = *fLookUpDyCorr[iPhi];
    TMatrixF &mDzCorr   = *fLookUpDzCorr[iPhi];

    for (Int_t ir=0; ir<fNR; ++ir){
      const Double_t r=fLimitsR(ir);
      x[0]=r;

      for (Int_t iz=0; iz<fNZ; ++iz){
        const Double_t z=fLimitsZ(iz);
        x[2]=z;

        const Double_t drphi = resDist.GetBinContent(resDist.GetBin(x));
        Double_t dx[3]={0.,drphi,0.};

        // transform rphi distortions (local y, so dy') to a global distortion
        // assume no radial distortion (dx' = 0)
        // assume no residual distortion in z for the moment
        Double_t cs=TMath::Cos(phi), sn=TMath::Sin(phi), lx=dx[0];
        dx[0]=lx*cs - dx[1]*sn; dx[1]=lx*sn + dx[1]*cs;

        mDxDist(ir,iz)=dx[0];
        mDyDist(ir,iz)=dx[1];
        mDzDist(ir,iz)=dx[2];

        mDxCorr(ir,iz)=-dx[0];
        mDyCorr(ir,iz)=-dx[1];
        mDzCorr(ir,iz)=-dx[2];
      }
    }
  }
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::CreateResidual(AliTPCCorrection *distortion, AliTPCCorrection* correction)
{
  /// create lookup table from residual distortions calculated from distorted - correction

  ResetTables();
  InitTables();

  Float_t x[3]={0.,0.,0.};

  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    const Double_t phi=fLimitsPhi(iPhi);
    //
    TMatrixF &mDxDist   = *fLookUpDxDist[iPhi];
    TMatrixF &mDyDist   = *fLookUpDyDist[iPhi];
    TMatrixF &mDzDist   = *fLookUpDzDist[iPhi];
    //
    TMatrixF &mDxCorr   = *fLookUpDxCorr[iPhi];
    TMatrixF &mDyCorr   = *fLookUpDyCorr[iPhi];
    TMatrixF &mDzCorr   = *fLookUpDzCorr[iPhi];

    for (Int_t ir=0; ir<fNR; ++ir){
      const Double_t r=fLimitsR(ir);
      x[0]=r * TMath::Cos(phi);
      x[1]=r * TMath::Sin(phi);

      for (Int_t iz=0; iz<fNZ; ++iz){
        const Double_t z=fLimitsZ(iz);
        x[2]=z;

        //original point
        Float_t xdc[3]={x[0], x[1], x[2]};

        Int_t roc=TMath::Nint(phi*TMath::RadToDeg()/20.)%18;
        if (r>133.) roc+=36;
        if (z<0)    roc+=18;

        //get residual distortion
        distortion->DistortPoint(xdc, roc);
        correction->CorrectPoint(xdc, roc);
        Float_t dx[3]={xdc[0]-x[0], xdc[1]-x[1], xdc[2]-x[2]};

        mDxDist(ir,iz)=dx[0];
        mDyDist(ir,iz)=dx[1];
        mDzDist(ir,iz)=dx[2];

        mDxCorr(ir,iz)=-dx[0];
        mDyCorr(ir,iz)=-dx[1];
        mDzCorr(ir,iz)=-dx[2];
      }
    }
  }
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::InitTablesPhiBin(Int_t iPhi)
{
  ///

  // check if already initialised
  if (iPhi<0||iPhi>=fNPhi) return;
  if (fLookUpDxCorr[iPhi]) return;

  fLookUpDxCorr[iPhi] = new TMatrixF(fNR,fNZ);
  fLookUpDyCorr[iPhi] = new TMatrixF(fNR,fNZ);
  fLookUpDzCorr[iPhi] = new TMatrixF(fNR,fNZ);

  fLookUpDxDist[iPhi] = new TMatrixF(fNR,fNZ);
  fLookUpDyDist[iPhi] = new TMatrixF(fNR,fNZ);
  fLookUpDzDist[iPhi] = new TMatrixF(fNR,fNZ);

}
//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::InitTableArrays()
{
  ///

  // needs to be before the table creation to set the limits
  SetupDefaultLimits();

  fLookUpDxCorr = new TMatrixF*[fNPhi];
  fLookUpDyCorr = new TMatrixF*[fNPhi];
  fLookUpDzCorr = new TMatrixF*[fNPhi];

  fLookUpDxDist = new TMatrixF*[fNPhi];
  fLookUpDyDist = new TMatrixF*[fNPhi];
  fLookUpDzDist = new TMatrixF*[fNPhi];

  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    fLookUpDxCorr[iPhi] = 0x0;
    fLookUpDyCorr[iPhi] = 0x0;
    fLookUpDzCorr[iPhi] = 0x0;

    fLookUpDxDist[iPhi] = 0x0;
    fLookUpDyDist[iPhi] = 0x0;
    fLookUpDzDist[iPhi] = 0x0;
  }
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::ResetTables()
{
  /// Reset the lookup tables

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
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::SetupDefaultLimits()
{
  /// Set default limits for tables

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

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::ResetLimits()
{
  fNR   = 0;
  fNPhi = 0;
  fNZ   = 0;

  fLimitsR.  ResizeTo(1);
  fLimitsPhi.ResizeTo(1);
  fLimitsZ.  ResizeTo(1);
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::MergePhiTables(const char* files)
{
  /// merge all lookup tables stored in 'files' with this one
  /// assume that each lookup table in each file has only one phi bin

  ResetTables();
  ResetLimits(); // use limits from the first file assuming they are all the same

  TString sfiles=gSystem->GetFromPipe(Form("ls %s",files));
  TObjArray *arrFiles=sfiles.Tokenize("\n");

  for (Int_t i=0; i<arrFiles->GetEntriesFast();++i){
    TFile f(arrFiles->At(i)->GetName());
    if (!f.IsOpen() || f.IsZombie()) continue;
    AliTPCCorrectionLookupTable *lt=dynamic_cast<AliTPCCorrectionLookupTable*>(f.Get(f.GetListOfKeys()->First()->GetName()));
    if (!lt) {
      AliError(Form("first object in file '%s' is not of type AliTPCCorrectionLookupTable!",f.GetName()));
      continue;
    }

    if (!fNR) {
      fNR        = lt->fNR;
      fNPhi      = lt->fNPhi;
      fNZ        = lt->fNZ;
      fLimitsR   = lt->fLimitsR;
      fLimitsZ   = lt->fLimitsZ;
      fLimitsPhi = lt->fLimitsPhi;
      InitTableArrays();
    } else {
      if (fNR!=lt->fNR || fNPhi!=lt->fNPhi || fNZ!=lt->fNZ ){
        AliError(Form("Current limits don't macht the ones in file '%s'",f.GetName()));
        continue;
      }
    }

    for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi) {
      if (!lt->fLookUpDxCorr[iPhi]) continue;

      AliInfo(Form("Adding phi bin '%d' from file '%s'",iPhi,f.GetName()));

      InitTablesPhiBin(iPhi);

      *fLookUpDxDist[iPhi] = *lt->fLookUpDxDist[iPhi];
      *fLookUpDyDist[iPhi] = *lt->fLookUpDyDist[iPhi];
      *fLookUpDzDist[iPhi] = *lt->fLookUpDzDist[iPhi];
      //
      *fLookUpDxCorr[iPhi] = *lt->fLookUpDxCorr[iPhi];
      *fLookUpDyCorr[iPhi] = *lt->fLookUpDyCorr[iPhi];
      *fLookUpDzCorr[iPhi] = *lt->fLookUpDzCorr[iPhi];
      break;
    }
  }

  //check of all phi bins are initialised
  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi) {
    if (!fLookUpDxCorr[iPhi]) {
      AliFatal(Form("Phi bin '%d' not initialised from files!",iPhi));
    }
  }

  delete arrFiles;
}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::BuildExactInverse()
{
  /// this method build the exact inverse of the standard distortion map
  /// for the the distortion man first needs to be calculated
  /// then the correction map will be overwritten

  Float_t x[3]    = {0.,0.,0.};
  Float_t x2[3]   = {0.,0.,0.};
  Float_t xref[3] = {0.,0.,0.};
  Float_t xd[3]   = {0.,0.,0.};
  Float_t dx[3]   = {0.,0.,0.};

  // reset correction matrices
  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    TMatrixF &mDxCorr   = *fLookUpDxCorr[iPhi];
    TMatrixF &mDyCorr   = *fLookUpDyCorr[iPhi];
    TMatrixF &mDzCorr   = *fLookUpDzCorr[iPhi];

    for (Int_t ir=0; ir<fNR; ++ir){
      for (Int_t iz=0; iz<fNZ; ++iz){
        mDxCorr(ir,iz) = -1000.;
        mDyCorr(ir,iz) = -1000.;
        mDzCorr(ir,iz) = -1000.;
      }
    }
  }

  // get interplolated corrections on standard grid
  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    Double_t phi=fLimitsPhi(iPhi);
    TMatrixF &mDxDist   = *fLookUpDxDist[iPhi];
    TMatrixF &mDyDist   = *fLookUpDyDist[iPhi];
    TMatrixF &mDzDist   = *fLookUpDzDist[iPhi];

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

        dx[0] = mDxDist(ir,iz);
        dx[1] = mDyDist(ir,iz);
        dx[2] = mDzDist(ir,iz);

        xd[0] = x[0]+dx[0];
        xd[1] = x[1]+dx[1];
        xd[2] = x[2]+dx[2];

        const Double_t phid = TVector2::Phi_0_2pi(TMath::ATan2(xd[1],xd[0]));
        const Double_t rd   = TMath::Sqrt(xd[0]*xd[0] + xd[1]*xd[1]);
        const Double_t zd   = xd[2];

        Int_t ilow = 0, jlow = 0, klow = 0 ;

        Search( fLimitsR.GetNrows(),   fLimitsR.GetMatrixArray(),   rd,   ilow   ) ;
        Search( fLimitsZ.GetNrows(),   fLimitsZ.GetMatrixArray(),   zd,   jlow   ) ;
        Search( fLimitsPhi.GetNrows(), fLimitsPhi.GetMatrixArray(), phid, klow   ) ;

        if ( ilow < 0 ) ilow = 0 ;   // check if out of range
        if ( jlow < 0 ) jlow = 0 ;
        if ( klow < 0 ) klow = 0 ;
        if ( ilow >= fLimitsR.GetNrows())   ilow = fLimitsR.GetNrows() - 1;
        if ( jlow >= fLimitsZ.GetNrows())   jlow = fLimitsZ.GetNrows() - 1;
        if ( klow >= fLimitsPhi.GetNrows()) klow = fLimitsPhi.GetNrows() - 1;

        const Double_t phiRef = fLimitsPhi[klow];
        const Double_t rRef   = fLimitsR[ilow];
        const Double_t zRef   = fLimitsZ[jlow];

        TMatrixF &mDxCorr   = *fLookUpDxCorr[klow];
        if ( mDxCorr(ilow, jlow) > -1000. ) continue;
        TMatrixF &mDyCorr   = *fLookUpDyCorr[klow];
        TMatrixF &mDzCorr   = *fLookUpDzCorr[klow];

        xref[0]= rRef * TMath::Cos(phiRef);
        xref[1]= rRef * TMath::Sin(phiRef);
        xref[2]= zRef;

        FindClosestPosition(ir,iz,iPhi, xref, x2);

        GetDistortion(x2,roc,dx);

        mDxCorr(ilow, jlow) = -dx[0];
        mDyCorr(ilow, jlow) = -dx[1];
        mDzCorr(ilow, jlow) = -dx[2];

//         printf("%3d %3d %3d\n",iPhi, ir, iz);
//         printf("%3d %3d %3d\n",klow, ilow, jlow);
//         printf("x2:   %.5f %.5f %.5f\n", x2[0], x2[1], x2[2]);
//         printf("x2d:  %.5f %.5f %.5f\n", x2[0]+dx[0], x2[1]+dx[1], x2[2]+dx[2]);
//         printf("xref: %.5f %.5f %.5f\n", xref[0], xref[1], xref[2]);
//         printf("xrd:  %.5f %.5f %.5f\n", x2[0]+dx[0]-xref[0], x2[1]+dx[1]-xref[1], x2[2]+dx[2]-xref[2]);
//         printf("phid: %.5f %.5f %.5f\n", phid,rd,zd);
//         printf("phir: %.5f %.5f %.5f\n", phiRef,rRef,zRef);
//         printf("\n");
      }
    }
  }

  // fill remaining empty bins
  // The last ein first phi bin entries must be identical, fill those first
  {
  TMatrixF &mDxCorr   = *fLookUpDxCorr[0];
  TMatrixF &mDyCorr   = *fLookUpDyCorr[0];
  TMatrixF &mDzCorr   = *fLookUpDzCorr[0];

  TMatrixF &mDxCorr2  = *fLookUpDxCorr[fNPhi-1];
  TMatrixF &mDyCorr2  = *fLookUpDyCorr[fNPhi-1];
  TMatrixF &mDzCorr2  = *fLookUpDzCorr[fNPhi-1];

  for (Int_t ir=0; ir<fNR; ++ir){
    for (Int_t iz=0; iz<fNZ; ++iz){
      mDxCorr2(ir,iz) = mDxCorr(ir,iz);
      mDyCorr2(ir,iz) = mDyCorr(ir,iz);
      mDzCorr2(ir,iz) = mDzCorr(ir,iz);
    }
  }
  }

  for (Int_t iPhi=0; iPhi<fNPhi; ++iPhi){
    TMatrixF &mDxCorr   = *fLookUpDxCorr[iPhi];
    TMatrixF &mDyCorr   = *fLookUpDyCorr[iPhi];
    TMatrixF &mDzCorr   = *fLookUpDzCorr[iPhi];

    Double_t phi=fLimitsPhi(iPhi);
    for (Int_t ir=0; ir<fNR; ++ir){
      Double_t r=fLimitsR(ir);
      x[0]=r * TMath::Cos(phi);
      x[1]=r * TMath::Sin(phi);

      for (Int_t iz=0; iz<fNZ; ++iz){
        if (mDxCorr(ir,iz) > -999.) continue;

        Double_t z=fLimitsZ(iz);
        x[2]=z;

        //TODO: change hardcoded value for r>133.?
        Int_t roc=TMath::Nint(phi*TMath::RadToDeg()/20.)%18;
        if (r>133.) roc+=36;
        if (z<0)    roc+=18;

        // get last point
        dx[0] = mDxCorr(ir,iz-1);
        dx[1] = mDyCorr(ir,iz-1);
        dx[2] = mDzCorr(ir,iz-1);

        xd[0] = x[0]+dx[0];
        xd[1] = x[1]+dx[1];
        xd[2] = x[2]+dx[2];

        // get distorted point
        const Double_t phid = TVector2::Phi_0_2pi(TMath::ATan2(xd[1],xd[0]));
        const Double_t rd   = TMath::Sqrt(xd[0]*xd[0] + xd[1]*xd[1]);
        const Double_t zd   = xd[2];

        Int_t ilow = 0, jlow = 0, klow = 0 ;

        Search( fLimitsR.GetNrows(),   fLimitsR.GetMatrixArray(),   rd,   ilow   ) ;
        Search( fLimitsZ.GetNrows(),   fLimitsZ.GetMatrixArray(),   zd,   jlow   ) ;
        Search( fLimitsPhi.GetNrows(), fLimitsPhi.GetMatrixArray(), phid, klow   ) ;

        if ( ilow < 0 ) ilow = 0 ;   // check if out of range
        if ( jlow < 0 ) jlow = 0 ;
        if ( klow < 0 ) klow = 0 ;
        if ( ilow >= fLimitsR.GetNrows())   ilow = fLimitsR.GetNrows() - 1;
        if ( jlow >= fLimitsZ.GetNrows())   jlow = fLimitsZ.GetNrows() - 1;
        if ( klow >= fLimitsPhi.GetNrows()) klow = fLimitsPhi.GetNrows() - 1;

        FindClosestPosition(ilow,jlow,klow, x, x2);

        GetDistortion(x2,roc,dx);

        mDxCorr(ir, iz) = -dx[0];
        mDyCorr(ir, iz) = -dx[1];
        mDzCorr(ir, iz) = -dx[2];
      }
    }
  }

}

//_________________________________________________________________________________________
void AliTPCCorrectionLookupTable::FindClosestPosition(const Int_t binR, const Int_t binZ, const Int_t binPhi,
                                                      const Float_t xref[3], Float_t xret[3])
{
  ///

//   static TLinearFitter fitx(2,"pol2");
//   static TLinearFitter fity(2,"pol2");
//   static TLinearFitter fitz(2,"pol2");
  static TLinearFitter fitx(4,"hyp3");
  static TLinearFitter fity(4,"hyp3");
  static TLinearFitter fitz(4,"hyp3");
  fitx.ClearPoints();
  fity.ClearPoints();
  fitz.ClearPoints();

  const Int_t nPoints=3;
  Int_t counter=0;
  Int_t rMin=binR;
  Int_t zMin=binZ;
  Int_t phiMin=binPhi;

  counter=nPoints/2;
  while (rMin>0 && counter--) --rMin;
  counter=nPoints/2;
  while (zMin>0 && counter--) --zMin;
  counter=nPoints/2;
  while (phiMin>0 && counter--) --phiMin;

  Int_t rMax   = rMin  +nPoints;
  Int_t zMax   = zMin  +nPoints;
  Int_t phiMax = phiMin+nPoints;

  while (rMax>=fNR) {--rMin; --rMax;}
  while (zMax>=fNZ) {--zMin; --zMax;}
  while (phiMax>=fNPhi) {--phiMin; --phiMax;}

  Float_t  x[3]    = {0.,0.,0.};
  Double_t xd[3]   = {0.,0.,0.};
  Float_t  dx[3]   = {0.,0.,0.};

  for (Int_t iPhi=phiMin; iPhi<phiMax; ++iPhi) {
    TMatrixF &mDxDist   = *fLookUpDxDist[iPhi];
    TMatrixF &mDyDist   = *fLookUpDyDist[iPhi];
    TMatrixF &mDzDist   = *fLookUpDzDist[iPhi];

    Double_t phi=fLimitsPhi(iPhi);
    for (Int_t ir=rMin; ir<rMax; ++ir){
      Double_t r=fLimitsR(ir);
      x[0]=r * TMath::Cos(phi);
      x[1]=r * TMath::Sin(phi);

      for (Int_t iz=zMin; iz<zMax; ++iz){
        Double_t z=fLimitsZ(iz);
        x[2]=z;

        dx[0] = mDxDist(ir,iz);
        dx[1] = mDyDist(ir,iz);
        dx[2] = mDzDist(ir,iz);

        xd[0] = x[0]+dx[0];
        xd[1] = x[1]+dx[1];
        xd[2] = x[2]+dx[2];

        fitx.AddPoint(xd,   x[0]);
        fity.AddPoint(xd, x[1]);
        fitz.AddPoint(xd, x[2]);
      }
    }
  }

  fitx.Eval();
  fity.Eval();
  fitz.Eval();
  xret[0] = fitx.GetParameter(0) + fitx.GetParameter(1)*xref[0]
                                 + fitx.GetParameter(2)*xref[1]
                                 + fitx.GetParameter(3)*xref[2];
  xret[1] = fity.GetParameter(0) + fity.GetParameter(1)*xref[0]
                                 + fity.GetParameter(2)*xref[1]
                                 + fity.GetParameter(3)*xref[2];
  xret[2] = fitz.GetParameter(0) + fitz.GetParameter(1)*xref[0]
                                 + fitz.GetParameter(2)*xref[1]
                                 + fitz.GetParameter(3)*xref[2];
//   xret[0] = fitx.GetParameter(0) + fitx.GetParameter(1)*xref[0] + fitx.GetParameter(2)*xref[0]*xref[0];
//   xret[1] = fity.GetParameter(0) + fity.GetParameter(1)*xref[1] + fity.GetParameter(2)*xref[1]*xref[1];
//   xret[2] = fitz.GetParameter(0) + fitz.GetParameter(1)*xref[2] + fitz.GetParameter(2)*xref[2]*xref[2];
}


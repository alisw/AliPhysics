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

#include <Riostream.h>
#include <TMath.h>
#include <TVector2.h>

#include "AliESDtrack.h"
#include "AliTRDgeometry.h" 
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "AliTRDtracklet.h"

ClassImp(AliTRDtrack)

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Represents a reconstructed TRD track                                     //
//  Local TRD Kalman track                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack()
  :AliKalmanTrack()
  ,fSeedLab(-1)
  ,fdEdx(0)
  ,fdEdxT(0)
  ,fDE(0)
  ,fAlpha(0)
  ,fX(0)
  ,fStopped(kFALSE)
  ,fY(0)
  ,fZ(0)
  ,fE(0)
  ,fT(0)
  ,fC(0)
  ,fCyy(1e10)
  ,fCzy(0)
  ,fCzz(1e10)
  ,fCey(0)
  ,fCez(0)
  ,fCee(1e10)
  ,fCty(0)
  ,fCtz(0)
  ,fCte(0)
  ,fCtt(1e10)
  ,fCcy(0)
  ,fCcz(0)
  ,fCce(0)
  ,fCct(0)
  ,fCcc(1e10)
  ,fLhElectron(0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(1e10)
  ,fBackupTrack(0x0)
{
  //
  // AliTRDtrack default constructor
  //

  Int_t  i = 0;
  Int_t  j = 0;
  UInt_t k = 0;

  for (i = 0; i < kNplane; i++) {
    for (j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = 0;
    }
    fTimBinPlane[i] = -1;
  }

  for (k = 0; k < kMAXCLUSTERSPERTRACK; k++) {
    fIndex[k]       = 0;
    fIndexBackup[k] = 0;
    fdQdl[k]        = 0;
  }

  for (i = 0; i < 3; i++) {
    fBudget[i]      = 0;
  }

}

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliTRDcluster *c, UInt_t index, 
                         const Double_t xx[5], const Double_t cc[15], 
                         Double_t xref, Double_t alpha) 
  :AliKalmanTrack() 
  ,fSeedLab(-1)
  ,fdEdx(0.0)
  ,fdEdxT(0.0)
  ,fDE(0.0)
  ,fAlpha(alpha)
  ,fX(xref)
  ,fStopped(kFALSE)
  ,fY(xx[0])
  ,fZ(xx[1])
  ,fE(xx[2])
  ,fT(xx[3])
  ,fC(xx[4])
  ,fCyy(cc[0])
  ,fCzy(cc[1])
  ,fCzz(cc[2])
  ,fCey(cc[3])  
  ,fCez(cc[4])  
  ,fCee(cc[5])
  ,fCty(cc[6])  
  ,fCtz(cc[7])  
  ,fCte(cc[8])  
  ,fCtt(cc[9])
  ,fCcy(cc[10]) 
  ,fCcz(cc[11]) 
  ,fCce(cc[12]) 
  ,fCct(cc[13]) 
  ,fCcc(cc[14])  
  ,fLhElectron(0.0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(1e10)
  ,fBackupTrack(0x0)
{
  //
  // AliTRDtrack main constructor
  //

  Int_t  i = 0;
  Int_t  j = 0;
  UInt_t k = 0;

  if (fAlpha <  -TMath::Pi()) fAlpha += 2.0 * TMath::Pi();
  if (fAlpha >=  TMath::Pi()) fAlpha -= 2.0 * TMath::Pi();   
 
  SaveLocalConvConst();

  fIndex[0] = index;
  SetNumberOfClusters(1);

  for (i = 0; i < kNplane; i++) {
    for (j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = 0;
    }
    fTimBinPlane[i] = -1;
  }

  Double_t q = TMath::Abs(c->GetQ());
  Double_t s = fX * fC - fE;
  Double_t t = fT;
  if (s*s < 1.0) {
    q *= TMath::Sqrt((1-s*s)/(1+t*t));
  }
  fdQdl[0] = q;
  
  for (k = 1; k < kMAXCLUSTERSPERTRACK; k++) {
    fdQdl[k]        = 0;
    fIndex[k]       = 0;
    fIndexBackup[k] = 0; 
  }

  for (i = 0; i < 3; i++) { 
    fBudget[i]      = 0;
  }

}                              
           
//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliTRDtrack &t) 
  :AliKalmanTrack(t) 
  ,fSeedLab(t.fSeedLab)
  ,fdEdx(t.fdEdx)
  ,fdEdxT(t.fdEdxT)
  ,fDE(t.fDE)
  ,fAlpha(t.fAlpha)
  ,fX(t.fX)
  ,fStopped(t.fStopped)
  ,fY(t.fY)
  ,fZ(t.fZ)
  ,fE(t.fE)
  ,fT(t.fT)
  ,fC(t.fC)
  ,fCyy(t.fCyy)
  ,fCzy(t.fCzy)
  ,fCzz(t.fCzz)
  ,fCey(t.fCey)  
  ,fCez(t.fCez)  
  ,fCee(t.fCee)
  ,fCty(t.fCty)  
  ,fCtz(t.fCtz)  
  ,fCte(t.fCte)  
  ,fCtt(t.fCtt)
  ,fCcy(t.fCcy) 
  ,fCcz(t.fCcz) 
  ,fCce(t.fCce) 
  ,fCct(t.fCct)
  ,fCcc(t.fCcc)  
  ,fLhElectron(0.0)
  ,fNWrong(t.fNWrong)
  ,fNRotate(t.fNRotate)
  ,fNCross(t.fNCross)
  ,fNExpected(t.fNExpected)
  ,fNLast(t.fNLast)
  ,fNExpectedLast(t.fNExpectedLast)
  ,fNdedx(t.fNdedx)
  ,fChi2Last(t.fChi2Last)
  ,fBackupTrack(0x0)
{
  //
  // Copy constructor.
  //

  Int_t  i = 0;
  Int_t  j = 0;
  UInt_t k = 0;

  for (i = 0; i < kNplane; i++) {
    for (j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = t.fdEdxPlane[i][j];
    }
    fTimBinPlane[i] = t.fTimBinPlane[i];
    fTracklets[i]   = t.fTracklets[i];
  }

  Int_t n = t.GetNumberOfClusters(); 
  for (i = 0; i < n; i++) {
    fIndex[i]       = t.fIndex[i];
    fIndexBackup[i] = t.fIndex[i];
    fdQdl[i]        = t.fdQdl[i];
  }
  for (k = n; k < kMAXCLUSTERSPERTRACK; k++) {
    fdQdl[k]        = 0;
    fIndex[k]       = 0;
    fIndexBackup[k] = 0; 
  }

  for (i = 0; i < 6; i++) {
    fTracklets[i] = t.fTracklets[i];
  }

  for (i = 0; i < 3; i++) { 
    fBudget[i]    = t.fBudget[i];
  }

}                                

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliKalmanTrack &t, Double_t alpha) 
  :AliKalmanTrack(t) 
  ,fSeedLab(-1)
  ,fdEdx(t.GetPIDsignal())
  ,fdEdxT(0)
  ,fDE(0)
  ,fAlpha(alpha)
  ,fX(0)
  ,fStopped(kFALSE)
  ,fY(0)
  ,fZ(0)
  ,fE(0)
  ,fT(0)
  ,fC(0)
  ,fCyy(0)
  ,fCzy(0)
  ,fCzz(0)
  ,fCey(0)
  ,fCez(0)
  ,fCee(0)
  ,fCty(0)
  ,fCtz(0)
  ,fCte(0)
  ,fCtt(0)
  ,fCcy(0)
  ,fCcz(0)
  ,fCce(0)
  ,fCct(0)
  ,fCcc(0)
  ,fLhElectron(0.0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(0.0)
  ,fBackupTrack(0x0)
{
  //
  // Constructor from AliTPCtrack or AliITStrack .
  //

  Int_t  i = 0;
  Int_t  j = 0;
  UInt_t k = 0;

  SetChi2(0.0);
  SetNumberOfClusters(0);

  for (i = 0; i < kNplane; i++) {
    for (j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = 0.0;
    }
    fTimBinPlane[i] = -1;
  }

  if      (fAlpha < -TMath::Pi()) {
    fAlpha += 2.0 * TMath::Pi();
  }
  else if (fAlpha >= TMath::Pi()) {
    fAlpha -= 2.0 * TMath::Pi();
  }

  Double_t x;
  Double_t p[5]; 
  t.GetExternalParameters(x,p);
  fX = x;
  fY = p[0];
  fZ = p[1];
  fT = p[3]; 
  x  = GetLocalConvConst();
  fC = p[4] / x;
  fE = fC * fX - p[2];   

  // Conversion of the covariance matrix
  Double_t c[15]; 
  t.GetExternalCovariance(c);
  c[10] /= x; 
  c[11] /= x; 
  c[12] /= x; 
  c[13] /= x; 
  c[14] /= x*x;

  Double_t c22 = fX*fX * c[14] - 2.0*fX*c[12] + c[ 5];
  Double_t c32 = fX    * c[13] -        c[ 8];
  Double_t c20 = fX    * c[10] -        c[ 3]; 
  Double_t c21 = fX    * c[11] -        c[ 4];
  Double_t c42 = fX    * c[14] -        c[12];

  fCyy = c[ 0];
  fCzy = c[ 1];  fCzz = c[ 2];
  fCey = c20;    fCez = c21;    fCee = c22;
  fCty = c[ 6];  fCtz = c[ 7];  fCte = c32;  fCtt = c[ 9];
  fCcy = c[10];  fCcz = c[11];  fCce = c42;  fCct = c[13];  fCcc = c[14];  

  for (k = 0; k < kMAXCLUSTERSPERTRACK; k++) {
    fdQdl[k]        = 0;
    fIndex[k]       = 0;
    fIndexBackup[k] = 0;
  }
  
  for (i = 0; i < 3; i++) { 
    fBudget[i]      = 0;
  }

}              

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliESDtrack &t) 
  :AliKalmanTrack() 
  ,fSeedLab(-1)
  ,fdEdx(t.GetTRDsignal())
  ,fdEdxT(0)
  ,fDE(0)
  ,fAlpha(t.GetAlpha())
  ,fX(0)
  ,fStopped(kFALSE)
  ,fY(0)
  ,fZ(0)
  ,fE(0)
  ,fT(0)
  ,fC(0)
  ,fCyy(1e10)
  ,fCzy(0)
  ,fCzz(1e10)
  ,fCey(0)
  ,fCez(0)
  ,fCee(1e10)
  ,fCty(0)
  ,fCtz(0)
  ,fCte(0)
  ,fCtt(1e10)
  ,fCcy(0)
  ,fCcz(0)
  ,fCce(0)
  ,fCct(0)
  ,fCcc(1e10)
  ,fLhElectron(0.0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(0.0)
  ,fBackupTrack(0x0)
{
  //
  // Constructor from AliESDtrack
  //

  Int_t  i = 0;
  Int_t  j = 0;
  UInt_t k = 0;

  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetMass(t.GetMass());
  SetNumberOfClusters(t.GetTRDclusters(fIndex)); 

  Int_t ncl = t.GetTRDclusters(fIndexBackup);
  for (k = ncl; k < kMAXCLUSTERSPERTRACK; k++) {
    fIndexBackup[k] = 0;
    fIndex[k]       = 0;
  }

  for (i = 0; i < kNplane; i++) {
    for (j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = t.GetTRDsignals(i,j);
    }
    fTimBinPlane[i] = t.GetTRDTimBin(i);
  }

  if      (fAlpha <  -TMath::Pi()) {
    fAlpha += 2.0 * TMath::Pi();
  }
  else if (fAlpha >=  TMath::Pi()) {
    fAlpha -= 2.0 * TMath::Pi();
  }

  // Conversion of the covariance matrix
  Double_t x;
  Double_t p[5]; 
  t.GetExternalParameters(x,p);
  Double_t c[15]; 
  t.GetExternalCovariance(c);
  if (t.GetStatus() & AliESDtrack::kTRDbackup) {
    t.GetOuterExternalParameters(fAlpha,x,p);
    t.GetOuterExternalCovariance(c);
    if      (fAlpha <  -TMath::Pi()) {
      fAlpha += 2.0 * TMath::Pi();
    }
    else if (fAlpha >=  TMath::Pi()) {
      fAlpha -= 2.0 * TMath::Pi();
    }
  }

  fX = x;
  fY = p[0];
  fZ = p[1]; 
  SaveLocalConvConst();
  fT = p[3]; 
  x  = GetLocalConvConst();
  fC = p[4] / x;
  fE = fC*fX - p[2];   

  c[10] /= x; 
  c[11] /= x; 
  c[12] /= x; 
  c[13] /= x;
  c[14] /= x*x;

  Double_t c22 = fX*fX * c[14] - 2.0*fX*c[12] + c[ 5];
  Double_t c32 = fX    * c[13] - c[ 8];
  Double_t c20 = fX    * c[10] - c[ 3];
  Double_t c21 = fX    * c[11] - c[ 4];
  Double_t c42 = fX    * c[14] - c[12];

  fCyy = c[ 0];
  fCzy = c[ 1];  fCzz = c[ 2];
  fCey = c20;    fCez = c21;    fCee = c22;
  fCty = c[ 6];  fCtz = c[ 7];  fCte = c32;  fCtt = c[ 9];
  fCcy = c[10];  fCcz = c[11];  fCce = c42;  fCct = c[13];  fCcc = c[14];  

  for (k = 0; k < kMAXCLUSTERSPERTRACK; k++) {
    fdQdl[k] = 0;
    //fIndex[k] = 0; //MI store indexes
  }

  for (i = 0; i < 3; i++) { 
    fBudget[i] = 0;
  }
  if ((t.GetStatus() & AliESDtrack::kTIME) == 0) {
    return;
  }

  StartTimeIntegral();
  Double_t times[10]; 
  t.GetIntegratedTimes(times); 
  SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());

}  

//____________________________________________________________________________
AliTRDtrack::~AliTRDtrack()
{
  //
  // Destructor
  //

  if (fBackupTrack) {
    delete fBackupTrack;
  }
  fBackupTrack = 0;

}

//____________________________________________________________________________
AliTRDtrack &AliTRDtrack::operator=(const AliTRDtrack &t)
{
  //
  // Assignment operator
  //

  fLhElectron    = 0.0;
  fNWrong        = 0;
  fStopped       = 0;
  fNRotate       = 0;
  fNExpected     = 0;
  fNExpectedLast = 0;
  fNdedx         = 0;
  fNCross        = 0;
  fNLast         = 0;
  fChi2Last      = 0;
  fBackupTrack   = 0;

  fAlpha = t.GetAlpha();
  if      (fAlpha <  -TMath::Pi()) {
    fAlpha += 2.0 * TMath::Pi();
  }
  else if (fAlpha >=  TMath::Pi()) {
    fAlpha -= 2.0 * TMath::Pi();
  }

  return *this;

}

//____________________________________________________________________________
Float_t AliTRDtrack::StatusForTOF()
{
  //
  // Defines the status of the TOF extrapolation
  //

  // Definition of res ????
  Float_t res = (0.2 + 0.8 * (fN / (fNExpected + 5.0)))
              * (0.4 + 0.6 * fTracklets[5].GetN() / 20.0);
  res *= (0.25 + 0.8 * 40.0 / (40.0 + fBudget[2]));
  return res;

  // This part of the function is never reached ????
  // What defines these parameters ????
  Int_t status = 0;
  if (GetNumberOfClusters()                <  20)  return 0;
  if ((fN                                  > 110) && 
      (fChi2/(Float_t(fN))                 < 3.0)) return 3; // Gold
  if ((fNLast                              >  30) && 
      (fChi2Last/(Float_t(fNLast))         < 3.0)) return 3; // Gold
  if ((fNLast                              >  20) && 
      (fChi2Last/(Float_t(fNLast))         < 2.0)) return 3; // Gold
  if ((fNLast/(fNExpectedLast+3.0)         > 0.8) && 
      (fChi2Last/Float_t(fNLast)           < 5.0) &&
      (fNLast                              >  20)) return 2; // Silver
  if ((fNLast                              >   5) &&  
      (((fNLast+1.0)/(fNExpectedLast+1.0)) > 0.8) &&
      (fChi2Last/(fNLast-5.0)              < 6.0)) return 1; 
  
  return status;

}
            
//_____________________________________________________________________________
void AliTRDtrack::GetExternalCovariance(Double_t cc[15]) const 
{
  //
  // This function returns external representation of the covriance matrix.
  //

  Double_t a   = GetLocalConvConst();

  Double_t c22 = fX*fX*fCcc  - 2.0*fX*fCce+fCee;
  Double_t c32 = fX*fCct-fCte;
  Double_t c20 = fX*fCcy-fCey;
  Double_t c21 = fX*fCcz-fCez;
  Double_t c42 = fX*fCcc-fCce;

  cc[ 0]=fCyy;
  cc[ 1]=fCzy;   cc[ 2]=fCzz;
  cc[ 3]=c20;    cc[ 4]=c21;    cc[ 5]=c22;
  cc[ 6]=fCty;   cc[ 7]=fCtz;   cc[ 8]=c32;   cc[ 9]=fCtt;
  cc[10]=fCcy*a; cc[11]=fCcz*a; cc[12]=c42*a; cc[13]=fCct*a; cc[14]=fCcc*a*a; 
  
}               
                       
//_____________________________________________________________________________
void AliTRDtrack::GetCovariance(Double_t cc[15]) const 
{
  //
  // Returns the track covariance matrix
  //

  cc[ 0]=fCyy;
  cc[ 1]=fCzy; cc[ 2]=fCzz;
  cc[ 3]=fCey; cc[ 4]=fCez; cc[ 5]=fCee;
  cc[ 6]=fCcy; cc[ 7]=fCcz; cc[ 8]=fCce; cc[ 9]=fCcc;
  cc[10]=fCty; cc[11]=fCtz; cc[12]=fCte; cc[13]=fCct; cc[14]=fCtt;
  
}    

//_____________________________________________________________________________
Int_t AliTRDtrack::Compare(const TObject *o) const 
{
  //
  // Compares tracks according to their Y2 or curvature
  //

  AliTRDtrack *t=(AliTRDtrack*)o;
  //  Double_t co=t->GetSigmaY2();
  //  Double_t c =GetSigmaY2();

  Double_t co=TMath::Abs(t->GetC());
  Double_t c =TMath::Abs(GetC());  

  if      (c>co) return 1;
  else if (c<co) return -1;
  return 0;

}                

//_____________________________________________________________________________
void AliTRDtrack::CookdEdx(Double_t low, Double_t up) 
{
  //
  // Calculates the dE/dX within the "low" and "up" cuts.
  //

  Int_t i = 0;

  // Number of clusters used for dedx
  Int_t nc = fNdedx; 
  // Require at least 10 clusters for a dedx measurement
  // Should fdEdxT not also be set here ????
  if (nc < 10) {
    SetdEdx(0);
    return;
  }

  // Get the charge deposits from each cluster
  Float_t sorted[kMAXCLUSTERSPERTRACK];
  for (i = 0; i < nc; i++) {
    sorted[i] = fdQdl[i];
  }

  // Lower and upper bound
  Int_t nl = Int_t(low * nc);
  Int_t nu = Int_t( up * nc);

  // Sum up the total energy deposit
  Float_t dedx = 0.0;
  for (i = 0; i < nc; i++) {
    dedx += sorted[i];
  }
  // Why is dedx divided by (nu - nl) ????
  // And why is it not saved ????
  if ((nu - nl) != 0) {
    dedx /= (nu - nl);
  }

  // Calculate the truncated mean
  // Can fdQdl be negative ???? Why is it nor checked above ????
  for (i = 0; i < nc; i++) {
    sorted[i] = TMath::Abs(fdQdl[i]);
  }

  // Sort the dedx values by amplitude
  Int_t *index = new Int_t[nc];
  TMath::Sort(nc,sorted,index,kFALSE);
  // Sum up the truncated charge between nl and nu  
  dedx = 0.0;
  for (i = nl; i <= nu; i++) {
    dedx += sorted[index[i]];
  }
  dedx /= (nu - nl + 1);

  // Why are fdEdx and fdEdxT set to the same value ????
  fdEdxT = dedx;
  SetdEdx(dedx);

  delete [] index;

}                     

//_____________________________________________________________________________
Int_t AliTRDtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho)
{
  //
  // Propagates a track of particle with mass=pm to a reference plane 
  // defined by x=xk through media of density=rho and radiationLength=x0
  //

  if (xk == fX) {
    return 1;
  }

  if (TMath::Abs(fC*xk - fE) >= 0.9) {
    return 0;
  }

  Double_t lcc = GetLocalConvConst();

  Double_t oldX = fX;
  Double_t oldY = fY;
  Double_t oldZ = fZ;  

  Double_t x1 = fX;
  Double_t x2 = x1 + (xk - x1);
  Double_t dx = x2 - x1;
  Double_t y1 = fY;
  Double_t z1 = fZ;
  Double_t c1 = fC*x1 - fE;
  if ((c1*c1) > 1) {
    return 0;
  }

  Double_t r1 = TMath::Sqrt(1.0 - c1*c1);
  Double_t c2 = fC*x2 - fE; 
  if ((c2*c2) > 1) {
    return 0;
  }

  Double_t r2 = TMath::Sqrt(1.0 - c2*c2);

  fY += dx*(c1+c2) / (r1+r2);
  fZ += dx*(c1+c2) / (c1*r2 + c2*r1) * fT;

  // f = F - 1
  Double_t rr  =  r1+r2;
  Double_t cc  =  c1+c2;
  Double_t xx  =  x1+x2;
  Double_t f02 = -dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  Double_t f04 =  dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  Double_t cr  =  c1*r2+c2*r1;
  Double_t f12 = -dx*fT*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  Double_t f13 =  dx*cc/cr;
  Double_t f14 =  dx*fT*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);

  // b = C*ft
  Double_t b00 = f02*fCey + f04*fCcy;
  Double_t b01 = f12*fCey + f14*fCcy + f13*fCty;
  Double_t b10 = f02*fCez + f04*fCcz;
  Double_t b11 = f12*fCez + f14*fCcz + f13*fCtz;
  Double_t b20 = f02*fCee + f04*fCce;
  Double_t b21 = f12*fCee + f14*fCce + f13*fCte;
  Double_t b30 = f02*fCte + f04*fCct;
  Double_t b31 = f12*fCte + f14*fCct + f13*fCtt;
  Double_t b40 = f02*fCce + f04*fCcc;
  Double_t b41 = f12*fCce + f14*fCcc + f13*fCct;

  // a = f*b = f*C*ft
  Double_t a00 = f02*b20 + f04*b40;
  Double_t a01 = f02*b21 + f04*b41;
  Double_t a11 = f12*b21 + f14*b41 + f13*b31;

  // F*C*Ft = C + (a + b + bt)
  fCyy += a00 + 2.0*b00;
  fCzy += a01 + b01 + b10;
  fCey += b20;
  fCty += b30;
  fCcy += b40;
  fCzz += a11 + 2.0*b11;
  fCez += b21;
  fCtz += b31;
  fCcz += b41;

  fX = x2;                                                     

  // Change of the magnetic field
  SaveLocalConvConst();
  cc =  fC;
  fC *= lcc / GetLocalConvConst();
  fE += fX  * (fC-cc);

  // Multiple scattering
  // What is 14.1 ????
  Double_t d      = TMath::Sqrt((x1-fX)*(x1-fX) + (y1-fY)*(y1-fY) + (z1-fZ)*(z1-fZ));
  Double_t p2     = (1.0 + GetTgl()*GetTgl()) / (Get1Pt()*Get1Pt());
  Double_t beta2  = p2 / (p2 + GetMass()*GetMass());
  Double_t theta2 = 14.1*14.1 / (beta2*p2*1e6) * d / x0 * rho;
  Double_t ey     = fC*fX - fE;
  Double_t ez     = fT;
  Double_t xz     = fC*ez;
  Double_t zz1    = ez*ez + 1.0;
  Double_t xy     = fE + ey;
  
  fCee += (2.0*ey*ez*ez*fE + 1.0 - ey*ey + ez*ez + fE*fE*ez*ez) * theta2;
  fCte += ez*zz1*xy*theta2;
  fCtt += zz1*zz1*theta2;
  fCce += xz*ez*xy*theta2;
  fCct += xz*zz1*theta2;
  fCcc += xz*xz*theta2;

  // Energy losses
  // What is 5940.0 ???? and 0.153e-3 ????
  if ((5940.0*beta2 / (1.0 - beta2 + 1e-10) - beta2) < 0.0) {
    return 0;
  }
  Double_t dE     = 0.153e-3/beta2 * (TMath::Log(5940.0*beta2 / (1.0 - beta2 + 1e-10)) - beta2) 
                  * d * rho;
  Float_t  budget = d * rho;
  fBudget[0] +=budget;
 
  // Suspicious part - think about it ????
  Double_t kinE =  TMath::Sqrt(p2);
  if (dE > 0.8 * kinE) {
    dE = 0.8 * kinE;
  }
  if (dE <        0.0) {
    dE = 0.0;       // Not valid region for Bethe bloch 
  }
  fDE += dE;
  if (x1 < x2) {
    dE = -dE;
  }
  cc  = fC;
  fC *= (1.0 - TMath::Sqrt(p2 + GetMass()*GetMass()) / p2 * dE);
  fE += fX * (fC - cc);    

  // Energy loss fluctuation 
  // Why 0.07 ????
  Double_t sigmade = 0.07 * TMath::Sqrt(TMath::Abs(dE));  
  Double_t sigmac2 = sigmade*sigmade * fC*fC * (p2 + GetMass()*GetMass()) / (p2*p2);
  fCcc += sigmac2;
  fCee += fX*fX * sigmac2;  

  // Track time measurement
  if (x1 < x2) {
    if (IsStartedTimeIntegral()) {
      Double_t l2 = TMath::Sqrt((fX-oldX)*(fX-oldX) 
                              + (fY-oldY)*(fY-oldY) 
                              + (fZ-oldZ)*(fZ-oldZ));
      if (TMath::Abs(l2*fC) > 0.0001){
        // <ake correction for curvature if neccesary
        l2 = 0.5 * TMath::Sqrt((fX-oldX)*(fX-oldX) 
                             + (fY-oldY)*(fY-oldY));
        l2 = 2.0 * TMath::ASin(l2 * fC) / fC;
        l2 = TMath::Sqrt(l2*l2 + (fZ-oldZ)*(fZ-oldZ));
      }
      AddTimeStep(l2);
    }
  }

  return 1;            

}     

//_____________________________________________________________________________
Int_t AliTRDtrack::Update(const AliTRDcluster *c, Double_t chisq
                        , UInt_t index, Double_t h01)
{
  //
  // Assignes a found cluster to the track and updates track information
  //

  Bool_t fNoTilt = kTRUE;
  // What is 0.003 ????
  if (TMath::Abs(h01) > 0.003) {
    fNoTilt = kFALSE;
  }

  // Add angular effect to the error contribution
  Float_t tangent2  = (fC*fX-fE) * (fC*fX-fE);
  if (tangent2 < 0.90000){
    tangent2 = tangent2 / (1.0 - tangent2);
  }
  // What is 0.04 ????
  Float_t errang    = tangent2 * 0.04;
  Float_t padlength = TMath::Sqrt(c->GetSigmaZ2() * 12.0);

  Double_t r00 = c->GetSigmaY2() + errang;
  Double_t r01 = 0.0;
  Double_t r11 = c->GetSigmaZ2() * 100.0;
  r00 += fCyy; 
  r01 += fCzy; 
  r11 += fCzz;
  Double_t det = r00*r11 - r01*r01;
  Double_t tmp = r00; 
  r00 =  r11 / det; 
  r11 =  tmp / det; 
  r01 = -r01 / det;

  Double_t k00 = fCyy*r00 + fCzy*r01;
  Double_t k01 = fCyy*r01 + fCzy*r11;
  Double_t k10 = fCzy*r00 + fCzz*r01;
  Double_t k11 = fCzy*r01 + fCzz*r11;
  Double_t k20 = fCey*r00 + fCez*r01; 
  Double_t k21 = fCey*r01 + fCez*r11;
  Double_t k30 = fCty*r00 + fCtz*r01;
  Double_t k31 = fCty*r01 + fCtz*r11;
  Double_t k40 = fCcy*r00 + fCcz*r01; 
  Double_t k41 = fCcy*r01 + fCcz*r11;

  Double_t dy  = c->GetY() - fY;
  Double_t dz  = c->GetZ() - fZ;
  Double_t cur = fC + k40*dy + k41*dz;
  Double_t eta = fE + k20*dy + k21*dz;

  if (fNoTilt) {

    if (TMath::Abs(cur*fX-eta) >= 0.9) {
      return 0;
    }
    fY += k00*dy + k01*dz;
    fZ += k10*dy + k11*dz;
    fE  = eta;
    fC  = cur;

  }
  else {

    // Empirical factor set by C.Xu in the first tilt version      
    // Is this factor still ok ???? 
    Double_t xuFactor = 100.0;  
    dy   = c->GetY() - fY; 
    dz   = c->GetZ() - fZ;     
    dy   = dy + h01*dz;

    Float_t add = 0.0;
    if (TMath::Abs(dz) > padlength/2.0) {
      Float_t dy2  = c->GetY() - fY;
      Float_t sign = (dz > 0.0) ? -1.0 : 1.0;
      dy2 += h01 * sign * padlength/2.0;	
      dy   = dy2;
      add  = 0.0;
    }
  
    r00  = c->GetSigmaY2() + errang + add;
    r01  = 0.0;
    r11  = c->GetSigmaZ2() * xuFactor; 
    r00 += (fCyy + 2.0*h01*fCzy + h01*h01*fCzz);
    r01 += (fCzy + h01*fCzz);
    r11 += fCzz;

    det  =  r00*r11 - r01*r01;
    tmp  =  r00; 
    r00  =  r11/det; 
    r11  =  tmp/det; 
    r01  = -r01/det;

    k00  = fCyy*r00 + fCzy * (r01 + h01*r00);
    k01  = fCyy*r01 + fCzy * (r11 + h01*r01);
    k10  = fCzy*r00 + fCzz * (r01 + h01*r00);
    k11  = fCzy*r01 + fCzz * (r11 + h01*r01);
    k20  = fCey*r00 + fCez * (r01 + h01*r00);
    k21  = fCey*r01 + fCez * (r11 + h01*r01);
    k30  = fCty*r00 + fCtz * (r01 + h01*r00);
    k31  = fCty*r01 + fCtz * (r11 + h01*r01);
    k40  = fCcy*r00 + fCcz * (r01 + h01*r00);
    k41  = fCcy*r01 + fCcz * (r11 + h01*r01);  

    cur  = fC + k40*dy + k41*dz; 
    eta  = fE + k20*dy + k21*dz;
    if (TMath::Abs(cur*fX - eta) >= 0.9) {
      return 0;
    }                           

    fY  += k00*dy + k01*dz;
    fZ  += k10*dy + k11*dz;
    fE   = eta;
    fT  += k30*dy + k31*dz;
    fC   = cur;
    
    k01 += h01*k00;
    k11 += h01*k10;
    k21 += h01*k20;
    k31 += h01*k30;
    k41 += h01*k40;  
    
  }

  Double_t c01 = fCzy;
  Double_t c02 = fCey;
  Double_t c03 = fCty;
  Double_t c04 = fCcy;
  Double_t c12 = fCez;
  Double_t c13 = fCtz;
  Double_t c14 = fCcz;

  fCyy -= k00*fCyy + k01*fCzy; 
  fCzy -= k00*c01  + k01*fCzz;
  fCey -= k00*c02  + k01*c12; 
  fCty -= k00*c03  + k01*c13;
  fCcy -= k00*c04  + k01*c14;
  
  fCzz -= k10*c01  + k11*fCzz;
  fCez -= k10*c02  + k11*c12;
  fCtz -= k10*c03  + k11*c13;
  fCcz -= k10*c04  + k11*c14;
  
  fCee -= k20*c02  + k21*c12;
  fCte -= k20*c03  + k21*c13;
  fCce -= k20*c04  + k21*c14;
  
  fCtt -= k30*c03  + k31*c13;
  fCct -= k40*c03  + k41*c13;  
  
  fCcc -= k40*c04  + k41*c14;                 

  Int_t n = GetNumberOfClusters();
  fIndex[n] = index;
  SetNumberOfClusters(n+1);

  SetChi2(GetChi2()+chisq);

  return 1;     

}                     

//_____________________________________________________________________________
Int_t AliTRDtrack::UpdateMI(const AliTRDcluster *c, Double_t chisq
                          , UInt_t index, Double_t h01 
			  , Int_t /*plane*/)
{
  //
  // Assignes found cluster to the track and updates track information
  //

  Bool_t fNoTilt = kTRUE;
  if (TMath::Abs(h01) > 0.003) {
    fNoTilt = kFALSE;
  }

  //
  // Add angular effect to the error contribution and make correction
  // Still needed ???? 
  // AliTRDclusterCorrection *corrector = AliTRDclusterCorrection::GetCorrection();
  // 

  Double_t tangent2 = (fC*fX-fE) * (fC*fX-fE);
  if (tangent2 < 0.9) {
    tangent2 = tangent2 / (1.0 - tangent2);
  }
  Double_t tangent  = TMath::Sqrt(tangent2);
  if ((fC*fX-fE) < 0.0) {
    tangent *= -1.0;
  }

  // Where are the parameters from ????
  Double_t errang = tangent2 * 0.04;
  Double_t errsys = 0.025*0.025 * 20.0;  // Systematic error part 

  Float_t  extend = 1.0;
  if (c->GetNPads() == 4) extend = 2.0;

  /////////////////////////////////////////////////////////////////////////////
  //
  // Is this still needed or will it be needed ????
  //
  //if (c->GetNPads() == 5) extend = 3.0;
  //if (c->GetNPads() == 6) extend = 3.0;
  //if (c->GetQ() < 15) {
  //  return 1;
  //}
  //
  // Will this be needed ????
  /*
  if (corrector !=0 ) {
  //if (0){
    correction = corrector->GetCorrection(plane,c->GetLocalTimeBin(),tangent);
    if (TMath::Abs(correction)>0){
      //if we have info 
      errang     = corrector->GetSigma(plane,c->GetLocalTimeBin(),tangent);
      errang    *= errang;      
      errang    += tangent2*0.04;
    }
  }
  */
  //  
  //  Double_t padlength = TMath::Sqrt(c->GetSigmaZ2()*12.);
  /////////////////////////////////////////////////////////////////////////////

  Double_t r00 = (c->GetSigmaY2() + errang + errsys) * extend;
  Double_t r01 = 0.0;
  Double_t r11 = c->GetSigmaZ2()*10000.0;
  r00 += fCyy; 
  r01 += fCzy; 
  r11 += fCzz;
  Double_t det =r00*r11 - r01*r01;
  Double_t tmp =r00;
  r00  =  r11 / det; 
  r11  =  tmp / det; 
  r01  = -r01 / det;

  Double_t k00 = fCyy*r00 + fCzy*r01;
  Double_t k01 = fCyy*r01 + fCzy*r11;
  Double_t k10 = fCzy*r00 + fCzz*r01;
  Double_t k11 = fCzy*r01 + fCzz*r11;
  Double_t k20 = fCey*r00 + fCez*r01;
  Double_t k21 = fCey*r01 + fCez*r11;
  Double_t k30 = fCty*r00 + fCtz*r01;
  Double_t k31 = fCty*r01 + fCtz*r11;
  Double_t k40 = fCcy*r00 + fCcz*r01;
  Double_t k41 = fCcy*r01 + fCcz*r11;

  Double_t dy  = c->GetY() - fY;
  Double_t dz  = c->GetZ() - fZ;
  Double_t cur = fC + k40*dy + k41*dz;
  Double_t eta = fE + k20*dy + k21*dz;

  if (fNoTilt) {

    if (TMath::Abs(cur*fX - eta) >= 0.9) {
      return 0;
    }

    fY += k00*dy + k01*dz;
    fZ += k10*dy + k11*dz;
    fE  = eta;
    fC  = cur;

  }
  else {

    Double_t padlength = TMath::Sqrt(c->GetSigmaZ2() * 12.0);
    // Empirical factor set by C.Xu in the first tilt version 
    Double_t xuFactor  = 1000.0;  

    dy = c->GetY() - fY; 
    dz = c->GetZ() - fZ;     
    //dy = dy + h01*dz + correction; // Still needed ????
    
    Double_t tiltdz = dz;
    if (TMath::Abs(tiltdz) > padlength/2.0) {
      tiltdz = TMath::Sign(padlength/2.0,dz);
    }
    dy = dy + h01*tiltdz;

    Double_t add = 0.0;
    if (TMath::Abs(dz) > padlength/2.0) {
      //Double_t dy2 = c->GetY() - fY;     // Still needed ????
      //Double_t sign = (dz>0) ? -1.: 1.;
      //dy2-=h01*sign*padlength/2.;	
      //dy = dy2;
      add = 1.0;
    }

    Double_t s00 = (c->GetSigmaY2() + errang) * extend + errsys + add;  // Error pad
    Double_t s11 = c->GetSigmaZ2() * xuFactor;                          // Error pad-row
    //
    r00  = fCyy + 2*fCzy*h01 + fCzz*h01*h01 + s00;
    r01  = fCzy + fCzz*h01;
    r11  = fCzz + s11;
    det  = r00*r11 - r01*r01;

    // Inverse matrix
    tmp  =  r00; 
    r00  =  r11 / det; 
    r11  =  tmp / det; 
    r01  = -r01 / det;

    // K matrix
    k00  = fCyy*r00 + fCzy * (r01 + h01*r00);
    k01  = fCyy*r01 + fCzy * (r11 + h01*r01);
    k10  = fCzy*r00 + fCzz * (r01 + h01*r00);
    k11  = fCzy*r01 + fCzz * (r11 + h01*r01);
    k20  = fCey*r00 + fCez * (r01 + h01*r00);
    k21  = fCey*r01 + fCez * (r11 + h01*r01);
    k30  = fCty*r00 + fCtz * (r01 + h01*r00);
    k31  = fCty*r01 + fCtz * (r11 + h01*r01);
    k40  = fCcy*r00 + fCcz * (r01 + h01*r00);
    k41  = fCcy*r01 + fCcz * (r11 + h01*r01);  
    
    // Update measurement
    cur  = fC + k40*dy + k41*dz; 
    eta  = fE + k20*dy + k21*dz;
    if (TMath::Abs(cur*fX - eta) >= 0.9) {
      return 0;
    }                           
    fY  += k00*dy + k01*dz;
    fZ  += k10*dy + k11*dz;
    fE   = eta;
    fT  += k30*dy + k31*dz;
    fC   = cur;
    
    k01 += h01*k00;
    k11 += h01*k10;
    k21 += h01*k20;
    k31 += h01*k30;
    k41 += h01*k40;  
    
  }

  // Update the covariance matrix
  Double_t oldyy = fCyy;
  Double_t oldzz = fCzz; 
  //Double_t oldee = fCee;
  //Double_t oldcc = fCcc;
  Double_t oldzy = fCzy;
  Double_t oldey = fCey;
  Double_t oldty = fCty;
  Double_t oldcy = fCcy;
  Double_t oldez = fCez;
  Double_t oldtz = fCtz;
  Double_t oldcz = fCcz;
  //Double_t oldte = fCte;
  //Double_t oldce = fCce;
  //Double_t oldct = fCct;

  fCyy -= k00*oldyy + k01*oldzy;   
  fCzy -= k10*oldyy + k11*oldzy;
  fCey -= k20*oldyy + k21*oldzy;   
  fCty -= k30*oldyy + k31*oldzy;
  fCcy -= k40*oldyy + k41*oldzy;  
  
  fCzz -= k10*oldzy + k11*oldzz;
  fCez -= k20*oldzy + k21*oldzz;   
  fCtz -= k30*oldzy + k31*oldzz;
  fCcz -= k40*oldzy + k41*oldzz;
  
  fCee -= k20*oldey + k21*oldez;   
  fCte -= k30*oldey + k31*oldez;
  fCce -= k40*oldey + k41*oldez;
  
  fCtt -= k30*oldty + k31*oldtz;
  fCct -= k40*oldty + k41*oldtz;
  
  fCcc -= k40*oldcy + k41*oldcz;                 

  Int_t n = GetNumberOfClusters();
  fIndex[n] = index;
  SetNumberOfClusters(n+1);

  SetChi2(GetChi2() + chisq);

  return 1;      

}                     

//_____________________________________________________________________________
Int_t AliTRDtrack::UpdateMI(const AliTRDtracklet &tracklet)
{
  //
  // Assignes found tracklet to the track and updates track information
  //
  
  Double_t r00 = (tracklet.GetTrackletSigma2());
  Double_t r01 = 0.0;
  Double_t r11 = 10000.0;
  r00 += fCyy; 
  r01 += fCzy; 
  r11 += fCzz;
  
  Double_t det = r00*r11 - r01*r01;
  Double_t tmp = r00; 
  r00 =  r11 / det; 
  r11 =  tmp / det; 
  r01 = -r01 / det;

  Double_t dy  = tracklet.GetY() - fY;
  Double_t dz  = tracklet.GetZ() - fZ;

  Double_t s00 = tracklet.GetTrackletSigma2();  // Error pad
  Double_t s11 = 100000.0;                      // Error pad-row
  Float_t  h01 = tracklet.GetTilt();
  
  r00 = fCyy + fCzz*h01*h01 + s00;
  r01 = fCzy;
  r11 = fCzz + s11;
  det = r00*r11 - r01*r01;

  // Inverse matrix
  tmp =  r00; 
  r00 =  r11 / det; 
  r11 =  tmp / det; 
  r01 = -r01 / det;

  // K matrix
  Double_t k00 = fCyy*r00 + fCzy*r01;
  Double_t k01 = fCyy*r01 + fCzy*r11;
  Double_t k10 = fCzy*r00 + fCzz*r01;
  Double_t k11 = fCzy*r01 + fCzz*r11;
  Double_t k20 = fCey*r00 + fCez*r01;
  Double_t k21 = fCey*r01 + fCez*r11;
  Double_t k30 = fCty*r00 + fCtz*r01;
  Double_t k31 = fCty*r01 + fCtz*r11;
  Double_t k40 = fCcy*r00 + fCcz*r01;
  Double_t k41 = fCcy*r01 + fCcz*r11;
  
  // Update measurement
  Double_t cur = fC + k40*dy + k41*dz;
  Double_t eta = fE + k20*dy + k21*dz;  
  if (TMath::Abs(cur*fX-eta) >= 0.9) {
    return 0;
  }                           

  fY += k00*dy + k01*dz;
  fZ += k10*dy + k11*dz;
  fE  = eta;
  fT += k30*dy + k31*dz;
  fC  = cur;
    
  // Update the covariance matrix
  Double_t oldyy = fCyy;
  Double_t oldzz = fCzz; 
  //Double_t oldee = fCee;
  //Double_t oldcc = fCcc;
  Double_t oldzy = fCzy;
  Double_t oldey = fCey;
  Double_t oldty = fCty;
  Double_t oldcy = fCcy;
  Double_t oldez = fCez;
  Double_t oldtz = fCtz;
  Double_t oldcz = fCcz;
  //Double_t oldte = fCte;
  //Double_t oldce = fCce;
  //Double_t oldct = fCct;

  fCyy -= k00*oldyy + k01*oldzy;   
  fCzy -= k10*oldyy + k11*oldzy;
  fCey -= k20*oldyy + k21*oldzy;   
  fCty -= k30*oldyy + k31*oldzy;
  fCcy -= k40*oldyy + k41*oldzy;  
  
  fCzz -= k10*oldzy + k11*oldzz;
  fCez -= k20*oldzy + k21*oldzz;   
  fCtz -= k30*oldzy + k31*oldzz;
  fCcz -= k40*oldzy + k41*oldzz;
  
  fCee -= k20*oldey + k21*oldez;   
  fCte -= k30*oldey + k31*oldez;
  fCce -= k40*oldey + k41*oldez;
  
  fCtt -= k30*oldty + k31*oldtz;
  fCct -= k40*oldty + k41*oldtz;
  
  fCcc -= k40*oldcy + k41*oldcz;                 

  return 1;      

}                     

//_____________________________________________________________________________
Int_t AliTRDtrack::Rotate(Double_t alpha, Bool_t absolute)
{
  //
  // Rotates the track parameters in the R*phi plane,
  // if the absolute rotation alpha is in global system.
  // Otherwise the rotation is relative to the current rotation angle
  //  

  if (absolute) {
    alpha -= fAlpha;
  }
  else{
    fNRotate++;
  }

  fAlpha += alpha;
  if (fAlpha <  -TMath::Pi()) {
    fAlpha += 2.0 * TMath::Pi();
  }
  if (fAlpha >=  TMath::Pi()) {
    fAlpha -= 2.0 * TMath::Pi();
  }

  Double_t x1 = fX;
  Double_t y1 = fY;
  Double_t ca = TMath::Cos(alpha);
  Double_t sa = TMath::Sin(alpha);
  Double_t r1 = fC*fX - fE;

  fX =  x1*ca + y1*sa;
  fY = -x1*sa + y1*ca;
  if ((r1*r1) > 1.0) {
    return 0;
  }
  fE = fE*ca + (fC*y1 + TMath::Sqrt(1.0 - r1*r1)) * sa;

  Double_t r2 = fC*fX - fE;
  if (TMath::Abs(r2) >= 0.9) {
    Int_t n = GetNumberOfClusters();
    if (n > 4) {
      AliError(Form("Rotation failed N = %d !\n",n));
    }
    return 0;
  }

  if ((r2*r2) > 1.0) {
    return 0;
  }
  Double_t y0 = fY + TMath::Sqrt(1.0 - r2*r2) / fC;
  if ((fY-y0)*fC >= 0.0) {
    Int_t n = GetNumberOfClusters();
    if (n > 4) {
      AliError(Form("Rotation failed N = %d !\n",n));
    }
    return 0;
  }

  // f = F - 1
  Double_t f00 = ca-1.0;
  Double_t f24 = (y1 - r1*x1/TMath::Sqrt(1.0 - r1*r1)) * sa;
  Double_t f20 = fC*sa; 
  Double_t f22 = (ca + sa*r1/TMath::Sqrt(1.0 - r1*r1)) - 1.0;

  // b = C*ft
  Double_t b00 = fCyy*f00;
  Double_t b02 = fCyy*f20 + fCcy*f24 + fCey*f22;
  Double_t b10 = fCzy*f00;
  Double_t b12 = fCzy*f20 + fCcz*f24 + fCez*f22;
  Double_t b20 = fCey*f00;
  Double_t b22 = fCey*f20 + fCce*f24 + fCee*f22;
  Double_t b30 = fCty*f00;
  Double_t b32 = fCty*f20 + fCct*f24 + fCte*f22;
  Double_t b40 = fCcy*f00;
  Double_t b42 = fCcy*f20 + fCcc*f24 + fCce*f22;

  // a = f*b = f*C*ft
  Double_t a00 = f00*b00;
  Double_t a02 = f00*b02;
  Double_t a22 = f20*b02  + f24*b42  + f22*b22;

  // F*C*Ft = C + (a + b + bt)
  fCyy += a00 + 2.0*b00;
  fCzy += b10;
  fCey += a02 +     b20 + b02;
  fCty += b30;
  fCcy += b40;
  fCez += b12;
  fCte += b32;
  fCee += a22 + 2.0*b22;
  fCce += b42;

  return 1;                            

}                         

//_____________________________________________________________________________
Double_t AliTRDtrack::GetPredictedChi2(const AliTRDcluster *c, Double_t h01) const
{
  //
  // Returns the track chi2
  //  

  Bool_t fNoTilt = kTRUE;
  if (TMath::Abs(h01) > 0.003) {
    fNoTilt = kFALSE;
  }

  Double_t chi2;
  Double_t dy;
  Double_t r00;
  Double_t r01;
  Double_t r11;

  if (fNoTilt) {

    dy   = c->GetY() - fY;
    r00  = c->GetSigmaY2();    
    chi2 = (dy*dy) / r00;    

  }
  else {

    Double_t padlength = TMath::Sqrt(c->GetSigmaZ2() * 12.0);

    r00  = c->GetSigmaY2(); 
    r01  = 0.0; 
    r11  = c->GetSigmaZ2();
    r00 += fCyy; 
    r01 += fCzy; 
    r11 += fCzz;

    Double_t det = r00*r11 - r01*r01;
    if (TMath::Abs(det) < 1.e-10) {
      Int_t n = GetNumberOfClusters(); 
      if (n > 4) {
        AliError(Form("Singular matrix N = %d!\n",n));
      }
      return 1e10;
    }

    Double_t tmp = r00; 
    r00  =  r11; 
    r11  =  tmp; 
    r01  = -r01;
    Double_t dy     = c->GetY() - fY;
    Double_t dz     = c->GetZ() - fZ;
    Double_t tiltdz = dz;
    if (TMath::Abs(tiltdz) > padlength/2.0) {
      tiltdz = TMath::Sign(padlength/2.0,dz);
    }
    dy = dy + h01*tiltdz;

    chi2 = (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz) / det; 

  }

  return chi2;

}      

//_________________________________________________________________________
void AliTRDtrack::GetPxPyPz(Double_t& px, Double_t& py, Double_t& pz) const
{
  //
  // Returns reconstructed track momentum in the global system.
  //

  Double_t pt = TMath::Abs(GetPt());
  Double_t r  = fC*fX - fE;

  Double_t y0; 
  if      (r >  1) { 
    py =  pt; 
    px = 0.0;
  }
  else if (r < -1) { 
    py = -pt; 
    px = 0.0;
  }
  else {
    y0 =  fY + TMath::Sqrt(1.0 - r*r) / fC;  
    px = -pt * (fY - y0) * fC; //cos(phi);
    py = -pt * (fE - fX*fC);   //sin(phi);
  }

  pz = pt*fT;
  Double_t tmp = px * TMath::Cos(fAlpha) 
               - py * TMath::Sin(fAlpha);
  py = px * TMath::Sin(fAlpha) 
     + py * TMath::Cos(fAlpha);
  px = tmp;            

}                                

//_________________________________________________________________________
void AliTRDtrack::GetGlobalXYZ(Double_t& x, Double_t& y, Double_t& z) const
{
  //
  // Returns reconstructed track coordinates in the global system.
  //

  x = fX; 
  y = fY; 
  z = fZ
; 
  Double_t tmp = x * TMath::Cos(fAlpha)
               - y * TMath::Sin(fAlpha);
  y = x * TMath::Sin(fAlpha) 
    + y * TMath::Cos(fAlpha);
  x = tmp;            

}                                

//_________________________________________________________________________
void AliTRDtrack::ResetCovariance() 
{
  //
  // Resets covariance matrix
  //

  fCyy *= 10.0;
  fCzy  =  0.0;  fCzz *= 10.0;
  fCey  =  0.0;  fCez  =  0.0;  fCee *= 10.0;
  fCty  =  0.0;  fCtz  =  0.0;  fCte  =  0.0;  fCtt *= 10.0;
  fCcy  =  0.0;  fCcz  =  0.0;  fCce  =  0.0;  fCct  =  0.0;  fCcc *= 10.0;  

}                                                         

//_____________________________________________________________________________
void AliTRDtrack::ResetCovariance(Float_t mult) 
{
  //
  // Resets covariance matrix
  //

  fCyy *= mult;
  fCzy *= 0.0;   fCzz *= 1.0;
  fCey *= 0.0;   fCez *= 0.0;   fCee *= mult;
  fCty *= 0.0;   fCtz *= 0.0;   fCte *= 0.0;   fCtt *= 1.0;
  fCcy *= 0.0;   fCcz *= 0.0;   fCce *= 0.0;   fCct *= 0.0;   fCcc *= mult;  

}                                                         

//_____________________________________________________________________________
void AliTRDtrack::MakeBackupTrack()
{
  //
  // Creates a backup track
  //

  if (fBackupTrack) {
    delete fBackupTrack;
  }

  fBackupTrack = new AliTRDtrack(*this);
  
}

//_____________________________________________________________________________
Int_t AliTRDtrack::GetProlongation(Double_t xk, Double_t &y, Double_t &z)
{
  //
  // Find prolongation at given x
  // Return 0 if it does not exist
  //
  
  Double_t c1 = fC*fX - fE;
  if (TMath::Abs(c1) > 1.0) {
    return 0;
  }

  Double_t r1 = TMath::Sqrt(1.0 - c1*c1);
  Double_t c2 = fC*xk - fE;
  if (TMath::Abs(c2) > 1.0) {
    return 0;  
  }

  Double_t r2 = TMath::Sqrt(1.0 - c2*c2);
  y = fY + (xk-fX)*(c1+c2)/(r1+r2);
  z = fZ + (xk-fX)*(c1+c2)/(c1*r2 + c2*r1)*fT;

  return 1;
  
}

//_____________________________________________________________________________
Int_t AliTRDtrack::PropagateToX(Double_t xr, Double_t step)
{
  //
  // Propagate track to a given x position 
  // Works inside of the 20 degree segmentation
  // (local cooordinate frame for TRD , TPC, TOF)
  // 
  // The material budget is taken from the geo manager
  //
 
  Double_t  xyz0[3];
  Double_t  xyz1[3];
  Double_t  y;
  Double_t  z;

  // Critical alpha  - cross sector indication
  const Double_t kAlphac  = TMath::Pi()/9.0;   
  const Double_t kTalphac = TMath::Tan(kAlphac*0.5);

  // Direction +-
  Double_t dir = (fX > xr) ? -1.0 : 1.0;

  for (Double_t x = fX + dir*step; dir*x < dir*xr; x += dir*step) {
    
    GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);	
    GetProlongation(x,y,z);
    xyz1[0] = x * TMath::Cos(fAlpha) + y * TMath::Sin(fAlpha); 
    xyz1[1] = x * TMath::Sin(fAlpha) - y * TMath::Cos(fAlpha);
    xyz1[2] = z;
    Double_t param[7];
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
    
    if ((param[0] > 0) && 
        (param[1] > 0)) {
      PropagateTo(x,param[1],param[0]);
    }
    if (fY >  fX*kTalphac) {
      Rotate(-kAlphac);
    }
    if (fY < -fX*kTalphac) {
      Rotate(kAlphac);
    }

  }

  PropagateTo(xr);

  return 0;

}

//_____________________________________________________________________________
Int_t AliTRDtrack::PropagateToR(Double_t r,Double_t step)
{
  //
  // Propagate a track to a given radial position
  // The rotation is always connected to the last track position
  //

  Double_t xyz0[3];
  Double_t xyz1[3];
  Double_t y;
  Double_t z; 

  // Direction +-
  Double_t radius = TMath::Sqrt(fX*fX + fY*fY);
  Double_t dir    = (radius > r) ? -1.0 : 1.0;   
  
  for (Double_t x = radius + dir*step; dir*x < dir*r; x += dir*step) {
    GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);	
    Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
    Rotate(alpha,kTRUE);
    GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);	
    GetProlongation(x,y,z);
    xyz1[0] = x * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
    xyz1[1] = x * TMath::Sin(alpha) - y * TMath::Cos(alpha);
    xyz1[2] = z;
    Double_t param[7];
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
    if (param[1] <= 0.0) {
      param[1] = 100000000.0;
    }
    PropagateTo(x,param[1],param[0]);
  } 

  GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);	
  Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
  Rotate(alpha,kTRUE);
  GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);	
  GetProlongation(r,y,z);
  xyz1[0] = r * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
  xyz1[1] = r * TMath::Sin(alpha) - y * TMath::Cos(alpha);
  xyz1[2] = z;
  Double_t param[7];
  AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
  
  if (param[1] <= 0.0) {
    param[1] = 100000000.0;
  }
  PropagateTo(r,param[1],param[0]);

  return 0;

}

//_____________________________________________________________________________
Int_t AliTRDtrack::GetSector() const
{
  //
  // Return the current sector
  //

  return Int_t(TVector2::Phi_0_2pi(fAlpha)
             / AliTRDgeometry::GetAlpha())
             % AliTRDgeometry::kNsect;

}

//_____________________________________________________________________________
Double_t AliTRDtrack::Get1Pt() const                       
{ 
  //
  // Returns the inverse Pt (1/GeV/c)
  // (or 1/"most probable pt", if the field is too weak)
  //

  if (TMath::Abs(GetLocalConvConst()) > kVeryBigConvConst) {
    return 1.0 / kMostProbableMomentum 
               / TMath::Sqrt(1.0 + GetTgl()*GetTgl());
  }

  return (TMath::Sign(1.0e-9,fC) + fC) * GetLocalConvConst();

}

//_____________________________________________________________________________
Double_t AliTRDtrack::GetP() const                         
{ 
  //
  // Returns the total momentum
  //

  return TMath::Abs(GetPt()) * TMath::Sqrt(1.0 + GetTgl()*GetTgl());  

}

//_____________________________________________________________________________
Double_t AliTRDtrack::GetYat(Double_t xk) const            
{     
  //
  // This function calculates the Y-coordinate of a track at 
  // the plane x = xk.
  // Needed for matching with the TOF (I.Belikov)
  //

  Double_t c1 = fC*fX - fE;
  Double_t r1 = TMath::Sqrt(1.0 - c1*c1);
  Double_t c2 = fC*xk - fE;
  Double_t r2 = TMath::Sqrt(1.0-  c2*c2);

  return fY + (xk-fX)*(c1+c2)/(r1+r2);

}

//_____________________________________________________________________________
void AliTRDtrack::SetSampledEdx(Float_t q, Int_t i)    
{
  //
  // The sampled energy loss
  //

  Double_t s = GetSnp();
  Double_t t = GetTgl();

  q *= TMath::Sqrt((1.0 - s*s) / (1.0 + t*t));
  fdQdl[i] = q;

}     

 //_____________________________________________________________________________
void AliTRDtrack::SetSampledEdx(Float_t q) 
{
  //
  // The sampled energy loss
  //

  Double_t s = GetSnp();
  Double_t t = GetTgl();

  q*= TMath::Sqrt((1.0 - s*s) / (1.0 + t*t));
  fdQdl[fNdedx] = q;
  fNdedx++;

}     

//_____________________________________________________________________________
void AliTRDtrack::GetXYZ(Float_t r[3]) const 
{
  //
  // Returns the position of the track in the global coord. system 
  //

  Double_t cs = TMath::Cos(fAlpha);
  Double_t sn = TMath::Sin(fAlpha);
  r[0] = fX*cs - fY*sn; 
  r[1] = fX*sn + fY*cs; 
  r[2] = fZ;

}

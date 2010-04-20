/**************************************************************************
 * Copyright(c) 2009-2011, ALICE Experiment at CERN, All rights reserved. *
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
///////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                           //
// The line is defined by equations (1)                                                      //
// a0*z+a1*x-a0*a1=0 and                                                                     //
// b0*z+b1*y-b0*b1=0                                                                         //
// where x,y,z are NOT the lab axes but z is the lab axis along which the track              //
// has the largest lever arm and x,y are the remaining 2 axis in                             //
// the order of fgkAxisID[z][0], fgkAxisID[z][1]                                             //
// The parameters are fParams[kA0,kB0,kA1,kB1] and the axis chosen as the independent        //
// var. is fParAxis (i.e. if fParAxis==kZ, then a0=ax,b0=bx, a1=ay,b1=by)                    //
//                                                                                           //
//                                                                                           //
// The helix is defined by the equations (2)                                                 //
// X(t) = (dr+R)*cos(phi0) - (R+sum{dRi})*cos(t+phi0) + sum{dRi*cos(phi0+ti)}                //
// Y(t) = (dr+R)*sin(phi0) - (R+sum{dRi})*sin(t+phi0) + sum{dRi*sin(phi0+ti)}                //
// Z(t) = dz - (R+sum{dRi})*t*tg(dip) + sum{dRi*ti}*tg(dip)                                  //
// where dRi is the change of the radius due to the ELoss at parameter ti                    //
//                                                                                           //
// Author: ruben.shahoyan@cern.ch                                                            //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "AliITSTPArrayFit.h"
#include "AliExternalTrackParam.h"
#include "AliSymMatrix.h"
#include "AliLog.h"
#include "AliParamSolver.h"
#include "AliGeomManager.h"
#include "AliITSgeomTGeo.h"
#include "AliTracker.h"
#include <TRandom.h>
#include <TArrayD.h>

ClassImp(AliITSTPArrayFit)

const Int_t  AliITSTPArrayFit::fgkAxisID[3][3] = { 
  {AliITSTPArrayFit::kY,AliITSTPArrayFit::kZ,AliITSTPArrayFit::kX},
  {AliITSTPArrayFit::kZ,AliITSTPArrayFit::kX,AliITSTPArrayFit::kY},
  {AliITSTPArrayFit::kX,AliITSTPArrayFit::kY,AliITSTPArrayFit::kZ} };

const Int_t  AliITSTPArrayFit::fgkAxisCID[3][6] = { 
  {AliITSTPArrayFit::kYY,AliITSTPArrayFit::kYZ,AliITSTPArrayFit::kXY,
   AliITSTPArrayFit::kZZ,AliITSTPArrayFit::kXZ,AliITSTPArrayFit::kXX},
  //
  {AliITSTPArrayFit::kZZ,AliITSTPArrayFit::kXZ,AliITSTPArrayFit::kYZ,
   AliITSTPArrayFit::kXX,AliITSTPArrayFit::kYX,AliITSTPArrayFit::kYY},
  //
  {AliITSTPArrayFit::kXX,AliITSTPArrayFit::kXY,AliITSTPArrayFit::kXZ,
   AliITSTPArrayFit::kYY,AliITSTPArrayFit::kYZ,AliITSTPArrayFit::kZZ}
};
//

const Double_t AliITSTPArrayFit::fgkAlmostZero = 1E-55;
const Double_t AliITSTPArrayFit::fgkCQConv = 0.299792458e-3;// R = PT/Bz/fgkCQConv with GeV,kGauss,cm
const Double_t AliITSTPArrayFit::fgkZSpanITS[AliITSTPArrayFit::kMaxLrITS] = {
  36. ,14.1,14.1,  38., 22.2,29.7, 51.   ,43.1,48.9};

const Double_t AliITSTPArrayFit::fgkRLayITS[AliITSTPArrayFit::kMaxLrITS] = {
  2.94, 3.9,7.6, 11.04, 15.0,23.9, 29.44 ,38.0,43.0};

const Int_t    AliITSTPArrayFit::fgkPassivLrITS[3] = 
  {AliITSTPArrayFit::kLrBeamPime,AliITSTPArrayFit::kLrShield1,AliITSTPArrayFit::kLrShield2};

const Int_t    AliITSTPArrayFit::fgkActiveLrITS[6] = 
  {AliITSTPArrayFit::kLrSPD1,AliITSTPArrayFit::kLrSPD2,
   AliITSTPArrayFit::kLrSDD1,AliITSTPArrayFit::kLrSDD2,
   AliITSTPArrayFit::kLrSSD1,AliITSTPArrayFit::kLrSSD2};

Double_t AliITSTPArrayFit::fgRhoLITS[AliITSTPArrayFit::kMaxLrITS] = {
  1.48e-01, 2.48e-01,2.57e-01, 1.34e-01, 3.34e-01,3.50e-01, 2.22e-01, 2.38e-01,2.25e-01};

//____________________________________________________
AliITSTPArrayFit::AliITSTPArrayFit() :
  fkPoints(0),fParSol(0),fBz(0),fCharge(0),fPntFirst(-1),
  fPntLast(-1),fNPBooked(0),fParAxis(-1),fCovI(0),fChi2NDF(0),
  fMaxIter(20),fIter(0),fEps(1e-6),fMass(0),fSwitch2Line(kFALSE),fMaxRforHelix(6.e5),
  fkAxID(0),fkAxCID(0),fCurT(0),
  fFirstPosT(0),fNElsPnt(0),fElsId(0),fElsDR(0)
{
  // default constructor
  for (int i=kMaxParam;i--;)   fParams[i] = 0;
  for (int i=kMaxParamSq;i--;) fParamsCov[i] = 0;
  SetMass();
}

//____________________________________________________
AliITSTPArrayFit::AliITSTPArrayFit(Int_t np) :
  fkPoints(0),fParSol(0),fBz(0),fCharge(0),fPntFirst(-1),
  fPntLast(-1),fNPBooked(np),fParAxis(-1),fCovI(0),fChi2NDF(0),
  fMaxIter(20),fIter(0),fEps(1e-6),fMass(0),fSwitch2Line(kFALSE),fMaxRforHelix(6.e5),
  fkAxID(0),fkAxCID(0),fCurT(0),fFirstPosT(0),fNElsPnt(0),fElsId(0),fElsDR(0)
{
  // constructor with booking of np points
  for (int i=kMaxParam;i--;)   fParams[i] = 0;
  for (int i=kMaxParamSq;i--;) fParamsCov[i] = 0;
  InitAux();
  SetEps();
  SetMass();
  SetMaxIterations();
}

//____________________________________________________
AliITSTPArrayFit::AliITSTPArrayFit(const AliITSTPArrayFit &src) : 
  TObject(src),fkPoints(src.fkPoints),fParSol(0),fBz(src.fBz),
  fCharge(src.fCharge),fPntFirst(src.fPntFirst),fPntLast(src.fPntLast),fNPBooked(src.fNPBooked),
  fParAxis(src.fParAxis),fCovI(0),fChi2NDF(0),fMaxIter(20),fIter(0),fEps(0),fMass(src.fMass),
  fSwitch2Line(src.fSwitch2Line),fMaxRforHelix(src.fMaxRforHelix),fkAxID(0),fkAxCID(0),fCurT(0),
  fFirstPosT(0),fNElsPnt(0),fElsId(0),fElsDR(0)
{
  // copy constructor
  InitAux();
  memcpy(fCovI,src.fCovI,fNPBooked*6*sizeof(Double_t));
  for (int i=kMaxParam;i--;)   fParams[i] = src.fParams[i];
  for (int i=kMaxParamSq;i--;) fParamsCov[i] = src.fParamsCov[i];
  memcpy(fCurT,src.fCurT,fNPBooked*sizeof(Double_t));
  SetEps(src.fEps);
  SetMaxIterations(src.fMaxIter);
  //  
}

//____________________________________________________
AliITSTPArrayFit &AliITSTPArrayFit::operator =(const AliITSTPArrayFit& src)
{
  // assignment operator
  if (this==&src) return *this;
  ((TObject*)this)->operator=(src);
  fkPoints   = src.fkPoints;
  if (!fParSol) fParSol = new AliParamSolver(*src.fParSol);
  else *fParSol = *src.fParSol; 
  fBz       = src.fBz; 
  fCharge   = src.fCharge;
  fNPBooked = src.fNPBooked;
  fPntFirst = src.fPntFirst;
  fPntLast  = src.fPntLast;
  InitAux();
  memcpy(fCovI,src.fCovI,fNPBooked*6*sizeof(Double_t));
  for (int i=kMaxParam;i--;)   fParams[i] = src.fParams[i];
  for (int i=kMaxParamSq;i--;) fParamsCov[i] = src.fParamsCov[i];
  SetParAxis(src.fParAxis);
  fNElsPnt   = src.fNElsPnt;
  fFirstPosT = src.fFirstPosT;
  memcpy(fCurT  ,src.fCurT  ,fNPBooked*sizeof(Double_t));
  memcpy(fElsId ,src.fElsId ,fNPBooked*sizeof(Int_t));
  memcpy(fElsDR ,src.fElsDR ,fNPBooked*sizeof(Double_t));
  memcpy(fCurT  ,src.fCurT  ,fNPBooked*sizeof(Double_t));
  SetEps(src.fEps);
  SetMaxIterations(src.fMaxIter);
  //
  return *this;
  //
}

//____________________________________________________
AliITSTPArrayFit::~AliITSTPArrayFit()
{
  // destructor
  delete   fParSol;
  delete[] fCovI;
  delete[] fCurT;
  delete[] fElsId;
  delete[] fElsDR;
}

//____________________________________________________
void AliITSTPArrayFit::Reset()
{
  // reset to process new track
  if (fParSol) fParSol->Clear();
  fkPoints=0; 
  fNElsPnt = 0;
  fFirstPosT = 0;
  //  fBz = 0;
  fCharge = 0;
  fIter = 0;
  fPntFirst=fPntLast=-1; 
  SetParAxis(-1);
  fSwitch2Line = kFALSE;
  ResetBit(kFitDoneBit|kCovInvBit);
}

//____________________________________________________
void AliITSTPArrayFit::AttachPoints(const AliTrackPointArray* points, Int_t pfirst,Int_t plast) 
{
  // create from piece of AliTrackPointArray
  Reset();
  fkPoints = points;
  int np = points->GetNPoints();
  if (fNPBooked<np) {
    fNPBooked = np;
    InitAux();
  }
  fPntFirst = pfirst<0 ? 0 : pfirst;
  fPntLast  = plast<fPntFirst ? np-1 : plast;
  //
  for (int i=kMaxParam;i--;)   fParams[i] = 0;
  for (int i=kMaxParamSq;i--;) fParamsCov[i] = 0;
  //
  InvertPointsCovMat();
  //
}

//____________________________________________________
Bool_t AliITSTPArrayFit::SetFirstLast(Int_t pfirst,Int_t plast) 
{
  // set first and last point to fit
  const AliTrackPointArray* pnts = fkPoints;
  if (!pnts) {AliError("TrackPointArray is not attached yet"); return kFALSE;}
  AttachPoints(pnts,pfirst,plast);
  return kTRUE;
  //
}

//____________________________________________________
Bool_t AliITSTPArrayFit::InvertPointsCovMat()
{
  // invert the cov.matrices of the points
  for (int i=fPntFirst;i<=fPntLast;i++) {
    //
    float *cov = (float*)fkPoints->GetCov() + i*6; // pointer on cov.matrix
    //
    Double_t t0 = cov[kYY]*cov[kZZ] - cov[kYZ]*cov[kYZ];
    Double_t t1 = cov[kXY]*cov[kZZ] - cov[kXZ]*cov[kYZ];
    Double_t t2 = cov[kXY]*cov[kYZ] - cov[kXZ]*cov[kYY];
    Double_t det = cov[kXX]*t0 - cov[kXY]*t1 + cov[kXZ]*t2;
    if (IsZero(det,1e-18)) { // one of errors is 0, fix this
      double norm[3];
      TGeoHMatrix hcov;
      TGeoRotation hrot,hrotI;
      GetNormal(norm,cov);
      double phi = TMath::ATan2(norm[1],norm[0]);
      hrot.SetAngles(-phi*TMath::RadToDeg(),0.,0.);
      hrotI.SetAngles(phi*TMath::RadToDeg(),0.,0.);
      //      
      Double_t hcovel[9];
      hcovel[0] = cov[kXX];
      hcovel[1] = cov[kXY];
      hcovel[2] = cov[kXZ];
      hcovel[3] = cov[kXY];
      hcovel[4] = cov[kYY];
      hcovel[5] = cov[kYZ];
      hcovel[6] = cov[kXZ];
      hcovel[7] = cov[kYZ];
      hcovel[8] = cov[kZZ];
      hcov.SetRotation(hcovel);
      //
      Double_t *hcovscl = hcov.GetRotationMatrix(); 
      //      printf(">> %+e %+e %+e\n   %+e %+e %+e\n   %+e %+e %+e\n\n",hcovscl[0],hcovscl[1],hcovscl[2],hcovscl[3],hcovscl[4],hcovscl[5],hcovscl[6],hcovscl[7],hcovscl[8]);
      //      printf("Rot by %+.e (%+.e %+.e %+.e)\n",phi, norm[0],norm[1],norm[2]);
      hcov.Multiply(&hrotI);
      hcov.MultiplyLeft(&hrot);
      //      printf("|| %+e %+e %+e\n   %+e %+e %+e\n   %+e %+e %+e\n\n",hcovscl[0],hcovscl[1],hcovscl[2],hcovscl[3],hcovscl[4],hcovscl[5],hcovscl[6],hcovscl[7],hcovscl[8]);
      if (hcovscl[0]<1e-8) hcovscl[0] = 1e-8;
      if (hcovscl[4]<1e-8) hcovscl[4] = 1e-8;
      if (hcovscl[8]<1e-8) hcovscl[8] = 1e-8;
      //      printf("** %+e %+e %+e\n   %+e %+e %+e\n   %+e %+e %+e\n\n",hcovscl[0],hcovscl[1],hcovscl[2],hcovscl[3],hcovscl[4],hcovscl[5],hcovscl[6],hcovscl[7],hcovscl[8]);
      hcov.Multiply(&hrot);
      hcov.MultiplyLeft(&hrotI);
      //      printf("^^ %+e %+e %+e\n   %+e %+e %+e\n   %+e %+e %+e\n\n",hcovscl[0],hcovscl[1],hcovscl[2],hcovscl[3],hcovscl[4],hcovscl[5],hcovscl[6],hcovscl[7],hcovscl[8]);
      cov[kXX] = hcovscl[0];
      cov[kXY] = hcovscl[1];
      cov[kXZ] = hcovscl[2];
      cov[kYY] = hcovscl[4];
      cov[kYZ] = hcovscl[5];
      cov[kZZ] = hcovscl[8];
    }
    t0 = cov[kYY]*cov[kZZ] - cov[kYZ]*cov[kYZ];
    t1 = cov[kXY]*cov[kZZ] - cov[kXZ]*cov[kYZ];
    t2 = cov[kXY]*cov[kYZ] - cov[kXZ]*cov[kYY];
    det = cov[kXX]*t0 - cov[kXY]*t1 + cov[kXZ]*t2;
    //
    AliDebug(2,Form("%+.4e %+.4e %+.4e -> %+.4e",t0,t1,t2,det));
    if (IsZero(det,fgkAlmostZero)) {
      AliInfo(Form("Cov.Matrix for point %d is singular",i));
      return kFALSE;
    }
    //
    Double_t *covI = GetCovI(i);
    covI[kXX] =  t0/det;
    covI[kXY] = -t1/det;
    covI[kXZ] =  t2/det;
    covI[kYY] =  (cov[kXX]*cov[kZZ] - cov[kXZ]*cov[kXZ])/det;
    covI[kYZ] =  (cov[kXY]*cov[kXZ] - cov[kXX]*cov[kYZ])/det;
    covI[kZZ] =  (cov[kXX]*cov[kYY] - cov[kXY]*cov[kXY])/det;
    //
  }
  SetCovInv();
  return kTRUE;
}

//____________________________________________________
void AliITSTPArrayFit::InitAux()
{
  // init auxiliary space
  if (fCovI) delete[] fCovI;
  if (fCurT) delete[] fCurT;
  //
  fCovI   = new Double_t[6*fNPBooked];
  fCurT   = new Double_t[fNPBooked+kMaxLrITS];
  fElsId  = new Int_t[fNPBooked+kMaxLrITS];
  fElsDR  = new Double_t[fNPBooked+kMaxLrITS];
  memset(fElsDR,0,(fNPBooked+kMaxLrITS)*sizeof(Double_t));
  memset(fCovI,0,fNPBooked*6*sizeof(Double_t));
  //
}

//____________________________________________________
Bool_t AliITSTPArrayFit::FitLineCrude()
{
  // perform linear fit w/o accounting the errors
  // fit is done in the parameterization
  // x = res[0] + res[1]*z
  // y = res[2] + res[3]*z
  // where x,y,z are NOT the lab axes but z is the lab axis along which the track 
  // has the largest lever arm and x,y are the remaining 2 axis in 
  // the order of fgkAxisID[z][0], fgkAxisID[z][1]
  //
  int np = fPntLast - fPntFirst + 1;
  if (np<2) {
    AliError("At least 2 points are needed for straight line fit");
    return kFALSE;
  }
  //
  if (fParAxis<0) SetParAxis(ChoseParAxis());
  Double_t sZ=0,sZZ=0,sY=0,sYZ=0,sX=0,sXZ=0,det=0;
  //
  const float *coord[3] = {fkPoints->GetX(),fkPoints->GetY(),fkPoints->GetZ()};
  const Float_t *varZ = coord[ fParAxis  ];
  const Float_t *varX = coord[ fkAxID[kX] ];
  const Float_t *varY = coord[ fkAxID[kY] ];
  //
  for (int i=fPntFirst;i<=fPntLast;i++) {
    sZ  += varZ[i];
    sZZ += varZ[i]*varZ[i];
    //
    sX  += varX[i];
    sXZ += varX[i]*varZ[i];
    //
    sY  += varY[i];
    sYZ += varY[i]*varZ[i];
  }
  det = sZZ*np-sZ*sZ;
  if (TMath::Abs(det)<fgkAlmostZero) return kFALSE;
  fParams[0] = (sX*sZZ-sZ*sXZ)/det;
  fParams[1] = (sXZ*np-sZ*sX)/det;
  //
  fParams[2] = (sY*sZZ-sZ*sYZ)/det;
  fParams[3] = (sYZ*np-sZ*sY)/det;
  //
  return kTRUE;
  //
}

//____________________________________________________
void AliITSTPArrayFit::SetParAxis(Int_t ax)
{
  // select the axis which will be used as a parameter for the line: longest baseline
  if (ax>kZ) {
    AliInfo(Form("Wrong axis choice: %d",ax));
    fParAxis = -1;
  }
  fParAxis = ax;
  if (ax>=0) {
    fkAxID  = fgkAxisID[ax];
    fkAxCID = fgkAxisCID[ax];
  }
  else {
    fkAxID = fkAxCID = 0;
  }
  //
}

//____________________________________________________
Int_t AliITSTPArrayFit::ChoseParAxis() const
{
  // select the variable with largest base as a parameter
  Double_t cmn[3]={1.e9,1.e9,1.e9},cmx[3]={-1.e9,-1.e9,-1.e9};
  //
  const float *coord[3] = {fkPoints->GetX(),fkPoints->GetY(),fkPoints->GetZ()};
  for (int i=fPntFirst;i<=fPntLast;i++) {
    for (int j=3;j--;) {
      Double_t val = coord[j][i];
      if (cmn[j]>val) cmn[j] = val;
      if (cmx[j]<val) cmx[j] = val;
    }
  }
  //
  int axis = kZ;
  if (cmx[axis]-cmn[axis] < cmx[kX]-cmn[kX]) axis = kX;
  if (cmx[axis]-cmn[axis] < cmx[kY]-cmn[kY]) axis = kY;
  return axis;
  //
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetPosition(Double_t *xyzPCA, const Double_t *xyz, const Double_t *covI) const
{
  // calculate the position of the track at PCA to xyz
  Double_t t = GetParPCA(xyz,covI);
  GetPosition(xyzPCA,t);
  return t;
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetPosition(Double_t *xyzPCA, const AliTrackPoint *pntCovInv, Bool_t useErr) const
{
  // calculate the position of the track at PCA to pntCovInv
  // NOTE: the covariance matrix of the point must be inverted
  double t;
  double xyz[3] = {pntCovInv->GetX(),pntCovInv->GetY(),pntCovInv->GetZ()};
  if (useErr) {
    Double_t covI[6];;
    for (int i=6;i--;) covI[i] = pntCovInv->GetCov()[i];
    t = GetParPCA(xyz,covI);
  }
  else t = GetParPCA(xyz);
  GetPosition(xyzPCA,t);
  return t;
}

//____________________________________________________
void AliITSTPArrayFit::GetResiduals(Double_t *resPCA, const AliTrackPoint *pntCovInv, Bool_t useErr) const
{
  // calculate the residuals  of the track at PCA to pntCovInv
  // NOTE: the covariance matrix of the point must be inverted
  GetPosition(resPCA,pntCovInv,useErr);
  resPCA[0] -= pntCovInv->GetX();
  resPCA[1] -= pntCovInv->GetY();
  resPCA[2] -= pntCovInv->GetZ();
}

//____________________________________________________
void AliITSTPArrayFit::GetResiduals(Double_t *resPCA, const Double_t *xyz, const Double_t *covI) const
{
  // calculate the residuals of the track at PCA to xyz
  GetPosition(resPCA,xyz,covI);
  resPCA[kX] -= xyz[kX];
  resPCA[kY] -= xyz[kY];
  resPCA[kZ] -= xyz[kZ];
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetParPCALine(const Double_t *xyz, const Double_t *covI) const
{
  // get parameter for the point with least weighted distance to the point
  //
  Double_t rhs,denom;
  Double_t dx = fParams[kA0]-xyz[ fkAxID[kX] ];
  Double_t dy = fParams[kA1]-xyz[ fkAxID[kY] ];
  Double_t dz =             -xyz[ fkAxID[kZ] ];
  //
  if (covI) {
    Double_t tx = fParams[kB0]*covI[ fkAxCID[kXX] ] + fParams[kB1]*covI[ fkAxCID[kXY] ] + covI[ fkAxCID[kXZ] ];
    Double_t ty = fParams[kB0]*covI[ fkAxCID[kXY] ] + fParams[kB1]*covI[ fkAxCID[kYY] ] + covI[ fkAxCID[kYZ] ];
    Double_t tz = fParams[kB0]*covI[ fkAxCID[kXZ] ] + fParams[kB1]*covI[ fkAxCID[kYZ] ] + covI[ fkAxCID[kZZ] ];
    rhs   = tx*dx + ty*dy + tz*dz;
    denom = -(fParams[kB0]*(covI[ fkAxCID[kXZ] ] + tx) + fParams[kB1]*(covI[ fkAxCID[kYZ] ] + ty) + covI[ fkAxCID[kZZ] ]);
  }
  else {
    rhs = fParams[kB0]*dx + fParams[kB1]*dy + dz;
    denom = -(fParams[kB0]*fParams[kB0] + fParams[kB1]*fParams[kB1] + 1);
  }
  //
  return rhs/denom;
  //
}

//____________________________________________________
void AliITSTPArrayFit::GetDResDPosLine(Double_t *dXYZdP, /*const Double_t *xyz,*/ const Double_t *covI) const
{
  // calculate detivative of the PCA residuals vs point position and fill in user provide
  // array in the format {dXdXp,dY/dXp,dZdXp, ... dXdZp,dYdZp,dZdZp}
  //
  Double_t dTdP[3];
  GetDtDPosLine(dTdP, /*xyz,*/ covI); // derivative of the t-param over point position
  //
  for (int i=3;i--;) {
    int var = fkAxID[i];
    Double_t *curd = dXYZdP + var*3;   // d/dCoord_i
    curd[ fkAxID[kX] ] = fParams[kB0]*dTdP[var];
    curd[ fkAxID[kY] ] = fParams[kB1]*dTdP[var];
    curd[ fkAxID[kZ] ] = dTdP[var];
    curd[    var     ]-= 1.;
  }
  //
}

//____________________________________________________
void AliITSTPArrayFit::GetDResDParamsLine(Double_t *dXYZdP, const Double_t *xyz, const Double_t *covI) const
{
  // calculate detivative of the PCA residuals vs line parameters and fill in user provide
  // array in the format {dXdP0,dYdP0,dZdP0, ... dXdPn,dYdPn,dZdPn}
  //
  Double_t dTdP[4];
  Double_t t = GetDtDParamsLine(dTdP, xyz, covI); // derivative of the t-param over line params
  //
  Double_t *curd = dXYZdP + kA0*3; // d/dA0
  curd[ fkAxID[kX] ] = fParams[kB0]*dTdP[kA0] + 1.;
  curd[ fkAxID[kY] ] = fParams[kB1]*dTdP[kA0];
  curd[ fkAxID[kZ] ] = dTdP[kA0];
  //
  curd = dXYZdP + kB0*3; // d/dB0
  curd[ fkAxID[kX] ] = fParams[kB0]*dTdP[kB0] + t;
  curd[ fkAxID[kY] ] = fParams[kB1]*dTdP[kB0];
  curd[ fkAxID[kZ] ] = dTdP[kB0];
  //
  curd = dXYZdP + kA1*3; // d/dA1
  curd[ fkAxID[kX] ] = fParams[kB0]*dTdP[kA1];
  curd[ fkAxID[kY] ] = fParams[kB1]*dTdP[kA1] + 1.;
  curd[ fkAxID[kZ] ] = dTdP[kA1];
  //
  curd = dXYZdP + kB1*3; // d/dB1
  curd[ fkAxID[kX] ] = fParams[kB0]*dTdP[kB1];
  curd[ fkAxID[kY] ] = fParams[kB1]*dTdP[kB1] + t;
  curd[ fkAxID[kZ] ] = dTdP[kB1];
  //
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetDtDParamsLine(Double_t *dtparam,const Double_t *xyz, const Double_t *covI) const
{
  // get t-param detivative over the parameters for the point with least weighted distance to the point
  //
  Double_t rhs,denom;
  Double_t dx = fParams[kA0]-xyz[ fkAxID[kX] ];
  Double_t dy = fParams[kA1]-xyz[ fkAxID[kY] ];
  Double_t dz =             -xyz[ fkAxID[kZ] ];
  Double_t rhsDA0,rhsDA1,rhsDB0,rhsDB1,denDB0,denDB1;
  //
  if (covI) {
    Double_t tx = fParams[kB0]*covI[ fkAxCID[kXX] ] + fParams[kB1]*covI[ fkAxCID[kXY] ] + covI[ fkAxCID[kXZ] ];
    Double_t ty = fParams[kB0]*covI[ fkAxCID[kXY] ] + fParams[kB1]*covI[ fkAxCID[kYY] ] + covI[ fkAxCID[kYZ] ];
    Double_t tz = fParams[kB0]*covI[ fkAxCID[kXZ] ] + fParams[kB1]*covI[ fkAxCID[kYZ] ] + covI[ fkAxCID[kZZ] ];
    rhs = tx*dx + ty*dy + tz*dz;
    denom = -(fParams[kB0]*(covI[ fkAxCID[kXZ] ] + tx) + fParams[kB1]*(covI[ fkAxCID[kYZ] ] + ty) + covI[ fkAxCID[kZZ] ]);
    //
    rhsDA0 = tx;
    rhsDA1 = ty;
    rhsDB0 = covI[ fkAxCID[kXX] ]*dx + covI[ fkAxCID[kXY] ]*dy + covI[ fkAxCID[kXZ] ]*dz;
    rhsDB1 = covI[ fkAxCID[kXY] ]*dx + covI[ fkAxCID[kYY] ]*dy + covI[ fkAxCID[kYZ] ]*dz;
    //
    denDB0 = -(tx + tx);
    denDB1 = -(ty + ty);
  }
  else {
    rhs = fParams[kB0]*dx + fParams[kB1]*dy + dz;
    denom = -(fParams[kB0]*fParams[kB0]	+ fParams[kB1]*fParams[kB1] + 1);
    //
    rhsDA0 = fParams[kB0];
    rhsDB0 = dx;
    rhsDA1 = fParams[kB1];
    rhsDB1 = dy;
    //
    denDB0 = -(fParams[kB0]+fParams[kB0]);
    denDB1 = -(fParams[kB1]+fParams[kB1]);
    //
  }
  //
  Double_t denom2 = denom*denom;
  dtparam[kA0] = rhsDA0/denom;    // denom does not depend on A0,A1
  dtparam[kA1] = rhsDA1/denom;
  dtparam[kB0] = rhsDB0/denom - rhs/denom2 * denDB0;
  dtparam[kB1] = rhsDB1/denom - rhs/denom2 * denDB1;
  //
  return rhs/denom;
}

//____________________________________________________
void AliITSTPArrayFit::GetDtDPosLine(Double_t *dtpos,/*const Double_t *xyz,*/ const Double_t *covI) const
{
  // get t-param detivative over the parameters for the point with least weighted distance to the point
  //
  //  Double_t rhs;
  //  Double_t dx = fParams[kA0]-xyz[ fkAxID[kX] ];
  //  Double_t dy = fParams[kA1]-xyz[ fkAxID[kY] ];
  //  Double_t dz =             -xyz[ fkAxID[kZ] ];
  Double_t denom;
  Double_t rhsDX,rhsDY,rhsDZ;
  //
  if (covI) {
    Double_t tx = fParams[kB0]*covI[ fkAxCID[kXX] ] + fParams[kB1]*covI[ fkAxCID[kXY] ] + covI[ fkAxCID[kXZ] ];
    Double_t ty = fParams[kB0]*covI[ fkAxCID[kXY] ] + fParams[kB1]*covI[ fkAxCID[kYY] ] + covI[ fkAxCID[kYZ] ];
    Double_t tz = fParams[kB0]*covI[ fkAxCID[kXZ] ] + fParams[kB1]*covI[ fkAxCID[kYZ] ] + covI[ fkAxCID[kZZ] ];
    // rhs = tx*dx + ty*dy + tz*dz;
    denom = -(fParams[kB0]*(covI[ fkAxCID[kXZ] ] + tx) + fParams[kB1]*(covI[ fkAxCID[kYZ] ] + ty) + covI[ fkAxCID[kZZ] ]);
    //
    rhsDX = -tx;
    rhsDY = -ty;
    rhsDZ = -tz;
  }
  else {
    // rhs = fParams[kB0]*dx + fParams[kB1]*dy + dz;
    denom = -(fParams[kB0]*fParams[kB0]	+ fParams[kB1]*fParams[kB1] + 1);
    //
    rhsDX = -fParams[kB0];
    rhsDY = -fParams[kB1];
    rhsDZ = -1;
    //
  }
  //
  dtpos[ fkAxID[kX] ] = rhsDX/denom;
  dtpos[ fkAxID[kY] ] = rhsDY/denom;
  dtpos[ fkAxID[kZ] ] = rhsDZ/denom;
  //
  //  return rhs/denom;
}

//____________________________________________________
void AliITSTPArrayFit::GetDResDParamsLine(Double_t *dXYZdP, Int_t ipnt) const
{
  // calculate detivative of the PCA residuals vs line parameters and fill in user provide
  // array in the format {dXdP0,dYdP0,dZdP0, ... dXdPn,dYdPn,dZdPn}
  //
  if (ipnt<fPntFirst || ipnt>fPntLast) {
    AliError(Form("Attempt to access the point %d not in the fitted points [%d:%d]",ipnt,fPntFirst,fPntLast));
    return;
  }
  GetDResDParamsLine(dXYZdP, GetPoint(ipnt) , IsCovIgnored() ? 0 : GetCovI(ipnt));
}

//____________________________________________________
void AliITSTPArrayFit::GetDResDPosLine(Double_t *dXYZdP, Int_t ipnt) const
{
  // calculate detivative of the PCA residuals vs point position and fill in user provide
  // array in the format {dXdXp,dY/dXp,dZdXp, ... dXdZp,dYdZp,dZdZp}
  //
  if (ipnt<fPntFirst || ipnt>fPntLast) {
    AliError(Form("Attempt to access the point %d not in the fitted points [%d:%d]",ipnt,fPntFirst,fPntLast));
    return;
  }
  GetDResDPosLine(dXYZdP,IsCovIgnored() ? 0 : GetCovI(ipnt));
}

//____________________________________________________
void AliITSTPArrayFit::GetDResDParams(Double_t *dXYZdP, Int_t ipnt)
{
  // calculate detivative of the PCA residuals vs track parameters and fill in user provide
  // array in the format {dXdP0,dYdP0,dZdP0, ... dXdPn,dYdPn,dZdPn}
  //
  if (ipnt<fPntFirst || ipnt>fPntLast) {
    AliError(Form("Attempt to access the point %d not in the fitted points [%d:%d]",ipnt,fPntFirst,fPntLast));
    return;
  }
  GetDResDParams(dXYZdP, GetPoint(ipnt) , IsCovIgnored() ? 0 : GetCovI(ipnt));
}

//____________________________________________________
void AliITSTPArrayFit::GetDResDPos(Double_t *dXYZdP, Int_t ipnt)
{
  // calculate detivative of the PCA residuals vs point position and fill in user provide
  // array in the format {dXdXp,dY/dXp,dZdXp, ... dXdZp,dYdZp,dZdZp} 
  //
  if (ipnt<fPntFirst || ipnt>fPntLast) {
    AliError(Form("Attempt to access the point %d not in the fitted points [%d:%d]",ipnt,fPntFirst,fPntLast));
    return;
  }
  GetDResDPos(dXYZdP, GetPoint(ipnt), IsCovIgnored() ? 0 : GetCovI(ipnt));
}

//____________________________________________________
void AliITSTPArrayFit::GetDResDParams(Double_t *dXYZdP, const Double_t *xyz, const Double_t *covI) 
{
  // get residual detivatives over the track parameters for the point with least weighted distance to the point
  //
  if (!IsHelix()) { // for the straight line calculate analytically
    GetDResDParamsLine(dXYZdP, xyz, covI);
    return;
  }
  //
  // calculate derivative numerically
  const Double_t delta = 0.01;
  Double_t xyzVar[4][3];
  //
  for (int ipar = 5;ipar--;) {
    double sav = fParams[ipar];
    fParams[ipar] -= delta;
    GetPosition(xyzVar[0],xyz,covI);
    fParams[ipar] += delta/2;
    GetPosition(xyzVar[1],xyz,covI);
    fParams[ipar] += delta;
    GetPosition(xyzVar[2],xyz,covI);
    fParams[ipar] += delta/2;
    GetPosition(xyzVar[3],xyz,covI);
    fParams[ipar] = sav;  // restore
    //
    double *curd = dXYZdP + 3*ipar;
    for (int i=3;i--;) curd[i] = (8.*(xyzVar[2][i]-xyzVar[1][i]) - (xyzVar[3][i]-xyzVar[0][i]))/6./delta;
  }
  //
}


//____________________________________________________
void AliITSTPArrayFit::GetDResDPos(Double_t *dXYZdP, const Double_t *xyz, const Double_t *covI) 
{
  // get residuals detivative over the point position for the point with least weighted distance to the point
  //

  if (!IsHelix()) { // for the straight line calculate analytically
    GetDResDPosLine(dXYZdP, /*xyz,*/ covI);
    return;
  }
  //
  // calculate derivative numerically
  const Double_t delta = 0.005;
  Double_t xyzVar[4][3];
  Double_t xyzv[3] = {xyz[0],xyz[1],xyz[2]};
  //
  for (int ipar = 3;ipar--;) {
    double sav = xyzv[ipar];
    xyzv[ipar] -= delta;
    GetPosition(xyzVar[0],xyzv,covI);
    xyzv[ipar] += delta/2;
    GetPosition(xyzVar[1],xyzv,covI);
    xyzv[ipar] += delta;
    GetPosition(xyzVar[2],xyzv,covI);
    xyzv[ipar] += delta/2;
    GetPosition(xyzVar[3],xyzv,covI);
    xyzv[ipar] = sav;  // restore
    //
    double *curd = dXYZdP + 3*ipar;
    for (int i=3;i--;) curd[i] = (8.*(xyzVar[2][i]-xyzVar[1][i]) - (xyzVar[3][i]-xyzVar[0][i]))/6./delta;
    curd[ipar] -= 1.;
  }
  //
}

//________________________________________________________________________________________________________
Double_t AliITSTPArrayFit::GetParPCAHelix(const Double_t* xyz, const Double_t* covI) const
{
  // find track parameter t (eq.2) corresponding to point of closest approach to xyz
  //
  Double_t phi  = GetParPCACircle(xyz[kX],xyz[kY]); 
  Double_t cs = TMath::Cos(fParams[kPhi0]);
  Double_t sn = TMath::Sin(fParams[kPhi0]);
  Double_t xc = (fParams[kD0]+fParams[kR0])*cs;
  Double_t yc = (fParams[kD0]+fParams[kR0])*sn;
  Double_t dchi2,ddchi2;
  //
  Double_t dzD  = -fParams[kR0]*fParams[kDip];
  Double_t dphi = 0;
  //
  double rEps = 1e-5/TMath::Abs(fParams[kR0]); // dphi corresponding to 0.1 micron
  if (rEps>fEps) rEps = fEps;
  //
  int it=0;
  do {
    cs = TMath::Cos(phi + fParams[kPhi0]);
    sn = TMath::Sin(phi + fParams[kPhi0]);
    //
    Double_t dxD  =  fParams[kR0]*sn;
    Double_t dyD  = -fParams[kR0]*cs;
    Double_t dxDD = -dyD;
    Double_t dyDD =  dxD;
    //
    Double_t dx   = xc - fParams[kR0]*cs - xyz[kX];
    Double_t dy   = yc - fParams[kR0]*sn - xyz[kY];
    Double_t dz   = fParams[kDZ] + dzD*phi- xyz[kZ];
    //
    if (covI) {
      Double_t tx = dx*covI[kXX] + dy*covI[kXY] + dz*covI[kXZ];
      Double_t ty = dx*covI[kXY] + dy*covI[kYY] + dz*covI[kYZ];
      Double_t tz = dx*covI[kXZ] + dy*covI[kYZ] + dz*covI[kZZ];
      //
      Double_t ttx = dxD*covI[kXX] + dyD*covI[kXY] + dzD*covI[kXZ];
      Double_t tty = dxD*covI[kXY] + dyD*covI[kYY] + dzD*covI[kYZ];
      Double_t ttz = dxD*covI[kXZ] + dyD*covI[kYZ] + dzD*covI[kZZ];
      //
      // chi2   = dx*tx + dy*ty + dz*tz;
      dchi2  = dxD*tx  + dyD*ty  + dzD*tz;
      ddchi2 = dxDD*tx + dyDD*ty           + dxD *ttx + dyD *tty + dzD *ttz;
      //
    }
    else {
      // chi2   = dx*dx + dy*dy + dz*dz;
      dchi2  = dxD*dx  + dyD*dy  + dzD*dz;
      ddchi2 = dxDD*dx + dyDD*dy +         + dxD*dxD + dyD*dyD + dzD*dzD;
    }
    //
    if (TMath::Abs(ddchi2)<fgkAlmostZero || TMath::Abs(dphi=dchi2/ddchi2)<rEps) break;
    phi -= dphi;    
  } while(++it<fMaxIter);

  //
  return phi;
}

//________________________________________________________________________________________________________
Double_t AliITSTPArrayFit::GetParPCACircle(Double_t x,Double_t y)  const
{
  // find track parameter t (eq.2) corresponding to point on the circle with closest approach to x,y
  //
  Double_t r = fParams[kD0]+fParams[kR0];
  Double_t t = TMath::ATan2( r*TMath::Sin(fParams[kPhi0])-y, r*TMath::Cos(fParams[kPhi0])-x ) - fParams[kPhi0]; 
  if (fParams[kR0] < 0) t += TMath::Pi();
  if (t > TMath::Pi())  t -= TMath::Pi()*2;
  if (t <-TMath::Pi())  t += TMath::Pi()*2;
  return t;
}

//________________________________________________________________________________________________________
Double_t AliITSTPArrayFit::GetHelixParAtR(Double_t r)  const
{
  // find helix parameter t (eq.2) corresponding to point on the circle of radius t
  //
  double gam = 1. - (r-fParams[kD0])*(r+fParams[kD0])/fParams[kR0]/(fParams[kD0]+fParams[kR0])/2.;
  return (TMath::Abs(gam)>1) ?  -1e9 : TMath::ACos(gam);
}

//________________________________________________________________________________________________________
Double_t AliITSTPArrayFit::CalcChi2NDF() const
{
  // calculate fit chi2/ndf
  Double_t chi2 = 0;
  Double_t dr[3]; // residuals
  //if (!IsFitDone()) return -1;
  for (int ipnt=fPntFirst;ipnt<=fPntLast;ipnt++) {
    GetResiduals(dr,ipnt);
    Double_t* covI = GetCovI(ipnt);
    chi2 += dr[kX]*(dr[kX]*covI[ kXX ]+dr[kY]*covI[ kXY ]+dr[kZ]*covI[ kXZ ])
      +     dr[kY]*(dr[kX]*covI[ kXY ]+dr[kY]*covI[ kYY ]+dr[kZ]*covI[ kYZ ])
      +     dr[kZ]*(dr[kX]*covI[ kXZ ]+dr[kY]*covI[ kYZ ]+dr[kZ]*covI[ kZZ ]);
  }
  int ndf = (fPntLast-fPntFirst+1)*3 - GetNParams();
  chi2 /= ndf;
  return chi2;
}

//________________________________________________________________________________________________________
void AliITSTPArrayFit::GetResiduals(Double_t *res,Int_t ipnt) const
{
  // calculate residuals at point
  if (ipnt<fPntFirst || ipnt>fPntLast) {
    AliError(Form("Attempt to access the point %d not in the fitted points [%d:%d]",ipnt,fPntFirst,fPntLast));
    return;
  }
  GetPosition(res,fCurT[ipnt]);
  res[kX] -= fkPoints->GetX()[ipnt];
  res[kY] -= fkPoints->GetY()[ipnt];
  res[kZ] -= fkPoints->GetZ()[ipnt];
}

//________________________________________________________________________________________________________
void AliITSTPArrayFit::GetPosition(Double_t *xyz, Double_t t) const
{
  // calculate track position for parameter value t
  if (IsHelix()) {
    //
    Double_t rrho = fParams[kD0]+fParams[kR0];
    Double_t xc = rrho*TMath::Cos(fParams[kPhi0]);
    Double_t yc = rrho*TMath::Sin(fParams[kPhi0]);
    Double_t r  = fParams[kR0];
    Double_t ze = 0;
    //
    if (IsELossON()) {
      if (t>0) {
	for (int i=fFirstPosT;i<fNElsPnt;i++) { // along the track direction
	  int indE = fElsId[i];
	  if ( t<fCurT[indE] ) break;       // does not reach this layer on  its way to t 
	  xc += fElsDR[indE] * TMath::Cos(fParams[kPhi0] + fCurT[indE]);
	  yc += fElsDR[indE] * TMath::Sin(fParams[kPhi0] + fCurT[indE]);
	  ze += fElsDR[indE] * fCurT[indE];
	  r  += fElsDR[indE];
	  //printf("ELoss@ %+.2e r:%+.3e got %+.3e\n",fCurT[indE],r,fElsDR[indE]);
	}
      }	else {
	for (int i=fFirstPosT;i--;) { // against the track direction
	  int indE = fElsId[i];
	  if ( t>=fCurT[indE] ) break;       // does not reach this layer on  its way to t 
	  xc += fElsDR[indE] * TMath::Cos(fParams[kPhi0] + fCurT[indE]);
	  yc += fElsDR[indE] * TMath::Sin(fParams[kPhi0] + fCurT[indE]);
	  ze += fElsDR[indE] * fCurT[indE];
	  r  += fElsDR[indE];
	  //printf("ELoss@ %+.2e r:%+.3e got %+.3e\n",fCurT[indE],r,fElsDR[indE]);
	}     
      }
    }
    //
    xyz[kZ] = fParams[kDZ] - fParams[kDip]*(t*r - ze);
    //
    t += fParams[kPhi0];    
    xyz[kX] = xc - r*TMath::Cos(t);
    xyz[kY] = yc - r*TMath::Sin(t);
    //    printf("t: %+.3e xyz:%+.2e %+.2e %+.2e | R %+.6e -> %+.6e | sign %d\n",t-fParams[kPhi0],xyz[0],xyz[1],xyz[2],fParams[kR0],r,GetSignQB());
  }
  else {
    xyz[ fkAxID[kX] ] = fParams[kA0] + fParams[kB0]*t;
    xyz[ fkAxID[kY] ] = fParams[kA1] + fParams[kB1]*t;
    xyz[ fParAxis   ] = t;
  }
}

//________________________________________________________________________________________________________
void AliITSTPArrayFit::GetDirCos(Double_t *dircos, Double_t t) const
{
  // calculate track direction cosines for parameter value t
  if (IsHelix()) {
    dircos[kZ] = -fParams[kDip];
    t += fParams[kPhi0];    
    dircos[kX] = TMath::Sin(t);
    dircos[kY] =-TMath::Cos(t);
    double gam = TMath::Sign(1/TMath::Sqrt(dircos[kZ]*dircos[kZ]+dircos[kY]*dircos[kY]+dircos[kX]*dircos[kX]),fParams[kR0]);
    for (int i=3;i--;) dircos[i] *= gam;
    if (GetSignQB()>0) for (int i=3;i--;) dircos[i] = -dircos[i]; // positive tracks move along decreasing t
  }
  else {
    double gam = 1/TMath::Sqrt( fParams[kB0]*fParams[kB0] + fParams[kB1]*fParams[kB1] + 1.);
    dircos[ fkAxID[kX] ] = fParams[kB0]*gam;
    dircos[ fkAxID[kY] ] = fParams[kB1]*gam;
    dircos[ fParAxis   ] = gam;
    // decide direction
    if (IsTypeCollision()) {
      static double xyzF[3],xyzL[3];
      GetPosition(xyzF,fPntFirst);
      GetPosition(xyzL,fPntLast);
      double dif = fCurT[fPntLast] - fCurT[fPntFirst];
      double dr = (xyzL[kX]-xyzF[kX])*(xyzL[kX]+xyzF[kX]) + (xyzL[kY]-xyzF[kY])*(xyzL[kY]+xyzF[kY]);
      if (dr*dif<0) for (int i=3;i--;) dircos[i] = -dircos[i]; // with increasing t the tracks comes closer to origin
    }
    else if (dircos[kY]>0) for (int i=3;i--;) dircos[i] = -dircos[i];  // cosmic tracks have negative angle to Y axis
  }
  //
}

//________________________________________________________________________________________________________
Double_t AliITSTPArrayFit::GetMachinePrec()
{
  // estimate machine precision
  Double_t eps=1.0,a;
  do { a = 1. + (eps=eps/2.0); } while(a>1.);
  return TMath::Abs(2.*eps);
}

//________________________________________________________________________________________________________
Bool_t AliITSTPArrayFit::FitHelixCrude(Int_t extQ)
{
  // crude estimate of helix parameters, w/o errors and Eloss.
  // Fast Riemann fit: Comp.Phy.Comm.131 (2000) 95
  //
  // if charge is not imposed (extQ==0) then it will be determined from the collision type
  //
  static TArrayD arrU,arrV,arrW;
  double *parrW,*parrU,*parrV;
  Bool_t eloss = IsELossON();
  //
  int np = fPntLast - fPntFirst + 1;
  if (np<3) { AliError("At least 3 points are needed for helix fit"); return kFALSE; }
  //
  const float *x=fkPoints->GetX(),*y=fkPoints->GetY(),*z=fkPoints->GetZ(),*cov=fkPoints->GetCov();
  //
  if (fPntLast>arrU.GetSize()) {
    arrU.Set(2*fPntLast);
    arrV.Set(2*fPntLast);
    arrW.Set(2*fPntLast);
  }
  parrU = arrU.GetArray();
  parrV = arrV.GetArray();
  parrW = arrW.GetArray();
  //
  double uav=0,vav=0,wav=0,muu=0,muv=0,muw=0,mvv=0,mvw=0,mww=0;
  int minRId = fPntFirst;
  //  
  // get points span
  double xmn=1e9,xmx=-1e9, ymn=1e9,ymx=-1e9;
  for (int i=fPntFirst;i<=fPntLast;i++) {
    parrW[i] = x[i]*x[i]+y[i]*y[i];
    if (parrW[i]<parrW[minRId]) minRId = i; // point closest to origin
    if (xmn>x[i]) xmn = x[i];
    if (xmx<x[i]) xmx = x[i];
    if (ymn>y[i]) ymn = y[i];
    if (ymx<y[i]) ymx = y[i];
  }
  int minRId1 = minRId>fPntFirst ? fPntFirst:fPntFirst+1;
  for (int i=fPntFirst;i<=fPntLast;i++) if (parrW[i]<parrW[minRId1] && i!=minRId) minRId1 = i; 
  //
  double xshift = (xmx+xmn)/2 + 10*(ymx-ymn); // shift origin to have uniform weights
  double yshift = (ymx+ymn)/2 - 10*(xmx-xmn);
  //  printf("X: %+e %+e Y: %+e %+e | shift: %+e %+e\n",xmn,xmx,ymn,ymx,xshift,yshift);
  //
  for (int i=fPntFirst;i<=fPntLast;i++) {
    double xs = x[i] - xshift;
    double ys = y[i] - yshift;
    double w = xs*xs + ys*ys;
    double scl = 1./(1.+w);
    int i0 = i-fPntFirst;
    wav += parrW[i0] = w*scl;
    uav += parrU[i0] = xs*scl;
    vav += parrV[i0] = ys*scl;
  }
  uav /= np;    vav /= np;   wav /= np;
  //
  for (int i=fPntFirst;i<=fPntLast;i++) {
    //
    // point next to closest
    int i0 = i-fPntFirst;
    if (parrW[i0]<parrW[minRId1-fPntFirst] && i!=minRId) minRId1 = i; 
    double u = parrU[i0] - uav;
    double v = parrV[i0] - vav;
    double w = parrW[i0] - wav;
    muu += u*u;
    muv += u*v;
    muw += u*w;
    mvv += v*v;
    mvw += v*w;
    mww += w*w;
  } 
  muu/=np; muv/=np; muw/=np; mvv/=np; mvw/=np; mww/=np;
  //
  // find eigenvalues:
  double trace3 = (muu + mvv + mww)/3.;
  double muut = muu-trace3;
  double mvvt = mvv-trace3;
  double mwwt = mww-trace3;
  double q = (muut*(mvvt*mwwt-mvw*mvw) - muv*(muv*mwwt-mvw*muw) + muw*(muv*mvw-mvvt*muw))/2;
  double p = (muut*muut+mvvt*mvvt+mwwt*mwwt+2*(muv*muv+muw*muw+mvw*mvw))/6;
  double dfpp = p*p*p-q*q;
  dfpp = dfpp>0 ? TMath::Sqrt(dfpp)/q : 0;
  double ph = TMath::ATan( dfpp )/3.;
  if (ph<0) ph += TMath::Pi()/3;
  p = p>0 ? TMath::Sqrt(p) : 0;
  const double kSqrt3 = 1.73205080;
  double snp = TMath::Sin(ph);
  double csp = TMath::Cos(ph);
  //  double eg1 = trace3 + 2*p*csp;
  double eg2 = trace3 - p*(csp+kSqrt3*snp); // smallest one
  //  double eg3 = trace3 - p*(csp-kSqrt3*snp);
  // eigenvector for min.eigenvalue
  muut = muu-eg2;
  mvvt = mvv-eg2;
  mwwt = mww-eg2;
  double n0 = muv*mvw-muw*mvvt;
  double n1 = muv*muw-mvw*muut;
  double n2 = muut*mvvt-muv*muv;
  // normalize to largest one
  double nrm = TMath::Abs(n0);
  if (nrm<TMath::Abs(n1)) nrm = TMath::Abs(n1);
  if (nrm<TMath::Abs(n2)) nrm = TMath::Abs(n2);
  n0/=nrm; n1/=nrm; n2/=nrm;
  //
  double cpar = -(uav*n0 + vav*n1 + wav*n2);
  double xc = -n0/(cpar+n2)/2 + xshift;
  double yc = -n1/(cpar+n2)/2 + yshift;
  double rad = TMath::Sqrt(n0*n0+n1*n1-4*cpar*(cpar+n2))/2./TMath::Abs(cpar+n2);
  //
  //  printf("Rad: %+e xc: %+e yc: %+e | X0: %+e Y0: %+e | X1: %+e Y1: %+e\n",rad,xc,yc, x[minRId],y[minRId],x[minRId1],y[minRId1]);

  // linear circle fit --------------------------------------------------- <<<
  //
  // decide sign(Q*B) and fill cicrle parameters ------------------------- >>>
  int sqb;
  if (extQ) {
    SetCharge(extQ); 
    sqb = fBz<0 ? -GetCharge():GetCharge();
  }
  else { 
    // determine the charge from the collision type and field sign
    // the negative Q*B will have positive Vc x dir product Z component
    // with Vc={-xc,-yc} : vector from circle center to the origin
    // and V0 - track direction vector (take {0,-1,1} for cosmics)
    // If Bz is not provided, assume positive Bz
    if ( IsTypeCosmics() ) sqb = xc>0 ? -1:1;
    else {
      // track direction vector as a - diference between the closest and the next to closest to origin points
      double v0X = x[minRId1] - x[minRId];
      double v0Y = y[minRId1] - y[minRId];
      sqb = (yc*v0X - xc*v0Y)<0 ? -1:1;
    }
    SetCharge( fBz<0 ? -sqb : sqb);
  }
  //
  Double_t phi = TMath::ATan2(yc,xc);
  if (sqb<0) phi += TMath::Pi();
  if      (phi > TMath::Pi()) phi -= 2.*TMath::Pi();
  else if (phi <-TMath::Pi()) phi += 2.*TMath::Pi();
  fParams[kPhi0] = phi;  
  fParams[kR0]   = sqb<0 ? -rad:rad;  
  fParams[kD0] = xc*TMath::Cos(phi) + yc*TMath::Sin(phi) - fParams[kR0];
  //
  // decide sign(Q*B) and fill cicrle parameters ------------------------- <<<
  //
  // find z-offset and dip + the parameter t of closest approach to hits - >>>
  //
  UInt_t hitLrPos=0;  // pattern of hit layers at pos
  UInt_t hitLrNeg=0;  // and negative t's

  Double_t ss=0,st=0,sz=0,stt=0,szt=0;
  for (int i=fPntFirst;i<=fPntLast;i++) {
    //
    Double_t ze2 = cov[i*6 + kZZ];
    Double_t t = TMath::ATan2(yc-y[i],xc-x[i]) - fParams[kPhi0]; // angle at measured z
    if (fParams[kR0]<0)  t += TMath::Pi();
    if      (t > TMath::Pi()) t -= TMath::Pi()*2;
    else if (t <-TMath::Pi()) t += TMath::Pi()*2;
    if (ze2<fgkAlmostZero) ze2 = 1E-8;
    ze2 = 1./ze2;
    ss += ze2;
    st += t*ze2;
    stt+= t*t*ze2;
    sz += z[i]*ze2;
    szt+= z[i]*t*ze2;
    //
    fCurT[i] = t; // parameter of the closest approach to the point
    //    printf("%d %+e %+e %+e %+e\n",i,x[i],y[i],z[i],t);
    if (eloss) {
      double r = TMath::Sqrt(x[i]*x[i]+y[i]*y[i]);
      int lr;
      for (lr=kMaxLrITS;lr--;) if ( IsZero(r-fgkRLayITS[ lr ],1.) ) break;
      if (lr<kMaxLrITS) {
	if (t>0) hitLrPos |= (1<<lr);  // set bit of the layer
	else     hitLrNeg |= (1<<lr);  // set bit of the layer
      }
    }
  }
  double det = ss*stt - st*st;
  if (TMath::Abs(det)<fgkAlmostZero) { // no Z dependence
    fParams[kDZ]  = sz/ss;
    fParams[kDip] = 0;
  }
  else {
    fParams[kDZ]  =  (sz*stt-st*szt)/det;
    fParams[kDip] = -(ss*szt-st*sz)/det/fParams[kR0];
  }
  //
  // find z-offset and dip + the parameter t of closest approach to hits - <<<
  //
  // fill info needed to account for ELoss ------------------------------- >>>
  if (eloss) {
    fNElsPnt = fPntLast - fPntFirst + 1;
    //
    // to account for the energy loss in the passive volumes, calculate the relevant t-parameters 
    double* tcur = fCurT + fPntFirst;
    double* ecur = fElsDR+ fPntFirst;
    //
    for (int ilp=3;ilp--;) {
      int id = fgkPassivLrITS[ilp];
      double tp = GetHelixParAtR( fgkRLayITS[ id ] );
      if (tp<0) continue; // does not hit this radius
      //
      tcur[fNElsPnt] = GetSignQB()>0 ? -tp : tp;
      ecur[fNElsPnt] = fgRhoLITS[ id ];
      fNElsPnt++;
      //      printf("Passive  on lr %d  %+e\n",ilp,tcur[fNElsPnt-1]);
      //
      if (IsTypeCosmics() && !IsZero(tp)) { // 2 crossings for cosmics
	tcur[fNElsPnt] = -tcur[fNElsPnt-1];
	ecur[fNElsPnt] =  ecur[fNElsPnt-1];
	fNElsPnt++;
	//printf("Passive* on lr %d  %+e\n",ilp,-tcur[fNElsPnt-1]);
      }
      //
    }
    // check if some active layers did not miss the hit, treat them as passive
    for (int ilp=6;ilp--;) {
      int id = fgkActiveLrITS[ilp];
      double tp = GetHelixParAtR( fgkRLayITS[ id ] );
      if (tp<0) continue; // does not hit this radius
      //
      if ( (GetSignQB()>0||IsTypeCosmics()) && !(hitLrNeg & (1<<id)) ) {
	tcur[fNElsPnt] = -tp;
	ecur[fNElsPnt] = fgRhoLITS[ id ];
	fNElsPnt++;
	//printf("Missed  on lr %d  %+e\n",ilp,-tp);
      }
      //
      if ( (GetSignQB()<0||IsTypeCosmics()) && !(hitLrPos & (1<<id)) ) {
	tcur[fNElsPnt] = tp;
	ecur[fNElsPnt] = fgRhoLITS[ id ];
	fNElsPnt++;
	//printf("Missed* on lr %d  %e\n",ilp,tp);
      }
    }
    //
    TMath::Sort(fNElsPnt,fCurT+fPntFirst,fElsId,kFALSE);    // index e-loss points in increasing order
    // find the position of smallest positive t-param
    for (fFirstPosT=0;fFirstPosT<fNElsPnt;fFirstPosT++) if (fCurT[ fElsId[ fFirstPosT ] ]>0) break;
    //
    Double_t cdip = 1./TMath::Sqrt(1.+fParams[kDip]*fParams[kDip]);
    Double_t ptot = TMath::Abs(fParams[kR0]*fgkCQConv*fBz/cdip); // momentum and energy
    Double_t etot = TMath::Sqrt(ptot*ptot + fMass*fMass);      // in the point of closest approach to beam
    Double_t normS[3];
    //
    // Positive t-params: along the track direction for negative track, against for positive
    Double_t pcur = ptot, ecurr = etot;
    for (int ip=fFirstPosT;ip<fNElsPnt;ip++) {
      int tID = fElsId[ip];
      Double_t t = fCurT[ tID ];
      //
      if (tID>fPntLast) { // this is not a hit layer but passive layer
	double php = TMath::ATan2(yc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t),
				  xc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t));
	normS[0] = -TMath::Cos(php);  // normal to the cylinder at intersection point
	normS[1] = -TMath::Sin(php);
	normS[2] = 0;
      }
      else GetNormal(normS,fkPoints->GetCov()+tID*6);   // vector normal to hit module
      fElsDR[tID] = GetDRofELoss(t,cdip,fElsDR[tID],normS,ptot,etot);
    }
    //
    // negaive t-params: against the track direction for negative track, along for positive
    pcur  = ptot;
    ecurr = etot;
    for (int ip=fFirstPosT;ip--;) {
      int tID = fElsId[ip];
      Double_t t = fCurT[ tID ];
      //
      if (tID>fPntLast) { // this is not a hit layer but passive layer
	double php = TMath::ATan2(yc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t),
				  xc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t));
	normS[0] = -TMath::Cos(php);  // normal to the cylinder at intersection point
	normS[1] = -TMath::Sin(php);
	normS[2] = 0;	
      }
      else GetNormal(normS,fkPoints->GetCov()+tID*6);   // vector normal to hit module
      //
      fElsDR[tID] = GetDRofELoss(t,cdip,fElsDR[tID],normS,ptot,etot);
    }
  }
  // fill info needed to account for ELoss ------------------------------- <<<
  //
  return kTRUE;
}

/*
//________________________________________________________________________________________________________
Bool_t AliITSTPArrayFit::FitHelixCrude(Int_t extQ)
{
  // crude estimate of helix parameters, w/o errors and Eloss.
  // 1st fit the circle (R,xc,yc) by minimizing
  // chi2 = sum{ (bx*xi + by*yi + xi^2+yi^2 + rho)^2 } vs bx,by,rho
  // with bx = -2*xc, by = -2*yc , rho = xc^2+yc^2 - R2
  //
  // if charge is not imposed (extQ==0) then it will be determined from the collision type
  //
  Bool_t eloss = IsELossON();
  //
  int np = fPntLast - fPntFirst + 1;
  if (np<3) { AliError("At least 3 points are needed for helix fit"); return kFALSE; }
  //
  const float *x=fkPoints->GetX(),*y=fkPoints->GetY(),*z=fkPoints->GetZ(),*cov=fkPoints->GetCov();
  //
  // linear circle fit --------------------------------------------------- >>>
  Double_t sxx=0,sxy=0,syy=0,sx=0,sy=0,rhs0=0,rhs1=0,rhs2=0,minR=1E9;
  int minRId = 0;
  for (int i=fPntFirst;i<=fPntLast;i++) {
    Double_t xx = x[i]*x[i];
    Double_t yy = y[i]*y[i];
    Double_t xy = x[i]*y[i];
    Double_t xxyy = xx + yy;
    //
    sxx += xx;
    sxy += xy;
    syy += yy;
    sx  += x[i];
    sy  += y[i];
    //
    rhs0 -= xxyy*x[i];
    rhs1 -= xxyy*y[i];
    rhs2 -= xxyy;
    // 
    // remember Id of the point closest to origin, to determine the charge  
    if (xxyy<minR) { minR   = xxyy; minRId = i; }
    //
    if (eloss) { // find layer id
      int lrid,volid = fkPoints->GetVolumeID()[i];
      if (volid>0) lrid = fgkActiveLrITS[AliGeomManager::VolUIDToLayer(fkPoints->GetVolumeID()[i])-1];
      else { // missing layer info, find from radius
	double r = TMath::Sqrt(xxyy);
	for (lrid=kMaxLrITS;lrid--;) if ( IsZero(r-fgkRLayITS[ lrid ],1.) ) break;
      }
      fElsDR[i] = (lrid>=0 && lrid<kMaxLrITS) ? fgRhoLITS[ lrid ] : 0;   // eloss for normal track
    }
    //
  }
  //
  Double_t mn00 = syy*np-sy*sy;
  Double_t mn01 = sxy*np-sy*sx;
  Double_t mn02 = sxy*sy-syy*sx;
  Double_t det  = sxx*mn00 - sxy*mn01 + sx*mn02; 
  if (TMath::Abs(det)<fgkAlmostZero) return kFALSE;
  //
  Double_t mn11 = sxx*np-sx*sx;
  Double_t mn12 = sxx*sy-sxy*sx;
  Double_t mn22 = sxx*syy-sxy*sxy;
  //
  Double_t mi00 =  mn00/det;
  Double_t mi01 = -mn01/det;
  Double_t mi02 =  mn02/det;
  Double_t mi11 =  mn11/det;
  Double_t mi12 = -mn12/det;
  Double_t mi22 =  mn22/det;
  //
  Double_t xc   = -(rhs0*mi00 + rhs1*mi01 + rhs2*mi02)/2;
  Double_t yc   = -(rhs0*mi01 + rhs1*mi11 + rhs2*mi12)/2;
  Double_t rho2 =  (rhs0*mi02 + rhs1*mi12 + rhs2*mi22);

  //
  // check
  double bx = -2*xc;
  double by = -2*yc;
  double sm0=0,sm1=0;
  for (int i=fPntFirst;i<=fPntLast;i++) {
    double dst = bx*x[i]+by*y[i]+x[i]*x[i]+y[i]*y[i]+rho2;
    sm0 += dst;
    sm1 += dst*dst;
    printf("Point %d: %+e %+e %+e\n",i,dst,sm0,sm1);
  }

  //
  Double_t rad = xc*xc + yc*yc - rho2;
  rad = (rad>fgkAlmostZero) ? (TMath::Sqrt(rad)):fgkAlmostZero;
  //
  //  printf("Rad: %+e xc: %+e yc: %+e\n",rad,xc,yc);

  // linear circle fit --------------------------------------------------- <<<
  //
  // decide sign(Q*B) and fill cicrle parameters ------------------------- >>>
  int sqb;
  if (extQ) {
    SetCharge(extQ); 
    sqb = fBz<0 ? -GetCharge():GetCharge();
  }
  else { 
    // determine the charge from the collision type and field sign
    // the negative Q*B will have positive Vc x V0 product Z component
    // with Vc={-xc,-yc} : vector from circle center to the origin
    // and V0 - track direction vector (take {0,-1,1} for cosmics)
    // If Bz is not provided, assume positive Bz
    sqb = ( IsTypeCosmics() ? xc:(yc*x[minRId]-xc*y[minRId]) ) > 0 ? -1:1;
    SetCharge( fBz<0 ? -sqb : sqb);
  }
  //
  Double_t phi = TMath::ATan2(yc,xc);
  if (sqb<0) phi += TMath::Pi();
  if      (phi > TMath::Pi()) phi -= 2.*TMath::Pi();
  else if (phi <-TMath::Pi()) phi += 2.*TMath::Pi();
  fParams[kPhi0] = phi;  
  fParams[kR0]   = sqb<0 ? -rad:rad;  
  fParams[kD0] = xc*TMath::Cos(phi) + yc*TMath::Sin(phi) - fParams[kR0];
  //
  // decide sign(Q*B) and fill cicrle parameters ------------------------- <<<
  //
  // find z-offset and dip + the parameter t of closest approach to hits - >>>
  //
  UInt_t hitLrPos=0;  // pattern of hit layers at pos
  UInt_t hitLrNeg=0;  // and negative t's

  Double_t ss=0,st=0,sz=0,stt=0,szt=0;
  for (int i=fPntFirst;i<=fPntLast;i++) {
    //
    Double_t ze2 = cov[i*6 + kZZ];
    Double_t t = TMath::ATan2(yc-y[i],xc-x[i]) - fParams[kPhi0]; // angle at measured z
    if (fParams[kR0]<0)  t += TMath::Pi();
    if      (t > TMath::Pi()) t -= TMath::Pi()*2;
    else if (t <-TMath::Pi()) t += TMath::Pi()*2;
    if (ze2<fgkAlmostZero) ze2 = 1E-8;
    ze2 = 1./ze2;
    ss += ze2;
    st += t*ze2;
    stt+= t*t*ze2;
    sz += z[i]*ze2;
    szt+= z[i]*t*ze2;
    //
    fCurT[i] = t; // parameter of the closest approach to the point
    //    printf("%d %+e %+e %+e %+e\n",i,x[i],y[i],z[i],t);
    if (eloss) {
      double r = TMath::Sqrt(x[i]*x[i]+y[i]*y[i]);
      int lr;
      for (lr=kMaxLrITS;lr--;) if ( IsZero(r-fgkRLayITS[ lr ],1.) ) break;
      if (lr<kMaxLrITS) {
	if (t>0) hitLrPos |= (1<<lr);  // set bit of the layer
	else     hitLrNeg |= (1<<lr);  // set bit of the layer
      }
    }
  }
  det = ss*stt - st*st;
  if (TMath::Abs(det)<fgkAlmostZero) { // no Z dependence
    fParams[kDZ]  = sz/ss;
    fParams[kDip] = 0;
  }
  else {
    fParams[kDZ]  =  (sz*stt-st*szt)/det;
    fParams[kDip] = -(ss*szt-st*sz)/det/fParams[kR0];
  }
  //
  // find z-offset and dip + the parameter t of closest approach to hits - <<<
  //
  // fill info needed to account for ELoss ------------------------------- >>>
  if (eloss) {
    fNElsPnt = fPntLast - fPntFirst + 1;
    //
    // to account for the energy loss in the passive volumes, calculate the relevant t-parameters 
    double* tcur = fCurT + fPntFirst;
    double* ecur = fElsDR+ fPntFirst;
    //
    for (int ilp=3;ilp--;) {
      int id = fgkPassivLrITS[ilp];
      double tp = GetHelixParAtR( fgkRLayITS[ id ] );
      if (tp<0) continue; // does not hit this radius
      //
      tcur[fNElsPnt] = GetSignQB()>0 ? -tp : tp;
      ecur[fNElsPnt] = fgRhoLITS[ id ];
      fNElsPnt++;
      //      printf("Passive  on lr %d  %+e\n",ilp,tcur[fNElsPnt-1]);
      //
      if (IsTypeCosmics() && !IsZero(tp)) { // 2 crossings for cosmics
	tcur[fNElsPnt] = -tcur[fNElsPnt-1];
	ecur[fNElsPnt] =  ecur[fNElsPnt-1];
	fNElsPnt++;
	//printf("Passive* on lr %d  %+e\n",ilp,-tcur[fNElsPnt-1]);
      }
      //
    }
    // check if some active layers did not miss the hit, treat them as passive
    for (int ilp=6;ilp--;) {
      int id = fgkActiveLrITS[ilp];
      double tp = GetHelixParAtR( fgkRLayITS[ id ] );
      if (tp<0) continue; // does not hit this radius
      //
      if ( (GetSignQB()>0||IsTypeCosmics()) && !(hitLrNeg & (1<<id)) ) {
	tcur[fNElsPnt] = -tp;
	ecur[fNElsPnt] = fgRhoLITS[ id ];
	fNElsPnt++;
	//printf("Missed  on lr %d  %+e\n",ilp,-tp);
      }
      //
      if ( (GetSignQB()<0||IsTypeCosmics()) && !(hitLrPos & (1<<id)) ) {
	tcur[fNElsPnt] = tp;
	ecur[fNElsPnt] = fgRhoLITS[ id ];
	fNElsPnt++;
	//printf("Missed* on lr %d  %e\n",ilp,tp);
      }
    }
    //
    TMath::Sort(fNElsPnt,fCurT+fPntFirst,fElsId,kFALSE);    // index e-loss points in increasing order
    // find the position of smallest positive t-param
    for (fFirstPosT=0;fFirstPosT<fNElsPnt;fFirstPosT++) if (fCurT[ fElsId[ fFirstPosT ] ]>0) break;
    //
    Double_t cdip = 1./TMath::Sqrt(1.+fParams[kDip]*fParams[kDip]);
    Double_t ptot = TMath::Abs(fParams[kR0]*fgkCQConv*fBz/cdip); // momentum and energy
    Double_t etot = TMath::Sqrt(ptot*ptot + fMass*fMass);      // in the point of closest approach to beam
    Double_t normS[3];
    //
    // Positive t-params: along the track direction for negative track, against for positive
    Double_t pcur = ptot, ecurr = etot;
    for (int ip=fFirstPosT;ip<fNElsPnt;ip++) {
      int tID = fElsId[ip];
      Double_t t = fCurT[ tID ];
      //
      if (tID>fPntLast) { // this is not a hit layer but passive layer
	double php = TMath::ATan2(yc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t),
				  xc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t));
	normS[0] = -TMath::Cos(php);  // normal to the cylinder at intersection point
	normS[1] = -TMath::Sin(php);
	normS[2] = 0;
      }
      else GetNormal(normS,fkPoints->GetCov()+tID*6);   // vector normal to hit module
      fElsDR[tID] = GetDRofELoss(t,cdip,fElsDR[tID],normS,ptot,etot);
    }
    //
    // negaive t-params: against the track direction for negative track, along for positive
    pcur  = ptot;
    ecurr = etot;
    for (int ip=fFirstPosT;ip--;) {
      int tID = fElsId[ip];
      Double_t t = fCurT[ tID ];
      //
      if (tID>fPntLast) { // this is not a hit layer but passive layer
	double php = TMath::ATan2(yc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t),
				  xc-fParams[kR0]*TMath::Cos(fParams[kPhi0]+t));
	normS[0] = -TMath::Cos(php);  // normal to the cylinder at intersection point
	normS[1] = -TMath::Sin(php);
	normS[2] = 0;	
      }
      else GetNormal(normS,fkPoints->GetCov()+tID*6);   // vector normal to hit module
      //
      fElsDR[tID] = GetDRofELoss(t,cdip,fElsDR[tID],normS,ptot,etot);
    }
  }
  // fill info needed to account for ELoss ------------------------------- <<<
  //
  return kTRUE;
}
*/
//____________________________________________________
Double_t AliITSTPArrayFit::FitHelix(Int_t extQ, Double_t extPT,Double_t extPTerr)
{
  // fit by helix accounting for the errors of all coordinates (and energy loss if requested)
  // 
  // If extQ is non-0, its sign is imposed as a charge of the track
  // If extPT>0 and extPTerr>=0, constrain to measured tr.momentum PT 
  // with corresponding error (err=0 -> rel.err=1e-6)
  //
  double chiprev=1e99;
  //const Double_t kMaxTEffect = 1E-6;
  Double_t dXYZdGlo[3*5],dXYZdLoc[3],xyzRes[3];
  //
  SetFitDone(kFALSE);
  fChi2NDF = -1;
  //
  if (!FitHelixCrude(extQ)) return -1; // get initial estimate, ignoring the errors
  //
  if (TMath::Abs(fParams[kR0])>fMaxRforHelix && extPT<=0) {
    fSwitch2Line = kTRUE;
    return FitLine();
  }
  //
  //
  if (!IsCovInv()) InvertPointsCovMat();    // prepare inverted errors
  if (!fParSol) fParSol = new AliParamSolver(5);
  fParSol->SetNGlobal(5);
  //
  //  printf("-1 | %+.2e %+.2e %+.2e %+.2e %+.2e | chi2: %+.4e\n",fParams[0],fParams[1],fParams[2],fParams[3],fParams[4],CalcChi2NDF());
  int iter = 0;
  fChi2NDF = 1e99;
  Bool_t converged = kFALSE;
  while(iter<fMaxIter) {
    chiprev = fChi2NDF;
    fParSol->Clear();
    for (int ip=fPntFirst;ip<=fPntLast;ip++) {
      //
      GetResiduals(xyzRes, ip); // current residuals at point ip
      Double_t rrho = fParams[kR0]+fParams[kD0];
      Double_t cs0  = TMath::Cos(fParams[kPhi0]);
      Double_t sn0  = TMath::Sin(fParams[kPhi0]);
      Double_t cst  = TMath::Cos(fCurT[ip]+fParams[kPhi0]);
      Double_t snt  = TMath::Sin(fCurT[ip]+fParams[kPhi0]);
      //
      int offs = kD0;                  // dXYZ/dD0
      dXYZdGlo[offs + kX] = cs0;
      dXYZdGlo[offs + kY] = sn0;
      dXYZdGlo[offs + kZ] = 0;
      //
      offs = kPhi0*3;                  // dXYZ/dPhi0
      dXYZdGlo[offs + kX] = -rrho*sn0 + fParams[kR0]*snt;
      dXYZdGlo[offs + kY] =  rrho*cs0 - fParams[kR0]*cst;
      dXYZdGlo[offs + kZ] = 0;
      //
      offs = kR0*3;                   // dXYZ/dR0
      if (TMath::Abs(fParams[kR0])<0.9*fMaxRforHelix) {
	dXYZdGlo[offs + kX] =  cs0 - cst;
	dXYZdGlo[offs + kY] =  sn0 - snt;
	dXYZdGlo[offs + kZ] =  -fParams[kDip]*fCurT[ip];
      }
      else {
	dXYZdGlo[offs + kX] = dXYZdGlo[offs + kY] = dXYZdGlo[offs + kZ] = 0;
	fParSol->AddConstraint(kR0,0,1.e2);
      }
      //
      offs = kDZ*3;                   // dXYZ/dDZ
      dXYZdGlo[offs + kX] =  0;
      dXYZdGlo[offs + kY] =  0;
      dXYZdGlo[offs + kZ] =  1.;
      //
      offs = kDip*3;                  // dXYZ/dDip
      dXYZdGlo[offs + kX] =  0;
      dXYZdGlo[offs + kY] =  0;
      dXYZdGlo[offs + kZ] = -fParams[kR0]*fCurT[ip];
      //
      //      /*
      dXYZdLoc[kX] =  fParams[kR0]*snt;
      dXYZdLoc[kY] = -fParams[kR0]*cst;
      dXYZdLoc[kZ] = -fParams[kR0]*fParams[kDip];
      //      */
      //      dXYZdLoc[0] = dXYZdLoc[1] = dXYZdLoc[2] = 0;
      //
      fParSol->AddEquation(dXYZdGlo,dXYZdLoc,xyzRes,GetCovI(ip));
    }
    //
    if (extPT>0) { // add constraint on pt
      if (extPTerr<fgkAlmostZero) extPTerr = 1e-6*extPT;
      Double_t cf = fBz*GetCharge()*fgkCQConv;
      Double_t err2i = extPTerr/cf;
      err2i = 1./err2i/err2i;
      //      printf("Constrain R to %+e\n",extPT/cf);
      fParSol->AddConstraint(kR0,-extPT/cf+fParams[kR0],err2i);
    }
    if (!fParSol->Solve()) { AliError("Failed to fit helix"); return -1; }
    Double_t *deltaG = fParSol->GetGlobals();
    //    Double_t *deltaT = fParSol->GetLocals();
    for (int ipar=5;ipar--;) fParams[ipar] -= deltaG[ipar];
    //
    if (TMath::Abs(fParams[kR0])>0.9*fMaxRforHelix) fParams[kR0] = TMath::Sign(0.9*fMaxRforHelix,fParams[kR0]);
    //
    for (int ip=fPntFirst;ip<=fPntLast;ip++) {
      fCurT[ip] = CalcParPCA(ip);
      //      printf("iter%d : delta%2d : expl: %+e | %+e | R=%+e p0=%+e\n",iter,ip,deltaT[ip-fPntFirst], fCurT[ip],fParams[kR0],fParams[kPhi0]);
      //      fCurT[ip] -= deltaT[ip-fPntFirst];
    }
    iter++;
    //
    fChi2NDF = CalcChi2NDF();
    //        printf("%2d | %+.2e %+.2e %+.2e %+.2e %+.2e | chi2: %+.4e %+.4e\n",iter,deltaG[0],deltaG[1],deltaG[2],deltaG[3],deltaG[4],fChi2NDF,fChi2NDF-chiprev);
    //        printf("->>  %+.2e %+.2e %+.2e %+.2e %+.2e | Chi2: %+.6e %+.6e\n",fParams[0],fParams[1],fParams[2],fParams[3],fParams[4],fChi2NDF,fChi2NDF-chiprev);
    double difchi2 = chiprev - fChi2NDF;
    if ( difchi2<fEps && TMath::Abs(difchi2)<1e-4) {converged = kTRUE; break;} 
    //    if (errT*TMath::Abs(fParams[kR0])<kMaxTEffect && errP<fEps) {converged = kTRUE; break;} 
  }
  //
  if (!converged) {
    AliDebug(2,Form("Max number of %d iteration reached, Current chi2:%.3e, chi2 change %+.3e",iter,
		    fChi2NDF,chiprev-fChi2NDF));
    for (int ip=fPntFirst;ip<=fPntLast;ip++)
      AliDebug(2,Form("P%2d| %+.3e %+.3e %+.3e\n",ip,fkPoints->GetX()[ip],fkPoints->GetY()[ip],fkPoints->GetZ()[ip]));

  }
  fIter = iter;
  SetCharge( fParams[kR0]>0 ? (fBz<0?-1:1):(fBz>0?-1:1) );
  SetFitDone();
  //  printf("F1>> %+.7e %+.7e %+.7e %+.7e %.7e\n",fParams[0],fParams[1],fParams[2],fParams[3],fParams[4]);
  //
  return fChi2NDF;
}

//____________________________________________________
Double_t AliITSTPArrayFit::FitLine()
{
  // fit by helix accounting for the errors of all coordinates (and energy loss if requested)
  // 
  double chiprev=1e99;
  //  const Double_t kMaxTEffect = 1.e-6;
  Double_t dXYZdGlo[3*4],dXYZdLoc[3],xyzRes[3];
  SetFitDone(kFALSE);
  fChi2NDF = -1;
  //
  if (fParAxis<0) SetParAxis(ChoseParAxis());
  //
  const float *xyzp[3]={fkPoints->GetX(),fkPoints->GetY(),fkPoints->GetZ()};
  if (!IsCovInv()) InvertPointsCovMat();
  if (!FitLineCrude()) return -1; // get initial estimate, ignoring the errors
  //
  if (!fParSol) fParSol = new AliParamSolver(5);
  fParSol->SetNGlobal(4);
  // initial set of parameters
  for (int ip=fPntFirst;ip<=fPntLast;ip++) fCurT[ip] = xyzp[fParAxis][ip]; // use measured param-coordinate
  //
  int iter = 0;
  Bool_t converged = kFALSE;
  fChi2NDF = 1e99;
  while(iter<fMaxIter) {
    chiprev = fChi2NDF;
    fParSol->Clear();
    for (int ip=fPntFirst;ip<=fPntLast;ip++) {
      //
      int offs;
      GetResiduals(xyzRes, ip); // current residuals at point ip
      //
      offs = kA0*3;                   // dXYZ/dA0
      dXYZdGlo[offs + fkAxID[kX]] = 1;
      dXYZdGlo[offs + fkAxID[kY]] = 0;
      dXYZdGlo[offs + fParAxis  ] = 0;
      //
      offs = kB0*3;                   // dXYZ/dB0
      dXYZdGlo[offs + fkAxID[kX]] = fCurT[ip];
      dXYZdGlo[offs + fkAxID[kY]] = 0;
      dXYZdGlo[offs + fParAxis  ] = 0;
      //
      offs = kA1*3;                   // dXYZ/dA1
      dXYZdGlo[offs + fkAxID[kX]] = 0;
      dXYZdGlo[offs + fkAxID[kY]] = 1;
      dXYZdGlo[offs + fParAxis  ] = 0;
      //
      offs = kB1*3;                   // dXYZ/dB1
      dXYZdGlo[offs + fkAxID[kX]] = 0;
      dXYZdGlo[offs + fkAxID[kY]] = fCurT[ip];
      dXYZdGlo[offs + fParAxis  ] = 0;
      //
      dXYZdLoc[ fkAxID[kX] ] =  fParams[kB0];  // dX/dt
      dXYZdLoc[ fkAxID[kY] ] =  fParams[kB1];  // dY/dt
      dXYZdLoc[ fParAxis   ] =  1;
      //
      fParSol->AddEquation(dXYZdGlo,dXYZdLoc,xyzRes,GetCovI(ip));
    }
    //
    if (!fParSol->Solve()) { AliError("Failed to fit line"); return -1; }
    Double_t *deltaG = fParSol->GetGlobals();
    Double_t *deltaT = fParSol->GetLocals();
    for (int ipar=4;ipar--;) fParams[ipar] -= deltaG[ipar];
    for (int ip=fPntFirst;ip<=fPntLast;ip++) fCurT[ip] -= deltaT[ip-fPntFirst];
    iter++;
    fChi2NDF = CalcChi2NDF();
    //    printf("%d %+e %+e | %+.2e %+.2e %+.2e %+.2e | chi2: %+.4e %+.4e\n",iter,errP,errT, deltaG[0],deltaG[1],deltaG[2],deltaG[3],fChi2NDF,fChi2NDF-chiprev);
    //    printf("->> %+.2e %+.2e %+.2e %+.2e %+.2e | Chi2: %+.6e %+.6e\n",fParams[0],fParams[1],fParams[2],fParams[3],fParams[4],fChi2NDF,fChi2NDF-chiprev);
    double difchi2 = chiprev - fChi2NDF;
    if ( difchi2<fEps && TMath::Abs(difchi2)<1e-4) {converged = kTRUE; break;} 
    chiprev = fChi2NDF;
    //    if (errT<kMaxTEffect && errP<fEps) {converged = kTRUE; break;} 
  }
  //
  if (!converged) {
    AliDebug(2,Form("Max number of %d iteration reached, Current chi2:%.3e, chi2 change %+.3e",iter,
		    fChi2NDF,chiprev-fChi2NDF));
    for (int ip=fPntFirst;ip<=fPntLast;ip++)
      AliDebug(2,Form("P%2d| %+.3e %+.3e %+.3e\n",ip,fkPoints->GetX()[ip],fkPoints->GetY()[ip],fkPoints->GetZ()[ip]));
  }
  fIter = iter;
  SetFitDone();
  //printf("F1>> %+.2e %+.2e %+.2e %+.2e\n",fParams[0],fParams[1],fParams[2],fParams[3]);
  return fChi2NDF;
  //
}

//____________________________________________________
void AliITSTPArrayFit::GetNormal(Double_t *norm, const Float_t *covMat) 
{
  // obtain the lab normal vector to the sensor from the covariance matrix
  // in such a way that when the local frame of the sensor coincides with 
  // the lab frame, the vector {0,1,0} is obtained
  Double_t tgxy = TMath::Tan(0.5*TMath::ATan2(2.*covMat[kXY],covMat[kYY]-covMat[kXX]));
  Double_t tgyz = TMath::Tan(0.5*TMath::ATan2(2.*covMat[kYZ],covMat[kZZ]-covMat[kYY]));
  norm[kY] = 1./TMath::Sqrt(1 + tgxy*tgxy + tgyz*tgyz);
  norm[kX] = norm[kY]*tgxy;
  norm[kZ] = norm[kY]*tgyz;
  //
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetDRofELoss(Double_t t,Double_t cdip,Double_t rhoL,const Double_t *normS, 
					Double_t &p,Double_t &e) const
{
  // Calculate energy loss of the particle at given t-param on the layer with rhoL (thickness*density) with
  // normal vector normS in the lab. The particle before eloss has energy "e" and momentum "p"
  // cdip = cosine of the dip angle = 1/sqrt(1+tgL^2)
  // Return the change DR of the radius due to the ELoss 
  //
  // NOTE: with B>0 the negative particles propagate along increasing t-param and positive 
  // particles - against.
  // t-param = 0 corresponds to the point of closest approach of the track to the beam.
  // Since the fitted helix parameters of the track are defined in this PCA point, when the correction
  // is applied upstream of the PCS, the energy must be increased (DR>0) rather than decreased (DR<0)
  //
  Double_t dirTr[3];
  dirTr[0] = -TMath::Sin(fParams[kPhi0]+t);
  dirTr[1] =  TMath::Cos(fParams[kPhi0]+t);
  dirTr[2] =  fParams[kDip];
  // cosine of the impact angle
  Double_t cosImp = cdip*TMath::Abs(dirTr[0]*normS[0]+dirTr[1]*normS[1]+dirTr[2]*normS[2]);
  //
  if (cosImp<0.3) cosImp = 0.3; //?
  Double_t dE = AliExternalTrackParam::BetheBlochSolid(p/fMass)*rhoL/cosImp;
  Double_t dP = e/p*dE;
  //
  if (t*GetSignQB() < 0) {
    dP =  -dP;
    dE =  -dE;
  }
  //
  if (p+dP<0) {
    AliInfo(Form("Estimated PLoss %.3f is larger than particle momentum %.3f. Skipping",dP,p));
    return 0;
  }
  //
  p += dP;
  e += dE;
  //
  return fCharge*dP*cdip/fBz/fgkCQConv;
}

//_____________________________________________________________
Double_t AliITSTPArrayFit::GetLineOffset(Int_t axis) const
{
  // return intercept of the parameterization coord = intercept + slope*t for given axis
  if (fParAxis<0) return -1E6; // no line fit
  if (axis==fParAxis) return 0;
  if (fParAxis==kX) return fParams[axis==kY ? kA0 : kA1 ];
  if (fParAxis==kY) return fParams[axis==kZ ? kA0 : kA1 ];
  return fParams[axis==kX ? kA0 : kA1 ];
}

//_____________________________________________________________
Double_t AliITSTPArrayFit::GetLineSlope(Int_t axis) const
{
  // return intercept of the parameterization coord = intercept + slope*t for given axis
  if (fParAxis<0) return -1E6; // no line fit
  if (axis==fParAxis) return 1.;
  if (fParAxis==kX) return fParams[axis==kY ? kB0 : kB1 ];
  if (fParAxis==kY) return fParams[axis==kZ ? kB0 : kB1 ];
  return fParams[axis==kX ? kB0 : kB1 ];
}

//_____________________________________________________________
void AliITSTPArrayFit::Print(Option_t *) const
{
  // print results of the fit
  //
  const char kCxyz[] = "XYZ";
  if (!fkPoints) return;
  //
  printf("Track of %3d points in Bz=%+.1f |Fit ",fPntLast-fPntFirst+1,fBz);
  if ( IsFitDone() ) {
    if (IsHelix())
      printf("Helix: Chi2: %5.1f | %+.2e %+.2e %+.2e %+.2e %+.2e\n",
	     fChi2NDF,fParams[kD0],fParams[kPhi0],fParams[kR0],fParams[kDZ],fParams[kDip]);
    else 
      printf("Line%c: Chi2: %5.1f | %+.2e %+.2e %+.2e %+.2e\n",
	     kCxyz[fParAxis],fChi2NDF,fParams[kA0],fParams[kB0],fParams[kA1],fParams[kB1]);
  }
  else printf("N/A\n");
}




//____________________________________________________
void AliITSTPArrayFit::BuildMaterialLUT(Int_t ntri) 
{
  // Fill a look-up table with mean material a la AliITSTrackerMI
  //
  if (!AliGeomManager::GetGeometry()) AliFatal("Geometry is not loaded");
  //
  // detector layer to check: dX,dZ,Ymin,Ymax
  const double kLayr[9][4] =  {{0.  ,60. , 2.80,3.00},  // beam pipe
			       {1.28,7.07,-0.20,0.22},  // SPD1
			       {1.28,7.07,-0.20,0.22},  // SPD2
			       {0.  ,76.0, 10.4,11.8},  // Shield1
			       {7.02,7.53,-1.00,4.50},  // SDD1
			       {7.02,7.53,-1.00,4.50},  // SDD2
			       {0.  ,102., 29.0,30.0},  // Shield2
			       {7.50,4.20,-0.15,4.50},  // SSD1
			       {7.50,4.20,-0.15,4.50}}; // SSD2
  //
  //
  // build <dens*L> for detectors (track hitting the sensor in normal direction)
  double pg1[3],pg2[3],res[7];
  //
  int sID = 0;
  int actLrID = 0;
  for (int lr=0;lr<9;lr++) {
    //
    Bool_t active = kFALSE;
    const double* tpars = kLayr[lr];
    //
    if (IsZero(tpars[0])) { // passive layer
      active = kFALSE;
      AliInfo(Form("Probing passive layer (total layer #%d)",lr));
    }  
    else {
      active = kTRUE;
      sID += AliGeomManager::LayerSize(++actLrID);
      AliInfo(Form("Probing sensors of active layer #%d (total layers #%d)",actLrID,lr));
    }
    double shift = TMath::Abs(tpars[2]-tpars[3])*1E-4;
    double rhol = 0;
    for (int i=ntri;i--;) {
      //
      if (active) {
	int ssID = sID -1 - AliGeomManager::LayerSize(actLrID)*gRandom->Rndm();
	pg1[0] = pg2[0] = (gRandom->Rndm()-0.5)*tpars[0] + shift; // local X
	pg2[0] -= 2*shift;
	pg1[1] = tpars[2];
	pg2[1] = tpars[3];
	pg1[2] = pg2[2] = (gRandom->Rndm()-0.5)*tpars[1] + shift; // local Z
	pg2[2] -= 2*shift;
	AliITSgeomTGeo::LocalToGlobal(ssID,pg1,pg1);    
	AliITSgeomTGeo::LocalToGlobal(ssID,pg2,pg2);	
      }
      else {
	double ang = gRandom->Rndm()*TMath::Pi()*2;
	pg1[0] = tpars[2]*TMath::Cos(ang)+shift;
	pg2[0] = tpars[3]*TMath::Cos(ang)-shift;
	pg1[1] = tpars[2]*TMath::Sin(ang);
	pg2[1] = tpars[3]*TMath::Sin(ang);
	pg1[2] = pg2[2] = (gRandom->Rndm()-0.5)*tpars[1]+shift; // local Z
	pg2[2] -= 2*shift;
      }

      //
      AliTracker::MeanMaterialBudget(pg1,pg2,res);
      rhol += res[0]*res[4];   // rho*L
    }
    fgRhoLITS[lr] = rhol/ntri;
    AliInfo(Form("Obtained <rho*L> = %e\n",fgRhoLITS[lr]));
  }
  //
  return;
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetPCA2PlaneInfo(Double_t *xyz, Double_t *dir, Int_t axis, Double_t axval) const
{
  // calculate the PCA to plane normal ti axis and crossing it at axval
  // fill the position and direction cosines at this point
  //
  double xyzp[3] = {0,0,0};                // create fake point
  xyzp[axis] = axval;
  double covI[6] = {1e-4,0,0,1e-4,0,1e-4}; // fake cov.matrix loose in all directions
  covI[4*axis - axis*(axis+1)/2] = 1e8;    // except axis
  //
  double t = GetPosition(xyz, xyzp, covI); // got pca
  //
  if (dir) GetDirCos(dir,t);
  return t;
}

//____________________________________________________
void AliITSTPArrayFit::GetT0Info(Double_t* xyz, Double_t *dir) const
{
  // get direction cosines for t = 0;
  GetPosition(xyz, 0);
  if (dir) GetDirCos(dir,0);
}

//____________________________________________________
Bool_t AliITSTPArrayFit::CalcErrorMatrix()
{
  // compute covariance matrix in lenear approximation of residuals on parameters
  static AliSymMatrix cv(5);
  static Double_t dres[5][3]; 
  cv.Reset();
  int npar = IsHelix() ? 5:4;
  //
  for (int ip=fPntFirst;ip<=fPntLast;ip++) {
    GetDResDParams(&dres[0][0],ip);
    Double_t* covI = GetCovI(ip);
    for (int i=npar;i--;) 
      for (int j=i+1;j--;)
	cv(i,j) += dres[i][kX]*(dres[j][kX]*covI[ kXX ] + dres[j][kY]*covI[ kXY ] + dres[j][kZ]*covI[ kXZ ])
	  +        dres[i][kY]*(dres[j][kX]*covI[ kXY ] + dres[j][kY]*covI[ kYY ] + dres[j][kZ]*covI[ kYZ ])
	  +        dres[i][kZ]*(dres[j][kX]*covI[ kXZ ] + dres[j][kY]*covI[ kYZ ] + dres[j][kZ]*covI[ kZZ ]);
  }
  cv.SetSizeUsed(npar);
  if (cv.InvertChol()) {
    //cv.Print("l");
    int cnt = 0;
    for (int i=npar;i--;) for (int j=i+1;j--;)fParamsCov[cnt++] = cv(i,j);
    return kTRUE;
  }
  //
  return kFALSE;
}

//____________________________________________________
Double_t AliITSTPArrayFit::CalcParPCA(Int_t ipnt) const
{
  // get parameter for the point with least weighted distance to the point
  const double *xyz  = GetPoint(ipnt);
  const double *covI = GetCovI(ipnt);
  if (IsHelix()) return GetParPCAHelix(xyz,covI);
  else             return GetParPCALine(xyz,covI);
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetPt() const 
{
  return IsFieldON()&&IsHelix() ? TMath::Abs(fParams[kR0]*fBz*fgkCQConv) : -1;
}

//____________________________________________________
Double_t AliITSTPArrayFit::GetP() const 
{
  if (!IsFieldON()) return -1;
  Double_t cdip = 1./TMath::Sqrt(1.+fParams[kDip]*fParams[kDip]);
  return TMath::Abs(fParams[kR0]*fgkCQConv*fBz/cdip); // momentum
}


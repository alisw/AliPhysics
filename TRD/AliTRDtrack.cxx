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

#include "AliTracker.h"
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

AliTRDtrack::AliTRDtrack():
  AliKalmanTrack(),
  fSeedLab(-1),
  fdEdx(0),
  fdEdxT(0),
  fDE(0),
  fStopped(kFALSE),
  fLhElectron(0),
  fNWrong(0),
  fNRotate(0),
  fNCross(0),
  fNExpected(0),
  fNLast(0),
  fNExpectedLast(0),
  fNdedx(0),
  fChi2Last(1e10),
  fBackupTrack(0x0)
{
  for (Int_t i=0; i<kNplane; i++) {
    for (Int_t j=0; j<kNslice; j++) {
      fdEdxPlane[i][j] = 0;
    }
    fTimBinPlane[i] = -1;
  }
  for (UInt_t i=0; i<kMAXCLUSTERSPERTRACK; i++) {
    fIndex[i] = 0;
    fIndexBackup[i] = 0;
    fdQdl[i] = 0;
  }
  for (Int_t i=0; i<3; i++) fBudget[i] = 0;
}

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliTRDcluster *c, Int_t index, 
                         const Double_t p[5], const Double_t cov[15], 
                         Double_t x, Double_t alpha) : 
  AliKalmanTrack(),
  fSeedLab(-1),
  fdEdx(0),
  fdEdxT(0),
  fDE(0),
  fStopped(kFALSE),
  fLhElectron(0),
  fNWrong(0),
  fNRotate(0),
  fNCross(0),
  fNExpected(0),
  fNLast(0),
  fNExpectedLast(0),
  fNdedx(0),
  fChi2Last(1e10),
  fBackupTrack(0x0) 
{
  //-----------------------------------------------------------------
  // This is the main track constructor.
  //-----------------------------------------------------------------
  Double_t cnv=1./(GetBz()*kB2C);

  Double_t pp[5]={
    p[0],
    p[1],
    x*p[4] - p[2],
    p[3],
    p[4]*cnv
  };

  Double_t c22 = x*x*cov[14] - 2*x*cov[12] + cov[5];
  Double_t c32 = x*cov[13] - cov[8];
  Double_t c20 = x*cov[10] - cov[3], 
           c21 = x*cov[11] - cov[4], c42 = x*cov[14] - cov[12];

  Double_t cc[15]={
    cov[0 ],
    cov[1 ],     cov[2 ],
    c20,         c21,         c22,
    cov[6 ],     cov[7 ],     c32,     cov[9 ],
    cov[10]*cnv, cov[11]*cnv, c42*cnv, cov[13]*cnv, cov[14]*cnv*cnv
  };

  Set(x,alpha,pp,cc);

  SetNumberOfClusters(1);
  fIndex[0]=index;

  for (Int_t i=0;i<kNplane;i++){
    for (Int_t j=0; j<kNslice; j++) {
      fdEdxPlane[i][j] = 0;
    }
    fTimBinPlane[i] = -1;
  }

  Double_t q = TMath::Abs(c->GetQ());
  Double_t s = GetSnp(), t=GetTgl();
  if(s*s < 1) q *= TMath::Sqrt((1-s*s)/(1+t*t));

  fdQdl[0] = q;  
  // initialisation [SR, GSI 18.02.2003] (i startd for 1)
  for(UInt_t i=1; i<kMAXCLUSTERSPERTRACK; i++) {
    fdQdl[i] = 0;
    fIndex[i] = 0;
    fIndexBackup[i] = 0;  //backup indexes MI    
  }
  for (Int_t i=0;i<3;i++) fBudget[i]=0;
}                              
           
//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliTRDtrack& t) : 
  AliKalmanTrack(t), 
  fSeedLab(t.GetSeedLabel()),
  fdEdx(t.fdEdx),
  fdEdxT(t.fdEdx),
  fDE(t.fDE),
  fStopped(t.fStopped),
  fLhElectron(0),
  fNWrong(t.fNWrong),
  fNRotate(t.fNRotate),
  fNCross(t.fNCross),
  fNExpected(t.fNExpected),
  fNLast(t.fNLast),
  fNExpectedLast(t.fNExpectedLast),
  fNdedx(t.fNdedx),
  fChi2Last(t.fChi2Last),
  fBackupTrack(0x0) 
{
  //
  // Copy constructor.
  //
  for (Int_t i=0;i<kNplane;i++){
    for (Int_t j=0; j<kNslice; j++) {
      fdEdxPlane[i][j] = t.fdEdxPlane[i][j];
    }
    fTimBinPlane[i] = t.fTimBinPlane[i];
    fTracklets[i]   = t.fTracklets[i];
  }

  Int_t n=t.GetNumberOfClusters(); 
  SetNumberOfClusters(n);
  for (Int_t i=0; i<n; i++) {
    fIndex[i]=t.fIndex[i];
    fIndexBackup[i]=t.fIndex[i];  // MI - backup indexes
    fdQdl[i]=t.fdQdl[i];
  }

  // initialisation (i starts from n) [SR, GSI, 18.02.2003]
  for(UInt_t i=n; i<kMAXCLUSTERSPERTRACK; i++) {
    fdQdl[i] = 0;
    fIndex[i] = 0;
    fIndexBackup[i] = 0;  //MI backup indexes
  }
  for (Int_t i=0;i<3;i++) fBudget[i]=t.fBudget[i];
}                                

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliKalmanTrack& t, Double_t /*alpha*/): 
  AliKalmanTrack(t), 
  fSeedLab(-1),
  fdEdx(t.GetPIDsignal()),
  fdEdxT(0),
  fDE(0),
  fStopped(kFALSE),
  fLhElectron(0.0),
  fNWrong(0),
  fNRotate(0),
  fNCross(0),
  fNExpected(0),
  fNLast(0),
  fNExpectedLast(0),
  fNdedx(0),
  fChi2Last(0.0),
  fBackupTrack(0x0)
{
  //
  // Constructor from AliTPCtrack or AliITStrack .
  //

  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetMass(t.GetMass());
  SetNumberOfClusters(0);

  for (Int_t i=0;i<kNplane;i++){
    for (Int_t j=0;j<kNslice;j++){
      fdEdxPlane[i][j] = 0.0;
    }
    fTimBinPlane[i] = -1;
  }

  // Initialization [SR, GSI, 18.02.2003]
  for(UInt_t i=0; i<kMAXCLUSTERSPERTRACK; i++) {
    fdQdl[i] = 0;
    fIndex[i] = 0;
    fIndexBackup[i] = 0;  // MI backup indexes    
  }
  
  for (Int_t i=0;i<3;i++) { fBudget[i]=0;};
}              

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliESDtrack &t):
  AliKalmanTrack(), 
  fSeedLab(-1),
  fdEdx(t.GetTRDsignal()),
  fdEdxT(0),
  fDE(0),
  fStopped(kFALSE),
  fLhElectron(0),
  fNWrong(0),
  fNRotate(0),
  fNCross(0),
  fNExpected(0),
  fNLast(0),
  fNExpectedLast(0),
  fNdedx(0),
  fChi2Last(1e10),
  fBackupTrack(0x0)
{
  //
  // Constructor from AliESDtrack
  //
  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetMass(t.GetMass());
  SetNumberOfClusters(t.GetTRDclusters(fIndex)); 

  Int_t ncl = t.GetTRDclusters(fIndexBackup);
  for (UInt_t i=ncl;i<kMAXCLUSTERSPERTRACK;i++) {
    fIndexBackup[i]=0;
    fIndex[i] = 0; //MI store indexes
  }
  for (Int_t i=0;i<kNplane;i++){
    for (Int_t j=0;j<kNslice;j++){
      fdEdxPlane[i][j] = t.GetTRDsignals(i,j);
    }
    fTimBinPlane[i] = t.GetTRDTimBin(i);
  }


  const AliExternalTrackParam *par=&t;
  if (t.GetStatus()&AliESDtrack::kTRDbackup) { 
    par=t.GetOuterParam();
    if (!par) {AliError("***** No backup info !!! ****"); par=&t;}
  }
  Set(par->GetX(),par->GetAlpha(),par->GetParameter(),par->GetCovariance());

  
  for (UInt_t i=0; i<kMAXCLUSTERSPERTRACK; i++) fdQdl[i] = 0;

  for (Int_t i=0;i<3;i++) fBudget[i]=0;

  if ((t.GetStatus()&AliESDtrack::kTIME) == 0) return;
  StartTimeIntegral();
  Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());

}  

//____________________________________________________________________________
AliTRDtrack::~AliTRDtrack()
{
  //
  // Destructor
  //

  if (fBackupTrack) delete fBackupTrack;
  fBackupTrack = 0;

}

//____________________________________________________________________________
Float_t AliTRDtrack::StatusForTOF()
{
  //
  // Defines the status of the TOF extrapolation
  //

  Float_t res = (0.2 + 0.8*(fN/(fNExpected+5.)))*(0.4+0.6*fTracklets[5].GetN()/20.);
  res *= (0.25+0.8*40./(40.+fBudget[2]));
  return res;

  Int_t status=0;
  if (GetNumberOfClusters()<20) return 0;   //
  if (fN>110&&fChi2/(Float_t(fN))<3) return 3;            //gold
  if (fNLast>30&&fChi2Last/(Float_t(fNLast))<3) return 3; //gold
  if (fNLast>20&&fChi2Last/(Float_t(fNLast))<2) return 3; //gold
  if (fNLast/(fNExpectedLast+3.)>0.8 && fChi2Last/Float_t(fNLast)<5&&fNLast>20) return 2; //silber
  if (fNLast>5 &&((fNLast+1.)/(fNExpectedLast+1.))>0.8&&fChi2Last/(fNLast-5.)<6)   return 1; 
  
  return status;

}

//_____________________________________________________________________________
Int_t AliTRDtrack::Compare(const TObject *o) const 
{
  //
  // Compares tracks according to their Y2 or curvature
  //

  AliTRDtrack *t = (AliTRDtrack *) o;
  //  Double_t co=t->GetSigmaY2();
  //  Double_t c =GetSigmaY2();

  Double_t co = TMath::Abs(t->GetC());
  Double_t c  = TMath::Abs(GetC());  

  if      (c > co) {
    return 1;
  }
  else if (c < co) {
    return -1;
  }
  return 0;

}                

//_____________________________________________________________________________
void AliTRDtrack::CookdEdx(Double_t low, Double_t up) 
{
  //
  // Calculates the truncated dE/dx within the "low" and "up" cuts.
  //

  Int_t   i  = 0;

  // Array to sort the dEdx values according to amplitude
  Float_t sorted[kMAXCLUSTERSPERTRACK];

  // Number of clusters used for dedx
  Int_t   nc = fNdedx; 

  // Require at least 10 clusters for a dedx measurement
  if (nc < 10) {
    SetdEdx(0);
    return;
  }

  // Lower and upper bound
  Int_t nl = Int_t(low * nc);
  Int_t nu = Int_t( up * nc);

  // Can fdQdl be negative ????
  for (i = 0; i < nc; i++) {
    sorted[i] = TMath::Abs(fdQdl[i]);
  }

  // Sort the dedx values by amplitude
  Int_t *index = new Int_t[nc];
  TMath::Sort(nc,sorted,index,kFALSE);

  // Sum up the truncated charge between nl and nu  
  Float_t dedx = 0.0;
  for (i = nl; i <= nu; i++) {
    dedx += sorted[index[i]];
  }
  dedx /= (nu - nl + 1.0);
  SetdEdx(dedx);

}                     

//_____________________________________________________________________________
Bool_t AliTRDtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho)
{
  // Propagates a track of particle with mass=pm to a reference plane 
  // defined by x=xk through media of density=rho and radiationLength=x0

  if (xk == GetX()) return kTRUE;

  Double_t oldX=GetX(), oldY=GetY(), oldZ=GetZ();

  Double_t bz=GetBz();
  if (!AliExternalTrackParam::PropagateTo(xk,bz)) return kFALSE;

  Double_t x=GetX(), y=GetY(), z=GetZ();
 Double_t d=TMath::Sqrt((x-oldX)*(x-oldX)+(y-oldY)*(y-oldY)+(z-oldZ)*(z-oldZ));
  if (oldX < xk)
  if (IsStartedTimeIntegral()) {
    Double_t l2=d;
    Double_t crv=GetC();
    if (TMath::Abs(l2*crv)>0.0001){
      // make correction for curvature if neccesary
      l2 = 0.5*TMath::Sqrt((x-oldX)*(x-oldX) + (y-oldY)*(y-oldY));
      l2 = 2*TMath::ASin(l2*crv)/crv;
      l2 = TMath::Sqrt(l2*l2+(z-oldZ)*(z-oldZ));
    }
    AddTimeStep(l2);
  }

  Double_t ll = (oldX < xk) ? -d : d;
  if (!AliExternalTrackParam::CorrectForMaterial(ll*rho/x0,x0,GetMass())) 
     return kFALSE;

  {//Energy losses************************
  Double_t p2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t beta2=p2/(p2 + GetMass()*GetMass());
  if ((5940*beta2/(1-beta2+1e-10) - beta2) < 0) return kFALSE;

  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2+1e-10)) - beta2)*d*rho;
  Float_t budget = d*rho;
  fBudget[0] += budget;
  /*
  // suspicious part - think about it ?
  Double_t kinE =  TMath::Sqrt(p2);
  if (dE>0.8*kinE) dE = 0.8*kinE;  //      
  if (dE<0)        dE = 0.0;       // not valid region for Bethe bloch 
  */
  //
  fDE+=dE;
  /*
  // Suspicious ! I.B.
  Double_t sigmade = 0.07*TMath::Sqrt(TMath::Abs(dE));   // energy loss fluctuation 
  Double_t sigmac2 = sigmade*sigmade*fC*fC*(p2+GetMass()*GetMass())/(p2*p2);
  fCcc += sigmac2;
  fCee += fX*fX*sigmac2;  
  */
  }

  return kTRUE;            

}     

//_____________________________________________________________________________
Bool_t AliTRDtrack::Update(const AliTRDcluster *c, Double_t chisq, Int_t index,
                          Double_t h01)
{
  // Assignes found cluster to the track and updates track information

  Bool_t fNoTilt = kTRUE;
  if(TMath::Abs(h01) > 0.003) fNoTilt = kFALSE;
  // add angular effect to the error contribution -  MI
  Float_t tangent2 = GetSnp()*GetSnp();
  if (tangent2 < 0.90000){
    tangent2 = tangent2/(1.-tangent2);
  }
  //Float_t errang = tangent2*0.04; //

  Double_t p[2]={c->GetY(), c->GetZ()};
  //Double_t cov[3]={c->GetSigmaY2()+errang, 0., c->GetSigmaZ2()*100.};
  Double_t sy2=c->GetSigmaY2()*4;
  Double_t sz2=c->GetSigmaZ2()*4;
  Double_t cov[3]={sy2 + h01*h01*sz2, h01*(sy2-sz2), sz2 + h01*h01*sy2};

  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;

  Int_t n=GetNumberOfClusters();
  fIndex[n]=index;
  SetNumberOfClusters(n+1);

  SetChi2(GetChi2()+chisq);

  return kTRUE;     
}                     

//_____________________________________________________________________________
Int_t AliTRDtrack::UpdateMI(const AliTRDcluster *c, Double_t chisq, Int_t index,                            Double_t h01, Int_t /*plane*/) {
  // Assignes found cluster to the track and updates track information
  Bool_t fNoTilt = kTRUE;
  if(TMath::Abs(h01) > 0.003) fNoTilt = kFALSE;
  // add angular effect to the error contribution and make correction  -  MI
  // 
  Double_t tangent2 = GetSnp()*GetSnp();
  if (tangent2 < 0.90000){
    tangent2 = tangent2/(1.-tangent2);
  }
  Double_t tangent = TMath::Sqrt(tangent2);
  if (GetSnp()<0) tangent*=-1;
  //  Double_t correction = 0*plane;
  /*
  Double_t errang = tangent2*0.04;  //
  Double_t errsys =0.025*0.025*20;  //systematic error part 

  Float_t extend =1;
  if (c->GetNPads()==4) extend=2;
  */
  //if (c->GetNPads()==5)  extend=3;
  //if (c->GetNPads()==6)  extend=3;
  //if (c->GetQ()<15) return 1;

  /*
  if (corrector!=0){
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
  //Double_t padlength = TMath::Sqrt(c->GetSigmaZ2()*12.);
  /*
    {
      Double_t dy=c->GetY() - GetY(), dz=c->GetZ() - GetZ();     
      printf("%e %e %e %e\n",dy,dz,padlength/2,h01);
    }
  */
  Double_t p[2]={c->GetY(), c->GetZ()};
  /*
  Double_t cov[3]={(c->GetSigmaY2()+errang+errsys)*extend, 0., 
  		   c->GetSigmaZ2()*10000.};
  */
  Double_t sy2=c->GetSigmaY2()*4;
  Double_t sz2=c->GetSigmaZ2()*4;
  Double_t cov[3]={sy2 + h01*h01*sz2, h01*(sy2-sz2), sz2 + h01*h01*sy2};

  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;

  Int_t n=GetNumberOfClusters();
  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chisq);

  return kTRUE;      
}                     
/*
//_____________________________________________________________________________
Int_t AliTRDtrack::UpdateMI(const AliTRDtracklet &tracklet)
{
  //
  // Assignes found tracklet to the track and updates track information
  //
  //
  Double_t r00=(tracklet.GetTrackletSigma2()), r01=0., r11= 10000.;
  r00+=fCyy; r01+=fCzy; r11+=fCzz;
  //
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;
  //

  Double_t dy=tracklet.GetY() - fY, dz=tracklet.GetZ() - fZ;

  
  Double_t s00 = tracklet.GetTrackletSigma2();  // error pad
  Double_t s11 = 100000;   // error pad-row
  Float_t  h01 = tracklet.GetTilt();
  //
  //  r00 = fCyy + 2*fCzy*h01 + fCzz*h01*h01+s00;
  r00 = fCyy + fCzz*h01*h01+s00;
  //  r01 = fCzy + fCzz*h01;
  r01 = fCzy ;
  r11 = fCzz + s11;
  det = r00*r11 - r01*r01;
  // inverse matrix
  tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

  Double_t k00=fCyy*r00+fCzy*r01, k01=fCyy*r01+fCzy*r11;
  Double_t k10=fCzy*r00+fCzz*r01, k11=fCzy*r01+fCzz*r11;
  Double_t k20=fCey*r00+fCez*r01, k21=fCey*r01+fCez*r11;
  Double_t k30=fCty*r00+fCtz*r01, k31=fCty*r01+fCtz*r11;
  Double_t k40=fCcy*r00+fCcz*r01, k41=fCcy*r01+fCcz*r11;
  
  // K matrix
//   k00=fCyy*r00+fCzy*(r01+h01*r00),k01=fCyy*r01+fCzy*(r11+h01*r01);
//   k10=fCzy*r00+fCzz*(r01+h01*r00),k11=fCzy*r01+fCzz*(r11+h01*r01);
//   k20=fCey*r00+fCez*(r01+h01*r00),k21=fCey*r01+fCez*(r11+h01*r01);
//   k30=fCty*r00+fCtz*(r01+h01*r00),k31=fCty*r01+fCtz*(r11+h01*r01);
//   k40=fCcy*r00+fCcz*(r01+h01*r00),k41=fCcy*r01+fCcz*(r11+h01*r01);  
  //
  //Update measurement
  Double_t cur=fC + k40*dy + k41*dz, eta=fE + k20*dy + k21*dz;  
  //  cur=fC + k40*dy + k41*dz; eta=fE + k20*dy + k21*dz;
  if (TMath::Abs(cur*fX-eta) >= 0.90000) {
    //Int_t n=GetNumberOfClusters();
    //      if (n>4) cerr<<n<<" AliTRDtrack warning: Filtering failed !\n";
    return 0;
  }                           
//   k01+=h01*k00;
//   k11+=h01*k10;
//   k21+=h01*k20;
//   k31+=h01*k30;
//   k41+=h01*k40;  


  fY += k00*dy + k01*dz;
  fZ += k10*dy + k11*dz;
  fE  = eta;
  fT += k30*dy + k31*dz;
  fC  = cur;
    
  
  //Update covariance
  //
  //
  Double_t oldyy = fCyy, oldzz = fCzz; //, oldee=fCee, oldcc =fCcc;
  Double_t oldzy = fCzy, oldey = fCey, oldty=fCty, oldcy =fCcy;
  Double_t oldez = fCez, oldtz = fCtz, oldcz=fCcz;
  //Double_t oldte = fCte, oldce = fCce;
  //Double_t oldct = fCct;

  fCyy-=k00*oldyy+k01*oldzy;   
  fCzy-=k10*oldyy+k11*oldzy;
  fCey-=k20*oldyy+k21*oldzy;   
  fCty-=k30*oldyy+k31*oldzy;
  fCcy-=k40*oldyy+k41*oldzy;  
  //
  fCzz-=k10*oldzy+k11*oldzz;
  fCez-=k20*oldzy+k21*oldzz;   
  fCtz-=k30*oldzy+k31*oldzz;
  fCcz-=k40*oldzy+k41*oldzz;
  //
  fCee-=k20*oldey+k21*oldez;   
  fCte-=k30*oldey+k31*oldez;
  fCce-=k40*oldey+k41*oldez;
  //
  fCtt-=k30*oldty+k31*oldtz;
  fCct-=k40*oldty+k41*oldtz;
  //
  fCcc-=k40*oldcy+k41*oldcz;                 
  //
  
  //Int_t n=GetNumberOfClusters();
  //fIndex[n]=index;
  //SetNumberOfClusters(n+1);

  //SetChi2(GetChi2()+chisq);
  //  cerr<<"in update: fIndex["<<fN<<"] = "<<index<<endl;

  return 1;      

}                     
*/

//_____________________________________________________________________________
Bool_t AliTRDtrack::Rotate(Double_t alpha, Bool_t absolute)
{
  // Rotates track parameters in R*phi plane
  // if absolute rotation alpha is in global system
  // otherwise alpha rotation is relative to the current rotation angle
  
  if (absolute) {
    alpha -= GetAlpha();
  }
  else{
    fNRotate++;
  }

  return AliExternalTrackParam::Rotate(GetAlpha()+alpha);
}                         

//_____________________________________________________________________________
Double_t AliTRDtrack::GetPredictedChi2(const AliTRDcluster *c, Double_t h01) const
{
  //
  // Returns the track chi2
  //  

  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t sy2=c->GetSigmaY2()*4;
  Double_t sz2=c->GetSigmaZ2()*4;
  Double_t cov[3]={sy2 + h01*h01*sz2, h01*(sy2-sz2), sz2 + h01*h01*sy2};

  return AliExternalTrackParam::GetPredictedChi2(p,cov);

  /*
  Bool_t fNoTilt = kTRUE;
  if(TMath::Abs(h01) > 0.003) fNoTilt = kFALSE;

  return (c->GetY() - GetY())*(c->GetY() - GetY())/c->GetSigmaY2();
  */

  /*
  Double_t chi2, dy, r00, r01, r11;

  if(fNoTilt) {
    dy=c->GetY() - fY;
    r00=c->GetSigmaY2();    
    chi2 = (dy*dy)/r00;    
  }
  else {
    Double_t padlength = TMath::Sqrt(c->GetSigmaZ2()*12);
    //
    r00=c->GetSigmaY2(); r01=0.; r11=c->GetSigmaZ2();
    r00+=fCyy; r01+=fCzy; r11+=fCzz;

    Double_t det=r00*r11 - r01*r01;
    if (TMath::Abs(det) < 1.e-10) {
      Int_t n=GetNumberOfClusters(); 
      if (n>4) cerr<<n<<" AliTRDtrack warning: Singular matrix !\n";
      return 1e10;
    }
    Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
    Double_t dy=c->GetY() - fY, dz=c->GetZ() - fZ;
    Double_t tiltdz = dz;
    if (TMath::Abs(tiltdz)>padlength/2.) {
      tiltdz = TMath::Sign(padlength/2,dz);
    }
    //    dy=dy+h01*dz;
    dy=dy+h01*tiltdz;

    chi2 = (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det; 
  }

  return chi2;
  */
}      

//_____________________________________________________________________________
void AliTRDtrack::MakeBackupTrack()
{
  //
  // Creates a backup track
  //

  if (fBackupTrack) delete fBackupTrack;
  fBackupTrack = new AliTRDtrack(*this);
  
}

//_____________________________________________________________________________
Int_t AliTRDtrack::GetProlongation(Double_t xk, Double_t &y, Double_t &z)
{
  //
  // Find prolongation at given x
  // return 0 if not exist
  
  Double_t bz=GetBz();

  if (!AliExternalTrackParam::GetYAt(xk,bz,y)) return 0;
  if (!AliExternalTrackParam::GetZAt(xk,bz,z)) return 0;

  return 1;  
}

//_____________________________________________________________________________
Int_t   AliTRDtrack::PropagateToX(Double_t xr, Double_t step)
{
  //
  // Propagate track to given x  position 
  // works inside of the 20 degree segmentation (local cooordinate frame for TRD , TPC, TOF)
  // 
  // material budget from geo manager
  // 
  Double_t  xyz0[3], xyz1[3],y,z;
  const Double_t kAlphac  = TMath::Pi()/9.;   
  const Double_t kTalphac = TMath::Tan(kAlphac*0.5);
  // critical alpha  - cross sector indication
  //
  Double_t dir = (GetX()>xr) ? -1.:1.;
  // direction +-
  for (Double_t x=GetX()+dir*step;dir*x<dir*xr;x+=dir*step){
    //
    GetXYZ(xyz0);	
    GetProlongation(x,y,z);
    xyz1[0] = x*TMath::Cos(GetAlpha())+y*TMath::Sin(GetAlpha()); 
    xyz1[1] = x*TMath::Sin(GetAlpha())-y*TMath::Cos(GetAlpha());
    xyz1[2] = z;
    Double_t param[7];
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
    //
    if (param[0]>0&&param[1]>0) PropagateTo(x,param[1],param[0]);
    if (GetY()>GetX()*kTalphac){
      Rotate(-kAlphac);
    }
    if (GetY()<-GetX()*kTalphac){
      Rotate(kAlphac);
    }
  }
  //
  PropagateTo(xr);

  return 0;

}

//_____________________________________________________________________________
Int_t   AliTRDtrack::PropagateToR(Double_t r,Double_t step)
{
  //
  // propagate track to the radial position
  // rotation always connected to the last track position
  //
  Double_t  xyz0[3], xyz1[3],y,z; 
  Double_t radius = TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  Double_t dir = (radius>r) ? -1.:1.;   // direction +-
  //
  for (Double_t x=radius+dir*step;dir*x<dir*r;x+=dir*step){
    GetXYZ(xyz0);	
    Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
    Rotate(alpha,kTRUE);
    GetXYZ(xyz0);	
    GetProlongation(x,y,z);
    xyz1[0] = x*TMath::Cos(alpha)+y*TMath::Sin(alpha); 
    xyz1[1] = x*TMath::Sin(alpha)-y*TMath::Cos(alpha);
    xyz1[2] = z;
    Double_t param[7];
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
    if (param[1]<=0) param[1] =100000000;
    PropagateTo(x,param[1],param[0]);
  } 
  GetXYZ(xyz0);	
  Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
  Rotate(alpha,kTRUE);
  GetXYZ(xyz0);	
  GetProlongation(r,y,z);
  xyz1[0] = r*TMath::Cos(alpha)+y*TMath::Sin(alpha); 
  xyz1[1] = r*TMath::Sin(alpha)-y*TMath::Cos(alpha);
  xyz1[2] = z;
  Double_t param[7];
  AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
  //
  if (param[1]<=0) param[1] =100000000;
  PropagateTo(r,param[1],param[0]);

  return 0;

}

//_____________________________________________________________________________
Int_t AliTRDtrack::GetSector() const
{
  //
  // Return the current sector
  //

  return Int_t(TVector2::Phi_0_2pi(GetAlpha())
             / AliTRDgeometry::GetAlpha())
             % AliTRDgeometry::kNsect;

}

//_____________________________________________________________________________
void AliTRDtrack::SetSampledEdx(Float_t q, Int_t i)    
{
  //
  // The sampled energy loss
  //

  Double_t s = GetSnp();
  Double_t t = GetTgl();
  q *= TMath::Sqrt((1-s*s)/(1+t*t));
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
  q*= TMath::Sqrt((1-s*s)/(1+t*t));
  fdQdl[fNdedx] = q;
  fNdedx++;

}     

Double_t AliTRDtrack::GetBz() const {
  //
  // returns Bz component of the magnetic field (kG)
  //
  if (AliTracker::UniformField()) return AliTracker::GetBz();
  Double_t r[3]; GetXYZ(r);
  return AliTracker::GetBz(r);
}



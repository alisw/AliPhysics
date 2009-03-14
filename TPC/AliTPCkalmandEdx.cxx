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

/* 
   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
   Not final version yet - fitting using histograms to be implemented
   Expected to be more robust.                     


   Model:
   TPC Kalman filter implementation for TPC dEdx alignment:
   The dEdx equalization for 3 types of the pad:
   Short (0.75 cm) Medium (1 cm) and Long(1.5 cm)

   Model of correction ci:
   corrected   raw    correction
   dEdxti    = dEdxri*ci = dEdxri*p01*(1+p1i*kY+p2i*kZ+p3i*dR+p4i/sqrt(1+kY^2+kZ^2)
   //
   Matching - update using 2 tracklets  ::UpdatedEdxPair
   dEdxti-dEdxtj=0 = dEdxri*ci-dEdxrj*cj
   //
   Matching - normalization of the signal ::UpdatedEdx 
   possible only for identified particles
   dEdxti = dN/dx * sqrt(1+kY^2+kZ^2)
*/

#include "TMath.h"
#include "TTreeStream.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "AliTPCROC.h"
#include "AliTPCkalmandEdx.h"


AliTPCkalmandEdx::AliTPCkalmandEdx():
  TNamed(),
  fState(0),
  fCovariance(0), 
  fMat1(0),                  
  fNpad(3),                  // number of pad types
  fNpar(5),                // number of parameters
  fNelem(3*5),                // number of elements
  fSampleSize(0),
  fInit(0)
{
  //
  // Default constructor
  //
}

AliTPCkalmandEdx::AliTPCkalmandEdx(const char* name, const char* title, Int_t sampleSize): 
  TNamed(name,title),
  fState(0),
  fCovariance(0),   
  fMat1(0),                  
  fNpad(3),                  // number of pad types
  fNpar(5),                // number of parameters
  fNelem(3*5),             // number of elements
  fSampleSize(sampleSize),
  fInit(0)
{
  //
  // Default constructor
  //
  Init();
}

AliTPCkalmandEdx::AliTPCkalmandEdx(const AliTPCkalmandEdx & kalman):
  TNamed(kalman),
  fState(0),
  fCovariance(0),   
  fMat1(0),                  
  fNpad(kalman.fNpad),                  // number of pad types
  fNpar(kalman.fNpar),                // number of parameters
  fNelem(kalman.fNelem),                // number of elements
  fSampleSize(kalman.fSampleSize)
{
  //
  // copy constructor
  //
  fState      = new TMatrixD(*(kalman.fState));
  fCovariance = new TMatrixD(*(kalman.fCovariance));
  fMat1       = new TMatrixD(*(kalman.fMat1));
} 

void AliTPCkalmandEdx::Init(){
  //
  // Default parameters setting
  //
  fState      = new TMatrixD(fNelem,1);
  fCovariance = new TMatrixD(fNelem,fNelem);
  fMat1       = new TMatrixD(fNelem,fNelem);

  fInit=0;
  for (Int_t i=0;i<3; i++) {
    fSample[i].ResizeTo(fSampleSize);
    fSampleStat[i].ResizeTo(2);
    fCounter[i]=0;
  }

  TMatrixD &vecXk=*fState;
  TMatrixD &covXk=*fCovariance;
  TMatrixD &mat1=*fMat1;
  //
  //
  for (Int_t i=0;i<fNelem;i++) 
    for(Int_t j=0;j<fNelem;j++) {covXk(i,j)=0;  mat1(i,j)=0;}
  //
  for (Int_t i=0;i<fNelem;i++)  {covXk(i,i)=1.; mat1(i,i)=1; vecXk(i,0)=0;}
  for (Int_t ipad=0;ipad<3;ipad++){
    vecXk(ipad*fNpar,0)=1;
  }
//   //
//   // set reference ipad=0
//   vecXk(1*fNpar+0,0)=1;
//   vecXk(1*fNpar+1,0)=0;
//   vecXk(1*fNpar+2,0)=0;
//   vecXk(1*fNpar+3,0)=0;
//   vecXk(1*fNpar+4,0)=0;
  

}

void AliTPCkalmandEdx::UpdatedEdxPair(Int_t ip0, Int_t ip1,
				      Double_t dedx0, Double_t dedx1, 
				      Double_t s0, Double_t s1, 
				      Double_t kY0, Double_t kY1,
				      Double_t kZ0, Double_t kZ1,
				      Double_t dR0, Double_t dR1, 
				      TTreeSRedirector *debug){
  //
  // update relative measurement
  // 0 = dEdxti-dEdxtj  = dEdxri*ci-dEdxrj*cj
  //
  // Model of correction ci:
  // dEdxti = dEdxri*ci = dEdxri*p01*(1+p1i*kY+p2i*kZ+p3i*dR+p4i/sqrt(1+kY^2+kZ^2)
  if (fInit<3) return;  // not initialized parameters
  const Double_t kchi2Cut = 3.*3.;

  Int_t nmeas = 1;
  TMatrixD vecXk=*fState;           // X vector
  TMatrixD covXk=*fCovariance;      // X covariance
  TMatrixD &mat1=*fMat1;            // update covariance matrix

  Double_t length0 = TMath::Sqrt(1+kY0*kY0+kZ0*kZ0);
  Double_t length1 = TMath::Sqrt(1+kY1*kY1+kZ1*kZ1);

  Double_t corr0 = vecXk(ip0*fNpar,0)*
    (1+vecXk(ip0*fNpar+1,0)*kY0+
     vecXk(ip0*fNpar+2,0)*kZ0+vecXk(ip0*fNpar+3,0)*dR0+vecXk(ip0*fNpar+4,0)/length0);
  Double_t corr1 = vecXk(ip1*fNpar,0)*
    (1+vecXk(ip1*fNpar+1,0)*kY1+
     vecXk(ip1*fNpar+2,0)*kZ1+vecXk(ip1*fNpar+3,0)*dR1+vecXk(ip1*fNpar+4,0)/length1);
  //
  TMatrixD vecZk(nmeas,nmeas);
  TMatrixD measR(nmeas,nmeas);
  TMatrixD matHk(nmeas,fNelem);     // vector to mesurement
  TMatrixD vecYk(nmeas,nmeas);          // Innovation or measurement residual
  TMatrixD matHkT(fNelem,nmeas);    // helper matrix Hk transpose
  TMatrixD matSk(nmeas,nmeas);      // Innovation (or residual) covariance
  TMatrixD matKk(fNelem,nmeas);     // Optimal Kalman gain
  TMatrixD covXk2(fNelem,fNelem);   // helper matrix
  TMatrixD covXk3(fNelem,fNelem);   // helper matrix

  //
  vecZk(0,0) =dedx1/dedx0;          // dedx ratio
  measR(0,0) =(s0*s0+s1*s1);
  //
  // reset matHk
  for (Int_t iel=0;iel<fNelem;iel++) 
    for (Int_t ip=0;ip<nmeas;ip++) matHk(ip,iel)=0; 
 
  //
  matHk(0, ip0*fNpar+0) = corr0/(corr1*vecXk(ip0*fNpar+0,0));
  matHk(0, ip0*fNpar+1) = (vecXk(ip0*fNpar+0,0)*kY0)/corr1;
  matHk(0, ip0*fNpar+2) = (vecXk(ip0*fNpar+0,0)*kZ0)/corr1;
  matHk(0, ip0*fNpar+3) = (vecXk(ip0*fNpar+0,0)*dR0)/corr1;
  matHk(0, ip0*fNpar+4) = (vecXk(ip0*fNpar+0,0)/length0)/corr1;
  //
  //  matHk(0, ip1*fNpar+0) = ;
  //matHk(0, ip1*fNpar+1) = ;
  //matHk(0, ip1*fNpar+2) = ;
  //matHk(0, ip1*fNpar+3) = ;
  //matHk(0, ip1*fNpar+4) = ;
  //
  //
  vecYk(0,0) = vecZk(0,0)- corr0/corr1;               // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk2= (mat1-(matKk*matHk));
  covXk3 =  covXk2*covXk;          
  covXk = covXk3;  
  Double_t chi2 = vecYk(0,0)* vecYk(0,0)*matSk(0,0);
  //
  //
  //
  if (debug){ // dump debug info
    char chname[1000];
    sprintf(chname,"updatePair_%s",GetName());
    (*debug)<<chname<<
      // input parameters					
      "ip0="<<ip0<<
      "ip1="<<ip1<<
      "dedx0="<<dedx0<<
      "dedx1="<<dedx1<<      
      "s0="<<s0<<
      "s1="<<s1<<
      "kY0="<<kY0<<
      "kY1="<<kY1<<
      "kZ0="<<kZ0<<
      "kZ1="<<kZ1<<
      "dR0="<<dR0<<
      "dR1="<<dR1<<
      "chi2="<<chi2<<
      "corr0="<<corr0<<
      "corr1="<<corr1<<
      "length0="<<length0<<
      "length1="<<length1<<
      //
      "vecYk.="<<&vecYk<<
      "vecZk.="<<&vecZk<<
      "matHk.="<<&matHk<<
      "matHkT.="<<&matHkT<<
      "matSk.="<<&matSk<<
      "matKk.="<<&matKk<<
      "covXk2.="<<&covXk2<<
      "covXk.="<<&covXk<<
      "cov.="<<fCovariance<<
      "vecXk.="<<&vecXk<<
      "vec.="<<fState<<
      "\n";
  }
  if (chi2<kchi2Cut){
    //    (*fCovariance)=covXk;
    //(*fState)=vecXk;
  }
  else{
    printf("chi2=%f\n",chi2);
  }
}

void AliTPCkalmandEdx::UpdatedEdx(Int_t ip0,
				  Double_t dedx0, 
				  Double_t dedxRef, 
				  Double_t s0,
				  Double_t kY0,
				  Double_t kZ0,
				  Double_t dR0,
				  TTreeSRedirector *debug){
  //
  // update relative measurement
  // dEdx  = dEdxti = dEdxri*ci
  //
  // Model of correction ci:
  // dEdxti = dEdxri*ci = dEdxri*p0i*(1+p1i*kY+p2i*kZ+p3i*dR+p4i/sqrt(1+kY^2+kZ^2)
  const Double_t kchi2Cut = 3.*3.;
  // removing of "big outliers
  if (fInit<3) return AdddEdx(ip0,dedx0,dedxRef);
  if (TMath::Abs(dedxRef/dedx0-fSampleStat[ip0][0])>4*fSampleStat[ip0][1]) return;
  //
  //
  Int_t nmeas = 1;
  TMatrixD vecXk=*fState;           // X vector
  TMatrixD covXk=*fCovariance;      // X covariance
  TMatrixD &mat1=*fMat1;            // update covariance matrix

  Double_t length0 = TMath::Sqrt(1+kY0*kY0+kZ0*kZ0);
  //
  Double_t corr0 = vecXk(ip0*fNpar,0)*
    (1+vecXk(ip0*fNpar+1,0)*kY0+
     vecXk(ip0*fNpar+2,0)*kZ0+vecXk(ip0*fNpar+3,0)*dR0+vecXk(ip0*fNpar+4,0)/length0);


  //
  TMatrixD vecZk(nmeas,nmeas);
  TMatrixD measR(nmeas,nmeas);
  TMatrixD matHk(nmeas,fNelem);     // vector to mesurement
  TMatrixD vecYk(nmeas,nmeas);          // Innovation or measurement residual
  TMatrixD matHkT(fNelem,nmeas);    // helper matrix Hk transpose
  TMatrixD matSk(nmeas,nmeas);      // Innovation (or residual) covariance
  TMatrixD matKk(fNelem,nmeas);     // Optimal Kalman gain
  TMatrixD covXk2(fNelem,fNelem);   // helper matrix
  TMatrixD covXk3(fNelem,fNelem);   // helper matrix

  vecZk(0,0) =dedxRef/dedx0;
  measR(0,0) =s0*s0*dedxRef*dedxRef/(dedx0*dedx0);
  //
  // reset matHk
  for (Int_t iel=0;iel<fNelem;iel++) 
    for (Int_t ip=0;ip<nmeas;ip++) matHk(ip,iel)=0; 
 
  //
  //
  matHk(0, ip0*fNpar+0) = corr0/vecXk(ip0*fNpar,0);
  matHk(0, ip0*fNpar+1) = kY0*vecXk(ip0*fNpar,0);
  matHk(0, ip0*fNpar+2) = kZ0*vecXk(ip0*fNpar,0);
  matHk(0, ip0*fNpar+3) = dR0*vecXk(ip0*fNpar,0);
  matHk(0, ip0*fNpar+4) = vecXk(ip0*fNpar,0)/length0;
  //
  //
  vecYk(0,0) = vecZk(0,0)-corr0;               // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk2= (mat1-(matKk*matHk));
  covXk3 =  covXk2*covXk;          
  covXk = covXk3;  
  Double_t chi2 = vecYk(0,0)*vecYk(0,0)*matSk(0,0);
  //
  //
  //
  if (debug){ // dump debug info
    char chname[1000];
    sprintf(chname,"update_%s",GetName());
    (*debug)<<chname<<
      //
      // input parameters					
      "ip0="<<ip0<<
      "dedxRef="<<dedxRef<<
      "dedx0="<<dedx0<<
      "s0="<<s0<<
      "kY0="<<kY0<<
      "kZ0="<<kZ0<<
      "dR0="<<dR0<<
      "chi2="<<chi2<<
      "corr0="<<corr0<<
      "length0="<<length0<<
      //
      "vecYk.="<<&vecYk<<
      "vecZk.="<<&vecZk<<
      "matHk.="<<&matHk<<
      "matHkT.="<<&matHkT<<
      "matSk.="<<&matSk<<
      "matKk.="<<&matKk<<
      "covXk2.="<<&covXk2<<
      "covXk.="<<&covXk<<
      "cov.="<<fCovariance<<
      "vecXk.="<<&vecXk<<
      "vec.="<<fState<<
      "\n";
  }
  if (chi2<kchi2Cut){
    (*fCovariance)=covXk;
    (*fState)=vecXk;
  }else{
    printf("chi2=%f\n",chi2);
  }
}

void AliTPCkalmandEdx::AdddEdx(Int_t ip0,Double_t dedx0, Double_t dedxRef){
  //
  //
  //
  if (fCounter[ip0]>=fSampleSize) return;
  fSample[ip0][fCounter[ip0]]=dedxRef/dedx0;
  fCounter[ip0]++;
  if (fCounter[ip0]==fSampleSize){
    fSampleStat[ip0].ResizeTo(2);
    fSampleStat[ip0][0] = TMath::Median(fSampleSize,fSample[ip0].GetMatrixArray());
    fSampleStat[ip0][1] = TMath::RMS(fSampleSize,fSample[ip0].GetMatrixArray());
    fInit++;
    (*fState)(ip0*fNpar,0)    = fSampleStat[ip0][0];
    (*fCovariance)(ip0*fNpar,ip0*fNpar) = fSampleStat[ip0][1]*fSampleStat[ip0][1]/fSampleSize;
    Double_t med2  = fSampleStat[ip0][0]*fSampleStat[ip0][0];
    //
    //
    (*fCovariance)(ip0*fNpar+1,ip0*fNpar+1)=0.01*med2;
    (*fCovariance)(ip0*fNpar+2,ip0*fNpar+2)=0.01*med2;
    (*fCovariance)(ip0*fNpar+3,ip0*fNpar+3)=0.01*med2;
    (*fCovariance)(ip0*fNpar+4,ip0*fNpar+4)=0.01*med2;
  }
}




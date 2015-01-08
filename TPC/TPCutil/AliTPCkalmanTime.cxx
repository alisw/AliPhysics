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
  TPC Kalman filter implementation for time dependent variables:


The drift velocity and the gas gain are changing in time.
The drift velocity and gas gain is a function of many parameters, but not all of 
them are known. We assume that the most important parameters are pressure and temperature
and the influence of other parameters (gas composition, and electric field) are only 
slowly varying in time and can be expressed by smooth function $x_{off}(t)$:
\begin{equation}
x(t) = x_{off}(t)+k_N\frac{\Delta{P/T}}{P/T}	
\end{equation}
where x(t) is the parameter which we observe.
\begin{equation}
\begin{split}
x(t)=\frac{\Delta{G}}{G_0}	\\
x(t)=\frac{\Delta{v_d}}{v_{d0}}	
\end{split}
\end{equation}

Kalman filter parameters are following:
\begin{itemize}
\item State vector  ($x_{off}(t)$, $k_N$) at given time
\item Covariance matrix
\end{itemize}

Kalman filter implent following functions:
\begin{itemize}
\item Prediction - adding covariance element $\sigma_{xoff}$
\item Update state vector with new measurement vector ($x_t,\frac{\Delta{P/T}}{P/T}$)
\end{itemize}


*/

#include "AliTPCkalmanTime.h"
#include "TTreeStream.h"
#include "TRandom.h"


AliTPCkalmanTime::AliTPCkalmanTime():
  TNamed(),
  fState(0),
  fCovariance(0),
  fTime(0)
{
  //
  // Default constructor
  //
}

AliTPCkalmanTime::AliTPCkalmanTime(Double_t time, Double_t xoff, Double_t k, Double_t sigmaxoff, Double_t sigmak): 
  TNamed(),
  fState(0),
  fCovariance(0),
  fTime(0)
{
  //
  // Default constructor
  //
  Init(time,xoff,k,sigmaxoff,sigmak);
}


void AliTPCkalmanTime::Init(Double_t time, Double_t xoff, Double_t k, Double_t sigmaxoff, Double_t sigmak){
  //
  // Default constructor
  //
  fState = new TMatrixD(2,1);
  fCovariance = new TMatrixD(2,2);
  (*fState)(0,0)= xoff;  // offset
  (*fState)(1,0)= k;     // slope of the taylor
  fTime=time;
  (*fCovariance)(0,0)=sigmaxoff*sigmaxoff;
  (*fCovariance)(1,1)=sigmak*sigmak;
  (*fCovariance)(0,1)=0;
  (*fCovariance)(1,0)=0;
}


void AliTPCkalmanTime::Propagate(Double_t time, Double_t sigma, TTreeSRedirector *debug){
  //
  // Propagate the Kalman
  //
  if (!fCovariance) return;
  if (!fState) return;
  Double_t deltaT  =time-fTime;  //delta time - param2 is the current time
  Double_t sigmaT2 =(deltaT*deltaT)*sigma*sigma; 
  if (debug){
    (*debug)<<"matP"<<
      "time="<<time<<
      "fTime="<<fTime<<
      "sigmaT2="<<sigmaT2<<
      "cov.="<<fCovariance<<
      "\n";
  }
  (*fCovariance)(0,0)+=sigmaT2;
  fTime=time;  
}

void  AliTPCkalmanTime::Update(Double_t x, Double_t xerr, Double_t ptratio,TTreeSRedirector *debug){
  //
  //
  //
  static TMatrixD matHk(1,2);    // vector to mesurement
  static TMatrixD measR(1,1);    // measurement error 
  static TMatrixD matQk(2,2);    // prediction noise vector
  static TMatrixD vecZk(1,1);    // measurement
  //
  static TMatrixD vecYk(1,1);    // Innovation or measurement residual
  static TMatrixD matHkT(2,1);
  static TMatrixD matSk(1,1);    // Innovation (or residual) covariance
  static TMatrixD matKk(2,1);    // Optimal Kalman gain
  static TMatrixD mat1(2,2);     // update covariance matrix
  static TMatrixD covXk2(2,2);
  //
  //
  TMatrixD covXk(*fCovariance);    // X covariance 
  TMatrixD vecXk(*fState);         // X vector
  //
  vecZk(0,0) = x;                // current mesurement
  measR(0,0) = xerr*xerr;        // measurement error                      
  matHk(0,0)=1;  matHk(0,1)= ptratio; 
  vecYk = vecZk-matHk*vecXk;     // Innovation or measurement residual
  //
  //
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  vecXk += matKk*vecYk;                    //  updated vector 
  mat1(0,0)=1; mat1(1,1)=1; 
  mat1(1,0)=0; mat1(0,1)=0;
  covXk2= (mat1-(matKk*matHk));
  //
  //
  //
  covXk =  covXk2*covXk;    
  //
  if (debug){
    (*debug)<<"matU"<<
      "time="<<fTime<<
      "x="<<x<<
      "xerr="<<xerr<<
      "pt="<<ptratio<<
      "matHk.="<<&matHk<<
      "matHkT.="<<&matHkT<<
      "matSk.="<<&matSk<<
      "matKk.="<<&matKk<<
      "covXk2.="<<&covXk2<<
      "covXk.="<<&covXk<<
      "cov.="<<fCovariance<<
      "vecYk.="<<&vecYk<<
      "vecXk.="<<&vecXk<<
      "vec.="<<fState<<
      "\n";
  }
  (*fCovariance)=covXk;
  (*fState)=vecXk;
}

void AliTPCkalmanTime::TestMC(const char * fname){
  //
  // Test of the class
  /*
    Usage:
    AliTPCkalmanTime::TestMC("testKalman.root");
    TFile f("testKalman.root");
    
   */
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector(fname);
  //
  const Double_t kp0=0;
  const Double_t kp1=1;
  const Int_t    klength=10*24*60*60;       // 10 days mesurement
  const Double_t ksigmap0=0.001/(24*60*60.); // 0.005 change in one day
  const Int_t    deltat=5*60;                // 5 minutes step
  const Double_t kmessError=0.0005;
  AliTPCkalmanTime testKalman;

  for (Int_t iter=0; iter<100;iter++){
    Double_t sp0      = kp0+gRandom->Gaus(0,0.01);    // variable to estimate -offset 
    Double_t sp1      = kp1+gRandom->Gaus(0,0.2);    // variable to estimate slope
    Double_t cp0      = sp0;    // variable to estimate 
    Double_t cp1      = sp1;
    
    //
    testKalman.Init(0,cp0+gRandom->Gaus(0,0.05),cp1+gRandom->Gaus(0,0.2),0.05,0.2);
    Double_t dptratio= 0;
    for (Int_t itime=0; itime<klength; itime+=deltat){
      dptratio+=gRandom->Gaus(0,0.0005);
      cp0     +=gRandom->Gaus(0,ksigmap0*deltat);
      //
      Double_t vdrift  = cp0+dptratio*cp1+gRandom->Gaus(0,kmessError);
      testKalman.Propagate(itime,ksigmap0,pcstream);
      Double_t fdrift = (*(testKalman.fState))(0,0) + dptratio*(*(testKalman.fState))(1,0);
      (*pcstream)<<"drift"<<
	"iter="<<iter<<
	"time="<<itime<<
	"vdrift="<<vdrift<<
	"fdrift="<<fdrift<<
	"pt="<<dptratio<<
	"sp0="<<sp0<<
	"sp1="<<sp1<<
	"cp0="<<cp0<<
	"cp1="<<cp1<<
	"vecXk.="<<testKalman.fState<<
	"covXk.="<<testKalman.fCovariance<<
	"\n";
      testKalman.Update(vdrift,kmessError,dptratio,pcstream);
    }
  }
  delete pcstream;
}


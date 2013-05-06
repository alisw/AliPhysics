/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
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

// -----------------------------------------------------------------------
//  Definitions the mathematical functions used in the DiHadronPID
//  analysis.
// -----------------------------------------------------------------------
//  Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include <iostream>

#include "AliExternalTrackParam.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"

#include "AliFunctionsDiHadronPID.h"

using namespace std;

// -----------------------------------------------------------------------
AliFunctionsDiHadronPID::AliFunctionsDiHadronPID()

{

	// Constructor.

} 

// -----------------------------------------------------------------------
AliFunctionsDiHadronPID::~AliFunctionsDiHadronPID()

{

	// Destructor.

} 

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian1D(const Double_t *x, const Double_t *par) {

	//  - Gaussian, I is the integral -
	// f(x) = I/(Sqrt(2*pi)*sigma) * exp{-(x-mu)^2/2sigma^2}
	// par = {I,mu,sigma}

	Double_t norm = par[0]/(TMath::Sqrt(2.*TMath::Pi())*par[2]);
	Double_t gaussian = TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*par[2]*par[2]));

	return (norm*gaussian);

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian1D(const Double_t xx, const Double_t integral, const Double_t mu, const Double_t sigma, const Double_t binwidth) {

	// The other implementation should make use of this one.
	Double_t norm = (binwidth*integral)/(TMath::Sqrt(2.*TMath::Pi())*sigma);
	Double_t gaussian = TMath::Exp(-(xx-mu)*(xx-mu)/(2.*sigma*sigma));

	return (norm*gaussian);

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian1DTail(const Double_t *x,const  Double_t *par) {

	// Gaussian with exponential tail on the right, I is the integral.
	// For function definition see: FitFunctions.nb

	Double_t integral = par[0];
	Double_t mu = par[1];		// Top of the gaussian.
	Double_t kappa = par[1] + par[2];	// Point at which the gaussian becomes an exponential (w.r.t. to mu).
	Double_t sigma_x = par[3];

	if (mu >= kappa) return 0.; 	// Function becomes ill-defined.

	Double_t beta = sigma_x*sigma_x/(kappa-mu);
	Double_t BB = TMath::Exp( (kappa*kappa-mu*mu)/(2.*sigma_x*sigma_x) );
	Double_t norm1 = beta*TMath::Exp( -(mu-kappa)*(mu-kappa)/(2.*sigma_x*sigma_x) );
	Double_t norm2 = TMath::Sqrt(TMath::Pi()/2.)*sigma_x*TMath::Erfc( (mu-kappa)/(TMath::Sqrt2()*sigma_x) );
	Double_t norm = norm1 + norm2;

	Double_t funcleft = (integral/norm)*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/(2.*sigma_x*sigma_x));
	Double_t funcright = (integral/norm)*BB*TMath::Exp(-x[0]/beta);

	if (x[0] <= kappa) return funcleft;
	else return funcright;

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian1DTail(const Double_t xx, const Double_t integral, const Double_t mu, const Double_t sigma, const Double_t tail, const Double_t binwidth) {

	// Gaussian with exponential tail on the right, I is the integral.
	// For function definition see: FitFunctions.nb

	Double_t kappa = mu + tail;

	if (mu >= kappa) return 0.; 	// Function becomes ill-defined.

	Double_t beta = sigma*sigma/(kappa-mu);
	Double_t BB = TMath::Exp( (kappa*kappa-mu*mu)/(2.*sigma*sigma) );
	Double_t norm1 = beta*TMath::Exp( -(mu-kappa)*(mu-kappa)/(2.*sigma*sigma) );
	Double_t norm2 = TMath::Sqrt(TMath::Pi()/2.)*sigma*TMath::Erfc( (mu-kappa)/(TMath::Sqrt2()*sigma) );
	Double_t norm = norm1 + norm2;

	Double_t funcleft = binwidth * (integral/norm)*TMath::Exp(-(xx-mu)*(xx-mu)/(2.*sigma*sigma));
	Double_t funcright = binwidth * (integral/norm)*BB*TMath::Exp(-xx/beta);

	if (xx <= kappa) return funcleft;
	else return funcright;

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian2D(const Double_t xx, const Double_t yy, const Double_t integral, 
	const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
	const Double_t binwidthx, const Double_t binwidthy) {

	// 2D Gaussian.
	Double_t GaussianX = Gaussian1D(xx, 1., mux, sigmax, binwidthx);
	Double_t GaussianY = Gaussian1D(yy, 1., muy, sigmay, binwidthy);

	return integral * GaussianX * GaussianY; 

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian2DTailX(const Double_t xx, const Double_t yy, const Double_t integral, 
	const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
	const Double_t tailx, const Double_t binwidthx, const Double_t binwidthy) {

	// 2D Gaussian with exponential tail in X direction.
	Double_t GaussianTailX = Gaussian1DTail(xx, 1., mux, sigmax, tailx, binwidthx);
	Double_t GaussianY = Gaussian1D(yy, 1., muy, sigmay, binwidthy);

	return integral * GaussianTailX * GaussianY;

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian2DTailY(const Double_t xx, const Double_t yy, const Double_t integral, 
	const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
	const Double_t taily, const Double_t binwidthx, const Double_t binwidthy) {

	// 2D Gaussian with exponential tail in Y direction.
	Double_t GaussianX = Gaussian1D(xx, 1., mux, sigmax, binwidthx);
	Double_t GaussianTailY = Gaussian1DTail(yy, 1., muy, sigmay, taily, binwidthy);

	return integral * GaussianX * GaussianTailY;

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Gaussian2DTailXY(const Double_t xx, const Double_t yy, const Double_t integral, 
	const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
	const Double_t tailx, const Double_t taily, const Double_t binwidthx, const Double_t binwidthy) {

	// 2D Gaussian with exponential tail in X- and Y direction.
	Double_t GaussianTailX = Gaussian1DTail(xx, 1., mux, sigmax, tailx, binwidthx);
	Double_t GaussianTailY = Gaussian1DTail(yy, 1., muy, sigmay, taily, binwidthy);

	return integral * GaussianTailX * GaussianTailY;

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::Exponent(const Double_t *x, const Double_t *par) {

	// f(x) = A * Exp(bx)
	// par = {A,b}

	return par[0]*TMath::Exp(par[1]*x[0]); 

}

// -----------------------------------------------------------------------
//  COMBINED FUNCTIONS
// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::SimpleTOFfit(const Double_t *x, const Double_t *par) {

	// Signal fitted with a Gaussian, mismatches by an exponent.
	return (Gaussian1D(x,&par[0]) + Exponent(x,&par[3]));

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::SimpleTOFfitWithTail(const Double_t *x, const Double_t *par) {

	// Signal fitted with a Gaussian with a tail, mismatches by an exponent.
	return (Gaussian1D(x,&par[0]) + Exponent(x,&par[4]));

}

// -----------------------------------------------------------------------
//  PENALTY FUNCTIONS
// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::PolyPenalty(const Double_t xx, const Double_t center, const Double_t flatwidth, const Int_t polyorder) {

	// Penalty function for a chi^2 fit. The function is defined as:
	// 1 											for |xx - center| < flatwidth,
	// (|xx - center| - flatwidth) ^ polyorder		for |xx - center| > flatwidth.

	Double_t fx = 1.;
	if (TMath::Abs(xx - center) > flatwidth) {
		fx = TMath::Power( (TMath::Abs(xx - center) - flatwidth), polyorder ) + 1.;
	}

	return fx;

}

// -----------------------------------------------------------------------
TCanvas* AliFunctionsDiHadronPID::TestPolyPenalty(const Double_t range, const Double_t center, const Double_t flatwidth, const Int_t polyorder) {

	// Creates an example of the TestPolyPenalty function.
	TF1* tf = new TF1("tf",Form("AliFunctionsDiHadronPID::PolyPenalty(x,[0],[1],%i)",polyorder),-range,range);
	tf->SetParameters(center,flatwidth);
	TCanvas* cvs = TCanvas::MakeDefCanvas();
	tf->Draw();

	return cvs;

}

// -----------------------------------------------------------------------
//  PID Expected signal functions.
// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::TOFExpTime(const Double_t pT, const Double_t eta, const Double_t mass) {

	// For description see ../Documents/TOFtime.tex

	Double_t AA = (2. * pT) / ( Charge() * BTPC() * GeVperkg() );
	Double_t BB = TMath::ASin( (Charge() * BTPC() * 0.01 * RTOF() * GeVperkg() ) / (2. * pT * C()) ); 
	Double_t CC = TMath::Sqrt( mass*mass/(pT*pT) + TMath::CosH(eta)*TMath::CosH(eta) );

	return (1.e12*AA*BB*CC);   // Time returned in ps.

}

// -----------------------------------------------------------------------
Double_t AliFunctionsDiHadronPID::TPCExpdEdX(const Double_t pT, const Double_t eta, const Double_t mass) {

	// Not so neat solution, however the easiest for now.

	// Prameters taken from the constructor of AliTPCPIDResponse:
	Double_t MIP = 50.;
	Double_t Kp[5] = {0.0283086, 2.63394e+01, 5.04114e-11, 2.12543, 4.88663};

	Double_t betaGamma = TMath::Abs( (pT * TMath::CosH(eta)) / mass );

	// Implementation as in AliTPCPIDResponse.
	return MIP * AliExternalTrackParam::BetheBlochAleph(betaGamma,Kp[0],Kp[1],Kp[2],Kp[3],Kp[4]);

}

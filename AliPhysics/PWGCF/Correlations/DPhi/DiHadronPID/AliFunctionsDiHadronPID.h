#ifndef ALIFUNCTIONSDIHADRONPID_H
#define ALIFUNCTIONSDIHADRONPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

#include "TCanvas.h"
#include "TMath.h"

//class TCanvas;

class AliFunctionsDiHadronPID {

public:
	AliFunctionsDiHadronPID();

protected:
	~AliFunctionsDiHadronPID();

public:
	// Math.
	static Int_t Power(Int_t base, Int_t power);

	// Natural constants.
	static Double_t Charge() {return 1.60217646e-19;}	// (C)
	static Double_t C() {return 2.99792458e+8;}			// (m/s)
	static Double_t Mpion() {return 0.13957018;}		// (GeV/c^2)
	static Double_t Mkaon() {return 0.493667;}			// (GeV/c^2)
	static Double_t Mproton() {return 0.938272046;}		// (GeV/c^2)
	static Double_t Mdeuteron() {return 2.01410178*GeVperu();}	// (GeV/c^2)
	static Double_t M(Int_t species) {
		switch(species) {
			case 0:return Mpion();
			case 1:return Mkaon();
			case 2:return Mproton();
			case 3:return Mdeuteron();
			default:return -999.; 
		}
	}

	// Conversions.
	static Double_t GeVperu() {return 0.931494061;}		// (GeV/c^2) per u
	static Double_t GeVperkg() {return 5.608524e+26;} 	// (GeV/c^2) per kg

	// Detector paramters.
	static Double_t RTOF() {return 385.;}				// Radius of TOF (cm).
	static Double_t BTPC() {return 0.5;}				// Magnetic field in TPC (T = kg C^-1 s^-1).

	// Fit Functions.
	static Double_t Gaussian1D(Double_t xx, Double_t integral, Double_t mu, Double_t sigma, Double_t binwidth = 1.);
	static Double_t Gaussian1DTail(Double_t xx, Double_t integral, Double_t mu, Double_t sigma, Double_t tail, Double_t binwidth = 1.);

	static Double_t Gaussian2D(Double_t xx, Double_t yy, Double_t integral, 
		Double_t mux, Double_t muy, Double_t sigmax, Double_t sigmay, 
		Double_t binwidthx = 1., Double_t binwidthy = 1.);

	static Double_t Gaussian2DTailX(Double_t xx, Double_t yy, Double_t integral, 
		Double_t mux, Double_t muy, Double_t sigmax, Double_t sigmay, 
		Double_t tailx, Double_t binwidthx = 1., Double_t binwidthy = 1.);

	static Double_t Gaussian2DTailY(Double_t xx, Double_t yy, Double_t integral, 
		Double_t mux, Double_t muy, Double_t sigmax, Double_t sigmay, 
		Double_t taily, Double_t binwidthx = 1., Double_t binwidthy = 1.);

	static Double_t Gaussian2DTailXY(Double_t xx, Double_t yy, Double_t integral, 
		Double_t mux, Double_t muy, Double_t sigmax, Double_t sigmay, 
		Double_t tailx, Double_t taily, Double_t binwidthx = 1., Double_t binwidthy = 1.);

	// Penalty Functions.
	static Double_t PolyPenalty(Double_t xx, Double_t center, Double_t flatwidth, Int_t polyorder);
	static TCanvas* TestPolyPenalty(Double_t range = 3., Double_t center = 1., Double_t flatwidth = 1., Int_t polyorder = 3);

	// PID Expected signal functions.
	static Double_t TOFExpTime(Double_t pT, Double_t eta, Double_t mass);
	static Double_t TPCExpdEdX(Double_t pT, Double_t eta, Double_t mass);

	// Standard Functions.
	static Double_t Exponent(Double_t xx, Int_t sign, Double_t p0, Double_t p1) {return (sign*TMath::Exp(p0 + xx*p1));}
	static Double_t Poly1(Double_t xx, Double_t p0, Double_t p1) {return (p0 + p1*xx);}
	static Double_t Poly2(Double_t xx, Double_t p0, Double_t p1, Double_t p2) {return (p0 + p1*xx + p2*xx*xx);}
	static Double_t Poly3(Double_t xx, Double_t p0, Double_t p1, Double_t p2, Double_t p3) {return (p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx);}
	static Double_t Poly4(Double_t xx, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4) {return (p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx + p4*xx*xx*xx*xx);}
	static Double_t Poly5(Double_t xx, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5) {return (p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx + p4*xx*xx*xx*xx + p5*xx*xx*xx*xx*xx);}				

};

#endif

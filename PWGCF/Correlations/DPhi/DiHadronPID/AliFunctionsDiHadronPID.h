#ifndef ALIFUNCTIONSDIHADRONPID_H
#define ALIFUNCTIONSDIHADRONPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

class TCanvas;

class AliFunctionsDiHadronPID {

public:
	AliFunctionsDiHadronPID();

protected:
	~AliFunctionsDiHadronPID();

public:
	// Natural constants.
	static Double_t Charge() {return 1.60217646e-19;}	// (C)
	static Double_t C() {return 2.99792458e+8;}			// (m/s)
	static Double_t Mpion() {return 0.13957018;}		// (GeV/c^2)
	static Double_t Mkaon() {return 0.493667;}			// (GeV/c^2)
	static Double_t Mproton() {return 0.938272046;}		// (GeV/c^2)
	static Double_t Mdeuteron() {return 2.01410178*GeVperu();}	// (GeV/c^2)
	static Double_t M(const Int_t species) {
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
	static Double_t Gaussian1D(const Double_t xx, const Double_t integral, const Double_t mu, const Double_t sigma, const Double_t binwidth = 1.);
	static Double_t Gaussian1DTail(const Double_t xx, const Double_t integral, const Double_t mu, const Double_t sigma, const Double_t tail, const Double_t binwidth = 1.);

	static Double_t Gaussian2D(const Double_t xx, const Double_t yy, const Double_t integral, 
		const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
		const Double_t binwidthx = 1., const Double_t binwidthy = 1.);

	static Double_t Gaussian2DTailX(const Double_t xx, const Double_t yy, const Double_t integral, 
		const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
		const Double_t tailx, const Double_t binwidthx = 1., const Double_t binwidthy = 1.);

	static Double_t Gaussian2DTailY(const Double_t xx, const Double_t yy, const Double_t integral, 
		const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
		const Double_t taily, const Double_t binwidthx = 1., const Double_t binwidthy = 1.);

	static Double_t Gaussian2DTailXY(const Double_t xx, const Double_t yy, const Double_t integral, 
		const Double_t mux, const Double_t muy, const Double_t sigmax, const Double_t sigmay, 
		const Double_t tailx, const Double_t taily, const Double_t binwidthx = 1., const Double_t binwidthy = 1.);

	// FUNCTIONS THAT ARE NOT BEING USED ANYMORE.
	static Double_t Gaussian1D(const Double_t *x, const Double_t *par);
	static Double_t Gaussian1DTail(const Double_t *x, const Double_t *par);	
	static Double_t Exponent(const Double_t *x, const Double_t *par);
	static Double_t SimpleTOFfit(const Double_t *x, const Double_t *par);
	static Double_t SimpleTOFfitWithTail(const Double_t *x, const Double_t *par);

	// Penalty Functions.
	static Double_t PolyPenalty(const Double_t xx, const Double_t center, const Double_t flatwidth, const Int_t polyorder);
	static TCanvas* TestPolyPenalty(const Double_t range = 3., const Double_t center = 1., const Double_t flatwidth = 1., const Int_t polyorder = 3);

	// PID Expected signal functions.
	static Double_t TOFExpTime(const Double_t pT, const Double_t eta, const Double_t mass);
	static Double_t TPCExpdEdX(const Double_t pT, const Double_t eta, const Double_t mass);

};

#endif

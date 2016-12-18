#ifndef ALIANACOULWAVE_H
#define ALIANACOULWAVE_H


//////////////////////////////////////////////////////////////////////////////////////
//
// CoulombWave class : Description
//
//  This class calculates the full coulomb correction factor, originally developed by
//  some NA44's collaborators in FORTARAN.
//  This function calculates the coulomb wave function and the probability of some
//  random points according to the gaussian source radii and q values.
//  The number of random points is set to 20 initially, but you can change it by using
//  Set_Gampt function. As the number increase, the calculation become accurate.
//  But it will take for long time to calculate.
//
//  The factor is obtained by Coul_Wave(Mass, Rinv, Rts, Rto, Rl, Qinv, Qts, Qto, Ql).
//  Mass               : Particle ideal mass.
//  Rinv, Rts, Rto, Rl : 1-d source radius, and ource radii in side, out, long directions.
//  Qinv, Qts, Qto, Ql : q invariant and q_side,out,long values.
//  *note: this function is revised as include Rinv. (01/26/03 A.E.)
//
//  This factor is assumed to be applied to actual pairs, so if you apply this factor
//  to background pairs, the factor should be 1/Coul_Wave().
//
//
//////////////////////////////////////////////////////////////////////////////////////

//#include <fstream>
#include <iostream>
#include <TRandom.h>
#include <TMath.h>
#include <complex>
//#include "math.h"

class AliAnaCoulWave {

	public:

		AliAnaCoulWave();				//default constructor
		virtual	~AliAnaCoulWave();

		void Set_Gampt(Int_t nopt) { fGam_pt=nopt; std::cout<<"GAM_PT = "<<fGam_pt<<std::endl; }
		void Set_Qtype(Int_t qtyp) { fQtype=qtyp; std::cout<<"Q Type = "<<fQtype<<std::endl; }
		Double_t Coul_Gamow(Double_t Mass, Double_t Qinv);
		Double_t Coul_Wave(const Double_t Mass, Double_t Rinv_i, Double_t Rts_i, Double_t Rto_i, Double_t Rl_i,
				Double_t Qinv, Double_t Qts, Double_t Qto, Double_t Ql);


		std::complex<float> CLog(std::complex<float> a_in);
		std::complex<float> Conjg(std::complex<float> a_in);
		std::complex<float> CPower(std::complex<float> Xin, std::complex<float> Yin);
		std::complex<float> CGamma(std::complex<float> a_in);
		std::complex<float> CHyper(std::complex<float> Ain,std::complex<float> Bin,std::complex<float> Zin,std::complex<float> f1in,std::complex<float> f2in,std::complex<float> hpin[20],std::complex<float> hp1in[4], std::complex<float> hp2in[4] );

	private:

		Int_t fGam_pt;
		Int_t	fQtype;
};


#endif /*__COULOMBWAVE_HH__*/

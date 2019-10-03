#include "Riostream.h"  //needed as include
#include "TMath.h"   //needed as include
#include "TProfile.h"   //needed as include
#include "TRandom.h"
#include "complex"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TList.h"
#include "TVector2.h"
#include "TSystem.h"

#include "AliAnaCoulWave.h"

ClassImp(AliAnaCoulWave)

	//----------------------------------------------------------------

	AliAnaCoulWave::AliAnaCoulWave():

		//fGam_pt = 20;
		//fQtype  = 0; //0: for QM analysis, 1: qinv based,  2: 3d-q based,  3: partialu 3d-q based
	fGam_pt(20)
,fQtype(0)

{
	//Constructor
}

AliAnaCoulWave::~AliAnaCoulWave()
{
	//Destructor
}

Double_t AliAnaCoulWave::Coul_Gamow(Double_t Mass, Double_t Qinv)
{

	Double_t eta = 2.0*3.141579*Mass/(137.036*Qinv);
	Double_t v = (exp(eta)-1)/eta;
	if(v >= 10.0) v = 10.0;
	return v;

}


Double_t AliAnaCoulWave::Coul_Wave(const Double_t Mass, Double_t Rinv_i, Double_t Rts_i, Double_t Rto_i, Double_t Rl_i,
		Double_t Qinv, Double_t Qts, Double_t Qto, Double_t Ql)
{

	Double_t Qtmp = 1.0;
	if( fQtype==0 || fQtype==1 || fQtype==3 ) { Qtmp = Qinv; }
	else { Qtmp = sqrt(fabs(Qts*Qts + Qto*Qto + Ql*Ql)); }


	Double_t seta = Mass/(137.036*Qtmp);
	std::complex<float> C_one(1.0,0.0);
	std::complex<float> Argtmp(1.0, seta);
	std::complex<float> Couly = CGamma(Argtmp);
	Couly = Couly*expf(-0.5*seta*3.141579);

	std::complex<float> A0(0.0, -seta);
	std::complex<float> B0(1.0, 0.0);
	std::complex<float> Ctop(1.0,0.0);
	std::complex<float> Cbot(1.0,0.0);
	std::complex<float> Chype[20] = { std::complex<float>( 0.0, 0.0 ) };
	std::complex<float> A = A0 - C_one;
	std::complex<float> B = B0 - C_one;
	std::complex<float> CJ( 0.0, 0.0 );

	for(Int_t ij=0; ij<20; ij++) {
		A = A + C_one;
		B = B + C_one;
		CJ = ij + 1.0;
		Ctop = Ctop*A;
		Cbot = Cbot*B*CJ;
		Chype[ij] = Ctop/Cbot;
		//    cout << "Input Chype[" << ij << "] = " << Chype[ij] << endl;
	}

	std::complex<float> C1top1(1.0,0.0);
	std::complex<float> C1top2(1.0,0.0);
	std::complex<float> C2top1(1.0,0.0);
	std::complex<float> C2top2(1.0,0.0);
	std::complex<float> Cbot0(1.0,0.0);
	std::complex<float> Chype1[4] = { std::complex<float>( 0.0, 0.0 ) };
	std::complex<float> Chype2[4] = { std::complex<float>( 0.0, 0.0 ) };

	for(Int_t ij=0; ij<4; ij++) {
		CJ = ij + 1.0;
		C1top1 = C1top1*(CJ+A0-C_one);
		C2top1 = C2top1*(CJ-A0);
		C1top2 = C1top2*(CJ-B0+A0);
		C2top2 = C2top2*(CJ+B0-A0-C_one);
		Cbot0 = Cbot0*CJ;
		Chype1[ij] = (C1top1*C1top2)/Cbot0;
		Chype2[ij] = (C2top1*C2top2)/Cbot0;
	}
	std::complex<float> Cfact1 = CGamma(B0)/CGamma(B0-A0);
	std::complex<float> Cfact2 = CGamma(B0)/CGamma(A0);

	// cout << "\nCfact1 = " << Cfact1 << "  Cfact2 = " << Cfact2 << endl;

	std::complex<float> Rctot(0.0,0.0);
	std::complex<float> Rcnon(0.0,0.0);

	for(Int_t ip=0; ip<fGam_pt; ip++) {

		Double_t R, ZK;
		if(fQtype==1) {
			Double_t Rinv = sqrt(2.0)*Rinv_i*gRandom->Gaus();
			R  = sqrt(fabs(Rinv*Rinv));
			ZK = Rinv*Qinv/Qtmp;
		} else {
			Double_t Rts = sqrt(2.0)*Rts_i*gRandom->Gaus();
			Double_t Rto = sqrt(2.0)*Rto_i*gRandom->Gaus();
			Double_t Rl  = sqrt(2.0)*Rl_i*gRandom->Gaus();
			R  = sqrt(fabs(Rts*Rts + Rto*Rto + Rl*Rl));
			if(fQtype==3) { ZK = (Rts*Qts + Rto*Qto + Rl*Ql)/sqrt(fabs(Qts*Qts + Qto*Qto + Ql*Ql)); }
			else         { ZK = (Rts*Qts + Rto*Qto + Rl*Ql)/Qtmp; }
		}

		std::complex<float> Z1(0.0,Qtmp/2*(R-ZK)/0.197327);
		std::complex<float> Z2(0.0,Qtmp/2*(R+ZK)/0.197327);
		std::complex<float> argtmp1(0.0,Qtmp/2*(ZK)/0.197327);
		std::complex<float> argtmp2(0.0,Qtmp/2*(-ZK)/0.197327);

		/*
			 if (ip<3) {
			 cout << "***************Nomber " << ip+1 << " *************************" << endl;
			 cout << "Rts = " << Rts << "  Rto = " << Rto << "  Rl = " << Rl << "  R = " << R << endl;
			 cout << "Qts = " << Qts << "  Qto = " << Qto << "  Ql = " << Ql << "  Qinv = " << Qinv << endl;
			 cout << "ZK = " << ZK << "  Z1 = " << Z1 << "  Z2 = " << Z2 << endl;
			 cout << "argtmp1 = " << argtmp1 << "  argtmp2 = " << argtmp2 << endl;
			 cout << "exp(argtmp1) = " << exp(argtmp1) << "  exp(argtmp2) = " << exp(argtmp2) << endl;
			 }
		 */

		std::complex<float> Cphi1 = Couly*exp(argtmp1)*CHyper(A0,B0,Z1,Cfact1,Cfact2,Chype,Chype1,Chype2);
		std::complex<float> Cphi2 = Couly*exp(argtmp2)*CHyper(A0,B0,Z2,Cfact1,Cfact2,Chype,Chype1,Chype2);
		std::complex<float> CphiS = (float)(1.0/sqrt(2.0))*(Cphi1+Cphi2);
		Rctot = Rctot + CphiS*Conjg(CphiS);
		std::complex<float> Cphi01 = exp(argtmp1);
		std::complex<float> Cphi02 = exp(argtmp2);
		std::complex<float> Cphi0S = (float)(1.0/sqrt(2.0))*(Cphi01+Cphi02);
		Rcnon = Rcnon + Cphi0S*Conjg(Cphi0S);

		/*
			 if (ip<3) {
			 cout << "Cphi1 = " << Cphi1 << "  Cphi2 = " << Cphi2 << endl;
			 cout << "CphiS = " << CphiS << "  Rctot = " << Rctot << endl;
			 cout << "Cphi01 = " << Cphi01 << "  Cphi02 = " << Cphi02 << endl;
			 cout << "Cphi0S = " << Cphi0S << "  Rcnon = " << Rcnon << endl;
			 }
		 */

	}

	Double_t v = Rcnon.real() / Rctot.real();

	if( v >= 10.0 ) v = 10.0;
	return v;
}


///////////////////////////////////////////////////////////
////                                                  /////
////  Member functions for mathmatical calcutlation   /////
////                                                  /////
///////////////////////////////////////////////////////////

std::complex<float> AliAnaCoulWave::CLog(std::complex<float> a_in)
{
	//  float out_imag = atan(a_in.imag()/a_in.real());
	//  float out_real = log(a_in.real()/cos(out_imag));
	Double_t out_imag = TMath::ATan2(a_in.imag(),a_in.real());
	std::complex<float> a_inabs = abs(a_in);
	Double_t out_real = log(a_inabs.real());
	std::complex<float> v(out_real, out_imag);

	return v;
}

std::complex<float> AliAnaCoulWave::Conjg(std::complex<float> a_in)
{
	std::complex<float> v(a_in.real(), -a_in.imag());

	return v;
}

std::complex<float> AliAnaCoulWave::CPower(std::complex<float> Xin,   std::complex<float> Yin)
{
	std::complex<float> v = exp(Yin*CLog(Xin));

	return v;
}

std::complex<float> AliAnaCoulWave::CGamma(std::complex<float> a_in)
{

	Float_t c[7] = { 2.5066282746310005, 76.18009172947146, -86.50532032941677
		,24.01409824083091,  -1.231739572450155, 0.1208650973866179e-2
			,-0.5395239384953e-5};

	std::complex<float> x( 0.0, 0.0 );
	std::complex<float> y( 0.0, 0.0 );
	std::complex<float> v( 0.0, 0.0 );
	if (a_in.real() < 1.0) {;
		x = a_in + (Float_t)1.0;
		y = a_in + (Float_t)1.0;
	}
	else {
		x = a_in;
		y = a_in;
	}

	std::complex<float> tmp = x + (float)5.5;
	tmp = tmp - (x+(float)0.5)*CLog(tmp);
	std::complex<float> ser(1.000000000190015,0.0);
	for (Int_t iq=1; iq<7; iq++) {
		y = y + (float)1.0;
		ser = ser + c[iq]/y;
	}

	if (a_in.real() < 1.0) {
		v = exp(-tmp + CLog(c[0]*ser/x))/a_in;
	}
	else {
		v = exp(-tmp + CLog(c[0]*ser/x));
	}

	return v;
}

std::complex<float> AliAnaCoulWave::CHyper(std::complex<float> Ain,   std::complex<float> Bin,  std::complex<float> Zin,
		std::complex<float> f1in,  std::complex<float> f2in, std::complex<float> hpin[20],
		std::complex<float> hp1in[4], std::complex<float> hp2in[4] )
{

	std::complex<float> v(0.0,0.0);
	std::complex<float> ZAbs = abs(Zin);
	//  cout << "ZAbs = " << ZAbs << endl;
	if( ZAbs.real() < 10.0 ) {
		std::complex<float> Cf1(1.0,0.0);
		std::complex<float> Czj(1.0,0.0);
		std::complex<float> CdAbs(0.0,0.0);
		std::complex<float> Cdelcf(0.0,0.0);
		for(int ij=0; ij<20; ij++) {
			Czj = Czj*Zin;
			Cdelcf = hpin[ij]*Czj;
			Cf1 = Cf1 + Cdelcf;
			CdAbs = abs(Cdelcf);
			if( CdAbs.real() < 0.001 ) break;
		}
		v = Cf1;
	}
	else {
		std::complex<float> Cf1(1.0,0.0);
		std::complex<float> Cf2(1.0,0.0);
		std::complex<float> Czj(1.0,0.0);
		std::complex<float> Czjm(1.0,0.0);
		std::complex<float> CW1(0.0,0.0);
		std::complex<float> CW2(0.0,0.0);
		for (int ij=0; ij<4; ij++) {
			Czj = Czj*Zin;
			Czjm = Czjm*(-Zin);
			Cf1 = Cf1 + hp1in[ij]/Czjm;
			Cf2 = Cf2 + hp2in[ij]/Czj;
			//      cout << "Cf1 & Cf2 = " << Cf1 << "  " << Cf2 << endl;
		}
		std::complex<float> Zminus = -Zin;
		std::complex<float> Aminus = -Ain;
		std::complex<float> ABminus = Ain - Bin;
		CW1 = f1in*(CPower(Zminus,Aminus))*Cf1;
		CW2 = f2in*(CPower(Zin,ABminus))*exp(Zin)*Cf2;
		/*
			 cout << "f1in & f2in = " << f1in << "  " << f2in << endl;
			 cout << "-Zin & -Ain = " << -Zin << "  " << -Ain << endl;
			 cout << "CP1 & CP2 = " << CPower(Zminus,Aminus) << "  " << CPower(Zin,ABminus) << endl;
			 cout << "CW1 & CW2 = " << CW1 << "  " << CW2 << endl;
		 */
		v = CW1+CW2;
	}

	return v;
}

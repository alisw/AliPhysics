/******************************************************************************
 *                      T H E R M I N A T O R                                 *
 *                   THERMal heavy-IoN generATOR                              *
 *                           version 1.0                                      *
 *                                                                            *
 * Authors of the model: Wojciech Broniowski, Wojciech.Broniowski@ifj.edu.pl, *
 *                       Wojciech Florkowski, Wojciech.Florkowski@ifj.edu.pl  *
 * Authors of the code:  Adam Kisiel, kisiel@if.pw.edu.pl                     *
 *                       Tomasz Taluc, ttaluc@if.pw.edu.pl                    *
 * Code designers: Adam Kisiel, Tomasz Taluc, Wojciech Broniowski,            *
 *                 Wojciech Florkowski                                        *
 *                                                                            *
 * For the detailed description of the program and furhter references         * 
 * to the description of the model plesase refer to: nucl-th/0504047,         *
 * accessibile at: http://www.arxiv.org/nucl-th/0504047                       *
 *                                                                            *
 * Homepage: http://hirg.if.pw.edu.pl/en/therminator/                         *
 *                                                                            *
 * This code can be freely used and redistributed. However if you decide to   *
 * make modifications to the code, please contact the authors, especially     *
 * if you plan to publish the results obtained with such modified code.       *
 * Any publication of results obtained using this code must include the       *
 * reference to nucl-th/0504047 and the published version of it, when         *
 * available.                                                                 *
 *                                                                            *
 *****************************************************************************/
#include <fstream>
#include "Integrator.h"
#include "ParticleType.h"
#include "Particle.h"
#include <TMath.h>
#include <TFile.h>

extern ReadPar *sRPInstance;
extern int      sTables;
extern int      sModel;

Integrator::Integrator(int aNpart)
{
  kFmToGev  = 0.197326960277;				/*MCH updated: kFmToGev  = 0.197;*/
  ReadParameters();

  char *tHash;
  tHash = ParameterHash();

  PRINT_MESSAGE("Hash for these parameters is: " << tHash);
  
  mNPart = aNpart;
  kTwoPi2 = TMath::Pi()*TMath::Pi()*2*2;		/*MCH*/
  kTwoPi3 = TMath::Pi()*TMath::Pi()*TMath::Pi()*2*2*2;
  mRandom = new TRandom2();
  //  mRandom->SetSeed2(41321, 8457);
  mRandom->SetSeed(41321);

  mFOHS = new Hypersurface(mFOHSlocation.Data());				/*MCH*/

  free (tHash);
}

double Integrator::CalcBE(double aX)
{
  if (aX>200.0) return 0.0;
  return 1.0/(TMath::Exp(aX)-1);
}

double Integrator::CalcFD(double aX)
{
  if (aX>200.0) return 0.0;
  return 1.0/(TMath::Exp(aX)+1);
}

double Integrator::Calka(double aMass, double aMiu, 
			 double *aRap, double *aPt, double *aPhiP, 
			 double *aAlfaP, double *aRho, double *aPhiS, double *aTime, double aSpin,
			 int aKeepPos)
{
  // Generate momentum components
  (*aRap) = mRandom->Rndm() * mRapRange - mRapRange/2.0;
  double tZet   = mRandom->Rndm();
  (*aPt) = tZet/(1-tZet);
  (*aPhiP) = mRandom->Rndm() * TMath::Pi() * 2.0;
  
  // Generate position components
  if (!aKeepPos) {
    (*aAlfaP) = mRandom->Rndm()*mAlfaRange - mAlfaRange/2.0;
    (*aRho) = mRandom->Rndm()*mRhoMax;
    (*aPhiS) = mRandom->Rndm()*TMath::Pi()*2;
    if (sModel == 8) {
      (*aTime) = mTau - mBWDelay*TMath::Log(mRandom->Rndm());
    }
  }
  
  double tFpod = 0.0;
  if ((sModel == 0) || (sModel == 3)) { // Single FreezeOut
    double tPU = (TMath::Sqrt(aMass*aMass+(*aPt)*(*aPt))*
		 cosh((*aAlfaP)-(*aRap))*
		 TMath::Sqrt(1+(((*aRho)*(*aRho))/(mTau*mTau)))
		 -(*aPt)*cos((*aPhiS)-(*aPhiP))*(*aRho)/mTau);
    if (fabs(floor(aSpin) - aSpin)< 0.01) 
      tFpod = 1/((1-tZet)*(1-tZet))*mTau*(*aPt)*(*aRho)*tPU*CalcBE((tPU- aMiu)/mTemp)/(kTwoPi3);
    else
      tFpod = 1/((1-tZet)*(1-tZet))*mTau*(*aPt)*(*aRho)*tPU*CalcFD((tPU- aMiu)/mTemp)/(kTwoPi3);
  }
/*MCH begin*/
  else if (sModel == 10) { // Lhyquid3D
    double tZeta  = mRandom->Rndm()*0.5*TMath::Pi();		// angle in rho-t plane
    double taHS   = mFOHS->fahs  ((*aPhiS),tZeta);		// velocity angle; 0.0=radial flow
    double tvHS   = mFOHS->fvhs  ((*aPhiS),tZeta);		// velocity
    double tdHS   = mFOHS->fdhs  ((*aPhiS),tZeta)/kFmToGev;	// distance in rho-t plane           [GeV^-1]
    double tDpdHS = mFOHS->fDpdhs((*aPhiS),tZeta)/kFmToGev;	// distance derivative over (*aPhiS) [GeV^-1]
    double tDzdHS = mFOHS->fDzdhs((*aPhiS),tZeta)/kFmToGev;	// distance derivative over tZeta    [GeV^-1]
    double tTemp  = mFOHS->TFO;					// freeze-out temparature
    double ttau0  = mFOHS->tau0/kFmToGev;			// tau 0
    double tMt    = TMath::Hypot(aMass, (*aPt));		// transverse mass
    (*aRho)       = tdHS*cos(tZeta);				// rho
    double tdPt   = 1.0/((1-tZet)*(1-tZet));			// dPt
    (*aTime)      = (ttau0+tdHS*sin(tZeta))*cosh(*aAlfaP);	// t

    double tPU    = 1.0/sqrt(1.0-tvHS*tvHS)*(tMt*cosh((*aRap)-(*aAlfaP))-(*aPt)*tvHS*cos((*aPhiS)-(*aPhiP)+taHS));
    double tFC    = tdHS*(ttau0+tdHS*sin(tZeta))*(
		      tdHS  *cos(tZeta)*( tMt*sin(tZeta)*cosh((*aRap)-(*aAlfaP))+(*aPt)*cos(tZeta)*cos((*aPhiS)-(*aPhiP)))+
		      tDzdHS*cos(tZeta)*(-tMt*cos(tZeta)*cosh((*aRap)-(*aAlfaP))+(*aPt)*sin(tZeta)*cos((*aPhiS)-(*aPhiP)))+
		      tDpdHS*(*aPt)*sin((*aPhiS)-(*aPhiP))
		    );
   if(tFC < 0.0) tFC = 0.0;
    if (fabs(floor(aSpin)-aSpin) < 0.01)
      tFpod =1.0/kTwoPi3 * tdPt*(*aPt) * tFC * CalcBE((tPU-aMiu)/tTemp);
    else
      tFpod =1.0/kTwoPi3 * tdPt*(*aPt) * tFC * CalcFD((tPU-aMiu)/tTemp);

  }
  else if (sModel == 11) { // Lhyquid2D
    (*aAlfaP)     = (*aRap);					// dirac delta (Y-ap)
    double tZeta  = mRandom->Rndm()*0.5*TMath::Pi();		// angle in rho-t plane
    double taHS   = mFOHS->fahs  ((*aPhiS),tZeta);		// velocity angle; 0.0=radial flow
    double tvHS   = mFOHS->fvhs  ((*aPhiS),tZeta);		// velocity
    double tdHS   = mFOHS->fdhs  ((*aPhiS),tZeta)/kFmToGev;	// distance in rho-t plane           [GeV^-1]
    double tDpdHS = mFOHS->fDpdhs((*aPhiS),tZeta)/kFmToGev;	// distance derivative over (*aPhiS) [GeV^-1]
    double tDzdHS = mFOHS->fDzdhs((*aPhiS),tZeta)/kFmToGev;	// distance derivative over tZeta    [GeV^-1]
    double tTemp  = mFOHS->TFO;					// freeze-out temparature
    double ttau0  = mFOHS->tau0/kFmToGev;			// tau 0
    double tMt    = TMath::Hypot(aMass, (*aPt));		// transverse mass
    (*aRho)       = tdHS*cos(tZeta);				// rho
    double tdPt   = 1.0/((1-tZet)*(1-tZet));			// dPt
    (*aTime)      = (ttau0+tdHS*sin(tZeta))*cosh(*aAlfaP);	// t

    double tPU    = 1.0/sqrt(1.0-tvHS*tvHS)*(tMt-(*aPt)*tvHS*cos((*aPhiS)-(*aPhiP)+taHS));
    double tFC    = tdHS*(
		      tdHS  *cos(tZeta)*( (*aPt)/tMt*cos(tZeta)*cos((*aPhiS)-(*aPhiP))+sin(tZeta) )+
		      tDzdHS*cos(tZeta)*( (*aPt)/tMt*sin(tZeta)*cos((*aPhiS)-(*aPhiP))-cos(tZeta) )+
		      tDpdHS*             (*aPt)/tMt*           sin((*aPhiS)-(*aPhiP))
		    );
    if(tFC < 0.0) tFC = 0.0;
    if (fabs(floor(aSpin)-aSpin) < 0.01)
      tFpod =1.0/kTwoPi2 * tdPt*(*aPt) * tFC * CalcBE((tPU-aMiu)/tTemp);
    else
      tFpod =1.0/kTwoPi2 * tdPt*(*aPt) * tFC * CalcFD((tPU-aMiu)/tTemp);
  }
/*MCH end*/  
  else if (sModel == 2) { // Blast-Wave with Vt
    double tMt  = TMath::Hypot(aMass, (*aPt));
    double tPre = (mTau + mBWA*(*aRho))*(*aPt)*(*aRho)/kTwoPi3;
    double tCHay = cosh((*aAlfaP)-(*aRap));
    double tCSsp = cos((*aPhiS) - (*aPhiP));
    double tBra = tMt*tCHay - mBWA*(*aPt)*tCSsp;
    double tGamma = 1.0/sqrt(1-mBWVt*mBWVt);
    double tPU  = tMt*tGamma*tCHay-(*aPt)*mBWVt*tGamma*tCSsp - aMiu;
    if (fabs(floor(aSpin) - aSpin)< 0.01) 
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcBE(tPU/mTemp);
    else
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcFD(tPU/mTemp);
  }
  else if (sModel == 6) { // Blast-Wave with Vt
    double tMt  = TMath::Hypot(aMass, (*aPt));
    double tPre = (mTau + mBWA*(*aRho))*(*aPt)*(*aRho)/kTwoPi3;
    double tCHay = cosh((*aAlfaP)-(*aRap));
    double tCSsp = cos((*aPhiS) - (*aPhiP));
    double tBra = tMt*tCHay - mBWA*(*aPt)*tCSsp;
    double tGamma = 1.0/sqrt(1-mBWVt*mBWVt);
    double tPU  = tMt*tGamma*tCHay-(*aPt)*mBWVt*tGamma*tCSsp - aMiu;
    if (fabs(floor(aSpin) - aSpin)< 0.01) 
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcBE(tPU/mTemp);
    else
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcFD(tPU/mTemp);
  }
  else if (sModel == 4) { // Blast-Wave with linear Vt profile
    double tMt  = TMath::Hypot(aMass, (*aPt));
    double tPre = (mTau + mBWA*(*aRho))*(*aPt)*(*aRho)/kTwoPi3;
    double tCHay = cosh((*aAlfaP)-(*aRap));
    double tCSsp = cos((*aPhiS) - (*aPhiP));
    double tBra = tMt*tCHay - mBWA*(*aPt)*tCSsp;
    double tVt = ((*aRho)/mRhoMax)/(mBWVt+((*aRho)/mRhoMax));
    double tGamma = 1.0/sqrt(1-tVt*tVt);
    double tPU  = tMt*tGamma*tCHay-(*aPt)*tVt*tGamma*tCSsp - aMiu;
    if (fabs(floor(aSpin) - aSpin)< 0.01) 
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcBE(tPU/mTemp);
    else
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcFD(tPU/mTemp);
  }
  else if (sModel == 7) { 
    // Blast-Wave with linear Vt profile
    // and delay in particle emission point - formation time
    double tMt  = TMath::Hypot(aMass, (*aPt));
    double tPre = (mTau + mBWA*(*aRho))*(*aPt)*(*aRho)/kTwoPi3;
    double tCHay = cosh((*aAlfaP)-(*aRap));
    double tCSsp = cos((*aPhiS) - (*aPhiP));
    double tBra = tMt*tCHay - mBWA*(*aPt)*tCSsp;
    double tVt = ((*aRho)/mRhoMax)/(mBWVt+((*aRho)/mRhoMax));
    double tGamma = 1.0/sqrt(1-tVt*tVt);
    double tPU  = tMt*tGamma*tCHay-(*aPt)*tVt*tGamma*tCSsp - aMiu;
    if (fabs(floor(aSpin) - aSpin)< 0.01) 
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcBE(tPU/mTemp);
    else
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcFD(tPU/mTemp);
  }
  else if (sModel == 8) { 
    // Blast-Wave with linear Vt profile
    // and delay in emission time - exponential decay
    double tTau   = (*aTime);
    double tMt    = TMath::Hypot(aMass, (*aPt));
    double tPre   = (tTau + mBWA*(*aRho))*(*aPt)*(*aRho)/kTwoPi3;
    double tCHay  = cosh((*aAlfaP)-(*aRap));
    double tCSsp  = cos((*aPhiS) - (*aPhiP));
    double tBra   = tMt*tCHay - mBWA*(*aPt)*tCSsp;
    double tVt    = ((*aRho)/mRhoMax)/(mBWVt+((*aRho)/mRhoMax));
    double tGamma = 1.0/sqrt(1-tVt*tVt);
    double tPU    = tMt*tGamma*tCHay-(*aPt)*tVt*tGamma*tCSsp - aMiu;
    if (fabs(floor(aSpin) - aSpin)< 0.01) 
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcBE(tPU/mTemp);
    else
      tFpod = 1/((1-tZet)*(1-tZet))*tPre*tBra*CalcFD(tPU/mTemp);
  }
/* -=[ Testing with Mathematica ]=-*/
/*if ((sModel == 10)||(sModel == 11)) {
  cout.precision(20);
  cout << endl;
  cout << "{phi  = " << (*aPhiS)  << ","  << endl;
  cout << " zeta = " <<   tZeta   << ","  << endl;
  cout << " z    = " <<   tZet    << ","  << endl;
  cout << " phip = " << (*aPhiP)  << ","  << endl;
  cout << " ap   = " << (*aAlfaP) << ","  << endl;
  cout << " Y    = " << (*aRap)   << ","  << endl;
  cout << " m    = " <<   aMass   << ","  << endl;
  cout << " s    = " <<   aSpin   << ","  << endl;
  cout << " mu   = " <<   aMiu    << "};" << endl;
  cout << endl;
  cout.precision(6);
  cout << "TF0  = " <<   tTemp  << endl;
  cout << "t    = " << (*aTime) << endl;
  cout << "rho  = " << (*aRho)  << endl;
  cout << "pT   = " << (*aPt)   << endl;
  cout << "dpT  = " <<   tdPt   << endl;
  cout << "mT   = " <<   tMt    << endl;
  cout << endl;
  cout << "aHS     = " << taHS     << endl;
  cout << "vHS     = " << tvHS     << endl;
  cout << "dHS     = " << tdHS     << endl;
  cout << "DpdHS   = " << tDpdHS   << endl;
  cout << "DzdHS   = " << tDzdHS   << endl;
  cout << "PU      = " << tPU      << endl;
  cout << "FC      = " << tFC      << endl;
  cout << "part1   = " << (tPU - aMiu)/tTemp << endl;
  cout << "part2BE = " << CalcBE( (tPU - aMiu)/tTemp ) << endl;
  cout << "part2FD = " << CalcFD( (tPU - aMiu)/tTemp ) << endl;
  cout << "Fpod    = " << tFpod    << endl;
  getchar();
}*/

  return tFpod;
}

double Integrator::GetMiu(double aIzo, double aBar, double aStr)
{
  return (mMiu_i*aIzo+mMiu_s*aStr+mMiu_b*aBar);
}

double Integrator::CalcFunPodCalk(double aMass, double aIzo, double aBar, double aStr, double aSpin)
{
  //  int proc = mNPart/100;
  double tFpod;
  double tMax = 0.0;
  double tMiu = GetMiu(aIzo, aBar, aStr);
  double tRap, tPhiP, tAlfaP, tRho, tPhiS, tTime; 
  double tPt;

  for (int tIpart=0; tIpart<mNPart; tIpart++)
    {
      tFpod = Calka(aMass,tMiu,&tRap,&tPt,&tPhiP,&tAlfaP,&tRho,&tPhiS,&tTime,aSpin);
      if (tFpod>tMax) tMax = tFpod;
    }
  return tMax;
}

double Integrator::Integrate(double aMass, double aIzo, double aBar, double aStr, double aSpin)
{
  double tFpod=0.0;
  double tFtest=0.0;
  double tMiu = GetMiu(aIzo, aBar, aStr);
  double tRap, tPhiP, tAlfaP, tRho, tPhiS, tTime; 
  double tPt;

  for (int tIpart=0; tIpart<mNPart; tIpart++)
    {
      tFpod = Calka(aMass,tMiu,&tRap,&tPt,&tPhiP,&tAlfaP,&tRho,&tPhiS,&tTime,aSpin);

      tFtest += tFpod;
    }

  double tCalka;
  if (sModel == 5)
    tCalka = mRapRange*(1.0)*mAlfaRange*4*TMath::Pi()*TMath::Pi()*(mTauf-mTau0)* tFtest / (mNPart*1.0);  
/*MCH begin*/
  else if (sModel == 10)
    tCalka = mRapRange*mAlfaRange*(1.0)*2.0*TMath::Pi()*2.0*TMath::Pi()*0.5*TMath::Pi()*tFtest / (mNPart*1.0);
  else if (sModel == 11)
    tCalka = mRapRange*(1.0)*2.0*TMath::Pi()*2.0*TMath::Pi()*0.5*TMath::Pi()*tFtest / (mNPart*1.0);
/*MCH end*/    
  else
    tCalka = mRapRange*(1.0)*mAlfaRange*4*TMath::Pi()*TMath::Pi()*mRhoMax* tFtest / (mNPart*1.0);  
  return tCalka;
}

void 
Integrator::Generate(ParticleType *aPartType, int aPartCount, Particle ***oParticles)
{
  int tIter = 0;
  double tFpod;
  double tFtest;
  double tMiu = GetMiu(1.0*aPartType->GetI3(), 1.0*aPartType->GetBarionN(), 1.0*aPartType->GetStrangeness());
  double tFMax;
  double tSpin = aPartType->GetSpin();
  double tRap, tPhiP, tAlfaP, tRho, tPhiS, tTime; 
  double tPt;

  PRINT_DEBUG_3("Gen for " << (aPartType->GetName()) << " i3 b s " << 1.0*aPartType->GetI3() << " " << 1.0*aPartType->GetBarionN() << " " << 1.0*aPartType->GetStrangeness() << " " << tMiu);
  
  tFMax = aPartType->GetFMax();

  (*oParticles) = (Particle **) malloc(sizeof(Particle *) * aPartCount);
  Particle *tBuf=0;
  
  while (tIter<aPartCount)
    {
      tFpod = Calka(aPartType->GetMass(),tMiu,&tRap,&tPt,&tPhiP,&tAlfaP,&tRho,&tPhiS,&tTime,tSpin);
      tFtest = mRandom->Rndm()*tFMax;
      if (tFtest<tFpod)
	{
	  if ((sModel == 0) || (sModel == 3)) // Single freeze-out
	    tBuf = new Particle(tRap, tPt, tPhiP, tAlfaP, tRho, tPhiS, TMath::Hypot(mTau, tRho), aPartType);
/*MCH begin*/
	  else if ((sModel == 10) || (sModel == 11)) {
	    tBuf = new Particle(tRap, tPt, tPhiP, tAlfaP, tRho, tPhiS, tTime/cosh(tAlfaP), aPartType);
	  }
/*MCH end*/
	  else if ((sModel == 1) || (sModel == 2) || (sModel == 4)) { // Blast-wave
	    double tTau = mTau + mBWA * tRho;
	    tBuf = new Particle(tRap, tPt, tPhiP, tAlfaP, tRho, tPhiS, tTau, aPartType);
	  }
	  else if (sModel == 7) { // Blast-wave
	    double tTau = mTau + mBWA * tRho;
	    double px = tPt*TMath::Cos(tPhiP);
	    double py = tPt*TMath::Sin(tPhiP);
	    double tMt = TMath::Hypot(aPartType->GetMass(),tPt);
	    double pz = tMt*TMath::SinH(tRap);
  
	    double rx = tRho*TMath::Cos(tPhiS);
	    double ry = tRho*TMath::Sin(tPhiS);

	    double rz = tTau*TMath::SinH(tAlfaP);
	    double rt = tTau*TMath::CosH(tAlfaP);
	    double dt = -mBWDelay * TMath::Log(mRandom->Rndm());
	    rt += dt;
	    double en = sqrt(tMt*tMt + pz*pz);
	    rx += dt * px/en;
	    ry += dt * py/en;
	    rz += dt * pz/en;

	    tBuf = new Particle(aPartType, px, py, pz, rx, ry, rz, rt);
	  }
	  else if (sModel == 8) { // Blast-wave
	    double tTau = tTime + mBWA * tRho;
	    tBuf = new Particle(tRap, tPt, tPhiP, tAlfaP, tRho, tPhiS, tTau, aPartType);
	  }
	  else if (sModel == 6) { // Blast-wave
	    double tTau = mTau + mBWA * tRho;
	    double px = tPt*TMath::Cos(tPhiP);
	    double py = tPt*TMath::Sin(tPhiP);
	    double tMt = TMath::Hypot(aPartType->GetMass(),tPt);
	    double pz = tMt*TMath::SinH(tRap);
  
	    double rx = tRho*TMath::Cos(tPhiS);
	    double ry = tRho*TMath::Sin(tPhiS);

	    double rz = tTau*TMath::SinH(tAlfaP);
	    double rt = tTau*TMath::CosH(tAlfaP);
	    rt += -mBWDelay * TMath::Log(mRandom->Rndm());

	    tBuf = new Particle(aPartType, px, py, pz, rx, ry, rz, rt);
	  }
	  else if (sModel == 5) {
	    double tTau = sqrt(tTime*tTime - tAlfaP*tAlfaP);
	    double tAlfa = 0.5*log((tTime+tAlfaP)/(tTime-tAlfaP));
	    tBuf = new Particle(tRap, tPt, tPhiP, tAlfa, tRho, tPhiS, tTau, aPartType);
	  }
	  (*oParticles)[tIter] = tBuf;
	  tIter++;
	}
    }
}

void  
Integrator::Randomize()
{
  TDatime tDat;

  //  mRandom->SetSeed2(tDat.Get(), (tDat.Get() % 11) * 7 + (tDat.Get() / 7));
  mRandom->SetSeed(tDat.Get());
}

void  
Integrator::ReadParameters()
{
  STR tModel;
  STR tTable;

  // First read the model parameters
  try {
    
    mRhoMax = atof(sRPInstance->getPar("RhoMax").Data()) / kFmToGev;
    mTau    = atof(sRPInstance->getPar("Tau").Data()) / kFmToGev;
    mTemp   = atof(sRPInstance->getPar("Temperature").Data());
    mMiu_i  = atof(sRPInstance->getPar("MiuI").Data());
    mMiu_s  = atof(sRPInstance->getPar("MiuS").Data());
    mMiu_b  = atof(sRPInstance->getPar("MiuB").Data());
    mFOHSlocation = sRPInstance->getPar("FOHSLocation");
  }
  catch (STR tError) {
    PRINT_DEBUG_1("Integrator::ReadParameters - Caught exception " << tError);
    PRINT_MESSAGE("Did not find one of the neccessary model parameters.");
    exit(0);
  }

  // Read additional parameters for BW
  if ((sModel == 1) || (sModel == 2) || (sModel == 4) || (sModel == 7) || (sModel == 8)) {
    try {
      mBWA  = atof(sRPInstance->getPar("BWA").Data());
    }
    catch (STR tError) {
      // Using default value of 0
      mBWA = 0.0;
    }
  }
  
  // Read additional parameters for MBWVt
  if ((sModel == 2) || (sModel == 4) || (sModel == 7) || (sModel == 8)) {
    try {
      mBWVt  = atof(sRPInstance->getPar("BWVt").Data());
    }
    catch (STR tError) {
      PRINT_DEBUG_1("Integrator::ReadParameters (BW Vt part) - Caught exception " << tError);
      PRINT_MESSAGE("Did not find one of the neccessary model parameters.");
      exit(0);
    }
  }

  // Read additional parameters for MBWVt
  if ((sModel ==6) || (sModel == 7) || (sModel == 8))
    try {
      mBWDelay  = atof(sRPInstance->getPar("BWDelay").Data()) / kFmToGev;
    }
    catch (STR tError) {
      PRINT_DEBUG_1("Integrator::ReadParameters (BW Vt Delay part) - Caught exception " << tError);
      PRINT_MESSAGE("Did not find one of the neccessary model parameters.");
      exit(0);
    }
  
  
  // Then the integration range
  try {
    mAlfaRange = atof(sRPInstance->getPar("AlphaRange").Data());
    mRapRange  = atof(sRPInstance->getPar("RapidityRange").Data());
  }
  catch (STR tError) {
    PRINT_DEBUG_1("Integrator::ReadParameters - Caught exception " << tError);
    PRINT_MESSAGE("Did not find one of the neccessary integration ranges.");
    exit(0);
  }

  // Read hydro-specific parameters
  if (sModel == 5) {
    try {
      mTauf   = atof(sRPInstance->getPar("TauF").Data());
      mTau0   = atof(sRPInstance->getPar("Tau0").Data());
      mLambda = atof(sRPInstance->getPar("Lambda").Data());
      mBN     = atof(sRPInstance->getPar("BN").Data());
      mAlfa   = atof(sRPInstance->getPar("Alfa").Data());
    }
    catch (STR tError) {
      PRINT_DEBUG_1("Integrator::ReadParameters - Caught exception " << tError);
      PRINT_MESSAGE("Did not find one of the neccessary hydro parameters.");
      exit(0);
    }
  }
}

char *  
Integrator::ParameterHash()
{
  char *tBuf;
  
  tBuf = (char *) malloc(sizeof(char *) * (15 * 3 + 3));
  
  if (sModel == 0) {
    sprintf(tBuf, "%1i%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f", sTables, mRhoMax*100, mTau*100, mTemp*1000, mMiu_i*10000, mMiu_s*100000, mMiu_b*10000, mAlfaRange*100, mRapRange*100);
  }
/*MCH begin*/
  else if ((sModel == 10) || (sModel == 11)) {
    sprintf(tBuf, "%1i%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f", sTables, mRhoMax*100, mTau*100, mTemp*1000, mMiu_i*10000, mMiu_s*100000, mMiu_b*10000, mAlfaRange*100, mRapRange*100);
  }
/*MCH end*/
  else if ((sModel == 2) || (sModel == 4))  {
    sprintf(tBuf, "%1i%1i%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f", sModel, sTables, mRhoMax*100, mTau*100, mTemp*1000, mMiu_i*10000, mMiu_s*100000, mMiu_b*10000, mBWA*100, mBWVt*1000, mAlfaRange*100, mRapRange*100);
  }
  else if ((sModel == 6) || (sModel == 7) || (sModel == 8)) {
    sprintf(tBuf, "%1i%1i%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f", sModel, sTables, mRhoMax*100, mTau*100, mTemp*1000, mMiu_i*10000, mMiu_s*100000, mMiu_b*10000, mBWA*100, mBWVt*1000, mAlfaRange*100, mRapRange*100, mBWDelay);
  }
  else if (sModel == 5)  {
    sprintf(tBuf, "%1i%1i%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f%03.0f", sModel, sTables, mRhoMax*100, mTemp*1000, mMiu_i*10000, mMiu_s*100000, mMiu_b*10000, mTau0*10, mTauf*10, mAlfa*1000, mBN*100, mLambda*100, mAlfaRange*100, mRapRange*100);
  }
  
  return tBuf;
}

void   
Integrator::ReadMultInteg(ParticleDB *aDB)
{
  // Make or read table with propabilities
  int tK;
  char *tHash;
  char  tIntName[1000];
  char  tMultName[1000];
  ifstream *tFileIn = NULL;
  ofstream *tFileOut = NULL;
  
  tHash = ParameterHash();
  
  if (mFOHSlocation != "") {
    strcpy(tIntName, mFOHSlocation.Data());
    strcat(tIntName, "/");
    strcat(tIntName, "fintegrandmax_");
    strcat(tIntName, tHash);
    strcat(tIntName, ".txt");
    tFileIn = new ifstream(tIntName);
  }
  else if (!((tFileIn) && (tFileIn->is_open()))) {
    strcpy(tIntName, "fintegrandmax_");
    strcat(tIntName, tHash);
    strcat(tIntName, ".txt");
    tFileIn = new ifstream(tIntName);
  }
  
  if ((tFileIn) && (tFileIn->is_open())) {
    PRINT_MESSAGE("Reading Max Integrand values from " << tIntName);

    char tPart[255];
    float tRFunMax;
    std::string tStrPart;
    
    while (!tFileIn->eof()) {
      (*tFileIn) >> tPart >> tRFunMax;
      
      tStrPart = tPart;
      aDB->GetParticleType(tStrPart)->SetFMax(tRFunMax);
    }
  }
  else {
    float tRFunMax;

    PRINT_MESSAGE("Max Integrand file " << tIntName << " not found. Generating...");
    
    tFileOut = new ofstream(tIntName);

    for(tK=0; tK<aDB->GetParticleTypeCount(); tK++)
      {
	tRFunMax = CalcFunPodCalk(aDB->GetParticleType(tK)->GetMass(),
			       aDB->GetParticleType(tK)->GetI3(),
			       aDB->GetParticleType(tK)->GetBarionN()*1.0,
			       aDB->GetParticleType(tK)->GetStrangeness()*1.0,
			       aDB->GetParticleType(tK)->GetSpin());
	
	(*tFileOut) << aDB->GetParticleType(tK)->GetName() << "	"
		<< tRFunMax <<endl;
	aDB->GetParticleType(tK)->SetFMax(tRFunMax);
	PRINT_DEBUG_1("particle "<<aDB->GetParticleType(tK)->GetNumber()<<":"
		      <<aDB->GetParticleType(tK)->GetName()
		      <<" done");
      }
    tFileOut->close();
  }
  
  // Calculate or read multiplicities
  if (mFOHSlocation != "") {
    strcpy(tMultName, mFOHSlocation.Data());
    strcat(tMultName, "/");
    strcat(tMultName, "fmultiplicity_");
    strcat(tMultName, tHash);
    strcat(tMultName, ".txt");
    tFileIn = new ifstream(tMultName);
  }
  else if (!((tFileIn) && (tFileIn->is_open()))) {
    strcpy(tMultName, "fmultiplicity_");
    strcat(tMultName, tHash);
    strcat(tMultName, ".txt");
    tFileIn = new ifstream(tMultName);
  }

  tFileIn = new ifstream(tMultName);
  if ((tFileIn) && (tFileIn->is_open())) {
    PRINT_MESSAGE("Reading Multiplicities from " << tMultName);
  }
  else {
    PRINT_MESSAGE("Multiplicities file " << tMultName << " not found. Generating...");
    
    tFileOut = new ofstream(tMultName);
      
    for (int tPart=0; tPart<aDB->GetParticleTypeCount(); tPart++) {
      double tMult = Integrate(aDB->GetParticleType(tPart)->GetMass(),
			       aDB->GetParticleType(tPart)->GetI3(),
			       aDB->GetParticleType(tPart)->GetBarionN()*1.0,
			       aDB->GetParticleType(tPart)->GetStrangeness()*1.0,
			       aDB->GetParticleType(tPart)->GetSpin());
      (*tFileOut) << (aDB->GetParticleType(tPart)->GetName()) << " " << tMult*TMath::Abs(aDB->GetParticleType(tPart)->GetSpin()*2+1) << endl;

      PRINT_DEBUG_1("particle "<<aDB->GetParticleType(tPart)->GetNumber()<<":"
		    <<aDB->GetParticleType(tPart)->GetName()
		    <<" done");
    }
    tFileOut->close();
  }
      
  free (tHash);
}


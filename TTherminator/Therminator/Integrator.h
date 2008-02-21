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
#ifndef _BFPW_INTEGRATOR_
#define _BFPW_INTEGRATOR_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <TRandom2.h>
#include "THGlobal.h"
#include "ReadPar.h"
#include "ParticleDB.h"
#include "Hypersurface.h"	/*MCH*/

class ParticleType;
class Particle;

class Integrator
{
 public:
  Integrator(int);
  ~Integrator(){}
	  
  double CalcBE(double);
  double CalcFD(double);
  double GetMiu(double,double,double);

  double Calka(double,double,double*,double*,double*,double*,double*,double*,double*,double,int keeppos=0);
  double CalcFunPodCalk(double,double,double,double,double);
  double Integrate(double,double,double,double,double);
  void   Generate(ParticleType *aPartType, int aPartCount, Particle*** oParticles);
  
  void   Randomize();
  void   ReadMultInteg(ParticleDB *aDB);
  char  *ParameterHash();

 private:
  void   ReadParameters();
  
  int    mNPart;
  double kFmToGev;
  double mRhoMax;
  double mTau;
  double mTemp;
  double mMiu_i;
  double mMiu_s;
  double mMiu_b;

  double mBWA;
  double mBWP;

  double mBWVt;
  
  double mBWDelay;

  double mAlfaRange;
  double mRapRange;
  
  double kTwoPi2;		/*MCH*/
  double kTwoPi3;

  TRandom2 *mRandom;

  // Hydro parameters
  double mTauf;
  double mTau0;
  double mLambda;
  double mBN;
  double mAlfa;

  Hypersurface *mFOHS;		/*MCH*/
  TString mFOHSlocation;
};

#endif


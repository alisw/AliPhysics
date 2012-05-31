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
#ifndef _BFPW_PARSER_
#define _BFPW_PARSER_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <curses.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "THGlobal.h"
#include "ReadPar.h"
#include "ParticleType.h"
#include "Integrator.h"
#include "ParticleDB.h"

class Parser
{
 private:
  char        mNameBuffer[369][20];	        //local table with particle names
  double      mMassBuffer[369];		//local table with particle masses
  double      mGammaBuffer[369];		//local table with particle gammas	
  double      mSpinBuffer[369];		//local table with particle spins

  int         mBarionBuffer[369];		//local table with particle barion number
  int         mStrangeBuffer[369];		//local table with partcle strangeless number
  double      mI3Buffer[369];       	//local table with particle izospin3
  int         mDecayBuffer[369][2];

  int         mParticleCount;		//number of particles
  int         mTypeCountBuffer[369];	//number of each particle
  
  ParticleDB *mDB;

  TString     mInputDir;
  
  double DeltaJ(double aJot1, double aJot2, double aJot);
  double ClebschGordan(double aJot, double aEm, double aJot1, double aEm1, double aJot2, double aEm2);
  void   ReadParameters();

 public:
  Parser(ParticleDB *aDB);
  ~Parser();
  void    ReadInput();	//read input1 "tables.m" and "i200STAR.m"
  void    ReadShare();	//read input from SHaRe format: particles.data & decays.data

  int     GetParticleCount(); 
	
  int     check(char *,char *); //check particules by name
	
  char   *GetParticleName(int);	
  double  GetParticleMass(int);
  double  GetParticleGamma(int);
  double  GetParticleSpin(int);	
  int     GetParticleBarionN(int);
  int     GetParticleStrangeN(int);
  double  GetParticleI3(int);
  int     GetParticleDecayChannelCount(int,int);
  int     GetParticleNumber(int);
};

#endif

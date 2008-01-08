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
#ifndef _BFPW_PARTICLE_TYPE_
#define _BFPW_PARTICLE_TYPE_

#include <iostream>
#include "THGlobal.h"

class DecayChannel;
class DecayTable;

class ParticleType
{
 private:

  int         mNumber;	        //particle number
  float       mMass;    	//mass
  int         mStrangeness;	//strangness
  int         mBarionN; 	//barion number
  int         mCharmN;  	//charm
  float       mSpin;    	//spin
  float       mI;               //izospin
  float       mI3;      	//izospin3
  float       mGamma;   	//szerokosc polowkowa
  string      mName;    	//particle name
  int         mPDGCode;
  int         mDecayChannelCount2;	//number of channels in this case
  int         mDecayChannelCount3;	//number of channels in this case
  DecayTable *mTable;
  float       mFMax;

 public:
  ParticleType();	//constructor
  ParticleType(const ParticleType& aParticleType);	//copying constructor
  ~ParticleType();	//destructor

  int    GetNumber() const;
  float  GetMass() const;
  int    GetStrangeness() const;
  int    GetBarionN() const;
  int    GetCharmN() const;
  float  GetSpin() const;
  float  GetI() const;
  float  GetI3() const;
  float  GetGamma() const;
  string GetName() const;
  int    GetDecayChannelCount2() const;
  int    GetDecayChannelCount3() const;
  int    GetCharge();
  int    GetPDGCode() const;
  float  GetFMax() const;

  void SetNumber(int aNumber);
  void SetMass(float aMass);
  void SetStrangeness(int aStrangeness);
  void SetBarionN(int aBarionN);
  void SetCharmN(int aCharmN);
  void SetSpin(float aSpin);
  void SetI(float aI);
  void SetI3(float aI3);
  void SetGamma(float aGamma);
  void SetName(char *aName);
  void SetDecayChannelCount2(int aDCCount2);
  void SetDecayChannelCount3(int aDCCount3);
  void SetPDGCode(int aCode);
  void SetFMax(float aFMax);

  DecayTable* GetTable() const;
  void        AddDecayChannel(DecayChannel aChannel);
  
  void WriteParticle(int);
  
  static int GetParticleTypeNumber(std::string aPartName);
};


#endif

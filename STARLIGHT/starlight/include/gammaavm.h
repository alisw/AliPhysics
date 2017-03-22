///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 276                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-09-13 20:54:42 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef GAMMAAVM_H
#define GAMMAAVM_H


#include <vector>

#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "eventchannel.h"
#include "upcevent.h"
#include "nBodyPhaseSpaceGen.h"


class Gammaavectormeson : public eventChannel
{
  
 public:
  Gammaavectormeson(const inputParameters& ipnut, beamBeamSystem& bbsystem);
  virtual ~Gammaavectormeson();
  starlightConstants::event produceEvent(int &ievent);
  
   upcEvent produceEvent();

  void pickwy(double &W, double &Y);
  void momenta(double W,double Y,double &E,double &px,double &py,double &pz,int &tcheck);
  double pTgamma(double E); 
  void vmpt(double W,double Y,double &E,double &px,double &py, double &pz,int &tcheck);
  void twoBodyDecay(starlightConstants::particleTypeEnum &ipid,double W,double px0,double py0,double pz0,double &px1,double &py1,double&pz1,double &px2,double &py2,double &pz2,int &iFbadevent);
  bool fourBodyDecay(starlightConstants::particleTypeEnum& ipid, const double E, const double W, const double* p, lorentzVector* decayMoms, int& iFbadevent);
  double getMass();
  double getWidth();
  virtual double getTheta(starlightConstants::particleTypeEnum ipid);
  double getSpin();
  double _VMbslope;
  virtual double getDaughterMass(starlightConstants::particleTypeEnum &ipid);                
  double pseudoRapidity(double px, double py, double pz);
  
 private:
  starlightConstants::particleTypeEnum _VMpidtest;
  int _VMnumw;
  int _VMnumy;
  int _VMinterferencemode;
  int _ProductionMode;
  int _TargetBeam; 
  int N0;
  int N1;
  int N2; 
  double _VMgamma_em;
  double _VMNPT;
  double _VMWmax;
  double _VMWmin;
  double _VMYmax;
  double _VMYmin;
  double _mass;
  double _width;
  double _VMptmax;
  double _VMdpt;
  int    _bslopeDef;
  double _bslopeVal;
  double _pEnergy;
  nBodyPhaseSpaceGen* _phaseSpaceGen;
  
};

class Gammaanarrowvm : public Gammaavectormeson
{
 public:
  Gammaanarrowvm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~Gammaanarrowvm();
};

class Gammaawidevm : public Gammaavectormeson
{  
 public:
  Gammaawidevm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~Gammaawidevm();
};

class Gammaaincoherentvm : public Gammaavectormeson
{  
 public:
  Gammaaincoherentvm(const inputParameters& input, beamBeamSystem& bbsystem);
  virtual ~Gammaaincoherentvm();
};

#endif  // GAMMAAVM_H

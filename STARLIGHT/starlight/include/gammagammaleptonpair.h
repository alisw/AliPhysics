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
// $Rev:: 293                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2017-11-11 15:46:05 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef GAMMAGAMMALEPTONPAIR_H
#define GAMMAGAMMALEPTONPAIR_H

#include <vector>

#include "starlightconstants.h"
#include "readinluminosity.h"
#include "starlightlimits.h"
#include "eventchannel.h"
#include "inputParameters.h"

class Gammagammaleptonpair : public eventChannel
{
 public:
  Gammagammaleptonpair(const inputParameters& input, randomGenerator* randy, beamBeamSystem& bbsystem);
  ~Gammagammaleptonpair();
  
  void twoLeptonCrossSection();
  void calculateTable();

  starlightConstants::event produceEvent(int &ievent);
  upcEvent produceEvent();

 private:
  double _sigmax[starlightLimits::MAXWBINS][starlightLimits::MAXYBINS];
  double _sigfint[starlightLimits::MAXWBINS];
  double _sigofw[starlightLimits::MAXWBINS];
  double _signormw;
  double _wdelta;  //Added 7/26/07 for passing sigmadelta to pickw
  double _remainwd;// "
  int _ivalwd;     // "
  double _dgammade[1000];
  double _tautolangle[101];
  
  double twoMuonCrossSection(double w);
  void pickw(double &w);
  void picky(double &y);
  
  void pairMomentum(double w,double y,double &E,double &px,double &py,double&pz);
  double pp_1(double E);
  double pp_2(double E);
  void twoBodyDecay(starlightConstants::particleTypeEnum &ipid,double W,double px0,double py0,double pz0,double &px1,double &py1,double&pz1,double &px2,double &py2,double &pz2,int &iFbadevent);
  double thetalep(double W,double theta);
  void tauDecay(double &px1,double &py1,double &pz1,double &E1,double &px2,double &py2,double &pz2,double &E2);
  
  double getMass();
  double getWidth();
  double getSpin();
  
  const starlightConstants::particleTypeEnum _GGlepInputpidtest;
  const int _GGlepInputnumw;
  const int _GGlepInputnumy;
  const double _GGlepInputGamma_em;
};


#endif  // GAMMAGAMMALEPTONPAIR_H

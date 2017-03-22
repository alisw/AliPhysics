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
// $Rev:: 255                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-04-06 15:01:46 +0200 #$: date of last commit
//
// Description:
//    Class needed for root output
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "eventchannel.h"


using namespace std;


//______________________________________________________________________________
eventChannel::eventChannel(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem)
	: readLuminosity(inputParametersInstance),
	  _bbs(bbsystem),
	  _nmbAttempts(0),
	  _nmbAccepted(0),
	  _totalChannelCrossSection(0)
{
  _ptCutEnabled  = inputParametersInstance.ptCutEnabled();
  _ptCutMin      = inputParametersInstance.ptCutMin();
  _ptCutMax      = inputParametersInstance.ptCutMax();
  _etaCutEnabled = inputParametersInstance.etaCutEnabled();
  _etaCutMin     = inputParametersInstance.etaCutMin();
  _etaCutMax     = inputParametersInstance.etaCutMax();
  _randy.SetSeed(inputParametersInstance.randomSeed());
}


//______________________________________________________________________________
eventChannel::~eventChannel()
{ }


//______________________________________________________________________________
void
eventChannel::transform(const double  betax,
                        const double  betay,
                        const double  betaz,
                        double&       E,
                        double&       px,
                        double&       py,
                        double&       pz,
                        int&          iFbadevent)
{
  // carries out a lorentz transform of the frame.  (Not a boost!)???
  const double E0  = E;
  const double px0 = px;
  const double py0 = py;
  const double pz0 = pz;

  const double beta = sqrt(betax * betax + betay * betay + betaz * betaz);
  if (beta >= 1)
	  iFbadevent = 1;
  const double gamma = 1. / sqrt(1. - beta * beta);
  const double gob   = (gamma - 1) / (beta * beta);

  E   = gamma * (E0 - betax * px0 - betay * py0 - betaz*  pz0);
  px  = -gamma * betax * E0 + (1. + gob * betax * betax) * px0
	  + gob * betax * betay * py0 + gob * betax * betaz * pz0;
  py  = -gamma * betay * E0 + gob * betay * betax * px0
	  + (1. + gob * betay * betay) * py0 + gob * betay * betaz  *pz0;
  pz  = -gamma * betaz * E0 + gob * betaz * betax * px0
	  + gob * betaz * betay * py0 + (1. + gob * betaz * betaz) * pz0;
}


//______________________________________________________________________________
double
eventChannel::pseudoRapidity(const double px,
                             const double py,
                             const double pz)
{
  const double pT= sqrt(px * px + py * py);
  const double p = sqrt(pz * pz + pT * pT);
  double eta = -99.9;  // instead of special value, std::numeric_limits<double>::quiet_NaN() should be used
  if ((p - pz) != 0)
	  eta = 0.5 * log((p + pz)/(p - pz));
  return eta;
}


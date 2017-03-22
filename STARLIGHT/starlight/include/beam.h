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
// $Rev:: 259                         $: revision of last commit
// $Author:: jseger                   $: author of last commit
// $Date:: 2016-04-19 02:58:25 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef BEAM_H
#define BEAM_H


//This calls inclues a single beam of nucleons
#include "nucleus.h"


class beam : public nucleus
{

public:

	beam(const int              Z,
	     const int              A,
	     const int		    productionMode,
	     const double	    beamLorentzGamma);
	
	~beam();

	double photonDensity(const double impactparameter,
	                     const double photonEnergy) const;  ///< calculates photon density (number of photons / (energy * area))

	double rapidity() const { return acosh(_beamLorentzGamma); }
	
	void setBeamLorentzGamma(double gamma) {_beamLorentzGamma = gamma;}
protected:

	double _beamLorentzGamma;  ///< Lorentz gamma factor of beams in collider frame (from inputParameters)

};


#endif  // BEAM_H

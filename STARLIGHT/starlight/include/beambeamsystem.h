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
// $Rev:: 213                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2015-08-15 23:08:02 +0200 #$: date of last commit
//
// Description:
//     this class covers a coliding beam system SK
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef BEAMBEAMSYSTEM_H
#define BEAMBEAMSYSTEM_H


#include "nucleus.h"
#include "beam.h"
#include <vector>
#include "inputParameters.h"

class beamBeamSystem
{

public:

	beamBeamSystem(const inputParameters& input,
		       const beam&            beam1,
	               const beam&            beam2);
	beamBeamSystem(const inputParameters& input);
	~beamBeamSystem();

	const beam& beam1() const { return _beam1; }  ///< returns beam particle 1
	const beam& beam2() const { return _beam2; }  ///< returns beam particle 2

	double probabilityOfBreakup(const double D) const;
	
	double cmsBoost() const { return _cmsBoost; }
	
	double beamLorentzGamma() const { return _beamLorentzGamma; }
	
	void init();

private:
	void generateBreakupProbabilities();
	double probabilityOfHadronBreakup(const double impactparameter);
	double probabilityOfPhotonBreakup(const double impactparameter, const int mode);

	double _pHadronBreakup;
	double _pPhotonBreakup;

	double _beamLorentzGamma;  ///< Lorentz gamma factor of beams in collider frame
        const double _beamLorentzGamma1;  ///< Lorentz gamma factor of beam1 in collider frame
        const double _beamLorentzGamma2;  ///< Lorentz gamma factor of beam2 in collider frame
	const int    _beamBreakupMode;   ///< \brief breakup mode for beam particles
	                           ///<
	                           ///< 1 = hard sphere nuclei (b > 2R),
	                           ///< 2 = both nuclei break up (XnXn),
	                           ///< 3 = a single neutron from each nucleus (1n1n),
	                           ///< 4 = neither nucleon breaks up (with b > 2R),
	                           ///< 5 = no hadronic break up (similar to option 1, but with the actual hadronic interaction)
	               
	beam   _beam1;             ///< beam particle 1
	beam   _beam2;             ///< beam particle 2

	double _cmsBoost;	   ///< Rapidity boost of the CMS wrt the lab system
	
	std::vector<double> _breakupProbabilities; ///< Vector containing breakup probabilities for impact parameters
	double _breakupImpactParameterStep; ///< Step size in the calculation of the breakup probs
	double _breakupCutOff;  ///< Cut off for minimum impact parameter probability
};


#endif  // BEAMBEAMSYSTEM_H

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
// $Rev:: 286                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2017-05-03 22:32:03 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef NARROWRESONANCECROSSSECTION_H
#define NARROWRESONANCECROSSSECTION_H


#include "photonNucleusCrossSection.h"


class narrowResonanceCrossSection : public photonNucleusCrossSection {

public:

	narrowResonanceCrossSection(const inputParameters& input, const beamBeamSystem&  bbsystem);
	~narrowResonanceCrossSection();

	void crossSectionCalculation(const double bwnormsave);

private:
	
	double _Ep;
	double _gamma1;
	double _gamma2;
	double _narrowYmax;
	double _narrowYmin;
	int    _narrowNumY;
	
};


#endif  // NARROWRESONANCECROSSSECTION_H

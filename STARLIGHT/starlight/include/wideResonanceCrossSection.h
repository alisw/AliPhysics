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
// $Rev:: 280                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-12-12 20:38:20 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef WIDERESONANCECROSSSECTION_H
#define WIDERESONANCECROSSSECTION_H


#include "photonNucleusCrossSection.h"
#include "inputParameters.h"

class wideResonanceCrossSection : public photonNucleusCrossSection {

public:

	wideResonanceCrossSection(const inputParameters& input, const beamBeamSystem& bbsystem);
	~wideResonanceCrossSection();

	void crossSectionCalculation(const double bwnormsave);

private:

	double _Ep;  // Proton Energy
	double _wideWmax;
	double _wideWmin;
	double _wideYmax;
	double _wideYmin;		
        int    _wideNumY;

};


#endif  // WIDERESONANCECROSSSECTION_H

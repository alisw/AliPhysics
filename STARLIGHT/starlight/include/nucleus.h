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
// $Rev:: 262                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-06-01 15:14:20 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef NUCLEUS_H
#define NUCLEUS_H


#include <cmath>


//This class holds the information for a target nucleus
class nucleus
{

public:
	nucleus();
	nucleus(const int    Z,
	        const int    A,
		const int     productionMode);
	~nucleus();
	
	void init();
 
	int    Z              () const { return _Z;                     }  ///< returns atomic number of nucleus
	int    A              () const { return _A;                     }  ///< returns nucleon number of nucleus
        int    productionMode () const { return _productionMode;        }

	double formFactor(const double t) const;
	// Calculates form factor for given squared 4-momentum transfer

	double dipoleFormFactor(const double t, const double t0) const;
	// Calculates dipole form factor with t0 as parameter 

	double thickness (const double b) const;
	// Calculates nuclear thickness function 

	double nuclearRadius() const { return _Radius; }
	double rho0() const { return _rho0; }
	
private:

	double woodSaxonSkinDepth() const { return 0.53; } // 0.53 fm skin depth
	double rws(const double r) const;

	int    _Z;                      ///< atomic number of nucleus
	int    _A;                      ///< nucleon number of nucleus
        int    _productionMode;

	double _r0;
	double _Radius;
	double _rho0;

};


#endif  // NUCLEUS_H

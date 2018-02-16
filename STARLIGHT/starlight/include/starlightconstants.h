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


#ifndef STARLIGHTCONSTANTS_H_INCLUDE
#define STARLIGHTCONSTANTS_H_INCLUDE


/*
 * Constants are set here
 */
namespace starlightConstants
{


	// constants
	static const double hbarc    = 0.1973269718;
	static const double hbarcmev = hbarc*1000.;
	static const double pi       = 3.141592654;
	static const double twoPi    = 2 * pi;
	static const double alpha    = 1/137.035999074;

	enum particleTypeEnum {
		UNKNOWN        = 0,
		ELECTRON       = 11,
		MUON           = 13,
		TAUON          = 15,
		TAUONDECAY     = 10015,
		PROTON         = 2212,
		PION           = 211,
		KAONCHARGE     = 321,
		KAONNEUTRAL    = 310,
		A2             = 115,
		ETA            = 221,
		F2             = 225,
		ETAPRIME       = 331,
		F2PRIME        = 335,
		ETAC           = 441,
		F0             = 9010221,
		ZOVERZ03       = 33,
		RHO            = 113,
		RHOZEUS        = 913,
		FOURPRONG      = 999,
		OMEGA          = 223,
		PHI            = 333,
		JPSI           = 443,
		JPSI_ee        = 443011,
		JPSI_mumu      = 443013,
		JPSI_ppbar     = 4432212,
		JPSI2S         = 444,
		JPSI2S_ee      = 444011,
		JPSI2S_mumu    = 444013,
		UPSILON        = 553,
		UPSILON_ee     = 553011,
		UPSILON_mumu   = 553013,
		UPSILON2S      = 554,
		UPSILON2S_ee   = 554011,
		UPSILON2S_mumu = 554013,
		UPSILON3S      = 555,
		UPSILON3S_ee   = 555011,
		UPSILON3S_mumu = 555013,
	        AXION          = 88,  //AXION HACK
        	PHOTON         = 22
	};

	enum decayTypeEnum {
		NOTKNOWN        = 0,
		NARROWVMDEFAULT = 1,
		WIDEVMDEFAULT   = 2,
		PSIFAMILY       = 3,
		LEPTONPAIR      = 4,
		SINGLEMESON     = 5
	};

	enum interactionTypeEnum {
		UNSPECIFIED         = 0,
		PHOTONPHOTON        = 1,
		PHOTONPOMERONNARROW = 2,
		PHOTONPOMERONWIDE   = 3,
                PHOTONPOMERONINCOHERENT = 4,
                PHOTONUCLEARSINGLE  = 5,
		PHOTONUCLEARDOUBLE  = 6,
		PHOTONUCLEARSINGLEPA = 7,
		PHOTONUCLEARSINGLEPAPY = 8
		
	};

        enum systemTypeEnum{
		NONSTANDARD         = 0,
		PP                  = 1,
		PA                  = 2,
		AA                  = 3
        };
		
        enum targetTypeEnum {
                NOTHADRON           = 0,
                NUCLEUS             = 1, //coherent gamma+A
                NUCLEON             = 2  //gamma+p or incoherent gamma+A
        };        
	
	//Structure for each event's set of tracks.
	struct event{
 
	public:

		int _numberOfTracks;
		double px[30],py[30],pz[30];
		int _fsParticle[30];
		int _charge[30];
		//To help track mothers and daughters produced through pythia.
		int _mother1[30];
		int _mother2[30];
		int _daughter1[30];
		int _daughter2[30];
		//Normally we just set vertices to 0
		//But for pythia, we decay additional states
		int _numberOfVertices;
		double _vertx[10],_verty[10],_vertz[10];	
	};


}  // starlightConstants


#endif  // STARLIGHTCONSTANTS_H_INCLUDE


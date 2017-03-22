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
// $Rev:: 283                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2017-03-07 18:17:50 +0100 #$: date of last commit
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

	// deuteron slope parameter
	static const double deuteronSlopePar = 9.5; // [(GeV/c)^{-2}]
	// particle masses
	static const double protonMass      = 0.938272046;   // [GeV/c^2]
	static const double pionChargedMass = 0.13957018;    // [GeV/c^2]
        static const double pionNeutralMass = 0.1349766;     // [GeV/c^2]
	static const double kaonChargedMass = 0.493677;      // [GeV/c^2]
	static const double mel      = 0.000510998928;       // [GeV/c^2]
	static const double muonMass = 0.1056583715;         // [GeV/c^2]
	static const double tauMass  = 1.77682;              // [GeV/c^2]
	
	static const double f0Mass = 0.990;                  // [GeV/c^2]
	static const double f0Width = 0.100;                 // [GeV/c^2]
	static const double f0BrPiPi = 1.0;                  // Branching ratio for pipi (set to 100%)
        static const double etaMass = 0.547862;              // [GeV/c^2]
        static const double etaWidth = 0.00000131;           // [GeV/c^2]
        static const double etaPrimeMass = 0.95778;          // [GeV/c^2]
        static const double etaPrimeWidth = 0.000198;        // [GeV/c^2]
        static const double etaCMass = 2.9836;               // [GeV/c^2]
        static const double etaCWidth = 0.0322;              // [GeV/c^2]
        static const double f2Mass = 1.2751;                 // [GeV/c^2]
        static const double f2Width = 0.1851;                // [GeV/c^2]
        static const double f2BrPiPi = 0.561;                // Branching ratio for pi+pi-
        static const double a2Mass = 1.3183;                 // [GeV/c^2]
        static const double a2Width = 0.105;                 // [GeV/c^2]
        static const double f2PrimeMass = 1.525;             // [GeV/c^2]
        static const double f2PrimeWidth = 0.073;            // [GeV/c^2]
        static const double f2PrimeBrKK = 0.887;             // Branching ratio for KKbar
        static const double zoverz03Mass = 1.540;            // [GeV/c^2]


        static const double f0PartialggWidth = 0.29E-6;      // [GeV/c^2]
        static const double etaPartialggWidth = 0.516E-6;    // [GeV/c^2]
        static const double etaPrimePartialggWidth = 4.35E-6;// [GeV/c^2]
        static const double etaCPartialggWidth = 5.0E-6;     // [GeV/c^2]
        static const double f2PartialggWidth = 3.03E-6;      // [GeV/c^2]
        static const double a2PartialggWidth = 1.0E-6;       // [GeV/c^2]
        static const double f2PrimePartialggWidth = 0.081E-6;// [GeV/c^2]
        static const double zoverz03PartialggWidth = 0.1E-6; // [GeV/c^2]

        static const double f0Spin = 0.0;
        static const double etaSpin = 0.0;
        static const double etaPrimeSpin = 0.0;
        static const double etaCSpin = 0.0;
        static const double f2Spin = 2.0;
        static const double a2Spin = 2.0;
        static const double f2PrimeSpin = 2.0;
        static const double zoverz03Spin = 2.0;
        static const double axionSpin = 0.0;  // AXION HACK

        static const double rho0Mass  = 0.769;               // [GeV/c^2]
	static const double rho0Width = 0.1517;              // [GeV/c^2]
	static const double rho0BrPiPi = 1.0;                // Branching ratio pi+pi-
        static const double rho0PrimeMass  = 1.540;          // [GeV/c^2]
	static const double rho0PrimeWidth = 0.570;          // [GeV/c^2]
	static const double rho0PrimeBrPiPi = 1.0;           // Branching ratio pi+pi- (set to 100%)
        static const double OmegaMass  = 0.78265;            // [GeV/c^2]
	static const double OmegaWidth = 0.00849;            // [GeV/c^2]
	static const double OmegaBrPiPi = 0.0153;            // Branching ratio pi+pi-
        static const double PhiMass  = 1.019461;             // [GeV/c^2]
	static const double PhiWidth = 0.004266;             // [GeV/c^2]
	static const double PhiBrKK = 0.489;                 // Branching ratio K+K-
        static const double JpsiMass = 3.096916;             // [GeV/c^2]
        static const double JpsiWidth = 0.0000929;           // [GeV/c^2]
	static const double JpsiBree = 0.05971;              // Branching ratio e+e-
	static const double JpsiBrmumu = 0.05961;            // Branching ratio mu+mu-
	static const double JpsiBrppbar = 0.002120;          // Branching ratio ppbar
        static const double Psi2SMass = 3.686109;            // [GeV/c^2]
        static const double Psi2SWidth = 0.000299;           // [GeV/c^2]
	static const double Psi2SBree = 0.00789;             // Branching ratio e+e-
	static const double Psi2SBrmumu = 0.0079;            // Branching ratio mu+mu-
        static const double Upsilon1SMass = 9.46030;         // [GeV/c^2]
        static const double Upsilon1SWidth = 0.00005402;     // [GeV/c^2]
	static const double Upsilon1SBree = 0.0238;          // Branching ratio e+e-
	static const double Upsilon1SBrmumu = 0.0248;        // Branching ratio mu+mu-
        static const double Upsilon2SMass = 10.02326;        // [GeV/c^2]
        static const double Upsilon2SWidth = 0.00003198;     // [GeV/c^2]
	static const double Upsilon2SBree = 0.0191;          // Branching ratio e+e-
	static const double Upsilon2SBrmumu = 0.0193;        // Branching ratio mu+mu-
        static const double Upsilon3SMass = 10.3552;         // [GeV/c^2]
        static const double Upsilon3SWidth = 0.00002032;     // [GeV/c^2]	
	static const double Upsilon3SBree = 0.0218;          // Branching ratio e+e- (set to same as mu+mu-)
	static const double Upsilon3SBrmumu = 0.0218;        // Branching ratio mu+mu-
	
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


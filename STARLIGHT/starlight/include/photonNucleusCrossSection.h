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
// $Rev:: 284                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2017-04-25 22:08:11 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef PHOTONNUCLEUSCROSSSECTION_H
#define PHOTONNUCLEUSCROSSSECTION_H


#include "starlightconstants.h"
#include "beambeamsystem.h"
#include "inputParameters.h"

class photonNucleusCrossSection {

public:

	photonNucleusCrossSection(const inputParameters& input, const beamBeamSystem&  bbsystem);
	~photonNucleusCrossSection();
  
	double         slopeParameter  () const { return _slopeParameter;   }  ///< returns slope of t-distribution [(GeV/c)^{-2}]
	double         getChannelMass  () const { return _channelMass;      }  ///< returns mass of the produced system [GeV/c^2]
	double         getBNORM        () const { return _BNORM;            }
	beamBeamSystem getbbs          () const { return _bbs;              }  ///< returns beamBeamSystem
	double         vmPhotonCoupling() const { return _vmPhotonCoupling; }  ///< vectormeson-photon coupling constant f_v / 4 pi (cf. Eq. 10 in KN PRC 60 (1999) 014903)
	double         getDefaultC     () const { return _defaultC;         }
	double         maxPhotonEnergy () const { return _maxPhotonEnergy;  }  ///< returns max photon energy in lab frame [GeV] (for vectormesons only)

	void crossSectionCalculation(const double bwnormsave);
	// Use the wide or narrow constructor to calculate sigma
	// wide/narrow will inherit this.
	double getcsgA(const double Egamma,
	               const double W,
                       const int beam);
	double photonFlux(const double Egamma,
                       const int beam);
	double sigmagp(const double Wgp);
	double sigma_A(const double sig_N, 
                       const int beam);
        double sigma_N(const double Wgp);
	double breitWigner(const double W,
	                   const double C);
	double nepoint(const double Egamma,
	               const double bmin);

	double getPhotonNucleusSigma () const {return _photonNucleusSigma;}
	void   setPhotonNucleusSigma (double sigma) {_photonNucleusSigma = sigma;}
	
protected:
	const unsigned int _nWbins;
	const unsigned int _nYbins;
	
	const double _wMin;
	const double _wMax;
	const double _yMax;
	
	const double _beamLorentzGamma;

	double _photonNucleusSigma; 

        int    _printDef; 
        int    _impulseSelected;
	
private:

	beamBeamSystem _bbs;
  
	// copied from inputParameters
	double                               _protonEnergy;
	starlightConstants::particleTypeEnum _particleType;
	int                                  _beamBreakupMode;     ///< breakup mode for beam particles
        int                                  _productionMode; 
	int                                  _sigmaNucleus; 

	// locally defined parameters
	double _slopeParameter;    ///< slope of t-distribution [(GeV/c)^{-2}]
	double _vmPhotonCoupling;  ///< vectormeson-photon coupling constant f_v / 4 pi (cf. Eq. 10 in KN PRC 60 (1999) 014903)
	double _ANORM;
	double _BNORM;
	double _defaultC;
	double _maxPhotonEnergy;  ///< max photon energy in lab frame [GeV] (for vectormesons only)
	double _width;            ///< width of the produced system  [GeV/c^2]
	double _channelMass;      ///< mass of the produced system  [GeV/c^2]
	
};


#endif  // PHOTONNUCLEUSCROSSSECTION_H

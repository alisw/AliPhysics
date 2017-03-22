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
 

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "starlightconfig.h"

#ifdef ENABLE_PYTHIA
#include "PythiaStarlight.h"
#endif

#ifdef ENABLE_DPMJET
#include "starlightdpmjet.h"
#endif

#ifdef ENABLE_PYTHIA6
#include "starlightpythia.h"
#endif

#include "reportingUtils.h"
#include "inputParameters.h"
#include "eventchannel.h"
#include "gammagammaleptonpair.h"
#include "gammagammasingle.h"
#include "gammaavm.h"
#include "twophotonluminosity.h"
#include "gammaaluminosity.h"
#include "incoherentPhotonNucleusLuminosity.h"
#include "upcevent.h"
#include "eventfilewriter.h"
#include "starlight.h"


using namespace std;
using namespace starlightConstants;


starlight::starlight() :
		_beamSystem            (0),
		_eventChannel          (0),
		_nmbEventsPerFile      (100),
		_isInitialised         (false),
		_inputParameters       (0)
{ }


starlight::~starlight()
{ }


bool
starlight::init()
{
	if(Starlight_VERSION_MAJOR == 9999)
	{
	  cout  << endl << "#########################################" << endl
	     	<< " Initialising Starlight version: trunk..." << endl
	     	<< "#########################################" << endl << endl;
	}
	else
	{
	  cout  << endl << "#########################################" << endl
	     	<< " Initialising Starlight version: " << Starlight_VERSION_MAJOR << "."
	     	<< Starlight_VERSION_MINOR << "." << Starlight_VERSION_MINOR_MINOR << "..." << endl
	        << "#########################################" << endl << endl;
	}

	_nmbEventsPerFile    = _inputParameters->nmbEvents();  // for now we write only one file...

	_beamSystem = new beamBeamSystem(*_inputParameters);
	
	// This sets the precision used by cout for printing
	cout.setf(ios_base::fixed,ios_base::floatfield);
	cout.precision(3);

        std::string _baseFileName;
        _baseFileName = _inputParameters->baseFileName();
       _lumLookUpTableFileName = _baseFileName + ".txt";

	const bool lumTableIsValid = luminosityTableIsValid();

	// Do some sanity checks of the input parameters here.
        if( _inputParameters->beam1Z() > _inputParameters->beam1A() ){
	  printErr << endl << "A must be >= Z; A beam1 = "<<_inputParameters->beam1A()<<", Z beam1 = "<<_inputParameters->beam1Z()<<". Terminating."<<endl ;
	  return false;
	}
        if( _inputParameters->beam2Z() > _inputParameters->beam2A() ){
	  printErr << endl << "A must be >= Z; A beam2 = "<<_inputParameters->beam2A()<<", Z beam2 = "<<_inputParameters->beam2Z()<<". Terminating."<<endl ;
	  return false;
	}
	if( _inputParameters->interactionType() == PHOTONPOMERONINCOHERENT && _inputParameters->beam1A() == 1 &&
	    _inputParameters->beam1Z() == 1 && _inputParameters->beam2A() == 1 && _inputParameters->beam2Z() ){
          printErr << endl << " Do not use PROD_MODE = 4 for pp collisions. Use PROD_MODE = 2 or 3 instead. Terminating."<<endl;
	  return false; 
	}

	bool createChannel = true;
	switch (_inputParameters->interactionType())	{
	case PHOTONPHOTON:
		if (!lumTableIsValid) {
			printInfo << "creating luminosity table for photon-photon channel" << endl;
			twoPhotonLuminosity(*_inputParameters, _beamSystem->beam1(), _beamSystem->beam2());
		}
		break;		
	case PHOTONPOMERONNARROW:  // narrow and wide resonances use
	case PHOTONPOMERONWIDE:    // the same luminosity function
		if (!lumTableIsValid) {
			printInfo << "creating luminosity table for coherent photon-Pomeron channel" <<endl;
			photonNucleusLuminosity lum(*_inputParameters, *_beamSystem);
		}
		break;
        case PHOTONPOMERONINCOHERENT:  // narrow and wide resonances use
                if (!lumTableIsValid) {
                        printInfo << "creating luminosity table for incoherent photon-Pomeron channel" << endl;
                        incoherentPhotonNucleusLuminosity lum(*_inputParameters, *_beamSystem);
                }
                break;
#ifdef ENABLE_DPMJET
	case PHOTONUCLEARSINGLE:
		createChannel = false;
		_eventChannel = new starlightDpmJet(*_inputParameters, *_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/DPMJET SINGLE" << std::endl;
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setSingleMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMinGammaEnergy(_inputParameters->minGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMaxGammaEnergy(_inputParameters->maxGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->init();
		break;
	case PHOTONUCLEARDOUBLE:
		createChannel = false;
		_eventChannel = new starlightDpmJet(*_inputParameters, *_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/DPMJET DOUBLE" << std::endl;
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setDoubleMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMinGammaEnergy(_inputParameters->minGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMaxGammaEnergy(_inputParameters->maxGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->init();
		break;
	case PHOTONUCLEARSINGLEPA:
		createChannel = false;
		_eventChannel = new starlightDpmJet(*_inputParameters, *_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/DPMJET SINGLE" << std::endl;
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setSingleMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setProtonMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMinGammaEnergy(_inputParameters->minGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMaxGammaEnergy(_inputParameters->maxGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->init();
		break;
#endif
#ifdef ENABLE_PYTHIA6
	case PHOTONUCLEARSINGLEPAPY:
		createChannel = false;
		_eventChannel = new starlightPythia(*_inputParameters, *_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/PYTHIA SINGLE" << std::endl;
		dynamic_cast<starlightPythia*>(_eventChannel)->setSingleMode();
		dynamic_cast<starlightPythia*>(_eventChannel)->setMinGammaEnergy(_inputParameters->minGammaEnergy());
		dynamic_cast<starlightPythia*>(_eventChannel)->setMaxGammaEnergy(_inputParameters->maxGammaEnergy());
		dynamic_cast<starlightPythia*>(_eventChannel)->init(_inputParameters->pythiaParams(), _inputParameters->pythiaFullEventRecord());
		break;
#endif
	default:
		{
			printWarn << "unknown interaction type '" << _inputParameters->interactionType() << "'."
			          << " cannot initialize starlight." << endl;
			return false;
		}
	}
	
	if(createChannel)
	{
	  if (!createEventChannel())
		  return false;
	}

	_isInitialised = true;
	return true;
}


upcEvent
starlight::produceEvent()
{
	if (!_isInitialised) {
		printErr << "trying to generate event but Starlight is not initialised. aborting." << endl;
		exit(-1);
	}
	return _eventChannel->produceEvent();
}


bool
starlight::luminosityTableIsValid() const
{
	printInfo << "using random seed = " << _inputParameters->randomSeed() << endl;

	ifstream lumLookUpTableFile(_lumLookUpTableFileName.c_str());
	lumLookUpTableFile.precision(15);
	if ((!lumLookUpTableFile) || (!lumLookUpTableFile.good())) {
		return false;
	}

	unsigned int beam1Z, beam1A, beam2Z, beam2A;
	double       beamLorentzGamma = 0;
	double       maxW = 0, minW = 0;
	unsigned int nmbWBins;
	double       maxRapidity = 0;
	unsigned int nmbRapidityBins;
	int          productionMode, beamBreakupMode;
	bool         interferenceEnabled = false;
	double       interferenceStrength = 0;
	bool         coherentProduction = false;
	double       incoherentFactor = 0, deuteronSlopePar = 0, maxPtInterference = 0;
	int          nmbPtBinsInterference;
	std::string  validationKey;
	if (!(lumLookUpTableFile
	      >> validationKey
	      >> beam1Z >> beam1A
	      >> beam2Z >> beam2A
	      >> beamLorentzGamma
	      >> maxW >> minW >> nmbWBins
	      >> maxRapidity >> nmbRapidityBins
	      >> productionMode
	      >> beamBreakupMode
	      >> interferenceEnabled >> interferenceStrength
	      >> maxPtInterference
	      >> nmbPtBinsInterference
	      >> coherentProduction >> incoherentFactor
	      >> deuteronSlopePar))
		// cannot read parameters from lookup table file
		return false;
        	
	std::string validationKeyEnd;
	while(!lumLookUpTableFile.eof())
	{
	  lumLookUpTableFile >> validationKeyEnd; 
	}
	lumLookUpTableFile.close();
	return (validationKey == _inputParameters->parameterValueKey() && validationKeyEnd == validationKey);
	return true;
}


bool
starlight::createEventChannel()
{
	switch (_inputParameters->prodParticleType()) {
	case ELECTRON:
	case MUON:
	case TAUON:
        case TAUONDECAY:
		{
			_eventChannel = new Gammagammaleptonpair(*_inputParameters, *_beamSystem);
			if (_eventChannel)
				return true;
			else {
				printWarn << "cannot construct Gammagammaleptonpair event channel." << endl;
				return false;
			}
		}
	case A2:        // jetset/pythia
	case ETA:       // jetset/pythia
	case ETAPRIME:  // jetset/pythia
	case ETAC:      // jetset/pythia
	case F0:        // jetset/pythia
		{
#ifdef ENABLE_PYTHIA
			// PythiaOutput = true;
 		        cout<<"Pythia is enabled!"<<endl;
// 			return true;
#else
			printWarn << "Starlight is not compiled against Pythia8; "
			          << "jetset event channel cannot be used." << endl;
 			return false;
#endif
		}
	case F2:
	case F2PRIME:
	case ZOVERZ03:
    case AXION: // AXION HACK
		{
		  //  #ifdef ENABLE_PYTHIA
	 	        cout<<" This is f2, f2prim, rho^0 rho^0, or axion "<<endl; 
			_eventChannel= new Gammagammasingle(*_inputParameters, *_beamSystem);
			if (_eventChannel)
				return true;
			else {
				printWarn << "cannot construct Gammagammasingle event channel." << endl;
				return false;
			}
		}
	case RHO:
	case RHOZEUS:
	case FOURPRONG:
	case OMEGA:  
	case PHI:
	case JPSI:
	case JPSI_ee:
	case JPSI_mumu:
	case JPSI_ppbar:
	case JPSI2S:
	case JPSI2S_ee:
	case JPSI2S_mumu:
	case UPSILON:
	case UPSILON_ee:
	case UPSILON_mumu:
	case UPSILON2S:
	case UPSILON2S_ee:
	case UPSILON2S_mumu:
	case UPSILON3S:
	case UPSILON3S_ee:
	case UPSILON3S_mumu:
		{
			if (_inputParameters->interactionType() == PHOTONPOMERONNARROW) {
				_eventChannel = new Gammaanarrowvm(*_inputParameters, *_beamSystem);
				if (_eventChannel)
					return true;
				else {
					printWarn << "cannot construct Gammaanarrowvm event channel." << endl;
					return false;
				}
			}

			if (_inputParameters->interactionType() == PHOTONPOMERONWIDE) {
				_eventChannel = new Gammaawidevm(*_inputParameters, *_beamSystem);
				if (_eventChannel)
					return true;
				else {
					printWarn << "cannot construct Gammaawidevm event channel." << endl;
					return false;
				}
			}

                        if (_inputParameters->interactionType() == PHOTONPOMERONINCOHERENT) {
                                _eventChannel = new Gammaaincoherentvm(*_inputParameters, *_beamSystem);
                                if (_eventChannel)
                                        return true;
                                else {
                                        printWarn << "cannot construct Gammaanarrowvm event channel." << endl;
                                        return false;
                                }
                        }

			printWarn << "interaction type '" << _inputParameters->interactionType() << "' "
			          << "cannot be used with particle type '" << _inputParameters->prodParticleType() << "'. "
			          << "cannot create event channel." << endl;
			return false;
		}
	default:
		{
			printWarn << "unknown event channel '" << _inputParameters->prodParticleType() << "'."
			          << " cannot create event channel." << endl;
			return false;
		}
	}
}

///////////////////////////////////////////////////////////////////////////
//
// Copyright 2010
//
// This file is part of starlight.
//
// starlight is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// starlight is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 293 $: revision of last commit // $Author:: butter  $: author of last commit
// $Date:: 2017-11-11 15:46:05 +0100 #$: date of last commit //
// Description:
//

#include <iostream>
#include <fstream>

#include "reportingUtils.h"
#include "starlightconstants.h"
#include "inputParameters.h"
#include "inputParser.h"
#include "starlightconfig.h"
#include <cmath>
#include <cstring>
#include "randomgenerator.h"

using namespace std;
using namespace starlightConstants;

#define REQUIRED true
#define NOT_REQUIRED false

//______________________________________________________________________________
inputParameters::inputParameters()
        : inputParametersBase(),
          _baseFileName          (this, "baseFileName","slight"),
 	  _beam1Z                (this, "BEAM_1_Z",0),
	  _beam1A                (this, "BEAM_1_A",0),
	  _beam2Z                (this, "BEAM_2_Z",0),
	  _beam2A                (this, "BEAM_2_A",0),
	  _beam1LorentzGamma     (this, "BEAM_1_GAMMA",0),
	  _beam2LorentzGamma     (this, "BEAM_2_GAMMA",0),
	  _maxW                  (this, "W_MAX",0),
	  _minW                  (this, "W_MIN",0),
	  _nmbWBins              (this, "W_N_BINS",0),
	  _maxRapidity           (this, "RAP_MAX",0),
	  _nmbRapidityBins       (this, "RAP_N_BINS",0),
	  _ptCutEnabled          (this, "CUT_PT",false, NOT_REQUIRED),
	  _ptCutMin              (this, "PT_MIN",0, NOT_REQUIRED),
	  _ptCutMax              (this, "PT_MAX",0, NOT_REQUIRED),
	  _etaCutEnabled         (this, "CUT_ETA",false, NOT_REQUIRED),
	  _etaCutMin             (this, "ETA_MIN",0, NOT_REQUIRED),
	  _etaCutMax             (this, "ETA_MAX",0, NOT_REQUIRED),
	  _productionMode        (this, "PROD_MODE",0),
	  _nmbEventsTot          (this, "N_EVENTS",0),
	  _prodParticleId        (this, "PROD_PID",0),
	  _randomSeed            (this, "RND_SEED",0),
	  _beamBreakupMode       (this, "BREAKUP_MODE",0),
	  _interferenceEnabled   (this, "INTERFERENCE",false),
	  _interferenceStrength  (this, "IF_STRENGTH",0),
	  _maxPtInterference     (this, "INT_PT_MAX",0),
	  _nmbPtBinsInterference (this, "INT_PT_N_BINS",0),
	  _ptBinWidthInterference(this, "INT_PT_WIDTH",0),
	  _protonEnergy          (this, "PROTON_ENERGY",0, NOT_REQUIRED),
	  _minGammaEnergy	 (this, "MIN_GAMMA_ENERGY",6.0, NOT_REQUIRED),
	  _maxGammaEnergy	 (this, "MAX_GAMMA_ENERGY",600000.0, NOT_REQUIRED),
	  _pythiaParams          (this, "PYTHIA_PARAMS","", NOT_REQUIRED),
	  _pythiaFullEventRecord (this, "PYTHIA_FULL_EVENTRECORD",false, NOT_REQUIRED),
	  _xsecCalcMethod	 (this, "XSEC_METHOD",0, NOT_REQUIRED),
          _axionMass             (this, "AXION_MASS",50, NOT_REQUIRED),  // AXION HACK
          _bslopeDefinition      (this, "BSLOPE_DEFINITION",0, NOT_REQUIRED),
	  _bslopeValue           (this, "BSLOPE_VALUE",4.0,NOT_REQUIRED),
          _bslope0               (this, "BSLOPE0",4.63,NOT_REQUIRED),
          _bslope_alphaprime     (this, "BSLOPE_ALPHAPRIME",0.164,NOT_REQUIRED),
	  _printVM               (this, "PRINT_VM",0,NOT_REQUIRED),
	  _impulseVM             (this, "SELECT_IMPULSE_VM",0,NOT_REQUIRED),
	  _quantumGlauber        (this, "QUANTUM_GLAUBER",0,NOT_REQUIRED),
	  _bmin                  (this, "BMIN",0,NOT_REQUIRED),
          _bmax                  (this, "BMAX",0,NOT_REQUIRED),

          _deuteronSlopePar      (this, "deuteronSlopePar"      , 9.5           , NOT_REQUIRED),
          _protonMass            (this, "protonMass"            , 0.938272046   , NOT_REQUIRED),
          _pionChargedMass       (this, "pionChargedMass"       , 0.13957018    , NOT_REQUIRED),
          _pionNeutralMass       (this, "pionNeutralMass"       , 0.1349766     , NOT_REQUIRED),
          _kaonChargedMass       (this, "kaonChargedMass"       , 0.493677      , NOT_REQUIRED),
          _mel                   (this, "mel"                   , 0.000510998928, NOT_REQUIRED),
          _muonMass              (this, "muonMass"              , 0.1056583715  , NOT_REQUIRED),
          _tauMass               (this, "tauMass"               , 1.77682       , NOT_REQUIRED),
          _f0Mass                (this, "f0Mass"                , 0.990         , NOT_REQUIRED),
          _f0Width               (this, "f0Width"               , 0.100         , NOT_REQUIRED),
          _f0BrPiPi              (this, "f0BrPiPi"              , 1.0           , NOT_REQUIRED),
          _etaMass               (this, "etaMass"               , 0.547862      , NOT_REQUIRED),
          _etaWidth              (this, "etaWidth"              , 0.00000131    , NOT_REQUIRED),
          _etaPrimeMass          (this, "etaPrimeMass"          , 0.95778       , NOT_REQUIRED),
          _etaPrimeWidth         (this, "etaPrimeWidth"         , 0.000198      , NOT_REQUIRED),
          _etaCMass              (this, "etaCMass"              , 2.9836        , NOT_REQUIRED),
          _etaCWidth             (this, "etaCWidth"             , 0.0322        , NOT_REQUIRED),
          _f2Mass                (this, "f2Mass"                , 1.2751        , NOT_REQUIRED),
          _f2Width               (this, "f2Width"               , 0.1851        , NOT_REQUIRED),
          _f2BrPiPi              (this, "f2BrPiPi"              , 0.561         , NOT_REQUIRED),
          _a2Mass                (this, "a2Mass"                , 1.3183        , NOT_REQUIRED),
          _a2Width               (this, "a2Width"               , 0.105         , NOT_REQUIRED),
          _f2PrimeMass           (this, "f2PrimeMass"           , 1.525         , NOT_REQUIRED),
          _f2PrimeWidth          (this, "f2PrimeWidth"          , 0.073         , NOT_REQUIRED),
          _f2PrimeBrKK           (this, "f2PrimeBrKK"           , 0.887         , NOT_REQUIRED),
          _zoverz03Mass          (this, "zoverz03Mass"          , 1.540         , NOT_REQUIRED),
          _f0PartialggWidth      (this, "f0PartialggWidth"      , 0.29E-6       , NOT_REQUIRED),
          _etaPartialggWidth     (this, "etaPartialggWidth"     , 0.516E-6      , NOT_REQUIRED),
          _etaPrimePartialggWidth(this, "etaPrimePartialggWidth", 4.35E-6       , NOT_REQUIRED),
          _etaCPartialggWidth    (this, "etaCPartialggWidth"    , 5.0E-6        , NOT_REQUIRED),
          _f2PartialggWidth      (this, "f2PartialggWidth"      , 3.03E-6       , NOT_REQUIRED),
          _a2PartialggWidth      (this, "a2PartialggWidth"      , 1.0E-6        , NOT_REQUIRED),
          _f2PrimePartialggWidth (this, "f2PrimePartialggWidth" , 0.081E-6      , NOT_REQUIRED),
          _zoverz03PartialggWidth(this, "zoverz03PartialggWidth", 0.1E-6        , NOT_REQUIRED),
          _f0Spin                (this, "f0Spin"                , 0.0           , NOT_REQUIRED),
          _etaSpin               (this, "etaSpin"               , 0.0           , NOT_REQUIRED),
          _etaPrimeSpin          (this, "etaPrimeSpin"          , 0.0           , NOT_REQUIRED),
          _etaCSpin              (this, "etaCSpin"              , 0.0           , NOT_REQUIRED),
          _f2Spin                (this, "f2Spin"                , 2.0           , NOT_REQUIRED),
          _a2Spin                (this, "a2Spin"                , 2.0           , NOT_REQUIRED),
          _f2PrimeSpin           (this, "f2PrimeSpin"           , 2.0           , NOT_REQUIRED),
          _zoverz03Spin          (this, "zoverz03Spin"          , 2.0           , NOT_REQUIRED),
          _axionSpin             (this, "axionSpin"             , 0.0           , NOT_REQUIRED),
          _rho0Mass              (this, "rho0Mass"              , 0.769         , NOT_REQUIRED),
          _rho0Width             (this, "rho0Width"             , 0.1517        , NOT_REQUIRED),
          _rho0BrPiPi            (this, "rho0BrPiPi"            , 1.0           , NOT_REQUIRED),
          _rho0PrimeMass         (this, "rho0PrimeMass"         , 1.540         , NOT_REQUIRED),
          _rho0PrimeWidth        (this, "rho0PrimeWidth"        , 0.570         , NOT_REQUIRED),
          _rho0PrimeBrPiPi       (this, "rho0PrimeBrPiPi"       , 1.0           , NOT_REQUIRED),
          _OmegaMass             (this, "OmegaMass"             , 0.78265       , NOT_REQUIRED),
          _OmegaWidth            (this, "OmegaWidth"            , 0.00849       , NOT_REQUIRED),
          _OmegaBrPiPi           (this, "OmegaBrPiPi"           , 0.0153        , NOT_REQUIRED),
          _PhiMass               (this, "PhiMass"               , 1.019461      , NOT_REQUIRED),
          _PhiWidth              (this, "PhiWidth"              , 0.004266      , NOT_REQUIRED),
          _PhiBrKK               (this, "PhiBrKK"               , 0.489         , NOT_REQUIRED),
          _JpsiMass              (this, "JpsiMass"              , 3.096916      , NOT_REQUIRED),
          _JpsiWidth             (this, "JpsiWidth"             , 0.0000929     , NOT_REQUIRED),
          _JpsiBree              (this, "JpsiBree"              , 0.05971       , NOT_REQUIRED),
          _JpsiBrmumu            (this, "JpsiBrmumu"            , 0.05961       , NOT_REQUIRED),
          _JpsiBrppbar           (this, "JpsiBrppbar"           , 0.002120      , NOT_REQUIRED),
          _Psi2SMass             (this, "Psi2SMass"             , 3.686109      , NOT_REQUIRED),
          _Psi2SWidth            (this, "Psi2SWidth"            , 0.000299      , NOT_REQUIRED),
          _Psi2SBree             (this, "Psi2SBree"             , 0.00789       , NOT_REQUIRED),
          _Psi2SBrmumu           (this, "Psi2SBrmumu"           , 0.0079        , NOT_REQUIRED),
          _Upsilon1SMass         (this, "Upsilon1SMass"         , 9.46030       , NOT_REQUIRED),
          _Upsilon1SWidth        (this, "Upsilon1SWidth"        , 0.00005402    , NOT_REQUIRED),
          _Upsilon1SBree         (this, "Upsilon1SBree"         , 0.0238        , NOT_REQUIRED),
          _Upsilon1SBrmumu       (this, "Upsilon1SBrmumu"       , 0.0248        , NOT_REQUIRED),
          _Upsilon2SMass         (this, "Upsilon2SMass"         , 10.02326      , NOT_REQUIRED),
          _Upsilon2SWidth        (this, "Upsilon2SWidth"        , 0.00003198    , NOT_REQUIRED),
          _Upsilon2SBree         (this, "Upsilon2SBree"         , 0.0191        , NOT_REQUIRED),
          _Upsilon2SBrmumu       (this, "Upsilon2SBrmumu"       , 0.0193        , NOT_REQUIRED),
          _Upsilon3SMass         (this, "Upsilon3SMass"         , 10.3552       , NOT_REQUIRED),
          _Upsilon3SWidth        (this, "Upsilon3SWidth"        , 0.00002032    , NOT_REQUIRED),
          _Upsilon3SBree         (this, "Upsilon3SBree"         , 0.0218        , NOT_REQUIRED),
          _Upsilon3SBrmumu       (this, "Upsilon3SBrmumu"       , 0.0218        , NOT_REQUIRED)
{
  // All parameters must be initialised in initialisation list!
  // If not: error: 'parameter<T, validate>::parameter() [with T = unsigned int, bool validate = true]' is private
  // or similar
}


//______________________________________________________________________________
inputParameters::~inputParameters()
{
}


//______________________________________________________________________________
bool
inputParameters::configureFromFile(const std::string &_configFileName)
{

	int nParameters = getInputParser().parseFile(_configFileName);

	if(nParameters == -1)
	{
		printWarn << "could not open file '" << _configFileName << "'" << endl;
	  return false;
	}


	if(getInputParser().validateParameters(cerr))
		printInfo << "successfully read input parameters from '" << _configFileName << "'" << endl;
	else {
 		printWarn << "problems reading input parameters from '" << _configFileName << "'" << endl
 		          << *this;
 		return false;
 	}

	/// temporary output test
	printWarn<< * this<<endl;
	return true;
}
 bool inputParameters::init()
 {

 	// Calculate beam gamma in CMS frame
 	double rap1 = acosh(beam1LorentzGamma());
	double rap2 = -acosh(beam2LorentzGamma());
	_beamLorentzGamma = cosh((rap1-rap2)/2);

	std::cout << "Rapidity beam 1: " << rap1 << ", rapidity beam 2: " << rap2 << ", rapidity CMS system: " << (rap1+rap2)/2 << ", beam gamma in CMS: " << _beamLorentzGamma<< std::endl;
	_ptBinWidthInterference = maxPtInterference() / nmbPtBinsInterference();
	_protonEnergy           = _beamLorentzGamma * protonMass();

	// check for deuteron or tritium - these must be the second beam
	if((beam1Z()==1) && (beam1A()==2)){
		if((beam2Z()==1) && (beam2A()==2)){
		printWarn << "deuteron-deuteron collisions are not supported" << endl;
		return false;}
		printWarn << "deuterium must be listed as the second nucleus" << endl;
		return false;}

	if( ((beam1Z()==1) && (beam1A()==3)) || ((beam2Z()==1) && (beam2A()==3)) ){
		printWarn << "tritium is not currently supported" << endl;
		return false;}
	// check that rho production uses wide resonance option
	if(_prodParticleId.value()==113 && _productionMode.value()==2){
		printWarn << endl<< "For rho meson production, you should choose the wide resonance option (production mode = 3)" << endl;
		return false;}

	// define interaction type
	switch (productionMode()) {
	case 1:
		_interactionType = PHOTONPHOTON;
		break;
	case 2:
		_interactionType = PHOTONPOMERONNARROW;
		break;
	case 3:
		_interactionType = PHOTONPOMERONWIDE;
		break;
        case 4:
                _interactionType = PHOTONPOMERONINCOHERENT;
                break;
	case 5:
		_interactionType = PHOTONUCLEARSINGLE;
		break;
	case 6:
		_interactionType = PHOTONUCLEARDOUBLE;
		break;
	case 7:
		_interactionType = PHOTONUCLEARSINGLEPA;
		break;
	case 8:
		_interactionType = PHOTONUCLEARSINGLEPAPY;
		break;
// 	case 9:
// 		_interactionType = PHOTONPHOTONKNIEHL;
// 		break;
//         case 10:
//                 _interactionType = PHOTONPHOTONKNIEHLMODIFIED;
//                 break;
	default:
		printWarn << "unknown production mode '" << _productionMode << "'" << endl;
		return false;
	}
	if( (_productionMode.value() ==4) && (_interferenceEnabled.value())) {
		printWarn << " cannot enable interference for incoherent production " << endl;
		return false;
	}

	double mass        = 0;
	double width       = 0;
	double defaultMinW = 0;  // default for _minW, unless it is defined later [GeV/c^2]
	double defaultMaxW = 0;  // default for _maxW, unless it is defined later [GeV/c^2]
	switch (prodParticleId()) {
	case 11:  // e+e- pair
		_particleType = ELECTRON;
		_decayType    = LEPTONPAIR;
		mass          = mel();
		defaultMinW   = 2*mass; // default is 0.01; up to 0.15 is safe for Summer 2000 triggering for e+e- pairs
                defaultMaxW     = sqrt(beam1LorentzGamma()*beam2LorentzGamma())*2*(starlightConstants::hbarc)/(1.2*pow(float(beam1A()),1./6.)*pow(float(beam2A()),1./6.)); // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
	case 13:  // mu+mu- pair
		_particleType = MUON;
		_decayType    = LEPTONPAIR;
		defaultMinW   = 2 * muonMass();
                defaultMaxW     = sqrt(beam1LorentzGamma()*beam2LorentzGamma())*2*(starlightConstants::hbarc)/(1.2*pow(float(beam1A()),1./6.)*pow(float(beam2A()),1./6.)); // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
	case 15:  // tau+tau- pair
		_particleType = TAUON;
		_decayType    = LEPTONPAIR;
		defaultMinW   = 2 * tauMass();
                defaultMaxW     = sqrt(beam1LorentzGamma()*beam2LorentzGamma())*2*(starlightConstants::hbarc)/(1.2*pow(float(beam1A()),1./6.)*pow(float(beam2A()),1./6.)); // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
	case 10015:  // tau+tau- pair
		_particleType = TAUONDECAY;
		_decayType    = LEPTONPAIR;
		defaultMinW   = 2 * tauMass();
                defaultMaxW   = sqrt(beam1LorentzGamma()*beam2LorentzGamma())*2*(starlightConstants::hbarc)/(1.2*pow(float(beam1A()),1./6.)*pow(float(beam2A()),1./6.)); // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;

// 	case 24:  // W+W- pair
// 		_particleType = W;
// 		_decayType    = WW;
// 		defaultMinW   = 2 * muonMass();
// 		break;
	case 115:  // a_2(1320)
		_particleType = A2;
		_decayType    = SINGLEMESON;
		mass          = a2Mass();
		width         = a2Width();
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW   = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
	case 221:  // eta
		_particleType = ETA;
		_decayType    = SINGLEMESON;
		mass          = etaMass();
		width         = etaWidth();
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW   = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
	case 225:  // f_2(1270)
		_particleType = F2;
		_decayType    = SINGLEMESON;
		mass          = f2Mass();
		width         = f2Width();
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = f2BrPiPi();
		break;
	case 331:  // eta'(958)
		_particleType = ETAPRIME;
		_decayType    = SINGLEMESON;
		mass          = etaPrimeMass();
		width         = etaPrimeWidth();
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
	case 335:  // f_2'(1525)
		_particleType = F2PRIME;
		_decayType    = SINGLEMESON;
		mass          = f2PrimeMass();
		width         = f2PrimeWidth();
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = f2PrimeBrKK();
		break;
	case 441:  // eta_c(1s)
		_particleType = ETAC;
		_decayType    = SINGLEMESON;
		mass          = etaCMass();
		width         = etaCWidth();
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
	case 9010221:  // f_0(980), was orginally called 10221? updated to standard number
		_particleType = F0;
		_decayType    = SINGLEMESON;
		mass          = f0Mass();
		width         = f0Width();
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
               _inputBranchingRatio = f0BrPiPi();
		break;
	case 33:  // Z"/Z0  This is the rho^0 rho^0 final state SRK
		_particleType = ZOVERZ03;
		_decayType    = SINGLEMESON;
                defaultMinW   = 4*pionChargedMass();
                defaultMaxW         = 1.6; // JES 6.17.2015 to avoid problems with no default
                _inputBranchingRatio = 1.0;
		break;
        case 88:  // axion// AXION HACK, till break statement
                _particleType = AXION;
                _decayType    = SINGLEMESON;
                mass          = _axionMass.value();
                width         = 1/(64*starlightConstants::pi)*mass*mass*mass/(1000*1000);//Fix Lambda=1000 GeV, rescaling is trivial.
                defaultMinW   = mass - 5*width; // JES 6.17.2015 to avoid problems with default of 0
                defaultMaxW         = mass + 5*width; // JES 6.17.2015 to avoid problems with no default
                break;    // AXION HACK, end
// 	case 25: // Higgs
// 		_particleType = HIGGS;
// 		_decayType    = SINGLEMESON;
// 		break;
	case 113:  // rho(770)
		_particleType = RHO;
		_decayType    = WIDEVMDEFAULT;
		mass          = rho0Mass();
		width         = rho0Width();
		defaultMinW   = 2 * pionChargedMass();
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = rho0BrPiPi();
		break;
	case 913:  // rho(770) with direct pi+pi- decay, interference given by ZEUS data
		_particleType = RHOZEUS;
		_decayType    = WIDEVMDEFAULT;
		mass          = rho0Mass();
		width         = rho0Width();
		defaultMinW   = 2 * pionChargedMass();
		defaultMaxW         = mass + 5 * width;  // use the same 1.5GeV max mass as ZEUS
		_inputBranchingRatio = rho0BrPiPi();
		break;
	case 999:  // pi+pi-pi+pi- phase space decay
		_particleType = FOURPRONG;
		_decayType    = WIDEVMDEFAULT;
		mass          = rho0PrimeMass();
		width         = rho0PrimeWidth();
		defaultMinW   = 4 * pionChargedMass();
                defaultMaxW     = sqrt(beam1LorentzGamma()*beam2LorentzGamma())*2*(starlightConstants::hbarc)/(1.2*pow(float(beam1A()),1./6.)*pow(float(beam2A()),1./6.)); // JES 6.17.2015 to avoid problems with no default
		_inputBranchingRatio = 1.0;
		break;
	case 223:  // omega(782)
		_particleType = OMEGA;
		_decayType    = NARROWVMDEFAULT;
		mass          = OmegaMass();
		width         = OmegaWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = OmegaBrPiPi();
		break;
	case 333:  // phi(1020)
		_particleType = PHI;
		_decayType    = NARROWVMDEFAULT;
		mass          = PhiMass();
		width         = PhiWidth();
		defaultMinW   = 2 * kaonChargedMass();
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = PhiBrKK();
		break;
	case 443:  // J/psi
		_particleType = JPSI;
		_decayType    = NARROWVMDEFAULT;
		mass          = JpsiMass();
		width         = JpsiWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = (JpsiBree() + JpsiBrmumu())/2.;
		break;
   	case 443011:  // J/psi
		_particleType = JPSI_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = JpsiMass();
		width         = JpsiWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = JpsiBree();
		break;
	case 443013:  // J/psi
	        cout<<"In inputParameters setting J/psi mass!"<<endl;
		_particleType = JPSI_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = JpsiMass();
		width         = JpsiWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = JpsiBrmumu();
		break;
	case 4432212:  // J/psi
	        cout<<"In inputParameters setting J/psi mass!"<<endl;
		_particleType = JPSI_ppbar;
		_decayType    = NARROWVMDEFAULT;
		mass          = JpsiMass();
		width         = JpsiWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = JpsiBrppbar();
		break;
	case 444:  // psi(2S)
		_particleType = JPSI2S;
		_decayType    = NARROWVMDEFAULT;
		mass          = Psi2SMass();
		width         = Psi2SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW         = mass + 5 * width;
		_inputBranchingRatio = (Psi2SBree() + Psi2SBrmumu())/2.;
		break;
	case 444011: // psi(2S)
		_particleType = JPSI2S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = Psi2SMass();
		width         = Psi2SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Psi2SBree();
		break;
	case 444013:  // psi(2S)
		_particleType = JPSI2S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = Psi2SMass();
		width         = Psi2SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Psi2SBrmumu();
		break;
	case 553:  // Upsilon(1S)
		_particleType = UPSILON;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon1SMass();
		width         = Upsilon1SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = (Upsilon1SBree() + Upsilon1SBrmumu())/2.;
		break;
	case 553011:  // Upsilon
		_particleType = UPSILON_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon1SMass();
		width         = Upsilon1SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Upsilon1SBree();
		break;
	case 553013:  // Upsilon
		_particleType = UPSILON_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon1SMass();
		width         = Upsilon1SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Upsilon1SBrmumu();
		break;
	case 554:  // Upsilon(2S)
		_particleType = UPSILON2S;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon2SMass();
		width         = Upsilon2SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = (Upsilon2SBree() + Upsilon2SBrmumu())/2.;
		break;
	case 554011:  // Upsilon(2S)
		_particleType = UPSILON2S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon2SMass();
		width         = Upsilon2SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Upsilon2SBree();
		break;
	case 554013:  // Upsilon(2S)
		_particleType = UPSILON2S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon2SMass();
		width         = Upsilon2SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Upsilon2SBrmumu();
		break;
	case 555:  // Upsilon(3S)
		_particleType = UPSILON3S;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon3SMass();
		width         = Upsilon3SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = (Upsilon3SBree() + Upsilon3SBrmumu())/2.;
		break;
	case 555011:  // Upsilon(3S)
		_particleType = UPSILON3S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon3SMass();
		width         = Upsilon3SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Upsilon3SBree();
		break;
	case 555013:  // Upsilon(3S)
		_particleType = UPSILON3S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = Upsilon3SMass();
		width         = Upsilon3SWidth();
		defaultMinW   = mass - 5 * width;
		defaultMaxW   = mass + 5 * width;
		_inputBranchingRatio = Upsilon3SBrmumu();
		break;
	default:
		printWarn << "unknown particle ID " << _prodParticleId << endl;
		return false;
	}  // _prodParticleId

	if (_minW.value() == -1)
		_minW = defaultMinW;
	if (_maxW.value() == -1)
		_maxW = defaultMaxW;
	if ( _maxW.value() <= _minW.value() ) {
		printWarn << "maxW must be greater than minW" << endl;
		return false;
	}

	printInfo << "using the following " << *this;

	return true;
}


//______________________________________________________________________________
ostream&
inputParameters::print(ostream& out) const
{
	out << "starlight parameters:" << endl
	    << "    base file name  ...................... '"  << _baseFileName.value() << "'" << endl
	    << "    beam 1 atomic number ................... " << _beam1Z.value() << endl
	    << "    beam 1 atomic mass number .............. " << _beam1A.value() << endl
	    << "    beam 2 atomic number ................... " << _beam2Z.value() << endl
	    << "    beam 2 atomic mass number .............. " << _beam2A.value() << endl
	    << "    Lorentz gamma of beams in CM frame ..... " << _beamLorentzGamma << endl
	    << "    mass W of produced hadronic system ..... " << _minW.value() << " < W < " << _maxW.value() << " GeV/c^2" << endl
	    << "    # of W bins ............................ " << _nmbWBins.value() << endl
	    << "    maximum absolute value for rapidity .... " << _maxRapidity.value() << endl
	    << "    # of rapidity bins ..................... " << _nmbRapidityBins.value() << endl
	    << "    cut in pT............................... " << yesNo(_ptCutEnabled.value()) << endl;
    if (_ptCutEnabled.value()) {
	out << "        minumum pT.......................... " << _ptCutMin.value() << " GeV/c" << endl
	    << "        maximum pT.......................... " << _ptCutMax.value() << " GeV/c" << endl;}
	out << "    cut in eta.............................. " << yesNo(_etaCutEnabled.value()) << endl;
    if (_etaCutEnabled.value()) {
	out << "        minumum eta......................... " << _etaCutMin.value() << endl
	    << "        maximum eta......................... " << _etaCutMax.value() << endl;}
        out << "    production mode ........................ " << _productionMode.value() << endl
	    << "    number of events to generate ........... " << _nmbEventsTot.value() << endl
	    << "    PDG ID of produced particle ............ " << _prodParticleId.value() << endl
	    << "    seed for random generator .............. " << _randomSeed.value() << endl
	    << "    breakup mode for beam particles ........ " << _beamBreakupMode.value() << endl
	    << "    interference enabled ................... " << yesNo(_interferenceEnabled.value()) << endl;
    if (_interferenceEnabled.value()) {
	out << "    interference strength .................. " << _interferenceStrength.value() << endl
	    << "    maximum p_T for interference calc. ..... " << _maxPtInterference.value() << " GeV/c" << endl
	    << "    # of p_T bins for interference calc. ... " << _nmbPtBinsInterference.value() << endl;}
    if (_productionMode.value()!=1) {
		if (_productionMode.value()==4) {
	 	  out  << "    coherent scattering off nucleus ........ no" << endl;}
		else {
	 	  out  << "    coherent scattering off nucleus ........ yes" << endl;
		}
	}
    out     <<"    Quantum Glauber parameter...............  "<<_quantumGlauber.value()<<endl;
    out     <<"    Impulse VM parameter....................  "<<_impulseVM.value()<<endl;
    //   if (_quantumGlauber.value()==1) {out << "    Quantum Glauber calculation being used"<< endl;}
    //if (_quantumGlauber.value()!=1) {out << "    Classical Glauber calculation being used"<< endl;}
    if (_beamBreakupMode.value()==8) {
      out <<"    Minimum impact parameter.................."<<_bmin.value()<<" fm"<<endl;
      out <<"    Maximum impact parameter.................."<<_bmax.value()<<" fm"<<endl;
    }

    // Add some checks here  SRK September, 2017
    if (_beamBreakupMode.value()==8 && _bmin.value() > _bmax.value()) {
      out <<"****Error bmin > bmax"<<endl;
      exit(-1);
    }
    if (_beamBreakupMode.value() == 8 && _productionMode.value() !=1) {
	out <<"****Error.  Cannot use beam breakup mode 8 (with bmin and bmax) with photonuclear interactions"<<endl;
	exit(-1);
    }
	return out;
}


//______________________________________________________________________________
ostream&
inputParameters::write(ostream& out) const
{

        out << "baseFileName"  << baseFileName         () <<endl
	    << "BEAM_1_Z"      << beam1Z               () <<endl
	    << "BEAM_2_Z"      << beam1A               () <<endl
	    << "BEAM_1_A"      << beam2Z               () <<endl
	    << "BEAM_2_A"      << beam2A               () <<endl
	    << "BEAM_GAMMA"    << beamLorentzGamma     () <<endl
	    << "W_MAX"         << maxW                 () <<endl
	    << "W_MIN"         << minW                 () <<endl
	    << "W_N_BINS"      << nmbWBins             () <<endl
	    << "RAP_MAX"       << maxRapidity          () <<endl
	    << "RAP_N_BINS"    << nmbRapidityBins      () <<endl
	    << "CUT_PT"        << ptCutEnabled         () <<endl
	    << "PT_MIN"        << ptCutMin             () <<endl
	    << "PT_MAX"        << ptCutMax             () <<endl
	    << "CUT_ETA"       << etaCutEnabled        () <<endl
	    << "ETA_MIN"       << etaCutMin            () <<endl
	    << "ETA_MAX"       << etaCutMax            () <<endl
	    << "PROD_MODE"     << productionMode       () <<endl
	    << "N_EVENTS"      << nmbEvents            () <<endl
	    << "PROD_PID"      << prodParticleId       () <<endl
	    << "RND_SEED"      << randomSeed           () <<endl
	    << "BREAKUP_MODE"  << beamBreakupMode      () <<endl
	    << "INTERFERENCE"  << interferenceEnabled  () <<endl
	    << "IF_STRENGTH"   << interferenceStrength () <<endl
	    << "INT_PT_MAX"    << maxPtInterference    () <<endl
	    << "INT_PT_N_BINS" << nmbPtBinsInterference() <<endl
	    << "QUANTUM_GLAUBER"<<quantumGlauber       () <<endl
	    << "BMIN"          <<bmin                  ()<<endl
	    << "BMAX"          <<bmax                  ()<<endl;

	return out;
}

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


#ifndef INPUTPARAMETERS_H
#define INPUTPARAMETERS_H


#include "starlightconstants.h"
#include "inputParser.h"
#include <string>
#include <ostream>
#include <vector>
#include <sstream>
#include <map>
#include <cstdlib>

class parameterbase {
public:
  virtual ~parameterbase() {}

  virtual std::string validationkey() const = 0;

protected:
private:
} ;

class inputParametersBase {
public:
  inputParametersBase()
    : _ip() {}
  virtual ~inputParametersBase() {}

  inputParser& getInputParser() {
    return _ip;
  }
  void addParPtr(const parameterbase* ptr) {
    _pv.push_back(ptr);
  }

  std::string parameterValueKey() const
  {
    std::string result;
    for (std::vector<const parameterbase* >::const_iterator i=_pv.begin(), end=_pv.end(); i!=end; ++i) {
      const parameterbase* p = *i;
      if (p)
	result += p->validationkey();
    }
    return result;
  }

protected:
private:
  inputParser _ip;  //! do not save
  std::vector<const parameterbase* > _pv; //! pointers to member variables in inputParameters
} ;


// The  parameter class
// validate parameter specifies if the parameter should be a part of the validity check of the current parameters
template<typename T, bool validate>
class parameter : public parameterbase
{
public:
    parameter(inputParametersBase *ip_base, const std::string &name, T value, bool required = true)
      : parameterbase()
      , _name(name)
      , _value(value)
      , _validate(validate)
      , _required(required) {
      if (ip_base) {
      	ip_base->getInputParser().addParameter(*this);
	ip_base->addParPtr(this);
      }
    }

    parameter &operator=(T v) { _value = v; return *this;}
    T* ptr() const {
        return const_cast<T*>(&_value);
    }

    std::string toString() const
    {
        std::stringstream s;
        s << _value;
        return s.str();
    }

    T value() const { return _value; }

    std::string name() const { return _name;}

    bool required() const { return _required; }

    void setValue(T v) { _value = v; }

    void setName(std::string name) { _name = name; }

    void setRequired(bool r) { _required = r; }

    // Validation key for this parameter
    std::string validationkey() const
    {
        return (_validate ? _name + ":" + toString() + "-" : std::string(""));
    }

    template<typename S, bool v>
    inline friend std::ostream& operator<<(std::ostream& os, const parameter<S,v>& par);

 private:
    std::string _name;

    T _value; // Value
    bool _validate; // true if a change in the parameter invalidates x-sec tables
    bool _required; // true if this is required option.
};

template<typename S, bool v>
inline std::ostream& operator<<(std::ostream& os, const parameter<S,v>& par)
{
    os << par._value;
    return os;
}

class inputParameters : public inputParametersBase {

public:
	inputParameters();
	~inputParameters();

	bool init();
	bool configureFromFile(const std::string &configFileName = "./config/slight.in");

        std::string  baseFileName          () const { return _baseFileName.value();           }
	unsigned int beam1Z                () const { return _beam1Z.value();                 }  ///< returns atomic number of beam particle 1
	unsigned int beam1A                () const { return _beam1A.value();                 }  ///< returns atomic mass number of beam particle 1
	unsigned int beam2Z                () const { return _beam2Z.value();                 }  ///< returns atomic number of beam particle 2
	unsigned int beam2A                () const { return _beam2A.value();                 }  ///< returns atomic mass number of beam particle 2
	double       beamLorentzGamma      () const { return _beamLorentzGamma;       	      }  ///< returns Lorentz gamma factor of both beams in beam CMS frame
	double       beam1LorentzGamma     () const { return _beam1LorentzGamma.value();      }  ///< returns Lorentz gamma factor of beam 1 in collider frame
	double       beam2LorentzGamma     () const { return _beam2LorentzGamma.value();      }  ///< returns Lorentz gamma factor of beam 2 in collider frame
	double       maxW                  () const { return _maxW.value();                   }  ///< returns maximum mass W of produced hadronic system [GeV/c^2]
	double       minW                  () const { return _minW.value();                   }  ///< returns minimum mass W of produced hadronic system [GeV/c^2]
	unsigned int nmbWBins              () const { return _nmbWBins.value();               }  ///< returns number of W bins in lookup table
	double       maxRapidity           () const { return _maxRapidity.value();            }  ///< returns maximum absolute value of rapidity
	unsigned int nmbRapidityBins       () const { return _nmbRapidityBins.value();        }  ///< returns number of rapidity bins in lookup table
	bool         ptCutEnabled          () const { return _ptCutEnabled.value();           }  ///< returns cut in pt
	double       ptCutMin              () const { return _ptCutMin.value();               }  ///< returns minimum pt
	double       ptCutMax              () const { return _ptCutMax.value();               }  ///< returns maximum pt
	bool         etaCutEnabled         () const { return _etaCutEnabled.value();          }  ///< returns cut in eta
	double       etaCutMin             () const { return _etaCutMin.value();              }  ///< returns minimum eta
	double       etaCutMax             () const { return _etaCutMax.value();              }  ///< returns maximum eta
	int          productionMode        () const { return _productionMode.value();         }  ///< returns production mode
	unsigned int nmbEvents             () const { return _nmbEventsTot.value();           }  ///< returns total number of events to generate
	int          prodParticleId        () const { return _prodParticleId.value();         }  ///< returns PDG particle ID of produced particle
	int          randomSeed            () const { return _randomSeed.value();             }  ///< returns seed for random number generator
	int          beamBreakupMode       () const { return _beamBreakupMode.value();        }  ///< returns breakup mode for beam particles
	bool         interferenceEnabled   () const { return _interferenceEnabled.value();    }  ///< returns whether interference is taken into account
	double       interferenceStrength  () const { return _interferenceStrength.value();   }  ///< returns percentage of interference
	double       maxPtInterference     () const { return _maxPtInterference.value();      }  ///< returns maximum p_T for interference calculation [GeV/c]
	int          nmbPtBinsInterference () const { return _nmbPtBinsInterference.value();  }  ///< returns number of p_T bins for interference calculation
	double       ptBinWidthInterference() const { return _ptBinWidthInterference.value(); }  ///< returns width of p_T bins for interference calculation [GeV/c]
	double 	     minGammaEnergy        () const { return _minGammaEnergy.value();         }  ///< returns minimum gamma energy in case of photo nuclear processes [GeV]
	double       maxGammaEnergy        () const { return _maxGammaEnergy.value();         }  ///< returns maximum gamma energy in case of photo nuclear processes [GeV]
	std::string  pythiaParams          () const { return _pythiaParams.value();           }  ///< returns parameters to be passed to pythia
	bool         pythiaFullEventRecord () const { return _pythiaFullEventRecord.value();  }  ///< returns if the full pythia event record should be printed
	int	     xsecCalcMethod        () const { return _xsecCalcMethod.value();         }  ///< returns the method used for the x-sec calculation
        double       axionMass             () const { return _axionMass.value();              }  ///< returns axion mass //AXION HACK
	int          bslopeDefinition      () const { return _bslopeDefinition.value();       }  ///< returns the definition of b-slope
	double       bslopeValue           () const { return _bslopeValue.value();            }  ///< returns the value of b-slope
        double       bslope0               () const { return _bslope0.value();                }  ///<
        double       bslope_alphaprime     () const { return _bslope_alphaprime.value();      }  ///<
	int          impulseVM             () const { return _impulseVM.value();              }  ///< returns the impulseVM value
	int          quantumGlauber        () const { return _quantumGlauber.value();         }  ///< returns the quantum glauber value
	double       bmin                  () const { return _bmin.value();                   }  // returns the minimum impact parameter for BREAKUP_MODE==6
	double       bmax                  () const { return _bmax.value();                   }  // returns the maximum impact parameter for BREAKUP_MODE==6

	starlightConstants::particleTypeEnum    prodParticleType     () const { return _particleType;    }  ///< returns type of produced particle
	starlightConstants::decayTypeEnum       prodParticleDecayType() const { return _decayType;       }  ///< returns decay type of produced particle
	starlightConstants::interactionTypeEnum interactionType      () const { return _interactionType; }  ///< returns interaction type
	double protonEnergy                () const { return _protonEnergy.value(); }
        double inputBranchingRatio         () const { return _inputBranchingRatio; }

        double deuteronSlopePar      () const {return _deuteronSlopePar      .value();}
        double protonMass            () const {return _protonMass            .value();}
        double pionChargedMass       () const {return _pionChargedMass       .value();}
        double pionNeutralMass       () const {return _pionNeutralMass       .value();}
        double kaonChargedMass       () const {return _kaonChargedMass       .value();}
        double mel                   () const {return _mel                   .value();}
        double muonMass              () const {return _muonMass              .value();}
        double tauMass               () const {return _tauMass               .value();}
        double f0Mass                () const {return _f0Mass                .value();}
        double f0Width               () const {return _f0Width               .value();}
        double f0BrPiPi              () const {return _f0BrPiPi              .value();}
        double etaMass               () const {return _etaMass               .value();}
        double etaWidth              () const {return _etaWidth              .value();}
        double etaPrimeMass          () const {return _etaPrimeMass          .value();}
        double etaPrimeWidth         () const {return _etaPrimeWidth         .value();}
        double etaCMass              () const {return _etaCMass              .value();}
        double etaCWidth             () const {return _etaCWidth             .value();}
        double f2Mass                () const {return _f2Mass                .value();}
        double f2Width               () const {return _f2Width               .value();}
        double f2BrPiPi              () const {return _f2BrPiPi              .value();}
        double a2Mass                () const {return _a2Mass                .value();}
        double a2Width               () const {return _a2Width               .value();}
        double f2PrimeMass           () const {return _f2PrimeMass           .value();}
        double f2PrimeWidth          () const {return _f2PrimeWidth          .value();}
        double f2PrimeBrKK           () const {return _f2PrimeBrKK           .value();}
        double zoverz03Mass          () const {return _zoverz03Mass          .value();}
        double f0PartialggWidth      () const {return _f0PartialggWidth      .value();}
        double etaPartialggWidth     () const {return _etaPartialggWidth     .value();}
        double etaPrimePartialggWidth() const {return _etaPrimePartialggWidth.value();}
        double etaCPartialggWidth    () const {return _etaCPartialggWidth    .value();}
        double f2PartialggWidth      () const {return _f2PartialggWidth      .value();}
        double a2PartialggWidth      () const {return _a2PartialggWidth      .value();}
        double f2PrimePartialggWidth () const {return _f2PrimePartialggWidth .value();}
        double zoverz03PartialggWidth() const {return _zoverz03PartialggWidth.value();}
        double f0Spin                () const {return _f0Spin                .value();}
        double etaSpin               () const {return _etaSpin               .value();}
        double etaPrimeSpin          () const {return _etaPrimeSpin          .value();}
        double etaCSpin              () const {return _etaCSpin              .value();}
        double f2Spin                () const {return _f2Spin                .value();}
        double a2Spin                () const {return _a2Spin                .value();}
        double f2PrimeSpin           () const {return _f2PrimeSpin           .value();}
        double zoverz03Spin          () const {return _zoverz03Spin          .value();}
        double axionSpin             () const {return _axionSpin             .value();}
        double rho0Mass              () const {return _rho0Mass              .value();}
        double rho0Width             () const {return _rho0Width             .value();}
        double rho0BrPiPi            () const {return _rho0BrPiPi            .value();}
        double rho0PrimeMass         () const {return _rho0PrimeMass         .value();}
        double rho0PrimeWidth        () const {return _rho0PrimeWidth        .value();}
        double rho0PrimeBrPiPi       () const {return _rho0PrimeBrPiPi       .value();}
        double OmegaMass             () const {return _OmegaMass             .value();}
        double OmegaWidth            () const {return _OmegaWidth            .value();}
        double OmegaBrPiPi           () const {return _OmegaBrPiPi           .value();}
        double PhiMass               () const {return _PhiMass               .value();}
        double PhiWidth              () const {return _PhiWidth              .value();}
        double PhiBrKK               () const {return _PhiBrKK               .value();}
        double JpsiMass              () const {return _JpsiMass              .value();}
        double JpsiWidth             () const {return _JpsiWidth             .value();}
        double JpsiBree              () const {return _JpsiBree              .value();}
        double JpsiBrmumu            () const {return _JpsiBrmumu            .value();}
        double JpsiBrppbar           () const {return _JpsiBrppbar           .value();}
        double Psi2SMass             () const {return _Psi2SMass             .value();}
        double Psi2SWidth            () const {return _Psi2SWidth            .value();}
        double Psi2SBree             () const {return _Psi2SBree             .value();}
        double Psi2SBrmumu           () const {return _Psi2SBrmumu           .value();}
        double Upsilon1SMass         () const {return _Upsilon1SMass         .value();}
        double Upsilon1SWidth        () const {return _Upsilon1SWidth        .value();}
        double Upsilon1SBree         () const {return _Upsilon1SBree         .value();}
        double Upsilon1SBrmumu       () const {return _Upsilon1SBrmumu       .value();}
        double Upsilon2SMass         () const {return _Upsilon2SMass         .value();}
        double Upsilon2SWidth        () const {return _Upsilon2SWidth        .value();}
        double Upsilon2SBree         () const {return _Upsilon2SBree         .value();}
        double Upsilon2SBrmumu       () const {return _Upsilon2SBrmumu       .value();}
        double Upsilon3SMass         () const {return _Upsilon3SMass         .value();}
        double Upsilon3SWidth        () const {return _Upsilon3SWidth        .value();}
        double Upsilon3SBree         () const {return _Upsilon3SBree         .value();}
        double Upsilon3SBrmumu       () const {return _Upsilon3SBrmumu       .value();}

        void setBaseFileName          (std::string v )  {  _baseFileName = v;     }
	void setBeam1Z                (unsigned int v)  {  _beam1Z = v;           }  ///< sets atomic number of beam particle 1
	void setBeam1A                (unsigned int v)  {  _beam1A = v;           }  ///< sets atomic mass number of beam particle 1
	void setBeam2Z                (unsigned int v)  {  _beam2Z = v;           }  ///< sets atomic number of beam particle 2
	void setBeam2A                (unsigned int v)  {  _beam2A = v;           }  ///< sets atomic mass number of beam particle 2
	void setBeamLorentzGamma      (double v)  {  _beamLorentzGamma = v;       }  ///< sets Lorentz gamma factor of both beams in beam CMS frame
	void setBeam1LorentzGamma     (double v)  {  _beam1LorentzGamma = v;      }  ///< sets Lorentz gamma factor of beam 1 in collider frame
	void setBeam2LorentzGamma     (double v)  {  _beam2LorentzGamma = v;      }  ///< sets Lorentz gamma factor of beam 2 in collider frame
	void setMaxW                  (double v)  {  _maxW = v;                   }  ///< sets maximum mass W of produced hadronic system [GeV/c^2]
	void setMinW                  (double v)  {  _minW = v;                   }  ///< sets minimum mass W of produced hadronic system [GeV/c^2]
	void setNmbWBins              (unsigned int v)  {  _nmbWBins = v;         }  ///< sets number of W bins in lookup table
	void setMaxRapidity           (double v)  {  _maxRapidity = v;            }  ///< sets maximum absolute value of rapidity
	void setNmbRapidityBins       (unsigned int v)  {  _nmbRapidityBins = v;  }  ///< sets number of rapidity bins in lookup table
	void setPtCutEnabled          (bool v)  {  _ptCutEnabled = v;             }  ///< sets cut in pt
	void setPtCutMin              (double v)  {  _ptCutMin = v;               }  ///< sets minimum pt
	void setPtCutMax              (double v)  {  _ptCutMax = v;               }  ///< sets maximum pt
	void setEtaCutEnabled         (bool v)  {  _etaCutEnabled = v;            }  ///< sets cut in eta
	void setEtaCutMin             (double v)  {  _etaCutMin = v;              }  ///< sets minimum eta
	void setEtaCutMax             (double v)  {  _etaCutMax = v;              }  ///< sets maximum eta
	void setProductionMode        (int v)  {  _productionMode = v;            }  ///< sets production mode
	void setNmbEvents             (unsigned int v)  {  _nmbEventsTot = v;     }  ///< sets total number of events to generate
	void setProdParticleId        (int v)  {  _prodParticleId = v;            }  ///< sets PDG particle ID of produced particle
	void setRandomSeed            (int v)  {  _randomSeed = v;                }  ///< sets seed for random number generator
	void setBeamBreakupMode       (int v)  {  _beamBreakupMode = v;           }  ///< sets breakup mode for beam particles
	void setInterferenceEnabled   (bool v)  {  _interferenceEnabled = v;      }  ///< sets whether interference is taken into account
	void setInterferenceStrength  (double v)  {  _interferenceStrength = v;   }  ///< sets percentage of interference
	void setMaxPtInterference     (double v)  {  _maxPtInterference = v;      }  ///< sets maximum p_T for voiderference calculation [GeV/c]
	void setNmbPtBinsInterference (int v)  {  _nmbPtBinsInterference = v;     }  ///< sets number of p_T bins for interference calculation
	void setPtBinWidthInterference(double v)  {  _ptBinWidthInterference = v; }  ///< sets width of p_T bins for voiderference calculation [GeV/c]
	void setMinGammaEnergy        (double v)  {  _minGammaEnergy = v;         }  ///< sets minimum gamma energy in case of photo nuclear processes [GeV]
	void setMaxGammaEnergy        (double v)  {  _maxGammaEnergy = v;         }  ///< sets maximum gamma energy in case of photo nuclear processes [GeV]
	void setPythiaParams          (std::string v)  {  _pythiaParams = v;      }  ///< sets parameters to be passed to pythia
	void setPythiaFullEventRecord (bool v)  {  _pythiaFullEventRecord = v;    }  ///< sets if the full pythia event record should be prvoided
	void setXsecCalcMethod        (int v)  {  _xsecCalcMethod = v;            }  ///< sets the method used for the x-sec calculation
	void setAxionMass        (double v)  {  _axionMass = v;                   }  ///< sets axion mass    //AXION HACK
	void setbslopeDefinition      (int v)  {  _bslopeDefinition = v;          }  ///< sets the definition of b slope
        void setbslopeValue           (double v)  {  _bslopeValue = v;            }  ///< sets the value of b slope
        void setbslope0               (double v)  {  _bslope0 = v;                }  ///< sets the value of the b-slope parameterization
        void setbslope_alphaprime     (double v)  {  _bslope_alphaprime = v;      }  ///< sets the value of the b-slope parameterization

	int  printVM                  () const { return _printVM.value();         }  ///< returns the printVM value
	void setprintVM               (int v)  {  _printVM = v;                   }  ///< sets the value of _printVM
        void setimpulseVM             (int v)  {  _impulseVM = v;                 }  ///< sets the value of _impulseVM
	void setquantumGlauber        (int v)  {  _quantumGlauber = v;            }  ///< sets the value of quantum_glauber
	void setbmin             (double v)    {  _bmin=v;                        }  ///< sets the minimum impact parameter (for BREAKUP_MODE==6
	void setbmax             (double v)    {  _bmax=v;                        }  ///< sets the minimum impact parameter (for BREAKUP_MODE==6

	void setProdParticleType      (starlightConstants::particleTypeEnum v)    { _particleType = v;    }  ///< sets type of produced particle
	void setProdParticleDecayType (starlightConstants::decayTypeEnum v)       { _decayType = v;       }  ///< sets decay type of produced particle
	void setInteractionType       (starlightConstants::interactionTypeEnum v) { _interactionType = v; }  ///< sets interaction type

	void setProtonEnergy        (double v)    { _protonEnergy = v;            }  ///< sets proton energy

	inline bool setParameter(std::string expression);

	std::ostream& print(std::ostream& out) const;  ///< prints parameter summary
	std::ostream& write(std::ostream& out) const;  ///< writes parameters back to an ostream
private:

// To indicate if the crossection table should be re-calculated if parameter changes
#define VALIDITY_CHECK true
#define NO_VALIDITY_CHECK false

	std::string _configFileName;  ///< path to configuration file (default = ./config/slight.in)

	// config file parameters
        parameter<std::string,NO_VALIDITY_CHECK>   _baseFileName;
	parameter<unsigned int,VALIDITY_CHECK>     _beam1Z;                  ///< atomic number of beam particle 1
	parameter<unsigned int,VALIDITY_CHECK>     _beam1A;                  ///< atomic mass number of beam particle 1
	parameter<unsigned int,VALIDITY_CHECK>     _beam2Z;                  ///< atomic number of beam particle 2
	parameter<unsigned int,VALIDITY_CHECK>     _beam2A;                  ///< atomic mass number of beam particle 2
	parameter<double, VALIDITY_CHECK>          _beam1LorentzGamma;       ///< Lorentz gamma factor of beam 1 in collider frame
	parameter<double, VALIDITY_CHECK>          _beam2LorentzGamma;       ///< Lorentz gamma factor of beam 2 in collider frame
	parameter<double, VALIDITY_CHECK>          _maxW;                    ///< maximum mass W of produced hadronic system [GeV/c^2]
	parameter<double, VALIDITY_CHECK>          _minW;                    ///< minimum mass W of produced hadronic system; if set to -1 default value is taken [GeV/c^2]
	parameter<unsigned int, VALIDITY_CHECK>    _nmbWBins;                ///< number of W bins in lookup table
	parameter<double, VALIDITY_CHECK>          _maxRapidity;             ///< maximum absolute value of rapidity
	parameter<unsigned int, VALIDITY_CHECK>    _nmbRapidityBins;         ///< number of rapidity bins in lookup table
	parameter<bool, VALIDITY_CHECK>            _ptCutEnabled;            ///< en/disables cut in pt
	parameter<double, VALIDITY_CHECK>          _ptCutMin;                ///< minimum pt, if cut is enabled
	parameter<double, VALIDITY_CHECK>          _ptCutMax;                ///< maximum pt, if cut is enabled
	parameter<bool, VALIDITY_CHECK>            _etaCutEnabled;           ///< en/disables cut in eta
	parameter<double, VALIDITY_CHECK>          _etaCutMin;               ///< minimum eta, if cut is enabled
	parameter<double, VALIDITY_CHECK>          _etaCutMax;               ///< maximum eta, if cut is enabled
	parameter<unsigned int, VALIDITY_CHECK>    _productionMode;          ///< \brief production mode
	                                                                     ///<
	                                                                     ///< 1 = photon-photon fusion,
	                                                                     ///< 2 = narrow vector meson resonance in photon-Pomeron fusion,
	                                                                     ///< 3 = Breit-Wigner vector meson resonance in photon-Pomeron fusion
	parameter<unsigned int, VALIDITY_CHECK>    _nmbEventsTot;            ///< total number of events to generate
	parameter<unsigned int, VALIDITY_CHECK>    _prodParticleId;          ///< PDG particle ID of produced particle
	parameter<unsigned int, VALIDITY_CHECK>    _randomSeed;              ///< seed for random number generator
	                                                                     ///<
	                                                                     ///< 1 = ASCII
	                                                                     ///< 2 = GSTARtext,
	                                                                     ///< 3 = PAW ntuple (not working)
	parameter<unsigned int, VALIDITY_CHECK>    _beamBreakupMode;         ///< \brief breakup mode for beam particles
	                                                                     ///<
	                                                                     ///< 1 = hard sphere nuclei (b > 2R),
	                                                                     ///< 2 = both nuclei break up (XnXn),
	                                                                     ///< 3 = a single neutron from each nucleus (1n1n),
	                                                                     ///< 4 = neither nucleon breaks up (with b > 2R),
	                                                                     ///< 5 = no hadronic break up (similar to option 1, but with the actual hadronic interaction)
	                                                                     ///  6 = set impact parameter range via bmax and bmin
	parameter<bool, VALIDITY_CHECK>            _interferenceEnabled;     ///< if VALIDITY_CHECK, interference is taken into account
	parameter<double, VALIDITY_CHECK>          _interferenceStrength;    ///< percentage of interference: from 0 = none to 1 = full
	parameter<double, VALIDITY_CHECK>          _maxPtInterference;       ///< maximum p_T for interference calculation [GeV/c]
	parameter<unsigned int, VALIDITY_CHECK>    _nmbPtBinsInterference;   ///< number of p_T bins for interference calculation
	parameter<double, VALIDITY_CHECK>          _ptBinWidthInterference;  ///< width of p_T bins for interference calculation [GeV/c]
	parameter<double, VALIDITY_CHECK>          _protonEnergy;
	parameter<double, VALIDITY_CHECK>          _minGammaEnergy;          ///< minimum gamma energy in case of photo nuclear processes [GeV]
	parameter<double, VALIDITY_CHECK>          _maxGammaEnergy;          ///< maximum gamma energy in case of photo nuclear processes [GeV]
	parameter<std::string,NO_VALIDITY_CHECK>   _pythiaParams;            ///< semi-colon separated parameters to pass to pythia, e.g. "mstj(1)=0;paru(13)=0.1"
	parameter<bool, NO_VALIDITY_CHECK>         _pythiaFullEventRecord;   ///< if the full pythia event record should be in the output
	parameter<unsigned int, VALIDITY_CHECK>    _xsecCalcMethod;	     ///< Select x-sec calc method. (0 is standard starlight method, 1 must be used for assym. collisions (e.g. p-A), but is slow)
        parameter<double, VALIDITY_CHECK>          _axionMass;               ///Axion mass//AXION HACK
        parameter<unsigned int, VALIDITY_CHECK>    _bslopeDefinition;        ///< Optional parameter to set different values of slope parameter
        parameter<double, VALIDITY_CHECK>          _bslopeValue;             ///< Value of slope parameter when _bslopeDefinition is set to 1
        parameter<double, VALIDITY_CHECK>          _bslope0;                 ///< Parameterization of slope parameter when _bslopeDefinition is set to 3
        parameter<double, VALIDITY_CHECK>          _bslope_alphaprime;       ///< Parameterization of slope parameter when _bslopeDefinition is set to 3
	parameter<unsigned int, VALIDITY_CHECK>    _printVM;                 ///< Optional parameter to set printing options for VM cross section
        parameter<unsigned int, VALIDITY_CHECK>    _impulseVM;               ///< Optional parameter to use impulse approximation (no nuclear effects) for VM cross section.
	parameter<unsigned int, VALIDITY_CHECK>    _quantumGlauber;          ///< Optional parameter.  Set = 1 to use Quantum Glauber calculation, rather than Classical Glauber
	parameter<double, VALIDITY_CHECK>          _bmin;                    ///< Optional parameter minimum impact parameter for b-range calculation
	parameter<double, VALIDITY_CHECK>          _bmax;                    /// < Optional parameter maximum impact parameter for b-range calculation

        parameter<double, VALIDITY_CHECK> _deuteronSlopePar      ;           ///< deuteron slope parameter (effective temperature) [(GeV/c)^-2]
        parameter<double, VALIDITY_CHECK> _protonMass            ;           ///< mass of the proton [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _pionChargedMass       ;           ///< mass of the pi^+/- [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _pionNeutralMass       ;           ///< mass of the pi^0 [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _kaonChargedMass       ;           ///< mass of the K^+/- [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _mel                   ;           ///< mass of the e^+/- [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _muonMass              ;           ///< mass of the mu^+/- [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _tauMass               ;           ///< mass of the tau^+/- [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f0Mass                ;           ///< mass of the f_0(980) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f0Width               ;           ///< width of the f_0(980) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f0BrPiPi              ;           ///< branching ratio f_0(980) -> pi^+ pi^- and pi^0 pi^0
        parameter<double, VALIDITY_CHECK> _etaMass               ;           ///< mass of the eta [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaWidth              ;           ///< width of the eta [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaPrimeMass          ;           ///< mass of the eta' [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaPrimeWidth         ;           ///< width of the eta' [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaCMass              ;           ///< mass of the eta_c [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaCWidth             ;           ///< width of the eta_c [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2Mass                ;           ///< mass of the f_2(1270) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2Width               ;           ///< width of the f_2(1270) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2BrPiPi              ;           ///< branching ratio f_2(1270) -> pi^+ pi^-
        parameter<double, VALIDITY_CHECK> _a2Mass                ;           ///< mass of the a_2(1320) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _a2Width               ;           ///< width of the a_2(1320) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2PrimeMass           ;           ///< mass of the f'_2(1525) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2PrimeWidth          ;           ///< width of the f'_2(1525) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2PrimeBrKK           ;           ///< branching ratio f'_2(1525) -> K^+ K^- and K^0 K^0bar
        parameter<double, VALIDITY_CHECK> _zoverz03Mass          ;           ///< mass of four-quark resonance (rho^0 pair production) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f0PartialggWidth      ;           ///< partial width f_0(980) -> g g [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaPartialggWidth     ;           ///< partial width eta -> g g [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaPrimePartialggWidth;           ///< partial width eta' -> g g [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _etaCPartialggWidth    ;           ///< partial width eta_c -> g g [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2PartialggWidth      ;           ///< partial width f_2(1270) -> g g [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _a2PartialggWidth      ;           ///< partial width a_2(1320) -> g g [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f2PrimePartialggWidth ;           ///< partial width f'_2(1525) -> g g [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _zoverz03PartialggWidth;           ///< partial width four-quark resonance -> g g (rho^0 pair production) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _f0Spin                ;           ///< spin of the f_0(980)
        parameter<double, VALIDITY_CHECK> _etaSpin               ;           ///< spin of the eta
        parameter<double, VALIDITY_CHECK> _etaPrimeSpin          ;           ///< spin of the eta'
        parameter<double, VALIDITY_CHECK> _etaCSpin              ;           ///< spin of the eta_c
        parameter<double, VALIDITY_CHECK> _f2Spin                ;           ///< spin of the f_2(1270)
        parameter<double, VALIDITY_CHECK> _a2Spin                ;           ///< spin of the a_2(1320)
        parameter<double, VALIDITY_CHECK> _f2PrimeSpin           ;           ///< spin of the f'_2(1525)
        parameter<double, VALIDITY_CHECK> _zoverz03Spin          ;           ///< spin of the four-quark resonance -> g g (rho^0 pair production)
        parameter<double, VALIDITY_CHECK> _axionSpin             ;           ///< spin of the axion
        parameter<double, VALIDITY_CHECK> _rho0Mass              ;           ///< mass of the rho^0 [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _rho0Width             ;           ///< width of the rho^0 [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _rho0BrPiPi            ;           ///< branching ratio rho^0 -> pi^+ pi^-
        parameter<double, VALIDITY_CHECK> _rho0PrimeMass         ;           ///< mass of the rho'^0 (4 pi^+/- final state) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _rho0PrimeWidth        ;           ///< width of the rho'^0 (4 pi^+/- final state) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _rho0PrimeBrPiPi       ;           ///< branching ratio rho'^0 -> pi^+ pi^-
        parameter<double, VALIDITY_CHECK> _OmegaMass             ;           ///< mass of the omega [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _OmegaWidth            ;           ///< width of the omega [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _OmegaBrPiPi           ;           ///< branching ratio omega -> pi^+ pi^-
        parameter<double, VALIDITY_CHECK> _PhiMass               ;           ///< mass of the phi [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _PhiWidth              ;           ///< width of the phi [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _PhiBrKK               ;           ///< branching ratio phi -> K^+ K^-
        parameter<double, VALIDITY_CHECK> _JpsiMass              ;           ///< mass of the J/psi [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _JpsiWidth             ;           ///< width of the J/psi [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _JpsiBree              ;           ///< branching ratio J/psi -> e^+ e^-
        parameter<double, VALIDITY_CHECK> _JpsiBrmumu            ;           ///< branching ratio J/psi -> mu^+ mu^-
        parameter<double, VALIDITY_CHECK> _JpsiBrppbar           ;           ///< branching ratio J/psi -> p pbar
        parameter<double, VALIDITY_CHECK> _Psi2SMass             ;           ///< mass of the psi(2S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Psi2SWidth            ;           ///< width of the psi(2S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Psi2SBree             ;           ///< branching ratio psi(2S) -> e^+ e^-
        parameter<double, VALIDITY_CHECK> _Psi2SBrmumu           ;           ///< branching ratio psi(2S) -> mu^+ mu^-
        parameter<double, VALIDITY_CHECK> _Upsilon1SMass         ;           ///< mass of the Upsilon(1S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Upsilon1SWidth        ;           ///< width of the Upsilon(1S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Upsilon1SBree         ;           ///< branching ratio Upsilon(1S) -> e^+ e^-
        parameter<double, VALIDITY_CHECK> _Upsilon1SBrmumu       ;           ///< branching ratio Upsilon(1S) -> mu^+ mu^-
        parameter<double, VALIDITY_CHECK> _Upsilon2SMass         ;           ///< mass of the Upsilon(2S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Upsilon2SWidth        ;           ///< width of the Upsilon(2S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Upsilon2SBree         ;           ///< branching ratio Upsilon(2S) -> e^+ e^-
        parameter<double, VALIDITY_CHECK> _Upsilon2SBrmumu       ;           ///< branching ratio Upsilon(2S) -> mu^+ mu^-
        parameter<double, VALIDITY_CHECK> _Upsilon3SMass         ;           ///< mass of the Upsilon(3S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Upsilon3SWidth        ;           ///< width of the Upsilon(3S) [GeV/c^2]
        parameter<double, VALIDITY_CHECK> _Upsilon3SBree         ;           ///< branching ratio Upsilon(3S) -> e^+ e^-
        parameter<double, VALIDITY_CHECK> _Upsilon3SBrmumu       ;           ///< branching ratio Upsilon(3S) -> mu^+ mu^-

	starlightConstants::particleTypeEnum       _particleType;
	starlightConstants::decayTypeEnum          _decayType;
	starlightConstants::interactionTypeEnum    _interactionType;

	double                         _beamLorentzGamma;         ///< Lorentz gamma factor of the beams in CMS frame, not an input parameter
	double                         _inputBranchingRatio;      ///< Branching ratio defined for each channel

};


bool inputParameters::setParameter(std::string expression)
{
    return getInputParser().parseString(expression);
}

inline
std::ostream&
operator <<(std::ostream&          out,
            const inputParameters& par)
{
	return par.print(out);
}


#endif  // INPUTPARAMETERS_H

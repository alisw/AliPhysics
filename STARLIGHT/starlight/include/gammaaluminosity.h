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
// $Rev:: 259                         $: revision of last commit
// $Author:: jseger                   $: author of last commit
// $Date:: 2016-04-19 02:58:25 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef GAMMAALUMINOSITY_H
#define GAMMAALUMINOSITY_H


#include "beambeamsystem.h"
#include "inputParameters.h"
#include "photonNucleusCrossSection.h"


class photonNucleusLuminosity : public photonNucleusCrossSection
{
 public:
  photonNucleusLuminosity(const inputParameters& input, beamBeamSystem& bbsystem);
  ~photonNucleusLuminosity();
  
 private:
  void photonNucleusDifferentialLuminosity();
  double *vmsigmapt(double W,double Egamma,double *SIGMAPT, int beam);
  void pttablegen();

  const double _ptBinWidthInterference;
  const double _interferenceStrength;
  const double _protonEnergy;
  const double _beamLorentzGamma;
  const std::string _baseFileName;
  const double _maxW;
  const double _minW;
  const unsigned int _nmbWBins;
  const double _maxRapidity;
  const unsigned int _nmbRapidityBins;
  const int _productionMode;
  const int _beamBreakupMode;
  const bool _interferenceEnabled;
  const double _maxPtInterference;
  const int _nmbPtBinsInterference;
};

#endif //GAMMAALUMINOSITY_H


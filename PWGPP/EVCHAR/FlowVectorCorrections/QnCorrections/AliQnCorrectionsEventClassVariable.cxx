/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
 
/// \file AliQnCorrectionsEventClassVariable.cxx
/// \brief Implementation of a variable used for defining an event class class

#include "AliQnCorrectionsEventClassVariable.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsEventClassVariable)
/// \endcond

/// Default constructor
AliQnCorrectionsEventClassVariable::AliQnCorrectionsEventClassVariable() :
TObject(),
fVarId(-1),
fNBins(0),
fNBinsPlusOne(0),
fBins(NULL),
fLabel("") {

}

/// Copy constructor
AliQnCorrectionsEventClassVariable::AliQnCorrectionsEventClassVariable(const AliQnCorrectionsEventClassVariable &ecv) :
TObject(ecv),
fVarId(ecv.fVarId),
fNBins(ecv.fNBins),
fNBinsPlusOne(ecv.fNBinsPlusOne),
fBins(NULL),
fLabel(ecv.fLabel) {

  fBins = new Double_t[ecv.fNBins + 1];
  for (Int_t i = 0; i < ecv.fNBins + 1; i++)
    fBins[i] = ecv.fBins[i];
}

/// Normal constructor
///
/// Allocates memory for the desired number of bins (plus one) and
/// build their lower and upper edges values
///
/// \param varId variable unique identity
/// \param varname variable name or label for a variable axis
/// \param nbins number of bins
/// \param min lower edge value for the first bin
/// \param max upper edge value for the last bin
AliQnCorrectionsEventClassVariable::AliQnCorrectionsEventClassVariable(Int_t varId, const char *varname, Int_t nbins, Double_t min, Double_t max) :
TObject(),
fVarId(varId),
fNBins(nbins),
fNBinsPlusOne(nbins+1),
fBins(NULL),
fLabel(varname) {

  fBins = new Double_t[fNBins + 1];
  Double_t low = min;
  Double_t width = (max - min) / fNBins;
  for (Int_t i = 0; i < fNBins + 1; i++) {
    fBins[i] = low;
    low += width;
  }
}

/// Normal constructor
///
/// Allocates memory for the desired number of bins (plus one) and
/// copies their passed lower and upper edges values
///
/// \param varId variable unique identity
/// \param varname variable name or label for a variable axis
/// \param nbins number of bins
/// \param bins array with bins lower edge value plus the upper of the last one
AliQnCorrectionsEventClassVariable::AliQnCorrectionsEventClassVariable(Int_t varId, const char *varname, Int_t nbins, Double_t *bins) :
TObject(),
fVarId(varId),
fNBins(nbins),
fNBinsPlusOne(nbins+1),
fBins(NULL),
fLabel(varname) {

  fBins = new Double_t[fNBins + 1];
  for (Int_t i = 0; i < fNBins + 1; i++) {
    fBins[i] = bins[i];
  }
}

/// Backwards compatible constructor
///
/// Allocates memory for the desired number of bins (plus one) and
/// build their lower and upper edges values
///
/// The passed array structure contains an array of pairs where the 1st
/// element of each pair is the lower edge of a coarse bin and the 2nd
/// element is the number of fine bins inside the coarse bin. The 2nd
/// element of the first pair is the total number of pairs
///
/// \param varId variable unique identity
/// \param varname variable name or label for a variable axis
/// \param binArray array with bin segments of different granularity
AliQnCorrectionsEventClassVariable::AliQnCorrectionsEventClassVariable(Int_t varId, const char *varname, Double_t binArray[][2]) :
TObject(),
fVarId(varId),
fNBins(0),
fNBinsPlusOne(0),
fBins(NULL),
fLabel(varname) {

  for(Int_t section = 1; section < (Int_t) binArray[0][1]; section++)
    fNBins += Int_t(binArray[section][1]);

  fNBinsPlusOne = fNBins+1;
  fBins = new Double_t [fNBins+1];

  Double_t low = binArray[0][0];

  fBins[0] = binArray[0][0];

  Int_t bin = 0;

  for(Int_t section = 1; section < binArray[0][1]; section++) {

    Double_t sectionWidth = (binArray[section][0]-binArray[section-1][0])/binArray[section][1];

    for(Int_t sectionBin = 0; sectionBin < binArray[section][1]; sectionBin++){
      fBins[bin] = low;
      low += sectionWidth;
      bin++;
    }
  }
  fBins[bin] = low;
}


/// Default destructor
///
/// Release heap memory if taken
AliQnCorrectionsEventClassVariable::~AliQnCorrectionsEventClassVariable() {

  if (fBins != NULL) {
    delete [] fBins;
    fBins = NULL;
  }
}


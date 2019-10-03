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
 
/// \file AliQnCorrectionsEventClassVariablesSet.cxx
/// \brief Implementation of the set of variables that define an event class class

#include "AliQnCorrectionsEventClassVariablesSet.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsEventClassVariablesSet)
/// \endcond

/// Gets the multidimensional configuration data
///
/// Fills the necessary information to construct a multidimensional
/// histogram
///
/// \param nbins storage for the number of bins of each variable
/// \param minvals storage for the lower values of each variable
/// \param maxvals storage for the upper values of each variable
void AliQnCorrectionsEventClassVariablesSet::GetMultidimensionalConfiguration(Int_t *nbins, Double_t *minvals, Double_t *maxvals) {
  for (Int_t var = 0; var < GetEntriesFast(); var++) {
    nbins[var] = At(var)->GetNBins();
    minvals[var] = At(var)->GetLowerEdge();
    maxvals[var] = At(var)->GetUpperEdge();
  }
}


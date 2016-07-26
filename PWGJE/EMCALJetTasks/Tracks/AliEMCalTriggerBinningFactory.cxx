/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"
#include <map>
#include <vector>
#include <TMath.h>
#include <TArrayD.h>

namespace EMCalTriggerPtAnalysis {

/**
 * Default constructor, nothing to do
 */
AliEMCalTriggerBinningFactory::AliEMCalTriggerBinningFactory() {
}

/**
 * Initialise binning component with default binning
 * \param[out] data the binning component to be initialised
 */
void AliEMCalTriggerBinningFactory::Create(AliEMCalTriggerBinningComponent* const data) {
  TArrayD binLimits;
  if(!data->GetBinning("pt")){
    CreateRAAPtBinning(binLimits);
    data->SetBinning("pt", binLimits);
  }
  if(!data->GetBinning("eta")){
    CreateDefaultEtaBinning(binLimits);
    data->SetBinning("eta", binLimits);
  }
  if(!data->GetBinning("phi")){
    CreateLinearBinning(binLimits, 100, 0, 2*TMath::Pi());
    data->SetBinning("phi", binLimits);
  }
  if(!data->GetBinning("zvertex")){
    CreateDefaultZVertexBinning(binLimits);
    data->SetBinning("zvertex", binLimits);
  }
  if(!data->GetBinning("centrality")){
    CreateLinearBinning(binLimits, 5, 0., 100.);
    data->SetBinning("centrality", binLimits);
  }
}

/**
 * Creating the default \f$ p_{t} \f$ binning.
 *
 * Definition used:
 * - from 0 to 2.5 GeV/c: 0.1 GeV/c bins
 * - from 2.5 to 7 GeV/c: 0.25 GeV/c bins
 * - from 7 to 10 GeV/c: 0.5 GeV/c bins
 * - from 10 to 15 GeV/c: 1 GeV/c bins
 * - from 15 to 20 GeV/c: 2.5 GeV/c bins
 * - from 20 to 30 GeV/c: 5 GeV/c bins
 * - from 30 to 100 GeV/c: 10 GeV/c bins
 * - from 100 to 200 GeV/c: 20 GeV/c bins
 *
 * \param[out] binning Array where to store the results.
 */
void AliEMCalTriggerBinningFactory::CreateMarkusPtBinning(TArrayD &binning) const{
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double,double>(2.5, 0.1));
  definitions.insert(std::pair<double,double>(7., 0.25));
  definitions.insert(std::pair<double,double>(10., 0.5));
  definitions.insert(std::pair<double,double>(15., 1.));
  definitions.insert(std::pair<double,double>(20., 2.5));
  definitions.insert(std::pair<double,double>(30., 5.));
  definitions.insert(std::pair<double,double>(100., 10.));
  definitions.insert(std::pair<double, double>(200., 20.));
  double currentval = 0;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Create \f$ p_{t} \f$ binning used in the \f$ R_{AA} \f$ analysis:
 *
 * Definitions are:
 * - from 0.15 to 1 GeV/c: 0.05 GeV/c bins
 * - from 1 to 2 GeV/c: 0.1 GeV/c bins
 * - from 2 to 4 GeV/c: 0.2 GeV/c bins
 * - from 4 to 7 GeV/c: 0.5 GeV/c bins
 * - from 7 to 16 GeV/c: 1 GeV/c bins
 * - from 16 to 36 GeV/c: 2 GeV/c bins
 * - from 36 to 40 GeV/c: 4 GeV/c bins
 * - from 40 to 50 GeV/c: 5 GeV/c bins
 * - from 50 to 100 GeV/c: 10 GeV/c bins
 *
 * \param[out] binning Array where to store the results
 */
void AliEMCalTriggerBinningFactory::CreateRAAPtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(36, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}


/**
 * Creating default z-Vertex binning. Bin size 5 cm.
 * \param[out] binning Array where to store the results.
 */
void AliEMCalTriggerBinningFactory::CreateDefaultZVertexBinning(TArrayD &binning) const {
  std::vector<double> mybinning;
  double currentval = -10;
  mybinning.push_back(currentval);
  while(currentval < 10.){
    currentval += 5.;
    mybinning.push_back(currentval);
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Creating default \f$ \eta \$f  binning. Bin size fixed at 0.1 units of rapidity.
 * \param[out] binning Array where to store the results.
 */
void AliEMCalTriggerBinningFactory::CreateDefaultEtaBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  double currentval = -0.8;
  mybinning.push_back(currentval);
  while(currentval < 0.8){
    currentval += 0.1;
    mybinning.push_back(currentval);
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Create any kind of linear binning from given ranges and stores it in the binning array.
 * \param[out] binning output array
 * \param[in] nbins Number of bins
 * \param[in] min lower range
 * \param[in] max upper range
 */
void AliEMCalTriggerBinningFactory::CreateLinearBinning(TArrayD& binning, int nbins, double min, double max){
  double binwidth = (max-min)/static_cast<double>(nbins);
  binning.Set(nbins+1);
  binning[0] = min;
  double currentlimit = min + binwidth;
  for(int ibin = 0; ibin < nbins; ibin++){
    binning[ibin+1] = currentlimit;
    currentlimit += binwidth;
  }
}

} /* namespace EMCalTriggerPtAnalysis */

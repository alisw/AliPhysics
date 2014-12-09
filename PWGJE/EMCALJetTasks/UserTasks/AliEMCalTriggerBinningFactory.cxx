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
/*
 * Create default binning and initalise the binning component with this. In case
 * users already defined a binning, do not overwrite this.
 *
 *   Author: Markus Fasel
 */
#include <map>
#include <vector>

#include "TArrayD.h"

#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerBinningFactory::AliEMCalTriggerBinningFactory() {
  /*
   * Default constructor, nothing to do
   */
}


//______________________________________________________________________________
void AliEMCalTriggerBinningFactory::Create(AliEMCalTriggerBinningComponent* const data) {
  /*
   * Initialise binning component with default binning
   *
   * @param data: the binning component to be initialised
   */
  AliEMCalTriggerBinningDimension *mybinning(NULL);
  TArrayD binLimits;
  if(!data->GetBinning("pt")){
    CreateDefaultPtBinning(binLimits);
    data->SetBinning("pt", binLimits);
  }
  if(!data->GetBinning("eta")){
    CreateDefaultEtaBinning(binLimits);
    data->SetBinning("eta", binLimits);
  }
  if(!data->GetBinning("zvertex")){
    CreateDefaultZVertexBinning(binLimits);
    data->SetBinning("zvertex", binLimits);
  }
}

//______________________________________________________________________________
void AliEMCalTriggerBinningFactory::CreateDefaultPtBinning(TArrayD &binning) const{
  /*
   * Creating the default pt binning.
   *
   * @param binning: Array where to store the results.
   */
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double,double>(2.5, 0.1));
  definitions.insert(std::pair<double,double>(7., 0.25));
  definitions.insert(std::pair<double,double>(15., 0.5));
  definitions.insert(std::pair<double,double>(25., 1.));
  definitions.insert(std::pair<double,double>(40., 2.5));
  definitions.insert(std::pair<double,double>(50., 5.));
  definitions.insert(std::pair<double,double>(100., 10.));
  double currentval = 0;
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

//______________________________________________________________________________
void AliEMCalTriggerBinningFactory::CreateDefaultZVertexBinning(TArrayD &binning) const {
  /*
   * Creating default z-Vertex binning.
   *
   * @param binning: Array where to store the results.
   */
  std::vector<double> mybinning;
  double currentval = -10;
  mybinning.push_back(currentval);
  while(currentval <= 10.){
    currentval += 5.;
    mybinning.push_back(currentval);
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

//______________________________________________________________________________
void AliEMCalTriggerBinningFactory::CreateDefaultEtaBinning(TArrayD& binning) const {
  /*
   * Creating default z-Vertex binning.
   *
   * @param binning: Array where to store the results.
   */
  std::vector<double> mybinning;
  double currentval = -0.8;
  mybinning.push_back(currentval);
  while(currentval <= 0.8){
    currentval += 0.1;
    mybinning.push_back(currentval);
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

} /* namespace EMCalTriggerPtAnalysis */

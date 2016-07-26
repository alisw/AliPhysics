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

AliEMCalTriggerBinningFactory::AliEMCalTriggerBinningFactory() {
  /*
   * See header file for details
   */
}

void AliEMCalTriggerBinningFactory::Create(AliEMCalTriggerBinningComponent* const data) {
  /*
   * See header file for details
   */
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

void AliEMCalTriggerBinningFactory::CreateMarkusPtBinning(TArrayD &binning) const{
  /*
   * See header file for details
   */
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

void AliEMCalTriggerBinningFactory::CreateRAAPtBinning(TArrayD& binning) const {
  /*
   * See header file for details
   */
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


void AliEMCalTriggerBinningFactory::CreateDefaultZVertexBinning(TArrayD &binning) const {
  /*
   * See header file for details
   */
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

void AliEMCalTriggerBinningFactory::CreateDefaultEtaBinning(TArrayD& binning) const {
  /*
   * See header file for details
   */
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

void AliEMCalTriggerBinningFactory::CreateLinearBinning(TArrayD& binning, int nbins, double min, double max){
  /*
   * See header file for details
   */
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

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <algorithm>
#include <cfloat>
#include <functional>
#include <vector>

#include <TArrayD.h>
#include <TMath.h>
#include <TCustomBinning.h>

ClassImp(TCustomBinning)

TCustomBinning::TCustomBinning() :
  TBinning(),
  fMinimum(-10000.),
  fMinimumSet(kFALSE),
  fSteps()
{
}

void TCustomBinning::AddStep(double maximum, double binwidth){
  if(fSteps.find(maximum) != fSteps.end()){
    printf("TCustomBinning - Maximum already defined\n");
    return;
  }
  fSteps.insert(std::pair<double, double>(maximum, binwidth));
}

void TCustomBinning::CreateBinEdges(TArrayD &binedges) const {
  if(!fMinimumSet){
    throw MinNotSetException();
  }

  // sort the ranges
  std::vector<double> keys;
  for(auto it : fSteps) keys.push_back(it.first);
  std::sort(keys.begin(), keys.end(), std::less<double>());
  std::map<double, double> sortedsteps;
  for(auto it : keys)
    sortedsteps.insert(std::pair<double, double>(it, fSteps.find(it)->second));


  std::vector<double> tmpedges;
  Double_t currentmin = fMinimum;
  tmpedges.push_back(currentmin);
  for(auto range : sortedsteps){
    while((currentmin < range.first) && (TMath::Abs(currentmin - range.first) > DBL_EPSILON)){
      currentmin += range.second;
      tmpedges.push_back(currentmin);
    }
  }

  binedges.Set(tmpedges.size());
  int bincounter = 0;
  for(auto binedge : tmpedges) binedges[bincounter++] = binedge;
}

TBinning *TCustomBinning::MakeCopy() const {
  return new TCustomBinning(*this);
}

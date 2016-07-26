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
#include <cfloat>
#include <vector>
#include <TArrayD.h>
#include <TLinearBinning.h>
#include <TMath.h>

TLinearBinning::TLinearBinning():
  TBinning(),
  fNbins(0),
  fMinimum(0),
  fMaximum(0),
  fLimitsSet(kFALSE)
{

}

TLinearBinning::TLinearBinning(Int_t nbins, Double_t min, Double_t max):
  TBinning(),
  fNbins(nbins),
  fMinimum(min),
  fMaximum(max),
  fLimitsSet(kTRUE)
{

}

void TLinearBinning::CreateBinEdges(TArrayD &binedges) const {
  if(!fLimitsSet){
    throw LimitsNotSetException();
  }
  double binwidth = (fMaximum - fMinimum) / fNbins;

  std::vector<double> tmpedges;
  double currentmin = fMinimum;
  tmpedges.push_back(currentmin);
  while((currentmin < fMaximum) && (TMath::Abs(currentmin - fMaximum) > DBL_EPSILON)){
    currentmin += binwidth;
    tmpedges.push_back(currentmin);
  }

  binedges.Set(tmpedges.size());
  int bincounter = 0;
  for(auto binedge : tmpedges) binedges[bincounter++] = binedge;
}

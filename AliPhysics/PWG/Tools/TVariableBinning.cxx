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
#include <TVariableBinning.h>

ClassImp(TVariableBinning)

TVariableBinning::TVariableBinning():
  TBinning(),
  fBinEdges()
{

}

TVariableBinning::TVariableBinning(Int_t nbins, const Double_t *binedges):
  TBinning(),
  fBinEdges()
{
  Set(nbins, binedges);
}

TVariableBinning::TVariableBinning(const TArrayD &binedges):
  TBinning(),
  fBinEdges(binedges)
{

}

void TVariableBinning::CreateBinEdges(TArrayD &binedges) const {
  if(!fBinEdges.GetSize()){
    throw LimitsNotSetException();
  }
  binedges = fBinEdges;
}

TBinning *TVariableBinning::MakeCopy() const {
  return new TVariableBinning(*this);
}

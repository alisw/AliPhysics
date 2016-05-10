#ifndef ALIEMCALLIST_H
#define ALIEMCALLIST_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * \class AliEmcalList
 * \brief Enhanced TList-derived class that implements correct merging for pt_hard binned production
 *
 * Must be activated using SetUseScaling(kTRUE). Otherwise the behavior is like a TList
 * Scaling is recursively applied also to all nested lists deriving from TCollection
 * fHistXsection and fHistTrials must be added directly to the list (not to a nested list)
 *
 * \author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * \date May 05, 2016
 */
//

class TList;

#include "TList.h"

class AliEmcalList : public TList {

public:
  AliEmcalList();
  ~AliEmcalList() {}
  Long64_t                    Merge(TCollection *hlist);
  void                        SetUseScaling(Bool_t val) {fUseScaling = val;}

private:
  // ####### Helper functions
  void                        ScaleAllHistograms(TCollection *hlist, Double_t scalingFactor);
  Double_t                    GetScalingFactor(TH1* xsection, TH1* ntrials);
  Bool_t                      IsLastMergeLevel(TCollection* collection);
  Int_t                       GetFilledBinNumber(TH1* hist);
  
  Bool_t                      fUseScaling;                    ///< if true, scaling will be done. if false AliEmcalList simplifies to TList

  /// \cond CLASSIMP
  ClassDef(AliEmcalList, 1);
  /// \endcond
};

#endif

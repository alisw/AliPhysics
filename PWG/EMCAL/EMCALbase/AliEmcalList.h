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
 * \ingroup EMCALCOREFW
 */
//

class TH1;

#include "TList.h"
#include "TString.h"

class AliEmcalList : public TList {

public:
  AliEmcalList();
  ~AliEmcalList() {}
  Long64_t                    Merge(TCollection *hlist);
  void                        SetUseScaling(Bool_t val) {fUseScaling = val;}
  Bool_t                      IsUseScaling() const { return fUseScaling; }
  void                        SetNameXsec(const char *name) { fNameXsec = name; }
  void                        SetNameTrials(const char *name) { fNameNTrials = name; }

private:
  // ####### Helper functions
  void                        ScaleAllHistograms(TCollection *hlist, Double_t scalingFactor);
  Double_t                    GetScalingFactor(const TH1* xsection, const TH1* ntrials) const;
  Bool_t                      IsLastMergeLevel(const TCollection* collection) const;
  Int_t                       GetFilledBinNumber(const TH1* hist) const;
  Bool_t                      IsScalingSupported(const TObject *scaleobject) const;
  
  Bool_t                      fUseScaling;                    ///< if true, scaling will be done. if false AliEmcalList simplifies to TList
  TString                     fNameXsec;                      ///< Name of the cross section histogram
  TString                     fNameNTrials;                   ///< Name of the histogram with the number of trials

  /// \cond CLASSIMP
  ClassDef(AliEmcalList, 2);
  /// \endcond
};

#endif

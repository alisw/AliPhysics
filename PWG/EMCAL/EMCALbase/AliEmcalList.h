/************************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIEMCALLIST_H
#define ALIEMCALLIST_H

class TH1;

#include "TList.h"
#include "TString.h"

/**
 * @class AliEmcalList
 * @brief Enhanced TList-derived class that implements correct merging for pt_hard binned production
 * @author Ruediger Haake <ruediger.haake@cern.ch>, CERN
 * @date May 05, 2016
 * @ingroup EMCALCOREFW
 *
 * Must be activated using SetUseScaling(kTRUE). Otherwise the behavior is like a TList
 * Scaling is recursively applied also to all nested lists deriving from TCollection
 * fHistXsection and fHistTrials must be added directly to the list (not to a nested list)
 */
class AliEmcalList : public TList {

public:

  /**
   * @brief Constructor
   */
  AliEmcalList();

  /**
   * @brief Destructor
   */
  ~AliEmcalList() {}
  /**
   * @brief Merge function including the reweighting of \f$p_{t}\f$-hard bins
   * @param hlist Collection of object to be merged
   * @return Number of entries to be merged
   * 
   * Overriding the standard Merge function in order to allow
   * for reweighting of the histograms based on the weight for 
   * the given \f$p_{t}\f$-hard bin: In case the last merging level
   * (merging of files from different \f$p_{t}\f$-hard bins) is reached,
   * all histograms in the list except from the usual scaling histograms
   * which are created by the AliAnalysisTaskEmcal or AliAnalysisTaskEmcalLight
   * are scaled by the cross section divided by the number of trials. The
   * merging level is determined from whether different bins in the scaling
   * histograms are filled for the various histograms to be merged.
   */
  Long64_t                    Merge(TCollection *hlist);

  /**
   * @brief Set the container to rescaing mode
   * @param val True for scaling mode, false for normal merging
   */
  void                        SetUseScaling(Bool_t val) {fUseScaling = val;}

  /**
   * @brief Check if the merging is in scaling mode
   * @return If true the merging is in scaling mode
   */
  Bool_t                      IsUseScaling() const { return fUseScaling; }

  /**
   * @brief Set the name of the cross section histogram used for the weight calculation
   * @param name Name of the cross section histogram
   */
  void                        SetNameXsec(const char *name) { fNameXsec = name; }
  
  /**
   * @brief Set the name of the histogram with the number of trials histogram used for the weight calculation
   * @param name Name of the histogram with the number of trials
   */
  void                        SetNameTrials(const char *name) { fNameNTrials = name; }

private:
  // ####### Helper functions

  /**
   * @brief Scaling all histograms in the histogram list by the weight of the \f$p_{t}\f$-hard bin
   * @param hlist List of histograms to be scaled
   * @param scalingFactor Scaling factor to be applied on the histograms
   * 
   * Scaling of all histograms in hlist recursively. Histograms which are created
   * by the AliAnalysisTaskEmcal or AliAnalysisTaskEmcalLight in order to obtain the
   * scale factors are not scaled.
   */
  void                        ScaleAllHistograms(TCollection *hlist, Double_t scalingFactor);
  
  /**
   * @brief Helper function scaling factor
   * @param xsection Histogram with the cross section for the different \f$p_{t}\f$-hard bins
   * @param ntrials  Histogram with the number of triaks for the different \f$p_{t}\f$-hard bins
   * @return Weight for the \f$p_{t}\f$-hard bin
   */
  Double_t                    GetScalingFactor(const TH1* xsection, const TH1* ntrials) const;

  /**
   * @brief Helper function to determine whether we are in last merge step
   * @param collection Collection of AliEmcalList objects
   */
  Bool_t                      IsLastMergeLevel(const TCollection* collection) const;

  /**
   * @brief Helper function that returns the bin in a TH1 that is filled
   * @param hist  Histogram
   * @return bin number that is filled. If no or more than one bin is filled, 0.
   */ 
  Int_t                       GetFilledBinNumber(const TH1* hist) const;

  /**
   * @brief Helper function checking whether type is supported for scaling
   * @param scaleobject Object for which to deterime scaling support
   * @return true if object can be scaled, false otherwise
   * 
   * Supported types are: 
   * - all TH1-derived / THnBase-derived histograms
   *  - objects inheriting from RooUnfoldResponse if AliPhysics is built with
   *    RooUnfold support
   */
  Bool_t                      IsScalingSupported(const TObject *scaleobject) const;
  
  Bool_t                      fUseScaling;                    ///< if true, scaling will be done. if false AliEmcalList simplifies to TList
  TString                     fNameXsec;                      ///< Name of the cross section histogram
  TString                     fNameNTrials;                   ///< Name of the histogram with the number of trials

  ClassDef(AliEmcalList, 2);
};

#endif

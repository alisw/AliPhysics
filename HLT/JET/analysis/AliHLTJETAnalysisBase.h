//-*- Mode: C++ -*-
#ifndef ALIHLTJETANALYSISBASE_H
#define ALIHLTJETANALYSISBASE_H
 
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETAnalysisBase.h
    @author Jochen Thaeder
    @date   
    @brief  Base functionality for HLT JET analysis package
*/

#include "TObject.h"

#include "AliHLTLogging.h"

/**
 * @class AliHLTJETAnalysisBase
 * This class contains basic constants for the Jet Analysis.
 *
 * @ingroup alihlt_jet
 * @ingroup alihlt_jet_analysis
 */

class AliHLTJETAnalysisBase : public TObject, public AliHLTLogging {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
 
  /** Standard constructor */
  AliHLTJETAnalysisBase();

  /** Destructor */
  ~AliHLTJETAnalysisBase();  
  
  /*
   * ---------------------------------------------------------------------------------
   *                          Histogram enum's - public
   * ---------------------------------------------------------------------------------
   */

  /** Plot type of Delta histograms */
  enum JetDeltaType_t {
    kDeltaAll,                /**< Delta of all jets */
    kDeltaLead,               /**< Delta of leading jets */
    kDeltaMatchedAll,         /**< Delta of all matched jets */
    kDeltaMatchedLead,        /**< Delta of all matched leading jets */
    kDeltaMax                 /**< Number of enum entries */
  };
 
  /** Plot type of Spectra histograms */
  enum JetSpectraType_t {
    kSpectraPythiaAll,        /**< Spectra of all pythia jets */
    kSpectraPythiaMatched,    /**< Spectra of matched pythia jets */
    kSpectraPythiaUnmatched,  /**< Spectra of unmatched pythia jets */
    kSpectraRecoAll,          /**< Spectra of all reco jets */
    kSpectraRecoMatched,      /**< Spectra of matched reco jets */
    kSpectraRecoUnmatched,    /**< Spectra of unmatched reco jets */
    kSpectraRecoLeadAll,      /**< Spectra of all leading reco jets */
    kSpectraRecoLeadMatched,  /**< Spectra of matched leading reco jets */
    kSpectraRecoLeadUnmatched,/**< Spectra of unmatched leading reco jets */
    kSpectraMax               /**< Number of enum entries */
  };
  
  /** Plot type of histograms */
  enum JetPlotType_t {
    kPlotAll,                 /**< All jets */
    kPlotLead,                /**< Leading jets */
    kPlotMax                  /**< Number of enum entries */
  };

  /** Type of histgrams */
  enum JetHistogramType_t {
    kHistDeltaEt,             /**< Delta E_t histogram */
    kHistDeltaEta,            /**< Delta Eta histogram */
    kHistDeltaPhi,            /**< Delta Phi histogram */
    kHistDeltaEtaDeltaPhi,    /**< Delta Eta Delta Phi histogram */
    kHistSpectraEt,           /**< Spectra E_t histogram */
    kHistSpectraEta,          /**< Spectra Eta histogram */
    kHistSpectraPhi,          /**< Spectra Phi histogram */
    kHistCorrelationsJetEt,   /**< Correlation Et histogram */
    kHistResolutionsJetEt,    /**< Resolution Et Jet histogram */
    kHistResolutionsDiJetEt,  /**< Resolution Et DiJet histogram */
    kHistMax                  /**< Number of enum entries */
  };

  /** Array of types of the Delta histograms */
  static const Char_t *fgkDeltaType[];        //! transient

  /** Array of types of the Spectra histograms */
  static const Char_t *fgkSpectraType[];      //! transient

  /** Array of types of histograms */
  static const Char_t *fgkPlotType[];         //! transient
 
  ///////////////////////////////////////////////////////////////////////////////////
  
 private:

  /** copy constructor prohibited */
  AliHLTJETAnalysisBase(const AliHLTJETAnalysisBase&);

  /** assignment operator prohibited */
  AliHLTJETAnalysisBase& operator=(const AliHLTJETAnalysisBase&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
    
  ClassDef(AliHLTJETAnalysisBase, 0)
};
#endif

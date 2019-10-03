//-*- Mode: C++ -*-
#ifndef ALIHLTJETANALYSISBASE_H
#define ALIHLTJETANALYSISBASE_H
 
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETAnalysisBase.h
    @author Jochen Thaeder <jochen@thaeder.de>
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
    kSpectraCmpAll,           /**< Spectra of all cmp jets */
    kSpectraCmpMatched,       /**< Spectra of matched cmp jets */
    kSpectraCmpUnmatched,     /**< Spectra of unmatched cmp jets */
    kSpectraCmpLeadAll,       /**< Spectra of all leading cmp jets */
    kSpectraCmpLeadMatched,   /**< Spectra of matched leading cmp jets */
    kSpectraCmpLeadUnmatched, /**< Spectra of unmatched leading cmp jets */
    kSpectraRecAll,           /**< Spectra of all rec jets */
    kSpectraRecMatched,       /**< Spectra of matched rec jets */
    kSpectraRecUnmatched,     /**< Spectra of unmatched rec jets */
    kSpectraRecLeadAll,       /**< Spectra of all leading rec jets */
    kSpectraRecLeadMatched,   /**< Spectra of matched leading rec jets */
    kSpectraRecLeadUnmatched, /**< Spectra of unmatched leading rec jets */
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

  /** Array of types of the MC Spectra histograms */
  static const Char_t *fgkSpectraTypeMC[];    //! transient

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

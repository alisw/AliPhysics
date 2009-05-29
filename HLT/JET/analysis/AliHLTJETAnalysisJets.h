//-*- Mode: C++ -*-

// $Id: AliHLTJETAnalysisJets.h  $

#ifndef ALIHLTJETANALYSISJETS_H
#define ALIHLTJETANALYSISJETS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETAnalysisJets.h
    @author Jochen Thaeder
    @date   
    @brief  Container holding analysis objects
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "TH1.h"
//#include "TClonesArray.h"

class TH1; 
class TClonesArray; 

class AliAODJet;
class AliMCEvent;

#include "AliHLTLogging.h"
#include "AliHLTMCEvent.h"

#include "AliHLTJETBase.h"
#include "AliHLTJETJets.h"

#include "AliHLTJETAnalysisBase.h"

/**
 * @class AliHLTJETAnalysisJets
 * This class is a container which holds TClonesArrys of 
 * histograms needed for the Jet Analysis.
 *
 * It need a ptr to MC information (pythia jets) and 
 * reconstructed jets in the form of AliHLTJETJets.
 *
 * @ingroup alihlt_jet
 * @ingroup alihlt_jet_analysis
 */

class AliHLTJETAnalysisJets : public TObject, public AliHLTLogging {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Constructor */
  AliHLTJETAnalysisJets();
  
  /** Destructor */
  ~AliHLTJETAnalysisJets();

  /*
   * ---------------------------------------------------------------------------------
   *                         Initialize / Setup / Reset - public
   * ---------------------------------------------------------------------------------
   */

  /** Setup analysis class 
   *  @return 0 on success, <0 on failure
   */
  Int_t Initialize();

  /** Reset input for current event */
  void ResetEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                                Setter - public
   * ---------------------------------------------------------------------------------
   */

  /** Set pythia jets in HLT environment
   *  @param mcEvent    Ptr to AliHLTMCEvent
   */
  void SetHLTMC( AliHLTMCEvent* mcEvent );

  /** Set pythia jets in Off-line environment
   *  @param mcEvent    Ptr to AliMCEvent
   */
  void SetMC( AliMCEvent* mcEvent );

  /** Set reconstructed jets
   *  @param jets       Ptr to AliHLTJETJets
   */
  void SetJets( AliHLTJETJets* jets );

  /*
   * ---------------------------------------------------------------------------------
   *                                 Getter - public
   * ---------------------------------------------------------------------------------
   */

  /** Returns histogram dependent of histogram type and plot type 
   *  @param histIdx    histogram type (Ptr to TClonesArray),
   *                    @see AliHLTJETAnalysisBase
   *  @param plotIdx    plot type (Entry Idx in TClonesArray),
   *                    @see AliHLTJETAnalysisBase
   *  @return           Ptr to TH1 on success, NULL on failure
   */
  TH1* GetHistogram ( Int_t histIdx, Int_t plotIdx );

  /*
   * ---------------------------------------------------------------------------------
   *                          Analysis - public
   * ---------------------------------------------------------------------------------
   */

  /** Anlayze Data Set 
   *  -Fill unmatched jets into histograms
   *  -Match jets
   *  -Fill matched jets into histograms
   *  @return 0 on success, <0 on failure
   */
  Int_t Analyze();

  ///////////////////////////////////////////////////////////////////////////////////
  
 private:
 
  /** copy constructor prohibited */
  AliHLTJETAnalysisJets(const AliHLTJETAnalysisJets&);
  
  /** assignment operator prohibited */
  AliHLTJETAnalysisJets& operator=(const AliHLTJETAnalysisJets&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Setup / Reset - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Setup Delta histograms */
  void SetupDeltaHistograms();

  /** Setup Spectra histograms */
  void SetupSpectraHistograms();

  /** Setup matched histograms */
  void SetupMatchedHistograms();

  /*
   * ---------------------------------------------------------------------------------
   *                             Analysis - private
   * ---------------------------------------------------------------------------------
   */

  /** Match Pythia and Reconstructed Jets  
   *  @return 0 on success, <0 on failure
   */
  Int_t MatchJets();
  
  /*
   * ---------------------------------------------------------------------------------
   *                                Fill - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Fill basic Spectra histogram */
  void FillBasicSpectraHistograms();

  /** Fill unmatched Delta histogram */
  void FillUnmatchedDeltaHistograms();

  /** Fill matched Delta histogram */
  void FillMatchedDeltaHistograms();

  /** Fill matched Spectra histogram */
  void FillMatchedSpectraHistograms();

  /** Fill matched histogram */
  void FillMatchedHistograms();
  
  /*
   * ---------------------------------------------------------------------------------
   *                               Helper - private
   * ---------------------------------------------------------------------------------
   */

  /** Setup histogram with common parameters 
   *  @param hist    Ptr to histogram 
   */
  void SetupHist( TH1* hist );

  /** Fill 1D histogram, in a TClonesArray
   *  @param array   Ptr to TClonesArray of histograms
   *  @param idx     Index in the TClonesArray
   *  @param valueX  x value
   */
  void FillHist( TClonesArray* array, Int_t idx, Float_t valueX  );
  
  /** Fill 2D histogram, in a TClonesArray
   *  @param array   Ptr to TClonesArray of histograms
   *  @param idx     Index in the TClonesArray
   *  @param valueX  x value
   *  @param valueY  y value
   */
  void FillHist( TClonesArray* array, Int_t idx, Float_t valueX,Float_t valueY );
  
  /** Get Distance^2 in eta phi space of 2 jets
   *  @param jet1    Ptr to jet 1
   *  @param jet2    Ptr to jet 2
   *  @return        Distance^2
   */
  Float_t GetDistance2( AliAODJet *jet1, AliAODJet *jet2); 

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++ Data Members 
  // ++
  // ++-> replaced every event
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // ---------------------------------------------------
  // -- Members filled per event
  // ---------------------------------------------------

  /** jets */
  AliHLTJETJets *fJets;                          //! transient      

  /** MC jets */
  AliHLTJETJets *fJetsMC;                        //! transient      

  /** Array of indices of the matched reconstructed jets */
  TArrayI       *fMatchedJets;                   //! transient      

  /** Array of indices of the matched MC jets */
  TArrayI       *fMatchedJetsMC;                 //! transient      

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++ Analysis Output 
  // ++
  // ++ -> Created once
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // ---------------------------------------------------
  // -- Difference in reconstruction
  // ---------------------------------------------------
  
  /** Matched Jets - delta Et */
  TClonesArray  *fDeltaEt;                       // see above

  /** Delta eta ( pythia - jetfinder ) */
  TClonesArray  *fDeltaEta;                      // see above
  
  /** Delta phi ( pythia - jetfinder ) */
  TClonesArray  *fDeltaPhi;                      // see above
  
  /** Delta eta, delta phi( pythia - jetfinder ) */
  TClonesArray  *fDeltaEtaDeltaPhi;              // see above
  
  // ---------------------------------------------------
  // -- Jet spectra
  // ---------------------------------------------------
  
  /** Jet spectra in Et */
  TClonesArray  *fSpectraEt;                     // see above

  /** Jet spectra in eta */
  TClonesArray  *fSpectraEta;                    // see above
  
  /** Jet spectra in phi */
  TClonesArray  *fSpectraPhi;                    // see above
           
  // ---------------------------------------------------
  // -- Correlations
  // ---------------------------------------------------

  /** Correleation pythia vs jet finder */
  TClonesArray  *fCorrelationsJetEt;             // see above

  // ---------------------------------------------------
  // -- Resolutions
  // ---------------------------------------------------
  
  /** Resolutions for Et for jetfinder - pythia fixed */
  TClonesArray  *fResolutionsJetEt;              // see above

  /** Resolutions for Et for jetfinder - nearside fixed */
  TClonesArray  *fResolutionsDiJetEt;            // see above

  ClassDef(AliHLTJETAnalysisJets, 2);
};
#endif

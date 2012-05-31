//-*- Mode: C++ -*-

// $Id: AliHLTMultiplicityCorrelations.h  $
#ifndef ALIHLTMULTIPLICITYCORRELATIONS_H
#define ALIHLTMULTIPLICITYCORRELATIONS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTMultiplicityCorrelations.h
    @author Jochen Thaeder
    @date   
    @brief  Correlation plots for multiplicity studies
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTLogging.h"

#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"

#include "TList.h"

class TH1;

/**
 * @class AliHLTMultiplicityCorrelations
 *
 * @ingroup alihlt_physics
 */

class AliHLTMultiplicityCorrelations : public TNamed, public AliHLTLogging {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Default Constructor */
  AliHLTMultiplicityCorrelations();
  
  /** Constructor */
  AliHLTMultiplicityCorrelations(Char_t* name, Char_t* title);
  
  /** Destructor */
  ~AliHLTMultiplicityCorrelations();

  /*
   * ---------------------------------------------------------------------------------
   *                         Initialize / Setup / Reset - public
   * ---------------------------------------------------------------------------------
   */

  /** Initialize class and members */
  Int_t Initialize();

  /** Initialize class and members */
  Int_t Initialize( const Char_t* listName );

  /*
   * ---------------------------------------------------------------------------------
   *                                Setter - public
   * ---------------------------------------------------------------------------------
   */

  /** Set ESD track cuts */
  void SetESDTrackCuts(AliESDtrackCuts *cuts) { fESDTrackCuts = cuts; }

  /** Set SPD clusters from inner and outer layer */
  void SetSPDClusters(Int_t inner, Int_t outer) { 
    fSpdNClustersInner = inner; fSpdNClustersOuter = outer; 
  }

  /** Set Binning of VZERO */
  void SetBinningVzero(Int_t i=1, Float_t f1=0., Float_t f2=1.) {
    fVzeroBinning = i; fVzeroBinningMin = f1; fVzeroBinningMax = f2;
  }
  
  /** Set Binning of TPC */
  void SetBinningTpc(Int_t i=1, Float_t f1=0., Float_t f2=1.) {
    fTpcBinning = i; fTpcBinningMin = f1; fTpcBinningMax = f2;
  }

  /** Set Binning of ZDC */
  void SetBinningZdc(Int_t i=1, Float_t f1=0., Float_t f2=1.) {
    fZdcBinning = i; fZdcBinningMin = f1; fZdcBinningMax = f2;
  }

  /** Set Binning of ZEM */
  void SetBinningZem(Int_t i=1, Float_t f1=0., Float_t f2=1.) {
    fZemBinning = i; fZemBinningMin = f1; fZemBinningMax = f2;
  }

  /** Set Binning of ZNP */
  void SetBinningZnp(Int_t /*i=1*/, Float_t /*f1=0.*/, Float_t /*f2=1.*/) {;}

  /** Set Binning of CALO */
  void SetBinningCalo(Int_t /*i=1*/, Float_t /*f1=0.*/, Float_t /*f2=1.*/) {;}

  /** Set Binning of SPD */
  void SetBinningSpd(Int_t i=1, Float_t f1=0., Float_t f2=1.) {
    fSpdBinning = i; fSpdBinningMin = f1; fSpdBinningMax = f2;
  }

  /** Set process PHOS */
  void SetProcessPhos(Bool_t /*b=kTRUE*/)  {;}
  /** Set process EMCAL */
  void SetProcessEmcal(Bool_t /*b=kTRUE*/) {;}

  /** Enable / Disable detectors */
  void SetProcessCALO(Bool_t /*b = kTRUE*/) {;}
  void SetProcessSPD(Bool_t b = kTRUE)  { fProcessSPD = b; }
  void SetProcessTPC(Bool_t b = kTRUE)  { fProcessTPC = b; }
  void SetProcessZDC(Bool_t b = kTRUE)  { fProcessZDC = b; }
  void SetProcessVZERO(Bool_t b = kTRUE){ fProcessVZERO = b; }
  void SetProcessCentrality(Bool_t b = kTRUE){ fProcessCentrality = b; }

  /** Centrality Estimator */
  void SetCentralityEstimator(TH1F* hist) { fCentHistV0Mpercentile = hist; }

  /*
   * ---------------------------------------------------------------------------------
   *                                 Getter - public
   * ---------------------------------------------------------------------------------
   */

  /** Get List of histograms */
  TList* GetHistList() const { return fHistList; }

  /*
   * ---------------------------------------------------------------------------------
   *                             Process - public
   * ---------------------------------------------------------------------------------
   */
  
  /** Process current event */
  Int_t ProcessEvent( AliESDEvent *esd, AliESDVZERO* esdVZERO, Int_t nSpdClusters);

  Int_t ProcessEvent( AliESDEvent *esd ) {
    return ProcessEvent(esd, NULL, 0);
  }


  ///////////////////////////////////////////////////////////////////////////////////
  
 private:
 
  /** copy constructor prohibited */
  AliHLTMultiplicityCorrelations(const AliHLTMultiplicityCorrelations&);
  
  /** assignment operator prohibited */
  AliHLTMultiplicityCorrelations& operator=(const AliHLTMultiplicityCorrelations&);

  /*
   * ---------------------------------------------------------------------------------
   *                         Initialize / Setup / Reset - private
   * ---------------------------------------------------------------------------------
   */

  /** Add esd object
   * param esd Ptr to AliESDEvent
   * return kTRUE if AliESDEvent and Vertex present
   */
  Bool_t AddESDEvent( AliESDEvent* esd );

  /** Setup histograms */
  Int_t SetupHistograms();

  /** Setup VZERO histograms */
  Int_t SetupVZERO();

  /** Setup ZDC histograms */
  Int_t SetupZDC();

  /** Setup TPC histograms */
  Int_t SetupTPC();

  /** Setup correlation histograms */
  Int_t SetupCorrelations();
  
  /** Setup SPD histograms */
  Int_t SetupSPD();

  /** Setup Centrality histograms */
  Int_t SetupCentrality();

  /*
   * ---------------------------------------------------------------------------------
   *                             Process - private
   * ---------------------------------------------------------------------------------
   */

  /** Process current event - TPC */
  Int_t ProcessTPC();

  /** Process current event - SPD */
  Int_t ProcessSPD();
  
  /** Process current event - VZERO */
  Int_t ProcessVZERO();

  /** Process current event - ZDC and correlations */
  Int_t ProcessZDC();

  /** Process current event - centrality */
  Int_t ProcessCentrality();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** List of histograms */ 
  TList           *fHistList;             // see above

  /** Ptr to AliESDEvent */
  AliESDEvent     *fESDEvent;             //! transient

  /** Ptr to ZDC object in AliESDEvent*/ 
  AliESDZDC       *fESDZDC;               //! transient

  /** Ptr to VZERO object in AliESDEvent*/ 
  AliESDVZERO     *fESDVZERO;             //! transient

  /** Ptr to AliESD track cuts */
  AliESDtrackCuts *fESDTrackCuts;         //! transient

  /** Centrality Estimator */
  TH1F* fCentHistV0Mpercentile;           //! transient

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  /** Process TPC information */
  Bool_t           fProcessTPC;           // see above

  /** Process SPD information */
  Bool_t           fProcessSPD;           // see above

  /** Process VZERO information */
  Bool_t           fProcessVZERO;         // see above

  /** Process ZDC information */
  Bool_t           fProcessZDC;           // see above

  /** Process Centrality information */
  Bool_t           fProcessCentrality;    // see above

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  /** N ESD tracks */
  Int_t fEsdTracks;                       // see above

  /** N ESD tracks accepted */
  Int_t fEsdTracksA;                      // see above

  /** N TPC tracks */
  Int_t fTpcTracks;                       // see above

  /** N TPC tracks accepted */
  Int_t fTpcTracksA;                      // see above

  /** N Reference TPC tracks */
  Int_t fTpcTracksRef;                    // see above

  /** VZERO mult */
  Float_t fVzeroMult;                     // see above

  /** VZERO trigger mult */
  Float_t fVzeroTriggerMult;              // see above

  /** Spd N clusters */
  Int_t   fSpdNClusters;                  // see above

  /** Spd N clusters inner layer*/
  Int_t   fSpdNClustersInner;             // see above

  /** Spd N clusters outer layer */
  Int_t   fSpdNClustersOuter;             // see above

  // -- -- -- 

  /** Binnning VZERO */
  Int_t   fVzeroBinning;                  // see above
  Float_t fVzeroBinningMin;               // see above
  Float_t fVzeroBinningMax;               // see above

  /** Binnning TPC */
  Int_t   fTpcBinning;                    // see above
  Float_t fTpcBinningMin;                 // see above
  Float_t fTpcBinningMax;                 // see above

  /** Binnning ZDC */
  Int_t   fZdcBinning;                    // see above
  Float_t fZdcBinningMin;                 // see above
  Float_t fZdcBinningMax;                 // see above

  /** Binnning ZEM */
  Int_t   fZemBinning;                    // see above
  Float_t fZemBinningMin;                 // see above
  Float_t fZemBinningMax;                 // see above

  /** Binnning SPD */
  Int_t   fSpdBinning;                    // see above
  Float_t fSpdBinningMin;                 // see above
  Float_t fSpdBinningMax;                 // see above
  
  ClassDef(AliHLTMultiplicityCorrelations, 4);
};
#endif

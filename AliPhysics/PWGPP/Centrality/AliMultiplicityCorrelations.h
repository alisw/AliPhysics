//-*- Mode: C++ -*-

// $Id: AliMultiplicityCorrelations.h  $
#ifndef ALIMULTIPLICITYCORRELATIONS_H
#define ALIMULTIPLICITYCORRELATIONS_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Correlation plots for multiplicity studies
// Authors: Jochen Thaeder <jochen@thaeder.de>

#include "AliLog.h"

#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "TList.h"

class TH1;

class AliMultiplicityCorrelations : public TNamed {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Default Constructor */
  AliMultiplicityCorrelations();
  
  /** Constructor */
  AliMultiplicityCorrelations(Char_t* name, Char_t* title);
  
  /** Destructor */
  ~AliMultiplicityCorrelations();

  /*
   * ---------------------------------------------------------------------------------
   *                         Initialize / Setup / Reset - public
   * ---------------------------------------------------------------------------------
   */

  /** Initialize class and members */
  Int_t Initialize() { return Initialize(""); }

  /** Initialize class and members */
  Int_t Initialize( const Char_t* listName );

  /*
   * ---------------------------------------------------------------------------------
   *                                Setter - public
   * ---------------------------------------------------------------------------------
   */

  void SetIsMC() { fIsMC = kTRUE; }

  /** Clean event sample */
  void SetCleanSample(Float_t min, Float_t max) { fCleanSample = kTRUE; fCleanMinKeep = min; fCleanMaxKeep = max; }

  /** Set ESD track cuts */
  void SetESDTrackCuts(AliESDtrackCuts *cuts) { fESDTrackCuts  = cuts; }
  void SetESDTrackCuts2(AliESDtrackCuts *cuts) { fESDTrackCuts2 = cuts; }

  /** Set SPD clusters from inner and outer layer */
  void SetSPDClusters(Float_t inner, Float_t outer) { fSpdNClustersInner = inner; fSpdNClustersOuter = outer; }

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

  /** Set Binning of SPD */
  void SetBinningSpd(Int_t i=1, Float_t f1=0., Float_t f2=1.) {
    fSpdBinning = i; fSpdBinningMin = f1; fSpdBinningMax = f2;
  }

  /** Enable / Disable detectors */
  void SetProcessSPD(Bool_t b = kTRUE)  { fProcessSPD = b; }
  void SetProcessTPC(Bool_t b = kTRUE)  { fProcessTPC = b; }
  void SetProcessZDC(Bool_t b = kTRUE)  { fProcessZDC = b; }
  void SetProcessVZERO(Bool_t b = kTRUE){ fProcessVZERO = b; }

  /*
   * ---------------------------------------------------------------------------------
   *                                 Getter - public
   * ---------------------------------------------------------------------------------
   */

  /** Get List of histograms */
  TList*  GetHistList()    const { return fHistList; }
  Int_t   GetNTracks()     const { return fEsdTracksA; }
  Int_t   GetNTracksTPC()  const { return fTpcTracksA; }
  Float_t GetVZEROA()      const { return fVzeroMultA; }
  Float_t GetVZEROC()      const { return fVzeroMultC; }
  Float_t GetVZEROCorr()   const { return fVzeroMult; }

  /*
   * ---------------------------------------------------------------------------------
   *                             Process - public
   * ---------------------------------------------------------------------------------
   */

  /** Process current event */
  Int_t ProcessEvent( AliESDEvent *esd );

  ///////////////////////////////////////////////////////////////////////////////////
  
  /** Corrected VZERO amplitude*/
  Float_t GetCorrVZERO(Float_t &v0CorrResc);

  /** Corrected SPD amplitude*/
  Float_t GetCorrSPD2(Float_t spd2raw,Float_t zv) const;

 private:
 
  /** copy constructor prohibited */
  AliMultiplicityCorrelations(const AliMultiplicityCorrelations&);
  
  /** assignment operator prohibited */
  AliMultiplicityCorrelations& operator=(const AliMultiplicityCorrelations&);

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
  
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  TList           *fHistList;             //  List of histograms

  Bool_t           fIsMC;                 //  If it is MC 

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  AliESDEvent     *fESDEvent;             //! Ptr to AliESDEvent
  AliESDZDC       *fESDZDC;               //! Ptr to ZDC object in AliESDEvent 
  AliESDVZERO     *fESDVZERO;             //! Ptr to VZERO object in AliESDEvent
  AliMultiplicity *fESDMultiplicity;      //! Ptr to AliMultiplicity in AliESDEvent
  AliESDtrackCuts *fESDTrackCuts;         //! Ptr to AliESDtrackCuts
  AliESDtrackCuts *fESDTrackCuts2;        //! Ptr to AliESDtrackCuts 2
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  Bool_t           fCleanSample;          //  Enables 'grass' cleaning
  Float_t          fCleanMinKeep;         //  Min of kept region
  Float_t          fCleanMaxKeep;         //  Max of kept region

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  Int_t            fRunNo;                //  RunNo for corrections
  Int_t            fCurrentRunNo;         //  Current RunNo

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  Bool_t           fProcessTPC;           //  Process TPC information
  Bool_t           fProcessSPD;           //  Process SPD information
  Bool_t           fProcessVZERO;         //  Process VZERO information
  Bool_t           fProcessZDC;           //  Process ZDC information

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  Int_t            fEsdTracks;            //  N ESD tracks
  Int_t            fEsdTracksA;           //  N ESD tracks accepted
  Int_t            fTpcTracks;            //  N TPC tracks
  Int_t            fTpcTracksA;           //  N TPC tracks accepted

  Float_t          fVzeroMult;            //  VZERO multiplicity
  Float_t          fVzeroMultA;           //  VZERO A multiplicity
  Float_t          fVzeroMultC;           //  VZERO C multiplicity

  Float_t          fSpdNClusters;         //  Spd N clusters
  Float_t          fSpdNClustersInner;    //  Spd N clusters Inner 
  Float_t          fSpdNClustersOuter;    //  Spd N clusters Outer

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  /** Binnning VZERO */
  Int_t   fVzeroBinning;                  // VZERO Binning nbin
  Float_t fVzeroBinningMin;               // VZERO Binning min  
  Float_t fVzeroBinningMax;               // VZERO Binning max

  /** Binnning TPC */
  Int_t   fTpcBinning;                    // TPC Binning nbin
  Float_t fTpcBinningMin;                 // TPC Binning min
  Float_t fTpcBinningMax;                 // TPC Binning max

  /** Binnning ZDC */
  Int_t   fZdcBinning;                    // ZDC Binning nbin
  Float_t fZdcBinningMin;                 // ZDC Binning min
  Float_t fZdcBinningMax;                 // ZDC Binning max

  /** Binnning ZEM */
  Int_t   fZemBinning;                    // ZEM Binning nbin
  Float_t fZemBinningMin;                 // ZEM Binning min
  Float_t fZemBinningMax;                 // ZEM Binning may

  /** Binnning SPD */
  Int_t   fSpdBinning;                    // SPD Binning nbin
  Float_t fSpdBinningMin;                 // SPD Binning min
  Float_t fSpdBinningMax;                 // SPD Binning max
  
  ClassDef(AliMultiplicityCorrelations, 3);
};
#endif

//-*- Mode: C++ -*-
#ifndef ALIHLTJETBASE_H
#define ALIHLTJETBASE_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETBase.h
    @author Jochen Thaeder
    @date   
    @brief  Base functionality for HLT JET package

*/

#include "AliHLTLogging.h"
#include "TObject.h"

/** Indices in grid position array */
enum GridIndex_t { 
  kIdxPrimary,    /**< 1D index for the main search area */
  kIdxOutter,     /**< 1D index for the outter search area */
  kIdxEtaPrimary, /**< 2D eta index for the main search area */
  kIdxPhiPrimary, /**< 2D phi index for the main search area */
  kIdxPhiOutter   /**< 2D phi index for the outter search area */
};

/**  Indices in array */
enum EtaPhiIndex_t { 
  kIdxEta, /**< Eta */
  kIdxPhi, /**< Phi */
  kIdxPt   /**< Pt */
};

/** Used track types */
enum TrackType_t { 
  kTrackMC,   /**< TParticle */
  kTrackESD,  /**< AliESDtrack */
  kTrackAOD   /**< AliAODtrack */
};

/**
 * @class AliHLTJETBase
 * This class contains a Seed for the JetFinder
 * 
 */

class AliHLTJETBase : public TObject, public AliHLTLogging {

 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
 
  /** Standard constructor */
  AliHLTJETBase();

  /** Destructor */
  ~AliHLTJETBase();  

  /** Types of jet algorithms */
  enum JetAlgorithmType_t { 
    kAntiKt,          /**< FastJet implementation of the Anti kt */
    kKt,              /**< FastJet implementation of the kt  */
    kFFSCSquareCell,  /**< Fast Fixed Seeded Cone, using a square cell */
    kFFSCRadiusCell,  /**< Fast Fixed Seeded Cone, using a radius cell */
    kJetAlgorithmMax  /**< Number of enum entries */
  };
    
  /** Array of types of the Jet Algorithms */
  static const Char_t *fgkJetAlgorithmType[];        //! transient


  ///////////////////////////////////////////////////////////////////////////////////
  
 private:

  /** copy constructor prohibited */
  AliHLTJETBase(const AliHLTJETBase&);

  /** assignment operator prohibited */
  AliHLTJETBase& operator=(const AliHLTJETBase&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
    
  ClassDef(AliHLTJETBase, 0)
};
#endif


#if 0


static Float_t GetDistance2( const Float_t eta1, const Float_t phi1, 
			     const Float_t eta2, const Float_t phi2);

static Int_t GetCellIndex( const Double_t* aEtaPhi, Int_t* aGridIndex );

static void XYZtoRPhiEta( const Double_t *xyz, Double_t *rpe );

static void XYZEtoRPhiEtaPt( const Double_t *xyze, Double_t *rpep );
static void XYZEtoRPhiEtaPt( const Float_t *xyze, Double_t *rpep );

static void XYZEtoEPhiEtaPt( const Double_t *xyze, Double_t *epep );
static void XYZEtoEPhiEtaPt( const Float_t *xyze, Double_t *epep );

static Double_t GetPtFromXYZ( const Double_t *pxpypz );
static Double_t GetPtFromXYZ( const Float_t *pxpypz );

static Double_t GetPhiFromXYZ( const Double_t *xyz );
static Double_t GetPhiFromXYZ( const Float_t *xyz );

static Double_t GetEtaFromXYZ( const Double_t *xyz );
static Double_t GetEtaFromXYZ( const Float_t *xyz );

#endif

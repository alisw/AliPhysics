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

enum GridIndex_t   { kIdxPrimary, kIdxOutter, kIdxEtaPrimary, kIdxPhiPrimary, kIdxPhiOutter };
enum EtaPhiIndex_t { kIdxEta, kIdxPhi, kIdxPt };
enum TrackType_t   { kTrackMC, kTrackESD, kTrackAOD };

/**
XXXX STILL TRUE???
 Different type of how to add tracks to a jet
    0 : check for radius compared to every track 
    in selected cell -- default
    1 : check for radius compared to center 
    of selected cell
    2 : take whole cell
*/
enum FinderType_t { kSquareCellRegion, kRadiusCellRegion, kRadius, kSquareArea };

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
   */

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

  static Double_t GetDistance2( const Double_t eta1, const Double_t phi1, 
				const Double_t eta2, const Double_t phi2);


#endif
